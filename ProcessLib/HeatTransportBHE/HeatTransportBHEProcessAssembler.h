/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include <vector>

#include "HeatTransportBHEProcessData.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
    namespace HeatTransportBHE
    {
        const unsigned NUM_NODAL_DOF = 1;

        class HeatTransportBHELocalAssemblerInterface
            : public ProcessLib::LocalAssemblerInterface,
            public NumLib::ExtrapolatableElement
        {
        public:
            virtual std::vector<double> const& getIntPtHeatFluxX(
                const double /*t*/,
                GlobalVector const& /*current_solution*/,
                NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                std::vector<double>& /*cache*/) const = 0;

            virtual std::vector<double> const& getIntPtHeatFluxY(
                const double /*t*/,
                GlobalVector const& /*current_solution*/,
                NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                std::vector<double>& /*cache*/) const = 0;

            virtual std::vector<double> const& getIntPtHeatFluxZ(
                const double /*t*/,
                GlobalVector const& /*current_solution*/,
                NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                std::vector<double>& /*cache*/) const = 0;
        };

        template <typename ShapeFunction, typename IntegrationMethod,
            unsigned GlobalDim>
        class LocalAssemblerData : public HeatTransportBHELocalAssemblerInterface
        {
            using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
            using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

            using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
                ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

            using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
            using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
            using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

        public:
            /// The thermal_conductivity factor is directly integrated into the local
            /// element matrix.
            LocalAssemblerData(MeshLib::Element const& element,
                std::size_t const local_matrix_size,
                bool is_axially_symmetric,
                unsigned const integration_order,
                HeatTransportBHEProcessData const& process_data)
                : _element(element),
                _process_data(process_data),
                _integration_method(integration_order),
                _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                    IntegrationMethod, GlobalDim>(
                        element, is_axially_symmetric, _integration_method)),
                _heat_fluxes(
                    GlobalDim,
                    std::vector<double>(_integration_method.getNumberOfPoints()))
            {
                // This assertion is valid only if all nodal d.o.f. use the same shape
                // matrices.
                assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
                (void)local_matrix_size;
            }

            void assemble(double const t, std::vector<double> const& local_x,
                std::vector<double>& local_M_data,
                std::vector<double>& local_K_data,
                std::vector<double>& /*local_b_data*/) override
            {
                auto const local_matrix_size = local_x.size();
                // This assertion is valid only if all nodal d.o.f. use the same shape
                // matrices.
                assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

                auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
                    local_M_data, local_matrix_size, local_matrix_size);
                auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
                    local_K_data, local_matrix_size, local_matrix_size);

                unsigned const n_integration_points =
                    _integration_method.getNumberOfPoints();

                SpatialPosition pos;
                pos.setElementID(_element.getID());

                for (unsigned ip = 0; ip < n_integration_points; ip++)
                {
                    pos.setIntegrationPoint(ip);
                    auto const& sm = _shape_matrices[ip];
                    auto const& wp = _integration_method.getWeightedPoint(ip);
                    auto const k_s = _process_data.thermal_conductivity_solid(t, pos)[0];
                    auto const k_f = _process_data.thermal_conductivity_fluid(t, pos)[0];
                    auto const k_g = _process_data.thermal_conductivity_gas(t, pos)[0];
                    auto const heat_capacity_s = _process_data.heat_capacity_solid(t, pos)[0];
                    auto const heat_capacity_f = _process_data.heat_capacity_fluid(t, pos)[0];
                    auto const heat_capacity_g = _process_data.heat_capacity_gas(t, pos)[0];
                    auto const density_s = _process_data.density_solid(t, pos)[0];
                    auto const density_f = _process_data.density_fluid(t, pos)[0];
                    auto const density_g = _process_data.density_gas(t, pos)[0];

                    // TODO: k_s, heat_capacity_s and density_s to be changed...
                    local_K.noalias() += sm.dNdx.transpose() * k_s * sm.dNdx * sm.detJ *
                        wp.getWeight() * sm.integralMeasure;
                    local_M.noalias() += sm.N.transpose() * density_s * heat_capacity_s *
                        sm.N * sm.detJ * wp.getWeight() *
                        sm.integralMeasure;
                }
            }

            void computeSecondaryVariableConcrete(
                const double t,
                std::vector<double> const& local_x) override
            {
                auto const local_matrix_size = local_x.size();
                // This assertion is valid only if all nodal d.o.f. use the same shape
                // matrices.
                assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

                unsigned const n_integration_points =
                    _integration_method.getNumberOfPoints();

                SpatialPosition pos;
                pos.setElementID(_element.getID());
                const auto local_x_vec =
                    MathLib::toVector<NodalVectorType>(local_x, local_matrix_size);

                for (unsigned ip = 0; ip < n_integration_points; ip++)
                {
                    pos.setIntegrationPoint(ip);
                    auto const& sm = _shape_matrices[ip];
                    auto const k_s = _process_data.thermal_conductivity_solid(t, pos)[0];
                    auto const k_f = _process_data.thermal_conductivity_fluid(t, pos)[0];
                    auto const k_g = _process_data.thermal_conductivity_gas(t, pos)[0];

                    // TODO: here do the average of all thermal conductivity values

                    // heat flux only computed for output.
                    GlobalDimVectorType const heat_flux = -k_s * sm.dNdx * local_x_vec; // TODO: k_s to be changed to average values

                    for (unsigned d = 0; d < GlobalDim; ++d)
                    {
                        _heat_fluxes[d][ip] = heat_flux[d];
                    }
                }
            }

            Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
                const unsigned integration_point) const override
            {
                auto const& N = _shape_matrices[integration_point].N;

                // assumes N is stored contiguously in memory
                return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
            }

            std::vector<double> const& getIntPtHeatFluxX(
                const double /*t*/,
                GlobalVector const& /*current_solution*/,
                NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                std::vector<double>& /*cache*/) const override
            {
                assert(!_heat_fluxes.empty());
                return _heat_fluxes[0];
            }

            std::vector<double> const& getIntPtHeatFluxY(
                const double /*t*/,
                GlobalVector const& /*current_solution*/,
                NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                std::vector<double>& /*cache*/) const override
            {
                assert(_heat_fluxes.size() > 1);
                return _heat_fluxes[1];
            }

            std::vector<double> const& getIntPtHeatFluxZ(
                const double /*t*/,
                GlobalVector const& /*current_solution*/,
                NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                std::vector<double>& /*cache*/) const override
            {
                assert(_heat_fluxes.size() > 2);
                return _heat_fluxes[2];
            }

        private:
            MeshLib::Element const& _element;
            HeatTransportBHEProcessData const& _process_data;

            IntegrationMethod const _integration_method;
            std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
                _shape_matrices;

            std::vector<std::vector<double>> _heat_fluxes;
        };

    }  // namespace HeatTransportBHE
}  // namespace ProcessLib
