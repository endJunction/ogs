/**
* \copyright
* Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include "HeatTransportBHELocalAssemblerSoilNearBHE.h"

#include <valarray>
#include <vector>

#include <Eigen/Eigen>

#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MathLib/Point3d.h"

#include "MeshLib/Node.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ProcessLib/LIE/Common/LevelSetFunction.h"
#include "ProcessLib/LIE/Common/Utils.h"

#include "HeatTransportBHEProcessAssemblerInterface.h"
#include "ProcessLib/HeatTransportBHE/LocalAssemblers/IntegrationPointDataSoilNearBHE.h"
#include "SecondaryData.h"

namespace ProcessLib
{
    namespace HeatTransportBHE
    {
        template <typename ShapeFunction,
            typename IntegrationMethod,
            int GlobalDim>
            HeatTransportBHELocalAssemblerSoilNearBHE<ShapeFunction,
            IntegrationMethod,
            GlobalDim>::
            HeatTransportBHELocalAssemblerSoilNearBHE(
                MeshLib::Element const& e,
                std::size_t const n_variables,
                std::size_t const /*local_matrix_size*/,
                std::vector<unsigned> const& dofIndex_to_localIndex,
                bool const is_axially_symmetric,
                unsigned const integration_order,
                HeatTransportBHEProcessData& process_data)
            : HeatTransportBHELocalAssemblerInterface(
                n_variables * ShapeFunction::NPOINTS * GlobalDim,
                dofIndex_to_localIndex),
            _process_data(process_data),
            _integration_method(integration_order),
            _element(e),
            _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                IntegrationMethod, GlobalDim>(e, is_axially_symmetric, _integration_method)),
            _is_axially_symmetric(is_axially_symmetric)
        {
            // need to make sure that the soil elements are 3-dimensional
            assert(_element.getDimension() == 3);

            unsigned const n_integration_points =
                _integration_method.getNumberOfPoints();
            const int nnodes = _element.getNumberOfNodes();

            _ip_data.reserve(n_integration_points);
            _secondary_data.N.resize(n_integration_points);

            auto mat_id = (*_process_data._mesh_prop_materialIDs)[e.getID()];
            auto BHE_id = _process_data._map_materialID_to_BHE_ID[mat_id];

            SpatialPosition x_position;
            x_position.setElementID(_element.getID());

            // ip data initialization
            for (unsigned ip = 0; ip < n_integration_points; ip++)
            {
                x_position.setIntegrationPoint(ip);

                IntegrationPointDataSoilNearBHE<ShapeMatricesType,
                                                BMatricesType, GlobalDim>
                    int_Point_Data_BHE(
                        *(_process_data._vec_BHE_property[BHE_id]));

                _ip_data.emplace_back(int_Point_Data_BHE);
                auto const& sm = _shape_matrices[ip];
                auto& ip_data = _ip_data[ip];
                ip_data.integration_weight =
                    _integration_method.getWeightedPoint(ip).getWeight() *
                    sm.integralMeasure * sm.detJ;

                _secondary_data.N[ip] = sm.N;
            }

            const int BHE_n_unknowns =
                _ip_data[0]._bhe_instance.get_n_unknowns();
        }

        template<typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
            void HeatTransportBHELocalAssemblerSoilNearBHE<ShapeFunction, IntegrationMethod, GlobalDim>::
                assemble(double const t, std::vector<double> const& local_x,
                    std::vector<double>& local_M_data,
                    std::vector<double>& local_K_data,
                    std::vector<double>& local_b_data)
            {
                assert(_element.getDimension() == GlobalDim);

                auto const local_matrix_size = local_x.size();

                const int BHE_n_unknowns =
                    _ip_data[0]._bhe_instance.get_n_unknowns();

                assert(local_matrix_size ==
                       ShapeFunction::NPOINTS *
                           (BHE_n_unknowns +
                            1) /*plus 1 becasue of soil temperature*/);

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

                    auto const k_f = _process_data.thermal_conductivity_fluid(t, pos)[0];
                    auto const k_g = _process_data.thermal_conductivity_gas(t, pos)[0];
                    auto const k_s = _process_data.thermal_conductivity_solid(t, pos)[0];

                    auto const heat_capacity_f = _process_data.heat_capacity_fluid(t, pos)[0];
                    auto const heat_capacity_g = _process_data.heat_capacity_gas(t, pos)[0];
                    auto const heat_capacity_s = _process_data.heat_capacity_solid(t, pos)[0];

                    auto const density_f = _process_data.density_fluid(t, pos)[0];
                    auto const density_g = _process_data.density_gas(t, pos)[0];
                    auto const density_s = _process_data.density_solid(t, pos)[0];

                    // for now only using the solid phase parameters

                    // assemble Conductance matrix, just the soil part
                    local_K
                        .block<ShapeFunction::NPOINTS, ShapeFunction::NPOINTS>(
                            0, 0)
                        .noalias() += sm.dNdx.transpose() * k_s * sm.dNdx *
                                      sm.detJ * wp.getWeight() *
                                      sm.integralMeasure;

                    // assemble Mass matrix, just the soil part
                    local_M
                        .block<ShapeFunction::NPOINTS, ShapeFunction::NPOINTS>(
                            0, 0)
                        .noalias() += sm.N.transpose() * density_s *
                                      heat_capacity_s * sm.N * sm.detJ *
                                      wp.getWeight() * sm.integralMeasure;
                }

                // debugging
                std::string sep = "\n----------------------------------------\n";
                Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
                std::cout << local_K.format(CleanFmt) << sep;
                std::cout << local_M.format(CleanFmt) << sep;
            }

        template <typename ShapeFunction,
            typename IntegrationMethod,
            int GlobalDim>
            void HeatTransportBHELocalAssemblerSoilNearBHE<ShapeFunction,
            IntegrationMethod,
            GlobalDim>::
            postTimestepConcrete(std::vector<double> const& /*local_x*/)
        {

        }

    }  // namespace HeatTransportBHE
}  // namespace ProcessLib
