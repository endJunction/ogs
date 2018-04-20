/**
* \copyright
* Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include "HeatTransportBHELocalAssemblerSoil.h"

#include <valarray>
#include <vector>

#include <Eigen/Eigen>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ProcessLib/HeatTransportBHE/HeatTransportBHEProcessData.h"

// #include "IntegrationPointDataMatrix.h"
#include "SecondaryData.h"
#include "HeatTransportBHEProcessAssemblerInterface.h"

namespace ProcessLib
{
    namespace HeatTransportBHE
    {
    template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
    HeatTransportBHELocalAssemblerSoil<ShapeFunction,
                                       IntegrationMethod,
                                       GlobalDim>::
        HeatTransportBHELocalAssemblerSoil(
            MeshLib::Element const& e,
            std::size_t const /*local_matrix_size*/,
            bool const is_axially_symmetric,
            unsigned const integration_order,
            HeatTransportBHEProcessData& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _shape_matrices(initShapeMatrices<ShapeFunction,
                                            ShapeMatricesType,
                                            IntegrationMethod,
                                            GlobalDim>(
              e, is_axially_symmetric, _integration_method)),
          _is_axially_symmetric(is_axially_symmetric)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        /*
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            
            _ip_data.emplace_back(*_process_data._material);
            auto& ip_data = _ip_data[ip];
            auto const& sm = shape_matrices[ip];
            ip_data.N = sm.N;
            ip_data.dNdx = sm.dNdx;
            ip_data.integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                sm.integralMeasure * sm.detJ;
            
            // Initialize current time step values
            static const int kelvin_vector_size =
                MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
            ip_data._sigma.setZero(kelvin_vector_size);
            ip_data._eps.setZero(kelvin_vector_size);

            // Previous time step values are not initialized and are set later.
            ip_data._sigma_prev.resize(kelvin_vector_size);
            ip_data._eps_prev.resize(kelvin_vector_size);

            ip_data._C.resize(kelvin_vector_size, kelvin_vector_size);

            _secondary_data.N[ip] = sm.N;
            
        }
        */
        }

        template<typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
            void HeatTransportBHELocalAssemblerSoil<ShapeFunction, IntegrationMethod, GlobalDim>::
                assemble(double const t,
                    std::vector<double> const& local_x,
                    std::vector<double>& local_M_data,
                    std::vector<double>& local_K_data,
                    std::vector<double>& local_b_data)
        {
            assert(_element.getDimension() == GlobalDim);

            auto const local_matrix_size = local_x.size();

            assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF_SOIL);

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

                // assemble Conductance matrix
                local_K.noalias() += sm.dNdx.transpose() * k_s * sm.dNdx * sm.detJ *
                    wp.getWeight() * sm.integralMeasure;

                // assemble Mass matrix
                local_M.noalias() += sm.N.transpose() * density_s * heat_capacity_s *
                    sm.N * sm.detJ * wp.getWeight() *
                    sm.integralMeasure;
            }
        }

        template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
            void HeatTransportBHELocalAssemblerSoil<ShapeFunction, IntegrationMethod, GlobalDim>::
            postTimestepConcrete(std::vector<double> const& /*local_x*/)
        {

        }

    }  // namespace HeatTransportBHE
}  // namespace ProcessLib
