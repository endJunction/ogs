/**
* \copyright
* Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include "HeatTransportBHELocalAssemblerBHE.h"

#include <Eigen/Eigen>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ProcessLib/HeatTransportBHE/LocalAssemblers/IntegrationPointDataBHE.h"

namespace ProcessLib
{
    namespace HeatTransportBHE
    {
        template <typename ShapeFunction, typename IntegrationMethod,
            int BHE_Dim>
            HeatTransportBHELocalAssemblerBHE<ShapeFunction, IntegrationMethod,
            BHE_Dim>::
            HeatTransportBHELocalAssemblerBHE(
                MeshLib::Element const& e,
                std::size_t const /*local_matrix_size*/,
                std::vector<unsigned> const& dofIndex_to_localIndex,
                bool const is_axially_symmetric,
                unsigned const integration_order,
                HeatTransportBHEProcessData& process_data)
            : HeatTransportBHELocalAssemblerInterface(
                ShapeFunction::NPOINTS * BHE_Dim,  // no intersection
                dofIndex_to_localIndex),
            _process_data(process_data),
            _integration_method(integration_order),
            _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                IntegrationMethod, BHE_Dim>(e, is_axially_symmetric, _integration_method)),
            _element(e)
        {
            // need to make sure that the BHE elements are one-dimensional
            assert(_element.getDimension() == 1);

            unsigned const n_integration_points =
                _integration_method.getNumberOfPoints();

            _ip_data.reserve(n_integration_points);
            _secondary_data.N.resize(n_integration_points);

            auto mat_id = (*_process_data._mesh_prop_materialIDs)[e.getID()];
            auto BHE_id = _process_data._map_materialID_to_BHE_ID[mat_id];
            // _BHE_instance = _process_data._vec_BHE_property[BHE_id].get();

            
            SpatialPosition x_position;
            x_position.setElementID(_element.getID());
            for (unsigned ip = 0; ip < n_integration_points; ip++)
            {
                x_position.setIntegrationPoint(ip);

                IntegrationPointDataBHE int_Point_Data_BHE(*(_process_data._vec_BHE_property[BHE_id]));
                _ip_data.emplace_back(int_Point_Data_BHE);
                auto const& sm = _shape_matrices[ip];
                auto& ip_data = _ip_data[ip];
                ip_data.integration_weight =
                    _integration_method.getWeightedPoint(ip).getWeight() *
                    sm.integralMeasure * sm.detJ;
                /*
                // ip_data._h_matrices.setZero(DisplacementDim,
                //     ShapeFunction::NPOINTS * DisplacementDim);

                computeHMatrix<DisplacementDim, ShapeFunction::NPOINTS,
                    typename ShapeMatricesType::NodalRowVectorType,
                    HMatrixType>(sm.N, ip_data._h_matrices);
                */

                _secondary_data.N[ip] = sm.N;
            }
            
        }

        template<typename ShapeFunction, typename IntegrationMethod, int BHE_Dim>
            void HeatTransportBHELocalAssemblerBHE<ShapeFunction, IntegrationMethod, BHE_Dim>::
                assemble(double const t, std::vector<double> const& local_x,
                         std::vector<double>& local_M_data,
                         std::vector<double>& local_K_data,
                         std::vector<double>& local_b_data)
        {
            unsigned const n_integration_points =
                _integration_method.getNumberOfPoints();

            SpatialPosition x_position;
            x_position.setElementID(_element.getID());

            /*
            auto const& nodal_jump = local_u;

            auto const& R = _fracture_property->R;

            // the index of a normal (normal to a fracture plane) component
            // in a displacement vector
            int const index_normal = DisplacementDim - 1;

            

            for (unsigned ip = 0; ip < n_integration_points; ip++)
            {
            x_position.setIntegrationPoint(ip);

            auto& ip_data = _ip_data[ip];
            auto const& integration_weight = ip_data.integration_weight;
            auto const& H = ip_data._h_matrices;
            auto& mat = ip_data._fracture_material;
            auto& sigma = ip_data._sigma;
            auto const& sigma_prev = ip_data._sigma_prev;
            auto& w = ip_data._w;
            auto const& w_prev = ip_data._w_prev;
            auto& C = ip_data._C;
            auto& state = *ip_data._material_state_variables;

            // displacement jumps
            w.noalias() = R * H * nodal_jump;

            // total aperture
            ip_data._aperture = ip_data._aperture0 + w[index_normal];

            // local C, local stress
            mat.computeConstitutiveRelation(
            t, x_position, ip_data._aperture0,
            Eigen::Matrix<double, DisplacementDim, 1>::Zero(),  // TODO (naumov)
            // Replace with
            // initial
            // stress values
            w_prev, w, sigma_prev, sigma, C, state);

            // r_[u] += H^T*Stress
            local_b.noalias() -=
            H.transpose() * R.transpose() * sigma * integration_weight;

            // J_[u][u] += H^T*C*H
            local_J.noalias() +=
            H.transpose() * R.transpose() * C * R * H * integration_weight;
            }
            */
        }

        template <typename ShapeFunction, typename IntegrationMethod,
            int BHE_Dim>
            void HeatTransportBHELocalAssemblerBHE<ShapeFunction, IntegrationMethod,
            BHE_Dim>::
            postTimestepConcrete(std::vector<double> const& /*local_x*/)
        {
            // double ele_b = 0;
            unsigned const n_integration_points =
                _integration_method.getNumberOfPoints();
            for (unsigned ip = 0; ip < n_integration_points; ip++)
            {
                // ele_b += _ip_data[ip]._aperture;
            }
            // ele_b /= n_integration_points;
            // (*_process_data._mesh_prop_b)[_element.getID()] = ele_b;
        }

    }  // namespace HeatTransportBHE
}  // namespace ProcessLib
