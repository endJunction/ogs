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
            // BHE must be assembled with one dimentional element
            assert(_element.getDimension() == 1);
            auto const local_matrix_size = local_x.size();
            const int BHE_n_unknowns =
                _ip_data[0]._bhe_instance.get_n_unknowns();
            assert(local_matrix_size ==
                   ShapeFunction::NPOINTS * BHE_n_unknowns);
            const int nnodes = _element.getNumberOfNodes();

            auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
                local_M_data, local_matrix_size, local_matrix_size);
            auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
                local_K_data, local_matrix_size, local_matrix_size);

            unsigned const n_integration_points =
                _integration_method.getNumberOfPoints();

            SpatialPosition x_position;
            x_position.setElementID(_element.getID());

            int shift = 0;

            for (unsigned ip = 0; ip < n_integration_points; ip++)
            {
                x_position.setIntegrationPoint(ip);
                auto& ip_data = _ip_data[ip];

                auto const& integration_weight = ip_data.integration_weight;
                auto const& sm = _shape_matrices[ip];
                auto const& wp = _integration_method.getWeightedPoint(ip);

                // looping over all unknowns.
                for (int idx_bhe_unknowns = 0;
                     idx_bhe_unknowns < BHE_n_unknowns; idx_bhe_unknowns++)
                {
                    // get coefficient of mass from corresponding BHE.
                    auto& mass_coeff =
                        ip_data._vec_mass_coefficients[idx_bhe_unknowns];
                    auto& laplace_mat =
                        ip_data._vec_mat_Laplace[idx_bhe_unknowns];
                    auto& advection_vec =
                        ip_data._vec_Advection_vectors[idx_bhe_unknowns];

                    // calculate shift.
                    shift = nnodes * idx_bhe_unknowns;
                    // local M
                    local_M.block(shift, shift, nnodes, nnodes).noalias() +=
                        sm.N.transpose() * mass_coeff * sm.N * sm.detJ *
                        wp.getWeight() * sm.integralMeasure;

                    // local K
                    // laplace part
                    local_K.block(shift, shift, nnodes, nnodes).noalias() +=
                        sm.dNdx.transpose() * laplace_mat * sm.dNdx * sm.detJ *
                        wp.getWeight() * sm.integralMeasure;
                    // advection part
                    local_K.block(shift, shift, nnodes, nnodes).noalias() +=
                        sm.N.transpose() * advection_vec.transpose() * sm.dNdx *
                        sm.detJ * wp.getWeight() * sm.integralMeasure;
                }
            }

            // debugging
            std::string sep = "\n----------------------------------------\n";
            Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
            std::cout << local_K.format(CleanFmt) << sep;
            std::cout << local_M.format(CleanFmt) << sep;
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
