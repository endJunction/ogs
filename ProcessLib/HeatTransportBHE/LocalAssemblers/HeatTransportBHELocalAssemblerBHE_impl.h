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
        using namespace BHE;

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

                IntegrationPointDataBHE int_Point_Data_BHE(*(_process_data._vec_BHE_property[BHE_id]));
                _ip_data.emplace_back(int_Point_Data_BHE);
                auto const& sm = _shape_matrices[ip];
                auto& ip_data = _ip_data[ip];
                ip_data.integration_weight =
                    _integration_method.getWeightedPoint(ip).getWeight() *
                    sm.integralMeasure * sm.detJ;

                _secondary_data.N[ip] = sm.N;
            }

            const int BHE_n_unknowns = _ip_data[0]._bhe_instance.get_n_unknowns();
            _R_matrix.setZero(nnodes * BHE_n_unknowns, nnodes * BHE_n_unknowns);
            _R_pi_s_matrix.setZero(nnodes * BHE_n_unknowns, nnodes);
            _R_s_matrix.setZero(nnodes, nnodes);
            // formulate the local BHE R matrix
            Eigen::MatrixXd matBHE_loc_R = Eigen::MatrixXd::Zero(nnodes, nnodes);
            for (int idx_bhe_unknowns = 0; idx_bhe_unknowns < BHE_n_unknowns; idx_bhe_unknowns++)
            {
                matBHE_loc_R.setZero();
                // Loop over Gauss points
                for (unsigned ip = 0; ip < n_integration_points; ip++)
                {
                    x_position.setIntegrationPoint(ip);
                    auto& ip_data = _ip_data[ip];

                    auto const& integration_weight = ip_data.integration_weight;
                    auto const& sm = _shape_matrices[ip];
                    auto const& wp = _integration_method.getWeightedPoint(ip);

                    // get coefficient of R matrix for corresponding BHE. 
                    auto R_coeff = _process_data._vec_BHE_property[BHE_id]->get_boundary_heat_exchange_coeff(idx_bhe_unknowns);

                    // calculate mass matrix for current unknown
                    matBHE_loc_R += sm.N.transpose() * R_coeff *
                        sm.N * sm.detJ * wp.getWeight() * sm.integralMeasure;
                }  // end of loop over integration point

                // The following assembly action is according to Diersch (2013) FEFLOW book
                // please refer to M.127 and M.128 on page 955 and 956
                switch (_process_data._vec_BHE_property[BHE_id]->get_type())
                {
                case BHE::BHE_TYPE::TYPE_1U:
                    switch (idx_bhe_unknowns)
                    {
                    case 0:  // PHI_fig
                        _R_matrix.block(0, 2 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(2 * nnodes, 0, nnodes, nnodes) += -1.0 * matBHE_loc_R;

                        _R_matrix.block(0, 0, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_i1
                        _R_matrix.block(2 * nnodes, 2 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_ig
                        break;
                    case 1:  // PHI_fog
                        _R_matrix.block(nnodes, 3 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(3 * nnodes, nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;

                        _R_matrix.block(nnodes, nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_o1
                        _R_matrix.block(3 * nnodes, 3 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_og
                        break;
                    case 2:  // PHI_gg
                        _R_matrix.block(2 * nnodes, 3 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(3 * nnodes, 2 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;

                        _R_matrix.block(2 * nnodes, 2 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R;  // K_ig  // notice we only have 1 PHI_gg term here. 
                        _R_matrix.block(3 * nnodes, 3 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R;  // K_og  // see Diersch 2013 FEFLOW book page 954 Table M.2
                        break;
                    case 3:  // PHI_gs
                        _R_s_matrix += 1.0 * matBHE_loc_R;

                        _R_pi_s_matrix.block(2 * nnodes, 0, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_pi_s_matrix.block(3 * nnodes, 0, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(2 * nnodes, 2 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R;  // K_ig
                        _R_matrix.block(3 * nnodes, 3 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R;  // K_og
                        break;
                    }
                    break;
                case BHE::BHE_TYPE::TYPE_2U:
                    switch (idx_bhe_unknowns)
                    {
                    case 0:  // R i1 i2
                        _R_matrix.block(0, 4 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(4 * nnodes, 0, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(nnodes, 5 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(5 * nnodes, nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;

                        _R_matrix.block(0, 0, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_i1
                        _R_matrix.block(nnodes, nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_i2
                        _R_matrix.block(4 * nnodes, 4 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_ig
                        _R_matrix.block(5 * nnodes, 5 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_ig
                        break;
                    case 1:  // R o1 o2
                        _R_matrix.block(2 * nnodes, 6 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(6 * nnodes, 2 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(3 * nnodes, 7 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(7 * nnodes, 3 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;

                        _R_matrix.block(2 * nnodes, 2 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_o1
                        _R_matrix.block(3 * nnodes, 3 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_o2
                        _R_matrix.block(6 * nnodes, 6 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_og
                        _R_matrix.block(7 * nnodes, 7 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_og
                        break;
                    case 2:  // R g1
                        _R_matrix.block(4 * nnodes, 6 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(6 * nnodes, 4 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(4 * nnodes, 7 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(7 * nnodes, 4 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(5 * nnodes, 6 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(6 * nnodes, 5 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(5 * nnodes, 7 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(7 * nnodes, 5 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;

                        _R_matrix.block(4 * nnodes, 4 * nnodes, nnodes, nnodes) += 2.0 * matBHE_loc_R; // K_ig
                        _R_matrix.block(5 * nnodes, 5 * nnodes, nnodes, nnodes) += 2.0 * matBHE_loc_R; // K_ig
                        _R_matrix.block(6 * nnodes, 6 * nnodes, nnodes, nnodes) += 2.0 * matBHE_loc_R; // K_og
                        _R_matrix.block(7 * nnodes, 7 * nnodes, nnodes, nnodes) += 2.0 * matBHE_loc_R; // K_og
                        break;
                    case 3:  // R g2
                        _R_matrix.block(6 * nnodes, 7 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(7 * nnodes, 6 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;

                        _R_matrix.block(4 * nnodes, 4 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_ig
                        _R_matrix.block(5 * nnodes, 5 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_ig
                        _R_matrix.block(6 * nnodes, 6 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_og
                        _R_matrix.block(7 * nnodes, 7 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_og
                        break;
                    case 4:  // R s
                        _R_s_matrix += 1.0 * matBHE_loc_R;

                        _R_pi_s_matrix.block(4 * nnodes, 0, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_pi_s_matrix.block(5 * nnodes, 0, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_pi_s_matrix.block(6 * nnodes, 0, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_pi_s_matrix.block(7 * nnodes, 0, nnodes, nnodes) += -1.0 * matBHE_loc_R;

                        _R_matrix.block(4 * nnodes, 4 * nnodes, nnodes, nnodes) += matBHE_loc_R; // K_gs
                        _R_matrix.block(5 * nnodes, 5 * nnodes, nnodes, nnodes) += matBHE_loc_R; // K_gs
                        _R_matrix.block(6 * nnodes, 6 * nnodes, nnodes, nnodes) += matBHE_loc_R; // K_gs
                        _R_matrix.block(7 * nnodes, 7 * nnodes, nnodes, nnodes) += matBHE_loc_R; // K_gs
                        break;
                    }
                    break;
                case BHE::BHE_TYPE::TYPE_CXA:
                    switch (idx_bhe_unknowns)
                    {
                    case 0:  // R i1
                        _R_matrix.block(0, 2 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(2 * nnodes, 0, nnodes, nnodes) += -1.0 * matBHE_loc_R;

                        _R_matrix.block(0, 0, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_i1
                        _R_matrix.block(2 * nnodes, 2 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_ig
                        break;
                    case 1:  // R io
                        _R_matrix.block(0, nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(nnodes, 0, nnodes, nnodes) += -1.0 * matBHE_loc_R;

                        _R_matrix.block(0, 0, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_i1
                        _R_matrix.block(nnodes, nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_o1
                        break;
                    case 2:  // R s
                        _R_s_matrix += matBHE_loc_R;

                        _R_pi_s_matrix.block(2 * nnodes, 0, nnodes, nnodes) += -1.0 * matBHE_loc_R;

                        _R_matrix.block(2 * nnodes, 2 * nnodes, nnodes, nnodes) += matBHE_loc_R; // K_gs
                        break;
                    }
                    break;
                case BHE::BHE_TYPE::TYPE_CXC:
                    switch (idx_bhe_unknowns)
                    {
                    case 0:  // R o1
                        _R_matrix.block(0, 2 * nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(2 * nnodes, 0, nnodes, nnodes) += -1.0 * matBHE_loc_R;

                        _R_matrix.block(nnodes, nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_o1
                        _R_matrix.block(2 * nnodes, 2 * nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_og
                        break;
                    case 1:  // R io
                        _R_matrix.block(0, nnodes, nnodes, nnodes) += -1.0 * matBHE_loc_R;
                        _R_matrix.block(nnodes, 0, nnodes, nnodes) += -1.0 * matBHE_loc_R;

                        _R_matrix.block(0, 0, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_i1
                        _R_matrix.block(nnodes, nnodes, nnodes, nnodes) += 1.0 * matBHE_loc_R; // K_o1
                        break;
                    case 2:  // R s
                        _R_s_matrix += matBHE_loc_R;

                        _R_pi_s_matrix.block(2 * nnodes, 0, nnodes, nnodes) += -1.0 * matBHE_loc_R;

                        _R_matrix.block(2 * nnodes, 2 * nnodes, nnodes, nnodes) += matBHE_loc_R; // K_gs
                        break;
                    }
                    break;
                }

            }  // end of loop over BHE unknowns

            // debugging
            // std::string sep = "\n----------------------------------------\n";
            // Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
            // std::cout << "_R_matrix: \n" << sep; 
            // std::cout << _R_matrix.format(CleanFmt) << sep;
            // std::cout << "_R_s_matrix: \n" << sep;
            // std::cout << _R_s_matrix.format(CleanFmt) << sep;
            // std::cout << "_R_pi_s_matrix: \n" << sep;
            // std::cout << _R_pi_s_matrix.format(CleanFmt) << sep;
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
                const int BHE_n_unknowns = _ip_data[0]._bhe_instance.get_n_unknowns();
                // plus one because the soil temperature is included in local_x
                assert(local_matrix_size ==
                       ShapeFunction::NPOINTS * (BHE_n_unknowns + 1));
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
                const int local_BHE_matrix_size =
                    ShapeFunction::NPOINTS * BHE_n_unknowns;
                const int shift_start =
                    local_matrix_size - local_BHE_matrix_size;

                // the mass and conductance matrix terms
                for (unsigned ip = 0; ip < n_integration_points; ip++)
                {
                    x_position.setIntegrationPoint(ip);
                    auto& ip_data = _ip_data[ip];

                    auto const& integration_weight = ip_data.integration_weight;
                    auto const& sm = _shape_matrices[ip];
                    auto const& wp = _integration_method.getWeightedPoint(ip);

                    // looping over all unknowns. 
                    for (int idx_bhe_unknowns = 0; idx_bhe_unknowns < BHE_n_unknowns; idx_bhe_unknowns++)
                    {
                        // get coefficient of mass from corresponding BHE. 
                        auto& mass_coeff = ip_data._vec_mass_coefficients[idx_bhe_unknowns];
                        auto& laplace_mat = ip_data._vec_mat_Laplace[idx_bhe_unknowns];
                        auto& advection_vec = ip_data._vec_Advection_vectors[idx_bhe_unknowns];

                        // calculate shift.
                        shift = shift_start + nnodes * idx_bhe_unknowns;
                        // local M
                        local_M.block(shift, shift, nnodes, nnodes).noalias() += sm.N.transpose() * mass_coeff *
                            sm.N * sm.detJ * wp.getWeight() * sm.integralMeasure;

                        // local K
                        // laplace part
                        local_K.block(shift, shift, nnodes, nnodes).noalias() += sm.dNdx.transpose() * laplace_mat * sm.dNdx * sm.detJ *
                            wp.getWeight() * sm.integralMeasure;
                        // advection part
                        local_K.block(shift, shift, nnodes, nnodes).noalias() += sm.N.transpose() * advection_vec.transpose() * sm.dNdx * sm.detJ *
                            wp.getWeight() * sm.integralMeasure;

                    }

                }

                // add the R matrix to local_K
                local_K.block(shift_start, shift_start, local_BHE_matrix_size,
                              local_BHE_matrix_size) += _R_matrix;

                // add the R_pi_s matrix to local_K
                local_K.block(shift_start, 0, local_BHE_matrix_size,
                              shift_start) += _R_pi_s_matrix;
                local_K.block(0, shift_start, shift_start,
                              local_BHE_matrix_size) +=
                    _R_pi_s_matrix.transpose();

                // add the R_s matrix to local_K
                local_K.block(0, 0, shift_start, shift_start) +=
                    2.0 * _R_s_matrix;

                // debugging
                // std::string sep =
                //     "\n----------------------------------------\n";
                // Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
                // std::cout << local_K.format(CleanFmt) << sep;
                // std::cout << local_M.format(CleanFmt) << sep;
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
