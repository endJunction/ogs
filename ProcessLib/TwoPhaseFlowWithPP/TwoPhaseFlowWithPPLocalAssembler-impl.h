/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**
* common nomenclature
* --------------primary variable----------------------
* pn_int_pt    pressure for nonwetting phase at each integration point
* pc_int_pt    capillary pressure at each integration point
* --------------secondary variable--------------------
* temperature              capillary pressure
* Sw wetting               phase saturation
* dSw_dpc                  derivative of wetting phase saturation with respect
* to capillary pressure
* rho_nonwet               density of nonwetting phase
* drhononwet_dpn           derivative of nonwetting phase density with respect
*to nonwetting phase pressure
* rho_wet                  density of wetting phase
* k_rel_nonwet             relative permeability of nonwetting phase
* mu_nonwet                viscosity of nonwetting phase
* lambda_nonwet            mobility of nonwetting phase
* k_rel_wet                relative permeability of wetting phase
* mu_wet                   viscosity of wetting phase
* lambda_wet               mobility of wetting phase
*/
#pragma once

#include "TwoPhaseFlowWithPPLocalAssembler.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "TwoPhaseFlowWithPPProcessData.h"
#include <iostream>

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
namespace MPL = MaterialPropertyLib;

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void TwoPhaseFlowWithPPLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assemble(double const t, std::vector<double> const& local_x,
                         std::vector<double>& local_M_data,
                         std::vector<double>& local_K_data,
                         std::vector<double>& local_b_data)
{
    auto const local_matrix_size = local_x.size();

    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);



    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<LocalVectorType>(
        local_b_data, local_matrix_size);

    auto Mgpg =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mgpc = local_M.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);

    auto Mlpg = local_M.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Mlpc = local_M.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    NodalMatrixType laplace_operator =
        NodalMatrixType::Zero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kgpg =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Klpg = local_K.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Klpc = local_K.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    auto Bg = local_b.template segment<nonwet_pressure_size>(
        nonwet_pressure_matrix_index);

    auto Bl =
        local_b.template segment<cap_pressure_size>(cap_pressure_matrix_index);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id =
        _process_data.material->getMaterialID(pos.getElementID().get());

//    const Eigen::MatrixXd& perm = _process_data.material->getPermeability(
//        material_id, t, pos, _element.getDimension());
//
//    assert(perm.rows() == _element.getDimension() || perm.rows() == 1);
//    GlobalDimMatrixType permeability = GlobalDimMatrixType::Zero(
//        _element.getDimension(), _element.getDimension());
//    if (perm.rows() == _element.getDimension())
//        permeability = perm;
//    else if (perm.rows() == 1)
//        permeability.diagonal().setConstant(perm(0, 0));

    GlobalDimMatrixType permeability =
            GlobalDimMatrixType::Zero(_element.getDimension(),
                    _element.getDimension());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
    	auto& ip_data = _ip_data[ip];

    	auto const& w = ip_data.integration_weight;
    	auto const& Np =  ip_data.N;
    	auto const& NpT = Np.transpose().eval();
    	auto const& gradNp = ip_data.dNdx;
    	auto const& gradNpT = gradNp.transpose().eval();

        double p_GR = 0.;
        double p_cap = 0.;
        NumLib::shapeFunctionInterpolate(local_x, Np, p_GR, p_cap);

        const double p_LR = p_GR - p_cap;

        const double temperature = _process_data.temperature(t, pos)[0];

        MPL::VariableArray variables;

        variables[MPL::Variables::capillary_pressure] = p_cap;
        variables[MPL::Variables::temperature] = temperature;

        auto const& solid_phase = medium->phase(0);
        auto const& liquid_phase = medium->phase(1);
        auto const& gas_phase = medium->phase(2);

        double const beta_p_SR = MPL::getScalar(solid_phase.property(
                MPL::PropertyEnum::compressibility), variables);

        variables[MPL::Variables::phase_pressure] = p_GR;
        double const rho_GR = MPL::getScalar(gas_phase.property(
                MPL::PropertyEnum::density), variables);
        double const beta_p_GR = MPL::getScalar(gas_phase.property(
                MPL::PropertyEnum::compressibility), variables);
        auto const mu_GR = MPL::getScalar(gas_phase.property(
                MPL::PropertyEnum::viscosity),variables);

        variables[MPL::Variables::phase_pressure] = p_LR;
        double const rho_LR = MPL::getScalar(liquid_phase.property(
                MPL::PropertyEnum::density), variables);
        double const beta_p_LR = MPL::getScalar(liquid_phase.property(
                MPL::PropertyEnum::compressibility), variables);
        auto const mu_LR = MPL::getScalar(liquid_phase.property(
                MPL::PropertyEnum::viscosity),variables);

        auto const s_L = MPL::getScalar(medium->property(
                MPL::PropertyEnum::saturation), variables);

        auto const s_G = 1. - s_L;
        auto const dsL_dp_cap = MPL::getScalarDerivative(medium->property(
                MPL::PropertyEnum::saturation), variables,
                MPL::Variables::capillary_pressure);

        auto const phi = MPL::getScalar(
            medium->property(MPL::PropertyEnum::porosity), variables);
        auto const biot = MPL::getScalar(
            medium->property(MPL::PropertyEnum::biot_coefficient), variables);

        variables[MPL::Variables::liquid_saturation] = s_L;

        auto const k_rel_LR = std::max(1.0e-9, MPL::getPair(medium->property(
                MPL::PropertyEnum::relative_permeability), variables)[0]);
        auto const k_rel_GR = std::max(1.0e-9, MPL::getPair(medium->property(
                MPL::PropertyEnum::relative_permeability), variables)[1]);


        _saturation[ip] = s_L;
        _pressure_wet[ip] = p_LR;
        _density_gas[ip] = rho_GR;
        _density_liquid[ip] = rho_LR;
        _viscosity_gas[ip] = mu_GR;
        _viscosity_liquid[ip] = mu_LR;
        _rel_perm_gas[ip] = k_rel_GR;
        _rel_perm_liquid[ip] = k_rel_LR;

        double const perm =
            MPL::getScalar(medium->property(MPL::permeability));

        permeability.diagonal().setConstant(perm);

        auto const lambda_GR = k_rel_GR / mu_GR;
        auto const lambda_LR = k_rel_LR / mu_LR;
        auto const SpS = (biot - phi) * beta_p_SR;

#define DBG_OUTPUT
#ifndef DBG_OUTPUT
        std::cout << std::setprecision(16);
        std::cout << "=====================================\n";
        std::cout << "               p_GR : " << p_GR  << "\n";
        std::cout << "              p_cap : " << p_cap  << "\n";
        std::cout << "                  T : " << temperature  << "\n";
        std::cout << "-------------------------------------\n";
        std::cout << "             rho_GR : " << rho_GR << "\n";
        std::cout << "             rho_LR : " << rho_LR << "\n";
        std::cout << "          beta_p_GR : " << beta_p_GR << "\n";
        std::cout << "          beta_p_LR : " << beta_p_LR << "\n";
        std::cout << "          beta_p_SR : " << beta_p_SR << "\n";
        std::cout << "              mu_GR : " << mu_GR << "\n";
        std::cout << "              mu_LR : " << mu_LR << "\n";
        std::cout << "                 Sw : " << s_L << "\n";
        std::cout << "                 Sg : " << s_G << "\n";
        std::cout << "            dSw_dpc : " << dsL_dp_cap << "\n";
        std::cout << "                phi : " << phi << "\n";
        std::cout << "           k_rel_LR : " << k_rel_LR << "\n";
        std::cout << "           k_rel_GR : " << k_rel_GR << "\n";
        std::cout << "          lambda_GR : " << lambda_GR << "\n";
        std::cout << "          lambda_LR : " << lambda_LR << "\n";
        std::cout << "       permeability :\n" << permeability  << "\n";
        std::cout << "=====================================\n";
#endif

        // Assemble M matrix
        // nonwetting

        laplace_operator.noalias() = gradNpT * permeability * gradNp * w;


        Mgpg.noalias() +=
        		NpT * rho_GR * s_G * (phi * beta_p_GR + SpS) * Np * w;
        Mgpc.noalias() -=
        		NpT * rho_GR * (phi*dsL_dp_cap +
        				s_G * (s_L + p_cap * dsL_dp_cap) * SpS) * Np * w;
        Mlpg.noalias() +=
        		NpT * rho_LR * s_L * (phi * beta_p_LR + SpS) * Np * w;

        Mlpc.noalias() +=
        		NpT * rho_LR * (phi * dsL_dp_cap -
        				s_L * (s_L + p_cap * dsL_dp_cap) * SpS) * Np * w;


        Kgpg.noalias() += rho_GR * lambda_GR * laplace_operator;

        Klpg.noalias() += rho_LR * lambda_LR * laplace_operator;

        Klpc.noalias() -= rho_LR * lambda_LR * laplace_operator;

        if (_process_data.has_gravity)
        {
            auto const& b = _process_data.specific_body_force;

            NodalVectorType gravity_operator = gradNpT * permeability * b * w;
            Bg.noalias() += rho_GR * rho_GR * lambda_GR * gravity_operator;
            Bl.noalias() += rho_LR * rho_LR * lambda_LR * gravity_operator;
        }  // end of has gravity
    }
    if (_process_data.has_mass_lumping)
    {
        for (unsigned row = 0; row < Mgpc.cols(); row++)
        {
            for (unsigned column = 0; column < Mgpc.cols(); column++)
            {
                if (row != column)
                {
                    Mgpc(row, row) += Mgpc(row, column);
                    Mgpc(row, column) = 0.0;
                    Mgpg(row, row) += Mgpg(row, column);
                    Mgpg(row, column) = 0.0;
                    Mlpc(row, row) += Mlpc(row, column);
                    Mlpc(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping
}

}  // end of namespace
} // end of namespace
