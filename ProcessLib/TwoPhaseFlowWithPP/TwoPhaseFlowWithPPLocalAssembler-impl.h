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

    auto Mgp =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mgpc = local_M.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);

    auto Mlpc = local_M.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    NodalMatrixType laplace_operator =
        NodalMatrixType::Zero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kgp =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Klp = local_K.template block<cap_pressure_size, nonwet_pressure_size>(
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
        double pc_int_pt = 0.;
        double pn_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, _ip_data[ip].N, pn_int_pt,
                                         pc_int_pt);

        const double p_GR = pn_int_pt;
        const double p_cap = pc_int_pt;

        const double p_LR = p_GR - p_cap;
        _pressure_wet[ip] = p_LR;

        const double temperature = _process_data.temperature(t, pos)[0];

        MPL::VariableArray variables;

        variables[MPL::Variables::capillary_pressure] = p_cap;
        variables[MPL::Variables::temperature] = temperature;

        auto const& solid_phase = medium->phase(0);
        auto const& liquid_phase = medium->phase(1);
        auto const& gas_phase = medium->phase(2);



//        double const rho_nonwet =
//            _process_data.material->getGasDensity(pn_int_pt, temperature);
//        double const rho_wet = _process_data.material->getLiquidDensity(
//            _pressure_wet[ip], temperature);

        variables[MPL::Variables::phase_pressure] = p_GR;
        double const rho_nonwet = MPL::getScalar(gas_phase.property(
                MPL::PropertyEnum::density), variables);
        auto const mu_nonwet = MPL::getScalar(gas_phase.property(
                MPL::PropertyEnum::viscosity),variables);

        variables[MPL::Variables::phase_pressure] = p_LR;
        double const rho_wet = MPL::getScalar(liquid_phase.property(
                MPL::PropertyEnum::density), variables);
        auto const mu_wet = MPL::getScalar(liquid_phase.property(
                MPL::PropertyEnum::viscosity),variables);

        auto const Sw = MPL::getScalar(medium->property(
                MPL::PropertyEnum::saturation), variables);
        auto const Sg = 1. - Sw;
        auto const dSw_dpc = MPL::getScalarDerivative(medium->property(
                MPL::PropertyEnum::saturation), variables,
                MPL::Variables::capillary_pressure);

        _saturation[ip] = Sw;
        _density_gas[ip] = rho_nonwet;
        _density_liquid[ip] = rho_wet;


        auto const porosity = MPL::getScalar(
            medium->property(MPL::PropertyEnum::porosity), variables);

        auto const drhononwet_dpn = MPL::getScalarDerivative(
                gas_phase.property(MPL::PropertyEnum::density), variables,
                MPL::Variables::phase_pressure);

        variables[MPL::Variables::liquid_saturation] = Sw;

        auto const k_rel_wet = MPL::getPair(medium->property(
                MPL::PropertyEnum::relative_permeability), variables)[0];
        auto const k_rel_nonwet = MPL::getPair(medium->property(
                MPL::PropertyEnum::relative_permeability), variables)[1];

        _rel_perm_gas[ip] = k_rel_nonwet;
        _rel_perm_liquid[ip] = k_rel_wet;

        double const perm =
            MPL::getScalar(medium->property(MPL::permeability));

        permeability.diagonal().setConstant(perm);

        double const lambda_nonwet = k_rel_nonwet / mu_nonwet;
        double const lambda_wet = k_rel_wet / mu_wet;

#define DBG_OUTPUT
#ifndef DBG_OUTPUT
        std::cout << std::setprecision(16);
        std::cout << "=====================================\n";
        std::cout << "               p_GR : " << p_GR  << "\n";
        std::cout << "              p_cap : " << p_cap  << "\n";
        std::cout << "                  T : " << temperature  << "\n";
        std::cout << "-------------------------------------\n";
        std::cout << "         rho_nonwet : " << rho_nonwet << "\n";
        std::cout << "            rho_wet : " << rho_wet << "\n";
        std::cout << "          mu_nonwet : " << mu_nonwet << "\n";
        std::cout << "             mu_wet : " << mu_wet << "\n";
        std::cout << "                 Sw : " << Sw << "\n";
        std::cout << "                 Sg : " << Sg << "\n";
        std::cout << "            dSw_dpc : " << dSw_dpc << "\n";
        std::cout << "           porosity : " << porosity << "\n";
        std::cout << "     drhononwet_dpn : " << drhononwet_dpn << "\n";
        std::cout << "          k_rel_wet : " << k_rel_wet << "\n";
        std::cout << "       k_rel_nonwet : " << k_rel_nonwet << "\n";
        std::cout << "      lambda_nonwet : " << lambda_nonwet << "\n";
        std::cout << "         lambda_wet : " << lambda_wet << "\n";
        std::cout << "       permeability :\n" << permeability  << "\n";
        std::cout << "=====================================\n";
#endif

        // Assemble M matrix
        // nonwetting

        Mgp.noalias() +=
            porosity * (1 - Sw) * drhononwet_dpn * _ip_data[ip].massOperator;
        Mgpc.noalias() +=
            -porosity * rho_nonwet * dSw_dpc * _ip_data[ip].massOperator;

        Mlpc.noalias() +=
            porosity * dSw_dpc * rho_wet * _ip_data[ip].massOperator;


        laplace_operator.noalias() = _ip_data[ip].dNdx.transpose() *
                                     permeability * _ip_data[ip].dNdx *
                                     _ip_data[ip].integration_weight;

        Kgp.noalias() += rho_nonwet * lambda_nonwet * laplace_operator;

        Klp.noalias() += rho_wet * lambda_wet * laplace_operator;
        Klpc.noalias() += -rho_wet * lambda_wet * laplace_operator;

        if (_process_data.has_gravity)
        {
            auto const& b = _process_data.specific_body_force;

            NodalVectorType gravity_operator = _ip_data[ip].dNdx.transpose() *
                                               permeability * b *
                                               _ip_data[ip].integration_weight;
            Bg.noalias() +=
                rho_nonwet * rho_nonwet * lambda_nonwet * gravity_operator;
            Bl.noalias() += rho_wet * rho_wet * lambda_wet * gravity_operator;
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
                    Mgp(row, row) += Mgp(row, column);
                    Mgp(row, column) = 0.0;
                    Mlpc(row, row) += Mlpc(row, column);
                    Mlpc(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping
}

}  // end of namespace
} // end of namespace
