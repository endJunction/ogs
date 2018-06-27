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

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
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


    auto p_GR_nodal = Eigen::Map<typename ShapeMatricesType::template VectorType<
        nonwet_pressure_size> const>(local_x.data() + nonwet_pressure_matrix_index, nonwet_pressure_size);

    auto p_cap_nodal = Eigen::Map<typename ShapeMatricesType::template VectorType<
        cap_pressure_size> const>(local_x.data() + cap_pressure_matrix_index, cap_pressure_size);

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

    const Eigen::MatrixXd& perm = _process_data.material->getPermeability(
        material_id, t, pos, _element.getDimension());
    assert(perm.rows() == _element.getDimension() || perm.rows() == 1);
    GlobalDimMatrixType permeability = GlobalDimMatrixType::Zero(
        _element.getDimension(), _element.getDimension());
    if (perm.rows() == _element.getDimension())
        permeability = perm;
    else if (perm.rows() == 1)
        permeability.diagonal().setConstant(perm(0, 0));

    _velocity_gas.clear();
    auto velocity_matrix_gas = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
                    _velocity_gas, GlobalDim, _ip_data.size());

    _velocity_liquid.clear();
    auto velocity_matrix_liquid = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
                    _velocity_liquid, GlobalDim, _ip_data.size());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        double pc_int_pt = 0.;
        double pn_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, _ip_data[ip].N, pn_int_pt,
                                         pc_int_pt);



        const double temperature = _process_data.temperature(t, pos)[0];
//        double const rho_nonwet =
//            _process_data.material->getGasDensity(pn_int_pt, temperature);
//        double const rho_wet = _process_data.material->getLiquidDensity(
//            _pressure_wet[ip], temperature);
//
//        double const Sw = _process_data.material->getSaturation(
//            material_id, t, pos, pn_int_pt, temperature, pc_int_pt);
//
//        _saturation[ip] = Sw;
//
//        double dSw_dpc = _process_data.material->getSaturationDerivative(
//            material_id, t, pos, pn_int_pt, temperature, Sw);
//
          double const porosity = _process_data.material->getPorosity(
              material_id, t, pos, pn_int_pt, temperature, 0);
//
//        // Assemble M matrix
//        // nonwetting
//        double const drhononwet_dpn =
//            _process_data.material->getGasDensityDerivative(pn_int_pt,
//                                                            temperature);
//



          // ---------------------------------------------------------------------

          double const parameter_a = _process_data.material->getLiquidDensity(
                  _pressure_wet[ip], temperature);

          double const parameter_n =
                  _process_data.material->getGasDensity(pn_int_pt, temperature);


          const double p_GR = pn_int_pt;
          const double p_cap = pc_int_pt;


         const double p_LR = p_GR - p_cap;


//         double const dRhoGRdpGR = 1.42857E-05;
//         double const dRhoLRdpLR = 3.57143E-07;

         double const rho_LR_0 = 995.;
         double const rho_GR_0 = 630.0;
         double const p_0 = 16000000;

 //        double const rho_LR = rho_LR_0 + dRhoLRdpLR * (p_LR - p_0);
 //        double const rho_GR = rho_GR_0 + dRhoGRdpGR * (p_GR - p_0);

         double const rho_LR = 1000.0;
         double const rho_GR = 1.0;

         double const dRhoGRdpGR = 0.0;
         double const dRhoLRdpLR = 0.0;

         const double s_L_res = 0.05;
         const double s_G_res = 0.05;
         const double s_L_max = 1 - s_G_res;
         const double l = 2.0;
         const double p_b = 5000.0;
         const double p_c_max = 1e9;
         const double e = 2.71828182845905;


         // Brooks-Corey
         const double s_e = std::pow(p_b/std::max(p_b, pc_int_pt), l);

//         // Fredlund
//         const double C = 1 - (std::log(1 + (pc_int_pt / p_b)) / std::log(1 + (p_c_max / p_b)));
//         const double s_e = C * 1 / (std::log(e + std::pow(pc_int_pt/parameter_a, parameter_n)));
//
         const double s_L = (1 - s_L_res - s_G_res)*s_e + s_L_res;
         const double dSEdPC = -l/pc_int_pt*s_L;
         const double dSLdPC = (1 - s_L_res - s_G_res) * dSEdPC;


         _density_gas[ip] = rho_GR;
         _density_liquid[ip] = rho_LR;

         _saturation[ip] = s_L;
         _pressure_wet[ip] = p_LR;



         auto const rho = 2000.0 * (1 - porosity) + porosity *
                 (s_L * rho_LR + (1. - s_L)* rho_GR);


         const double k_rel_GR = (1-s_e)*(1-s_e)*(1-std::pow(s_e,(2.0+l)/l));
         const double k_rel_LR = std::pow(s_e,(2.0+3.0*l)/l);

         const double mu_GR = 57.50623e-6;
         const double mu_LR = 0.51e-3;

         const double phi = porosity;
         const double s_G = 1.0 - s_L;


         const double beta_p_GR = 1./rho_GR*dRhoGRdpGR;
         const double beta_p_SR = 0.000001;
         const double beta_p_LR = 0.;

         const double lambda_GR = k_rel_GR / mu_GR;
         const double lambda_LR = k_rel_LR / mu_LR;


//         ---------------------------------------------------------------------



              Mgp.noalias() +=
                  porosity * (1 - s_L) * dRhoGRdpGR * _ip_data[ip].massOperator;
              Mgpc.noalias() +=
                  -porosity * rho_GR * dSLdPC * _ip_data[ip].massOperator;

              Mlpc.noalias() +=
                  porosity * dSLdPC * rho_LR * _ip_data[ip].massOperator;


              laplace_operator.noalias() = _ip_data[ip].dNdx.transpose() *
                                           permeability * _ip_data[ip].dNdx *
                                           _ip_data[ip].integration_weight;

              Kgp.noalias() += rho_GR * lambda_GR * laplace_operator;

              Klp.noalias() += rho_LR * lambda_LR * laplace_operator;
              Klpc.noalias() += -rho_LR * lambda_LR * laplace_operator;

//        ---------------------------------------------------------------------

//
//        Mgp.noalias() +=
//            porosity * (1 - Sw) * drhononwet_dpn * _ip_data[ip].massOperator;
//        Mgpc.noalias() +=
//            -porosity * rho_nonwet * dSw_dpc * _ip_data[ip].massOperator;
//
//        Mlpc.noalias() +=
//            porosity * dSw_dpc * rho_wet * _ip_data[ip].massOperator;
//
//        // nonwet
//        double const k_rel_nonwet =
//            _process_data.material->getNonwetRelativePermeability(
//                t, pos, pn_int_pt, temperature, Sw);
//        double const mu_nonwet =
//            _process_data.material->getGasViscosity(pn_int_pt, temperature);
//        double const lambda_nonwet = k_rel_nonwet / mu_nonwet;
//
//        // wet
//        double const k_rel_wet =
//            _process_data.material->getWetRelativePermeability(
//                t, pos, _pressure_wet[ip], temperature, Sw);
//        double const mu_wet = _process_data.material->getLiquidViscosity(
//            _pressure_wet[ip], temperature);
//        double const lambda_wet = k_rel_wet / mu_wet;
//
//        laplace_operator.noalias() = _ip_data[ip].dNdx.transpose() *
//                                     permeability * _ip_data[ip].dNdx *
//                                     _ip_data[ip].integration_weight;
//
//        Kgp.noalias() += rho_nonwet * lambda_nonwet * laplace_operator;
//
//        Klp.noalias() += rho_wet * lambda_wet * laplace_operator;
//        Klpc.noalias() += -rho_wet * lambda_wet * laplace_operator;

        if (_process_data.has_gravity)
       {

            //---------------------------------------------------------------------


             auto const& b = _process_data.specific_body_force;
             auto const& dNdx_T = _ip_data[ip].dNdx.transpose();
             auto const& K = permeability;
             auto const& w = _ip_data[ip].integration_weight;


             NodalVectorType gravity_operator = _ip_data[ip].dNdx.transpose() *
                                         permeability * b *
                                         _ip_data[ip].integration_weight;

             velocity_matrix_liquid.col(ip).noalias() =
                     -permeability*lambda_LR * _ip_data[ip].dNdx *
                     (p_GR_nodal-p_cap_nodal) + permeability*lambda_LR * rho_LR * b; // sign?

             velocity_matrix_gas.col(ip).noalias() =
                     -permeability*lambda_GR * _ip_data[ip].dNdx *
                     (p_GR_nodal) + permeability*lambda_GR * rho_GR * b; // sign?

             Bg.noalias() +=
                     rho_GR * rho_GR * lambda_GR * gravity_operator;
             Bl.noalias() += rho_LR * rho_LR * lambda_LR * gravity_operator;


             //---------------------------------------------------------------------


//
//            auto const& b = _process_data.specific_body_force;
//            auto const& dNdx_T = _ip_data[ip].dNdx.transpose();
//            auto const& K = permeability;
//            auto const& w = _ip_data[ip].integration_weight;
//
//            NodalVectorType gravity_operator = _ip_data[ip].dNdx.transpose() *
//                                        permeability * b *
//                                        _ip_data[ip].integration_weight;
//
//            Bg.noalias() +=
//                    rho_nonwet * rho_nonwet * lambda_nonwet * gravity_operator;
//            Bl.noalias() += rho_wet * rho_wet * lambda_wet * gravity_operator;
//



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
}  // end of namespace
