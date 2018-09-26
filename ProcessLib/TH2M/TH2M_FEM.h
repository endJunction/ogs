/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "IntegrationPointData.h"
#include "LocalAssemblerInterface.h"
#include "TH2MProcessData.h"

#include <iostream>

namespace ProcessLib
{
namespace TH2M
{
namespace MPL = MaterialPropertyLib;

/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N_u;
};

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
class TH2MLocalAssembler : public LocalAssemblerInterface
{
public:
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, DisplacementDim>;

    using GlobalDimMatrixType =
        typename ShapeMatricesTypePressure::GlobalDimMatrixType;

    TH2MLocalAssembler(TH2MLocalAssembler const&) = delete;
    TH2MLocalAssembler(TH2MLocalAssembler&&) = delete;

    TH2MLocalAssembler(MeshLib::Element const& e,
                       std::size_t const /*local_matrix_size*/,
                       bool const is_axially_symmetric,
                       unsigned const integration_order,
                       TH2MProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _is_axially_symmetric(is_axially_symmetric)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N_u.resize(n_integration_points);

        auto const shape_matrices_u =
            initShapeMatrices<ShapeFunctionDisplacement,
                              ShapeMatricesTypeDisplacement, IntegrationMethod,
                              DisplacementDim>(e, is_axially_symmetric,
                                               _integration_method);

        auto const shape_matrices_p =
            initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            // displacement (subscript u)
            _ip_data.emplace_back(*_process_data.material);
            auto& ip_data = _ip_data[ip];
            auto const& sm_u = shape_matrices_u[ip];
            ip_data.integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                sm_u.integralMeasure * sm_u.detJ;

            // Initialize current time step values
            ip_data.sigma_eff.setZero(kelvin_vector_size);
            ip_data.eps.setZero(kelvin_vector_size);

            // Previous time step values are not initialized and are set later.
            ip_data.eps_prev.resize(kelvin_vector_size);
            ip_data.sigma_eff_prev.resize(kelvin_vector_size);

            ip_data.N_u_op = ShapeMatricesTypeDisplacement::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

            for (int i = 0; i < DisplacementDim; ++i)
                ip_data.N_u_op
                    .template block<1, displacement_size / DisplacementDim>(
                        i, i * displacement_size / DisplacementDim)
                    .noalias() = sm_u.N;

            ip_data.N_u = sm_u.N;
            ip_data.dNdx_u = sm_u.dNdx;

            ip_data.N_p = shape_matrices_p[ip].N;
            ip_data.dNdx_p = shape_matrices_p[ip].dNdx;

            _secondary_data.N_u[ip] = shape_matrices_u[ip].N;
        }
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_rhs_data) override
    {
        assert(local_x.size() == gas_pressure_size + cap_pressure_size +
                                     temperature_size + displacement_size);

        const auto matrix_size = gas_pressure_size + cap_pressure_size +
                                 temperature_size + displacement_size;

        // primary variables
        auto gas_phase_pressure =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                          gas_pressure_size);

        auto capillary_pressure =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                cap_pressure_size> const>(local_x.data() + cap_pressure_index,
                                          cap_pressure_size);

        auto temperature =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                temperature_size> const>(local_x.data() + temperature_index,
                                         temperature_size);

        auto displacement =
            Eigen::Map<typename ShapeMatricesTypeDisplacement::
                           template VectorType<displacement_size> const>(
                local_x.data() + displacement_index, displacement_size);

        // create mass matrix
        auto local_M = MathLib::createZeroedMatrix<
            typename ShapeMatricesTypeDisplacement::template MatrixType<
                matrix_size, matrix_size>>(local_M_data, matrix_size,
                                           matrix_size);

        // create stiffness matrix:
        auto local_K = MathLib::createZeroedMatrix<
            typename ShapeMatricesTypeDisplacement::template MatrixType<
                matrix_size, matrix_size>>(local_K_data, matrix_size,
                                           matrix_size);

        // create rhs-vector:
        auto local_rhs =
            MathLib::createZeroedVector<typename ShapeMatricesTypeDisplacement::
                                            template VectorType<matrix_size>>(
                local_rhs_data, matrix_size);

        //         ********************************************************************
        //         ********************************************************************

        auto Mgpg =
            local_M.template block<gas_pressure_size, gas_pressure_size>(
                gas_pressure_index, gas_pressure_index);

        auto Mgpc =
            local_M.template block<gas_pressure_size, cap_pressure_size>(
                gas_pressure_index, cap_pressure_index);

        auto MgT = local_M.template block<gas_pressure_size, temperature_size>(
            gas_pressure_index, temperature_index);

        auto Mgus =
            local_M.template block<gas_pressure_size, displacement_size>(
                gas_pressure_index, displacement_index);

        auto Mlpg =
            local_M.template block<cap_pressure_size, gas_pressure_size>(
                cap_pressure_index, gas_pressure_index);

        auto Mlpc =
            local_M.template block<cap_pressure_size, cap_pressure_size>(
                cap_pressure_index, cap_pressure_index);

        auto MlT = local_M.template block<cap_pressure_size, temperature_size>(
            cap_pressure_index, temperature_index);

        auto Mlus =
            local_M.template block<cap_pressure_size, displacement_size>(
                cap_pressure_index, displacement_index);

        auto MeT = local_M.template block<temperature_size, temperature_size>(
            temperature_index, temperature_index);

        auto Mepg = local_M.template block<temperature_size, gas_pressure_size>(
            temperature_index, gas_pressure_index);

        auto Mepc = local_M.template block<temperature_size, cap_pressure_size>(
            temperature_index, cap_pressure_index);

        typename ShapeMatricesTypePressure::NodalMatrixType Laplace =
            ShapeMatricesTypePressure::NodalMatrixType::Zero(gas_pressure_size,
                                                             gas_pressure_size);

        auto Kgpg =
            local_K.template block<gas_pressure_size, gas_pressure_size>(
                gas_pressure_index, gas_pressure_index);

        auto Klpg =
            local_K.template block<cap_pressure_size, gas_pressure_size>(
                cap_pressure_index, gas_pressure_index);

        auto Klpc =
            local_K.template block<cap_pressure_size, cap_pressure_size>(
                cap_pressure_index, cap_pressure_index);

        auto Kepg = local_K.template block<gas_pressure_size, temperature_size>(
            gas_pressure_index, temperature_index);

        auto Kepc = local_K.template block<cap_pressure_size, temperature_size>(
            cap_pressure_index, temperature_index);

        auto KeT = local_K.template block<temperature_size, temperature_size>(
            temperature_index, temperature_index);

        typename ShapeMatricesTypePressure::NodalMatrixType Aepg =
            ShapeMatricesTypePressure::NodalMatrixType::Zero(temperature_size,
                                                             temperature_size);

        typename ShapeMatricesTypePressure::NodalMatrixType Aepc =
            ShapeMatricesTypePressure::NodalMatrixType::Zero(temperature_size,
                                                             temperature_size);

        typename ShapeMatricesTypePressure::NodalMatrixType AeT =
            ShapeMatricesTypePressure::NodalMatrixType::Zero(temperature_size,
                                                             temperature_size);

        typename ShapeMatricesTypePressure::NodalMatrixType LeT =
            ShapeMatricesTypePressure::NodalMatrixType::Zero(temperature_size,
                                                             temperature_size);

        auto Kupg =
            local_K.template block<displacement_size, gas_pressure_size>(
                displacement_index, gas_pressure_index);
        auto Kupc =
            local_K.template block<displacement_size, cap_pressure_size>(
                displacement_index, cap_pressure_index);
        auto Kuu=
            local_K.template block<displacement_size, displacement_size>(
                displacement_index, displacement_index);
        auto Bg =
            local_rhs.template segment<gas_pressure_size>(gas_pressure_index);

        auto Bl =
            local_rhs.template segment<cap_pressure_size>(cap_pressure_index);

        auto Bu =
            local_rhs.template segment<displacement_size>(displacement_index);
        //
        //          auto gravity_operator =
        //                  local_rhs.template
        //                  segment<gas_pressure_size>(gas_pressure_index);

        //         ********************************************************************

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        double const& dt = _process_data.dt;

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto& ip_data = _ip_data[ip];

            auto const& w = ip_data.integration_weight;
            auto const& N_u_op = ip_data.N_u_op;
            auto const& N_u = ip_data.N_u;
            auto const& dNdx_u = ip_data.dNdx_u;
            auto const& N_p = ip_data.N_p;
            auto const& dNdx_p = ip_data.dNdx_p;
            // auto const& mass_operator = N_p.transpose() * N_p * w;
            auto const x_coord =
                interpolateXCoordinate<ShapeFunctionDisplacement,
                                       ShapeMatricesTypeDisplacement>(_element,
                                                                      N_u);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunctionDisplacement::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx_u, N_u, x_coord,
                                                     _is_axially_symmetric);
            auto& eps = ip_data.eps;
            auto const& sigma_eff = ip_data.sigma_eff;

            //            double const S = _process_data.specific_storage(t,
            //            x_position)[0]; double const K_over_mu =
            //                _process_data.intrinsic_permeability(t,
            //                x_position)[0] / _process_data.fluid_viscosity(t,
            //                x_position)[0];

            //            auto const permeability =
            //            _process_data.intrinsic_permeability(t,
            //            x_position)[0];

            // get porous medium properties from process data
            auto& medium = _process_data.medium;

            auto const p_cap = capillary_pressure.dot(N_p);
            auto const p_GR = gas_phase_pressure.dot(N_p);
            auto const T = temperature.dot(N_p);
            // auto const u = N_u_op.transpose().dot(displacement);

            const double p_LR = p_GR - p_cap;

#define nDBG_OUTPUT

#ifdef DBG_OUTPUT
            std::cout << "= Shape functions: ===============\n";

            std::cout << " N_u_op:\n " << N_u_op << "\n";
            std::cout << "    N_u:\n " << N_u << "\n";
            std::cout << " dNdx_u:\n " << dNdx_u << "\n";
            std::cout << "    N_p:\n " << N_p << "\n";
            std::cout << " dNdx_p:\n " << dNdx_p << "\n";

            std::cout << "= Primary variables: =============\n";
            std::cout << " Integration Point:\n";
            std::cout << " p_cap: " << p_cap << "\n";
            std::cout << " p_GR: " << p_GR << "\n";
            std::cout << " T: " << T << "\n";
            // std::cout << " u: " << u << "\n";
            std::cout << "----------------------------------\n";
            std::cout << " Nodal values:\n";
            std::cout << " p_cap:\n " << capillary_pressure << "\n";
            std::cout << " p_GR:\n " << gas_phase_pressure << "\n";
            std::cout << " T:\n " << temperature << "\n";
             std::cout << " u:\n "  << "???" << "\n";

            std::cout << "= Mass Operator: ===============\n";
            std::cout << " N_u_op:\n " <<  N_p.transpose() * N_p * w << "\n";
            std::cout << "==================================\n";

#endif

            // insert all primary variables into one object
            MPL::VariableArray variables;
            variables[MPL::Variables::p_cap] = p_cap;
            variables[MPL::Variables::p_GR] = p_GR;
            variables[MPL::Variables::p_LR] = p_LR;
            variables[MPL::Variables::T] = T;
            // todo: displacement

            // get fluid phase properties
            auto const& solid_phase = medium.phase(0);
            auto const& liquid_phase = medium.phase(1);
            auto const& gas_phase = medium.phase(2);

            // intrinsic permeability
            double const permeability =
                MPL::getScalar(medium.property(MPL::permeability));
#ifdef DBG_OUTPUT
            std::cout << "   permeability: " << permeability << " \n";
#endif
            GlobalDimMatrixType permeability_tensor =
                GlobalDimMatrixType::Zero(DisplacementDim, DisplacementDim);
            permeability_tensor.diagonal().setConstant(permeability);

#ifdef DBG_OUTPUT
            std::cout << "   permeability_tensor: " << permeability_tensor
                      << " \n";
            std::cout << "==================================\n";
#endif
            auto const alpha_B = MPL::getScalar(
                    medium.property(MPL::PropertyEnum::biot_coefficient),
                    variables);

            auto const rho_SR =
                    MPL::getScalar(solid_phase.property(MPL::PropertyEnum::density),
                            variables);

            const double beta_p_SR = getScalar(
                solid_phase.property(MPL::PropertyEnum::compressibility));
            const double beta_p_PR =
                getScalar(medium.property(MPL::PropertyEnum::compressibility));
            const double beta_p_GR = getScalar(
                gas_phase.property(MPL::PropertyEnum::compressibility));
            const double beta_p_LR = getScalar(
                liquid_phase.property(MPL::PropertyEnum::compressibility));

#ifdef DBG_OUTPUT
            std::cout << "   beta_p_SR: " << beta_p_SR << " \n";
            std::cout << "   beta_p_PR: " << beta_p_PR << " \n";
            std::cout << "   beta_p_GR: " << beta_p_GR << " \n";
            std::cout << "   beta_p_LR: " << beta_p_LR << " \n";
            std::cout << "==================================\n";
#endif


                        auto const rho_LR = MPL::getScalar(
                                liquid_phase.property(MPL::PropertyEnum::density),
                                variables);

                        auto const rho_GR =
                                MPL::getScalar(gas_phase.property(MPL::PropertyEnum::density),
                                        variables);

            variables[MPL::Variables::liquid_density] = rho_LR;
            variables[MPL::Variables::gas_density] = rho_GR;

#ifdef DBG_OUTPUT
            std::cout << "   alpha_B: " << alpha_B << " \n";
            std::cout << "==================================\n";
            std::cout << "   rho_SR: " << rho_SR << " \n";
            std::cout << "   rho_LR: " << rho_LR << " \n";
            std::cout << "   rho_GR: " << rho_GR << " \n";
            std::cout << "==================================\n";
#endif
            auto const phi = MPL::getScalar(
                medium.property(MPL::PropertyEnum::porosity), variables);
            auto const phi_S = 1. - phi;

#ifdef DBG_OUTPUT
            std::cout << "  phi : " << phi << " \n";
            std::cout << "  phi_S : " << phi_S << " \n";
            std::cout << "==================================\n";
#endif


//            auto const s_L =
//                MPL::getScalar(medium.property(MPL::PropertyEnum::saturation),
//                               variables);


            auto const s_L_r = 0.2;
//            auto const s_G_r = 0.5;
            auto const s_a = -1.9722e-11;
            auto const s_b = 2.4279;


            auto const s_L = std::max(s_L_r,1. + s_a * std::pow (std::max(0.,p_cap), s_b));
//            auto const s_L = std::min(1.0-s_G_r,std::max(0.0, s_m*p_cap+1.0));

//            for (double test_p = -120000; test_p <= 120000; test_p +=5000)
//            {
//            double test_s_L = std::min(1.0-s_G_r,std::max(s_L_r,1. - 1.9722e-11 * std::pow (std::max(0.,test_p), 2.4279)));
//            std::cout << test_p << " " << test_s_L << "\n";
//            }
//
//            OGS_FATAL ("");

            auto const s_G = 1 - s_L;
            auto const s_e = (s_L - s_L_r) / (1.0 - s_L_r);



#ifdef DBG_OUTPUT
            std::cout << "  p_cap : " << p_cap << " \n";
            std::cout << "  s_L : " << s_L << " \n";
            std::cout << "  s_G: " << s_G << " \n";
#endif

            auto const mu_LR = MPL::getScalar(
                liquid_phase.property(MPL::PropertyEnum::viscosity),
                variables);
            auto const mu_GR =
                MPL::getScalar(gas_phase.property(MPL::PropertyEnum::viscosity),
                               variables);

#ifdef DBG_OUTPUT
            std::cout << "  mu_LR: " << mu_LR << " \n";
            std::cout << "  mu_GR: " << mu_GR << " \n";
            std::cout << "==================================\n";
#endif
            auto const rho =
                phi_S * rho_SR + phi * (s_L * rho_LR + s_G * rho_GR);
            auto const& b = _process_data.specific_body_force;

#ifdef DBG_OUTPUT
            std::cout << "  rho: " << rho << " \n";
            std::cout << "   Gravity vector: \n";
            std::cout << "   b: \n" << b << " \n";
            std::cout << "==================================\n";
#endif
            static int const KelvinVectorSize =
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value;

            auto const& identity2 =
                MathLib::KelvinVector::Invariants<KelvinVectorSize>::identity2;

            eps.noalias() = B * displacement;

            auto const e = // volume strain
                MathLib::KelvinVector::Invariants<KelvinVectorSize>::trace(eps);
#ifdef DBG_OUTPUT
            std::cout << "   e: " << e << " \n";
            std::cout << "==================================\n";
#endif

            auto C = ip_data.updateConstitutiveRelation(t, x_position, dt, displacement,
                    T, p_GR);


#ifdef DBG_OUTPUT
            // std::cout << "   M_g : " << molar_mass << " \n";
            std::cout << "   C : " << C << " \n";
            std::cout << "==================================\n";
#endif

            auto const drhoGRdpGR = MPL::getScalarDerivative(
                gas_phase.property(MPL::density), variables, MPL::p_GR);

            auto const drhoLRdpLR =
                MPL::getScalarDerivative(liquid_phase.property(MPL::density),
                                         variables, MPL::p_LR);

            /*
            auto const molar_mass =
                getScalar(gas_phase.property(MPL::PropertyEnum::molar_mass));
                */

#ifdef DBG_OUTPUT
            // std::cout << "   M_g : " << molar_mass << " \n";
//            std::cout << "==================================\n";
            std::cout << "   drho_gr_dp_gr : " << drhoGRdpGR << " \n";
            std::cout << "   drho_lr_dp_lr : " << drhoLRdpLR << " \n";
            std::cout << "==================================\n";
#endif

//            auto const dsLdpc = MPL::getScalarDerivative(
//                medium.property(MPL::PropertyEnum::saturation),
//                variables, MPL::Variables::p_cap);

//            auto const dsLdpc = -4.78830438E-11*std::pow(p_cap,1.4279);
            auto const dsLdpc = s_a*s_b*std::pow(std::max(0.,p_cap), s_b - 1.0);


#ifdef DBG_OUTPUT
            std::cout << "   dsLdpc : " << dsLdpc << " \n";
            std::cout << "==================================\n";
#endif

//            auto const k_rel_LR = MPL::getPair(
//                medium.property(MPL::PropertyEnum::relative_permeability),
//                variables)[0];
//            auto const k_rel_GR = MPL::getPair(
//                medium.property(MPL::PropertyEnum::relative_permeability),
//                variables)[1];

            auto const k_rel_LR = 1.0 - 2.207*std::pow((1.0 - s_L), 1.0121);
            auto const min_k_rel_GR = 0.0001;

            auto const k_rel_GR = (1.0 - s_e) * (1 - s_e)
                    		* (1.0 - std::pow(s_e, (5./3.))) + min_k_rel_GR;

    //        std::cout << "pCap: " << p_cap << "s_L: " << s_L << " s_G: " << s_G << " dsldpc: " << dsLdpc << " k_rel_LR: " << k_rel_LR << " k_rel_GR: " << k_rel_GR << "\n";

#ifdef DBG_OUTPUT
            std::cout << "    k_rel_LR: " << k_rel_LR << " \n";
            std::cout << "    k_rel_GR: " << k_rel_GR << " \n";
            std::cout << "==================================\n";
#endif

            // for output only
            ip_data.pressure_gas_linear = p_GR;
            ip_data.pressure_cap_linear = p_cap;
            ip_data.pressure_wet = p_LR;
            ip_data.density_gas = rho_GR;
            ip_data.density_liquid = rho_LR;
            ip_data.saturation = s_L;

            const double k_over_mu_GR = k_rel_GR / mu_GR;
            const double k_over_mu_LR = k_rel_LR / mu_LR;

#ifdef DBG_OUTPUT
            std::cout << "   k_over_mu_GR: " << k_over_mu_GR << " \n";
            std::cout << "   k_over_mu_LR: " << k_over_mu_LR << " \n";
            std::cout << "==================================\n";
#endif

            ip_data.velocity_liquid.noalias() =
                -permeability_tensor * k_over_mu_LR * dNdx_p *
                    (gas_phase_pressure - capillary_pressure) +
                permeability_tensor * k_over_mu_LR * rho_LR * b;

            ip_data.velocity_gas.noalias() =
                -permeability_tensor * k_over_mu_GR * dNdx_p *
                    gas_phase_pressure +
                permeability_tensor * k_over_mu_GR * rho_GR * b;

            const auto w_LS = ip_data.velocity_liquid;
            const auto w_GS = ip_data.velocity_gas;

#ifdef DBG_OUTPUT
            std::cout << "   Velocities: \n";
            std::cout << "   ----------------------------------\n";
            std::cout << "   w_GR: \n" << w_GS << " \n";
            std::cout << "   w_LR: \n" << w_LS << " \n";
            std::cout << "   Velocities : \n";
            std::cout << "==================================\n";
#endif

            Laplace.noalias() =
                dNdx_p.transpose() * permeability_tensor * dNdx_p * w;

#ifdef DBG_OUTPUT
            std::cout << "   Laplace: \n";
            std::cout << "   L:\n " << Laplace << " \n";
            std::cout << "==================================\n";
#endif
            Mgpg.noalias() += N_p.transpose() * rho_GR * s_G *
                              (phi * beta_p_GR + (alpha_B - phi) * beta_p_SR) *
                              N_p * w;

#ifdef DBG_OUTPUT
            std::cout << "   Mgpg:\n " << Mgpg << " \n";
            std::cout << "==================================\n";
#endif
            Mgpc.noalias() -=
                N_p.transpose() * rho_GR *
                ((phi + s_G * (alpha_B - phi) * beta_p_SR * p_cap) * dsLdpc +
                 s_G * s_L * (alpha_B - phi) * beta_p_SR) *
                N_p * w;

#ifdef DBG_OUTPUT
            std::cout << "   Mgpc:\n " << Mgpc << " \n";
            std::cout << "==================================\n";
#endif

            auto const drhoGRdT = MPL::getScalarDerivative(
                gas_phase.property(MPL::density), variables, MPL::T);

            const double beta_T_GR = rho_GR == 0. ? 0 : drhoGRdT / rho_GR;

            auto const drhoLRdT = MPL::getScalarDerivative(
                liquid_phase.property(MPL::density), variables, MPL::T);

            const double beta_T_LR = drhoLRdT / rho_LR;

            const double beta_T_SR = getScalar(
                solid_phase.property(MPL::PropertyEnum::thermal_expansivity));

#ifdef DBG_OUTPUT
            std::cout << " drhoGRdT: " << drhoGRdT << " \n";
            std::cout << " beta_T_GR: " << beta_T_GR << " \n";
            std::cout << " drhoLRdT: " << drhoLRdT << " \n";
            std::cout << " beta_T_LR: " << beta_T_LR << " \n";
            std::cout << " beta_T_SR: " << beta_T_SR << " \n";

            std::cout << "=================================\n";
#endif
            MgT.noalias() -= N_p.transpose() * s_G * rho_GR *
                             (phi * beta_T_GR + beta_T_SR * (alpha_B - phi)) *
                             N_p * w;

#ifdef DBG_OUTPUT
            std::cout << "   MgT:\n " << MgT << " \n";
            std::cout << "==================================\n";
#endif

            Mgus.noalias() += N_p.transpose() * s_G * rho_GR * alpha_B *
                              identity2.transpose() * B * w;
#ifdef DBG_OUTPUT
            std::cout << "   Mgus:\n " << Mgus << " \n";
            std::cout << "==================================\n";
#endif

            Mlpg.noalias() += N_p.transpose() * rho_LR * s_L *
                              (phi * beta_p_LR + (alpha_B - phi) * beta_p_SR) *
                              N_p * w;
#ifdef DBG_OUTPUT
            std::cout << "   Mlpg:\n " << Mlpg << " \n";
            std::cout << "==================================\n";
#endif

            const auto a_L = s_L * (alpha_B - phi) * beta_p_SR;

            Mlpc.noalias() +=
            		N_p.transpose() * rho_LR *
					(dsLdpc * (phi - a_L*p_cap) - s_L*(phi*beta_p_LR + a_L)) *  N_p * w;

#ifdef DBG_OUTPUT
            std::cout << "     a_L: " << a_L << " \n";
            std::cout << "==================================\n";
            std::cout << "   Mlpc:\n " << Mlpc << " \n";
            std::cout << "==================================\n";
#endif

            MlT.noalias() -= N_p.transpose() * s_L * rho_LR *
                             (phi * beta_T_LR + beta_T_SR * (alpha_B - phi)) *
                             N_p * w;

#ifdef DBG_OUTPUT
            std::cout << "   MlT:\n " << MlT << " \n";
            std::cout << "==================================\n";
#endif

            Mlus.noalias() += N_p.transpose() * rho_LR * s_L * alpha_B *
                              identity2.transpose() * B * w;

#ifdef DBG_OUTPUT
            std::cout << "   Mlus:\n " << Mlus << " \n";
            std::cout << "==================================\n";
#endif

            const double cp_G = getScalar(
                gas_phase.property(MPL::PropertyEnum::specific_heat_capacity));
            const double cp_L = getScalar(liquid_phase.property(
                MPL::PropertyEnum::specific_heat_capacity));
            const double cp_S = getScalar(solid_phase.property(
                MPL::PropertyEnum::specific_heat_capacity));

            const double phi_L = s_L * phi;
            const double phi_G = s_G * phi;

            const double rho_cp_eff = phi_G * rho_GR * cp_G +
                                  phi_L * rho_LR * cp_L + phi_S * rho_SR * cp_S;

#ifdef DBG_OUTPUT
            std::cout << "   cp_G: " << cp_G << " \n";
            std::cout << "   cp_L: " << cp_L << " \n";
            std::cout << "   cp_S: " << cp_S << " \n";

            std::cout << " cp_eff: " << rho_cp_eff << " \n";
            std::cout << "=================================\n";
#endif

            Mepg.noalias() -=
                N_p.transpose() *
                (phi_G*beta_T_GR + phi_L*beta_T_LR + phi_S*beta_T_SR) * T * N_p *
                w;

#ifdef DBG_OUTPUT
            std::cout << "   Mepg:\n " << Mepg << " \n";
            std::cout << "=================================\n";
#endif
            Mepc.noalias() +=
                N_p.transpose() *
                (phi_L*beta_T_LR - (s_L + p_cap*dsLdpc)*phi_S*beta_T_SR) * T * N_p *
                w;

#ifdef DBG_OUTPUT
            std::cout << "   Mepc:\n " << Mepc << " \n";
            std::cout << "=================================\n";
#endif

            MeT.noalias() += N_p.transpose() * rho_cp_eff * N_p * w;

#ifdef DBG_OUTPUT
            std::cout << "   MeT:\n " << MeT << " \n";
            std::cout << "=================================\n";
#endif

            Kgpg.noalias() += rho_GR * k_over_mu_GR * Laplace;

#ifdef DBG_OUTPUT
            std::cout << "   Kgpg:\n " << Kgpg << " \n";
            std::cout << "==================================\n";
#endif
            Klpg.noalias() += rho_LR * k_over_mu_LR * Laplace;
#ifdef DBG_OUTPUT
            std::cout << "   Klpg:\n " << Klpg << " \n";
            std::cout << "==================================\n";
#endif

            Klpc.noalias() -= rho_LR * k_over_mu_LR * Laplace;
#ifdef DBG_OUTPUT
            std::cout << "   Klpc:\n " << Klpc << " \n";
            std::cout << "==================================\n";
#endif

            Aepg.noalias() = - N_p.transpose() *
                             (phi_G*beta_T_GR*w_GS.transpose() +
                              phi_L*beta_T_LR*w_LS.transpose()) * T *
                             dNdx_p * w;

            Kepg.noalias() += Aepg;

#ifdef DBG_OUTPUT
            std::cout << "   Aepg:\n " << Aepg << " \n";
            std::cout << "   Kepg:\n " << Kepg << " \n";
            std::cout << "=================================\n";
            // Aepc(0, 0) = 3.4;
#endif

            Aepc.noalias() =
            		N_p.transpose() * (phi_L*beta_T_LR*T*w_LS.transpose())
					* dNdx_p * w;

            Kepc.noalias() += Aepc;

#ifdef DBG_OUTPUT
            std::cout << "   Aepc:\n " << Aepc << " \n";
            std::cout << "   Kepc:\n " << Kepc << " \n";
            std::cout << "=================================\n";
#endif

            AeT.noalias() = N_p.transpose() *
                            (s_G*rho_GR*cp_G*w_GS.transpose() +
                             s_L*rho_LR*cp_L*w_LS.transpose()) *
                            dNdx_p * w;

#ifdef DBG_OUTPUT
            std::cout << "   AeT:\n " << AeT << " \n";
            std::cout << "=================================\n";
#endif
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>
                conductivity_tensor;
            conductivity_tensor.setIdentity();

            const double lambda_G = getScalar(
                gas_phase.property(MPL::PropertyEnum::thermal_conductivity));
            const double lambda_L = getScalar(
                liquid_phase.property(MPL::PropertyEnum::thermal_conductivity));
            const double lambda_S = getScalar(
                solid_phase.property(MPL::PropertyEnum::thermal_conductivity));

            const double lambda_eff = phi_G * lambda_G +
                                      phi_L * lambda_L +
                                      phi_S * lambda_S;

            const auto conductivity_effective =
                lambda_eff * conductivity_tensor;

            LeT.noalias() =
                dNdx_p.transpose() * conductivity_effective * dNdx_p * w;

            KeT.noalias() += AeT + LeT;

#ifdef DBG_OUTPUT
            std::cout << "   LeT:\n " << LeT << " \n";
            std::cout << "   KeT:\n " << KeT << " \n";
            std::cout << "=================================\n";
#endif

            Kupg.noalias() -= B.transpose() * alpha_B * identity2 * N_p * w;

#ifdef DBG_OUTPUT
            std::cout << "   Kupg:\n " << Kupg << " \n";
            std::cout << "==================================\n";
#endif
            Kupc.noalias() +=
                B.transpose() * alpha_B * identity2 * s_L * N_p * w;
#ifdef DBG_OUTPUT
            std::cout << "   Kupc:\n " << Kupc << " \n";
            std::cout << "==================================\n";
#endif

            Kuu.noalias() +=
                B.transpose() * C * B * w;
#ifdef DBG_OUTPUT
            std::cout << "   Kuu:\n " << Kupc << " \n";
            std::cout << "==================================\n";
#endif

            auto const gravity_operator =
                (dNdx_p.transpose() * permeability_tensor * b * w).eval();

            Bg.noalias() += rho_GR * rho_GR * k_over_mu_GR * gravity_operator;

#ifdef DBG_OUTPUT
            std::cout << "   Bg:\n " << Bg << " \n";
            std::cout << "==================================\n";
#endif

            Bl.noalias() += rho_LR * rho_LR * k_over_mu_LR * gravity_operator;

#ifdef DBG_OUTPUT
            std::cout << "   Bl:\n " << Bl << " \n";
            std::cout << "==================================\n";
#endif

            Bu.noalias() += N_u_op.transpose() * rho * b * w;

#ifdef DBG_OUTPUT
            std::cout << "   Bu:\n " << Bu << " \n";
            std::cout << "==================================\n";
#endif
        }

#ifdef DBG_OUTPUT

                std::cout << "== Local M: ====\n";
                std::cout << local_M << "\n";
                std::cout << "================\n";
                std::cout << "== Local K: ====\n";
                std::cout << local_K << "\n";
                std::cout << "================\n";
                std::cout << "== Local f: ====\n";
                std::cout << local_rhs << "\n";
                std::cout << "================\n";

            //   OGS_FATAL ("##########################################");
#endif

//        for (unsigned row = 0; row < Mgpc.cols(); row++)
//        {
//            for (unsigned column = 0; column < Mgpc.cols(); column++)
//            {
//                if (row != column)
//                {
//                    Mgpc(row, row) += Mgpc(row, column);
//                    Mgpc(row, column) = 0.0;
//                    Mgpg(row, row) += Mgpg(row, column);
//                    Mgpg(row, column) = 0.0;
//                    Mlpc(row, row) += Mlpc(row, column);
//                    Mlpc(row, column) = 0.0;
//                }
//            }
  //      }
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override
    {

        assert(local_x.size() == gas_pressure_size + cap_pressure_size +
                                     temperature_size + displacement_size);

        const auto matrix_size = gas_pressure_size + cap_pressure_size +
                                 temperature_size + displacement_size;

        // primary variables
        auto gas_phase_pressure =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                          gas_pressure_size);

        auto gas_phase_pressure_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                gas_pressure_size> const>(local_xdot.data() + gas_pressure_index,
                                          gas_pressure_size);

        auto capillary_pressure =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                cap_pressure_size> const>(local_x.data() + cap_pressure_index,
                                          cap_pressure_size);

        auto capillary_pressure_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                cap_pressure_size> const>(local_xdot.data() + cap_pressure_index,
                                          cap_pressure_size);


        auto temperature =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                temperature_size> const>(local_x.data() + temperature_index,
                                         temperature_size);

        auto temperature_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                temperature_size> const>(local_xdot.data() + temperature_index,
                                         temperature_size);

        auto displacement =
            Eigen::Map<typename ShapeMatricesTypeDisplacement::
                           template VectorType<displacement_size> const>(
                local_x.data() + displacement_index, displacement_size);


        auto displacement_dot =
            Eigen::Map<typename ShapeMatricesTypeDisplacement::
                           template VectorType<displacement_size> const>(
                local_xdot.data() + displacement_index, displacement_size);


        // create Jacobian
        auto J = MathLib::createZeroedMatrix<
            typename ShapeMatricesTypeDisplacement::template MatrixType<
                matrix_size, matrix_size>>(
            local_Jac_data, matrix_size, matrix_size);

        // create residuum
        auto r = MathLib::createZeroedVector<
            typename ShapeMatricesTypeDisplacement::template VectorType<
                matrix_size>>(
            local_rhs_data, matrix_size);

        //         ********************************************************************
        //         ********************************************************************


        // Matrix templates

        using PPMatrix = typename ShapeMatricesTypeDisplacement::template
                MatrixType<Np_intPoints, Np_intPoints>;
        using PUMatrix = typename ShapeMatricesTypeDisplacement::template
                MatrixType<Np_intPoints, Nu_intPoints>;
        using UPMatrix = typename ShapeMatricesTypeDisplacement::template
                MatrixType<Nu_intPoints, Np_intPoints>;
        using UUMatrix = typename ShapeMatricesTypeDisplacement::template
                MatrixType<Nu_intPoints, Nu_intPoints>;


        using PVector = typename ShapeMatricesTypeDisplacement::template
        		VectorType<Np_intPoints>;

        using UVector = typename ShapeMatricesTypeDisplacement::template
        		VectorType<Nu_intPoints>;

        // mass matrices

        PPMatrix Mgpg =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix Mgpc =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix MgT =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PUMatrix Mgus =
        		PUMatrix::Zero(Np_intPoints, Nu_intPoints);

        PPMatrix Mlpg =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix Mlpc =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix MlT =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PUMatrix Mlus =
        		PUMatrix::Zero(Np_intPoints, Nu_intPoints);

        PPMatrix Mepg =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix Mepc =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix MeT =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);

        // mass matrix derivatives

        PPMatrix dMgpg_dpg =
                PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dMgpg_dpc =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix dMgpg_dT =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);

        PPMatrix dMgpc_dpg =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dMgpc_dpc =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix dMgpc_dT =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PUMatrix dMgpc_dus =
        		PUMatrix::Zero(Np_intPoints, Nu_intPoints);

        PPMatrix dMgT_dpg =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dMgT_dpc =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix dMgT_dT =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);

        PUMatrix dMgus_dpg =
        		PUMatrix::Zero(Np_intPoints, Nu_intPoints);
        PUMatrix dMgus_dpc =
        		PUMatrix::Zero(Np_intPoints,Nu_intPoints);
        PUMatrix dMgus_dT =
        		PUMatrix::Zero(Np_intPoints,Nu_intPoints);

        PPMatrix dMlpg_dpg =
                PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dMlpg_dpc =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix dMlpg_dT =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);


        PPMatrix dMlpc_dpg =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dMlpc_dpc =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix dMlpc_dT =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);

        PPMatrix dMlT_dpg =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dMlT_dpc =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix dMlT_dT =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);

        PUMatrix dMlus_dpg =
        		PUMatrix::Zero(Np_intPoints, Nu_intPoints);
        PUMatrix dMlus_dpc =
        		PUMatrix::Zero(Np_intPoints,Nu_intPoints);
        PUMatrix dMlus_dT =
        		PUMatrix::Zero(Np_intPoints,Nu_intPoints);


        PPMatrix dMepg_dpg =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dMepg_dpc =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dMepg_dT =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);

        PPMatrix dMepc_dpg =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dMepc_dpc =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dMepc_dT =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);

        PPMatrix dMeT_dpg =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dMeT_dpc =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dMeT_dT =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);


        // Laplace matrix coefficients
        PPMatrix Lgpg =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix Llpg =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix Llpc =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix LeT =
        		PPMatrix::Zero(Np_intPoints,Np_intPoints);

        // Laplace matrix derivatives
        PPMatrix dLgpg_dpg =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dLgpg_dpc =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dLgpg_dT =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);

        PPMatrix dLlpg_dpg =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dLlpg_dpc =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dLlpg_dT =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);

        PPMatrix dLlpc_dpg =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dLlpc_dpc =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);
        PPMatrix dLlpc_dT =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);

        PPMatrix dLeT_dpc =
        		PPMatrix::Zero(Np_intPoints, Np_intPoints);



        // advection matrices
        PPMatrix Aepg =
                PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix Aepc =
                PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix AeT =
                PPMatrix::Zero(Np_intPoints,Np_intPoints);

        PPMatrix dAepgPg_dpg =
                PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix dAepgPg_dpc =
                PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix dAepgPg_dT =
                PPMatrix::Zero(Np_intPoints,Np_intPoints);

        PPMatrix dAepcPc_dpg =
                PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix dAepcPc_dpc =
                PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix dAepcPc_dT =
                PPMatrix::Zero(Np_intPoints,Np_intPoints);

        PPMatrix dAeTT_dpg =
                PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix dAeTT_dpc =
                PPMatrix::Zero(Np_intPoints,Np_intPoints);
        PPMatrix dAeTT_dT =
                PPMatrix::Zero(Np_intPoints,Np_intPoints);



        // stiffness matrices
        UPMatrix Kupg =
        		UPMatrix::Zero(Nu_intPoints,Np_intPoints);
        UPMatrix Kupc =
        		UPMatrix::Zero(Nu_intPoints,Np_intPoints);

        UPMatrix dKupc_dpc =
        		UPMatrix::Zero(Nu_intPoints,Np_intPoints);

        UUMatrix dfu_du =
        		UUMatrix::Zero(Nu_intPoints,Nu_intPoints);

        // right hand side vectors
        PVector fg = PVector::Zero(Np_intPoints);
        PVector fl = PVector::Zero(Np_intPoints);
        PVector fe = PVector::Zero(Np_intPoints);
        UVector fu = UVector::Zero(Nu_intPoints);

        // right hand side derivatives
        PVector dfg_dpg = PVector::Zero(Np_intPoints);
        PVector dfg_dpc = PVector::Zero(Np_intPoints);
        PVector dfg_dT = PVector::Zero(Np_intPoints);
        PVector dfg_du = PVector::Zero(Np_intPoints);

        PVector dfl_dpg = PVector::Zero(Np_intPoints);
        PVector dfl_dpc = PVector::Zero(Np_intPoints);
        PVector dfl_dT = PVector::Zero(Np_intPoints);
        PVector dfl_du = PVector::Zero(Np_intPoints);

        PVector dfe_dpg = PVector::Zero(Np_intPoints);
        PVector dfe_dpc = PVector::Zero(Np_intPoints);
        PVector dfe_dT = PVector::Zero(Np_intPoints);
        PVector dfe_du = PVector::Zero(Np_intPoints);

        UVector dfu_dpg = UVector::Zero(Nu_intPoints);
        UVector dfu_dpc = UVector::Zero(Nu_intPoints);
        UVector dfu_dT = UVector::Zero(Nu_intPoints);


        typename ShapeMatricesTypePressure::NodalMatrixType Laplace =
        		ShapeMatricesTypePressure::NodalMatrixType::Zero(Np_intPoints,
        				Np_intPoints);

        // Jacobian blocks:
        auto drg_dpg =
        		J.template block<gas_pressure_size, gas_pressure_size>(
        				gas_pressure_index, gas_pressure_index);
        auto drg_dpc =
        		J.template block<gas_pressure_size, cap_pressure_size>(
        				gas_pressure_index, cap_pressure_index);
        auto drg_dT=
        		J.template block<gas_pressure_size, temperature_size>(
        				gas_pressure_index, temperature_index);
        auto drg_dus =
        		J.template block<gas_pressure_size, displacement_size>(
        				gas_pressure_index, displacement_index);
        auto drl_dpg =
        		J.template block<cap_pressure_size, gas_pressure_size>(
        				cap_pressure_index, gas_pressure_index);
        auto drl_dpc =
        		J.template block<cap_pressure_size, cap_pressure_size>(
        				cap_pressure_index, cap_pressure_index);
        auto drl_dT=
        		J.template block<cap_pressure_size, temperature_size>(
        				cap_pressure_index, temperature_index);
        auto drl_dus =
        		J.template block<cap_pressure_size, displacement_size>(
        				cap_pressure_index, displacement_index);
        auto dre_dpg =
        		J.template block<temperature_size, gas_pressure_size>(
        				temperature_index, gas_pressure_index);
        auto dre_dpc =
        		J.template block<temperature_size, cap_pressure_size>(
        				temperature_index, cap_pressure_index);
        auto dre_dT=
        		J.template block<temperature_size, temperature_size>(
        				temperature_index, temperature_index);
//        auto dre_dus =
//        		J.template block<temperature_size, displacement_size>(
//        				temperature_index, displacement_index);
        auto dru_dpg =
        		J.template block<displacement_size, gas_pressure_size>(
        				displacement_index, gas_pressure_index);
        auto dru_dpc =
        		J.template block<displacement_size, cap_pressure_size>(
        				displacement_index, cap_pressure_index);
        auto dru_dT =
        		J.template block<displacement_size, temperature_size>(
        				displacement_index, temperature_index);
        auto dru_dus =
        		J.template block<displacement_size, displacement_size>(
        				displacement_index, displacement_index);

        // resuduum segments
        auto rg =  r.template segment<gas_pressure_size>(gas_pressure_index);
        auto rl =  r.template segment<cap_pressure_size>(cap_pressure_index);
        auto re =  r.template segment<temperature_size>(temperature_index);
        auto ru =  r.template segment<displacement_size>(displacement_index);

         SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        double const& dt = _process_data.dt;

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto& ip_data = _ip_data[ip];

            auto const& w = ip_data.integration_weight;
            auto const& N_u_op = ip_data.N_u_op;
            auto const& N_u = ip_data.N_u;
            auto const& dNdx_u = ip_data.dNdx_u;
            auto const& N_p = ip_data.N_p;
            auto const& dNdx_p = ip_data.dNdx_p;

            auto const& Np =  ip_data.N_p;
            auto const& NpT = Np.transpose();

            auto const x_coord =
                interpolateXCoordinate<ShapeFunctionDisplacement,
                                       ShapeMatricesTypeDisplacement>(_element,
                                                                      N_u);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunctionDisplacement::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx_u, N_u, x_coord,
                                                     _is_axially_symmetric);

            static int const KelvinVectorSize =
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value;

            auto const& identity2 =
                MathLib::KelvinVector::Invariants<KelvinVectorSize>::identity2;

            auto const& mass_operator = N_p.transpose() * N_p * w;
            auto const& mass_operator2 = N_p.transpose() * identity2.transpose() * B * w;


            auto& eps = ip_data.eps;
            auto const& sigma_eff = ip_data.sigma_eff;

            // get medium properties
            auto& medium = _process_data.medium;

            // primary unknowns at integration point
            auto const p_GR = gas_phase_pressure.dot(N_p);
            auto const p_cap = capillary_pressure.dot(N_p);
            auto const T = temperature.dot(N_p);
            auto const u = N_u_op * displacement;

            // primary unknown derivatives at integration point
            auto const p_GR_dot = gas_phase_pressure_dot.dot(N_p);
            auto const p_cap_dot = capillary_pressure.dot(N_p);
            auto const T_dot = temperature_dot.dot(N_p);

            typename ShapeMatricesTypeDisplacement::GlobalDimVectorType u_ip(
            		DisplacementDim);
            for (int i = 0; i < u_ip.size(); ++i)
            {
            	NumLib::shapeFunctionInterpolate(
            			displacement.segment(i * ShapeFunctionDisplacement::NPOINTS,
            					ShapeFunctionDisplacement::NPOINTS),
								N_u, u_ip.coeffRef(i));
            }


            const double p_LR = p_GR - p_cap;

            eps.noalias() = B * displacement;

            auto const e = // volume strain
                MathLib::KelvinVector::Invariants<KelvinVectorSize>::trace(eps);

#define DBG_OUTPUT

#ifdef DBG_OUTPUT
            std::cout << "= Primary variables: =============\n";
            std::cout << " Integration Point:\n";
            std::cout << "       p_cap : " << p_cap << "\n";
            std::cout << "        p_GR : " << p_GR << "\n";
            std::cout << "           T : " << T << "\n";
            std::cout << "           u : " << u << "\n\n";
            std::cout << "----------------------------------\n";
            std::cout << "  Nodal values :\n";
            std::cout << "       p_cap :\n" << capillary_pressure << "\n\n";
            std::cout << "        p_GR :\n" << gas_phase_pressure << "\n\n";
            std::cout << "           T :\n" << temperature << "\n\n";
            std::cout << "           u :\n"  << displacement << "\n\n";
            std::cout << "----------------------------------\n";
            std::cout << "     x_coord :\n " << x_coord << "\n\n";
            std::cout << "==================================\n";

            std::cout << "= Shape functions: ===============\n";
            std::cout << "      N_u_op :\n" << N_u_op << "\n\n";
            std::cout << "         N_u :\n" << N_u << "\n\n";
            std::cout << "      dNdx_u :\n" << dNdx_u << "\n\n";
            std::cout << "         N_p :\n" << N_p << "\n\n";
            std::cout << "      dNdx_p :\n" << dNdx_p << "\n\n";
            std::cout << "           B :\n" << B << "\n\n";
            std::cout << "           m :\n" << identity2 << "\n\n";


            std::cout << "= Mass Operator: ===============\n";
            std::cout << " mass_op :\n" <<  N_p.transpose() * N_p * w << "\n\n";
            std::cout << " mass    :\n" <<  mass_operator << "\n";

            #endif

            // insert all primary variables into one object for GPML
            MPL::VariableArray variables;
            variables[MPL::Variables::p_cap] = p_cap;
            variables[MPL::Variables::p_GR] = p_GR;
            variables[MPL::Variables::p_LR] = p_LR;
            variables[MPL::Variables::T] = T;

            // get phase properties
            auto const& solid_phase = medium.phase(0);
            auto const& liquid_phase = medium.phase(1);
            auto const& gas_phase = medium.phase(2);


            /*
             * Porous medium properties
             */
            double const permeability =
                MPL::getScalar(medium.property(MPL::permeability));

            GlobalDimMatrixType permeability_tensor =
                GlobalDimMatrixType::Zero(DisplacementDim, DisplacementDim);
            permeability_tensor.diagonal().setConstant(permeability);

            auto const alpha_B = MPL::getScalar(
                    medium.property(MPL::PropertyEnum::biot_coefficient),
                    variables);

            auto const phi = MPL::getScalar(
                medium.property(MPL::PropertyEnum::porosity), variables);
            auto const phi_S = 1. - phi;

            // Saturation

//          auto const s_L =
//              MPL::getScalar(medium.property(MPL::PropertyEnum::saturation),
//                               variables);
            auto const s_L_r = 0.2;
            auto const s_a = -1.9722e-11;
            auto const s_b = 2.4279;

            auto const s_L = std::max(s_L_r,1. + s_a * std::pow (std::max(0.,p_cap), s_b));

            auto const s_G = 1 - s_L;
            auto const s_e = (s_L - s_L_r) / (1.0 - s_L_r);

            const double phi_L = s_L * phi;
            const double phi_G = s_G * phi;

//          auto const dsLdpc = MPL::getScalarDerivative(
//              medium.property(MPL::PropertyEnum::saturation),
//              variables, MPL::Variables::p_cap);

//          auto const dsLdpc = -4.78830438E-11*std::pow(p_cap,1.4279);
            auto const dsLdpc = s_a*s_b*std::pow(std::max(0.,p_cap), s_b - 1.0);
            auto const d2sLdpc2 = s_a*(s_b*(s_b - 1)) * std::pow(std::max(0.,p_cap), s_b - 2.0);


//          auto const k_rel_LR = MPL::getPair(
//              mediumte.property(MPL::PropertyEnum::relative_permeability),
//              variables)[0];
//          auto const k_rel_GR = MPL::getPair(
//              medium.property(MPL::PropertyEnum::relative_permeability),
//              variables)[1];

            auto const k_rel_LR = 1.0 - 2.207*std::pow((1.0 - s_L), 1.0121);
            auto const min_k_rel_GR = 0.0001;

            auto const k_rel_GR = (1.0 - s_e) * (1 - s_e)
                                                * (1.0 - std::pow(s_e, (5./3.))) + min_k_rel_GR;

            auto const dkrelGdsL = 0.0;
            auto const dkrelLdsL = 0.0;

            /*
             * Phase properties
             */
            // Densities
            auto const rho_SR =
                    MPL::getScalar(solid_phase.property(MPL::PropertyEnum::density),
                            variables);
            auto const rho_LR = MPL::getScalar(
                    liquid_phase.property(MPL::PropertyEnum::density),
                    variables);
            auto const rho_GR =
                    MPL::getScalar(gas_phase.property(MPL::PropertyEnum::density),
                            variables);

            auto const& b = _process_data.specific_body_force;

            variables[MPL::Variables::liquid_density] = rho_LR;
            variables[MPL::Variables::gas_density] = rho_GR;

//            auto const drhoGRdpGR = MPL::getScalarDerivative(
//                    gas_phase.property(MPL::density), variables, MPL::p_GR);
            auto const drhoGRdpGR = 0.;

//            auto const drhoLRdpLR =
//                                MPL::getScalarDerivative(liquid_phase.property(MPL::density),
//                                        variables, MPL::p_LR);
            auto const drhoLRdpLR = 0.;

            auto const drhoLRdpGR = drhoLRdpLR;
            auto const drhoLRdpCap = - drhoLRdpLR;

//            auto const drhoGRdT = MPL::getScalarDerivative(
//                    gas_phase.property(MPL::density), variables, MPL::T);
            auto const drhoGRdT = 0.;

//            auto const drhoLRdT = MPL::getScalarDerivative(
//                    liquid_phase.property(MPL::density), variables, MPL::T);
            auto const drhoLRdT = 0.;

            auto const drhoSRdpGR = 0.;
            auto const drhoSRdpCap = 0.;
            auto const drhoSRdT = 0.;

            // volume average density
            auto const rho =
            		phi_S * rho_SR + phi * (s_L * rho_LR + s_G * rho_GR);
            auto const drhodpGR = phi_G*drhoGRdpGR + phi_L*drhoLRdpGR + phi_S*drhoSRdpGR;
            auto const drhodpCap = phi*dsLdpc*(rho_LR-rho_GR) + phi_L*drhoLRdpCap + phi_S*drhoSRdpCap;
            auto const drhodT = phi_G*drhoGRdT + phi_L*drhoLRdT + phi_S*drhoSRdT;

            // Viscosities
            auto const mu_LR = MPL::getScalar(
                liquid_phase.property(MPL::PropertyEnum::viscosity),
                variables);
            auto const mu_GR =
                MPL::getScalar(gas_phase.property(MPL::PropertyEnum::viscosity),
                               variables);

            auto const dmuLRdpLR = 0.0;
            auto const dmuLRdT= 0.0;
            auto const dmuGRdpGR = 0.0;
            auto const dmuGRdT = 0.0;




            // heat capacities
            const double cp_G = getScalar(
                gas_phase.property(MPL::PropertyEnum::specific_heat_capacity));
            const double cp_L = getScalar(liquid_phase.property(
                MPL::PropertyEnum::specific_heat_capacity));
            const double cp_S = getScalar(solid_phase.property(
                MPL::PropertyEnum::specific_heat_capacity));

            // Elasticity matrix:
            auto C = ip_data.updateConstitutiveRelation(t, x_position, dt, displacement,
                    T, p_GR);

            // compressibilities
            const double beta_p_SR = getScalar(
                solid_phase.property(MPL::PropertyEnum::compressibility));
            const double beta_p_PR =
                getScalar(medium.property(MPL::PropertyEnum::compressibility));
            const double beta_p_GR = getScalar(
                gas_phase.property(MPL::PropertyEnum::compressibility));
            const double beta_p_LR = getScalar(
                liquid_phase.property(MPL::PropertyEnum::compressibility));


            // expansivities TODO
            const double beta_T_GR = rho_GR == 0. ? 0 : drhoGRdT / rho_GR;
            const double beta_T_LR = drhoLRdT / rho_LR;
            const double beta_T_SR = getScalar(
                solid_phase.property(MPL::PropertyEnum::thermal_expansivity));

            // heat conductivities
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>
                conductivity_tensor;
            conductivity_tensor.setIdentity();

            const double lambda_G = getScalar(
                gas_phase.property(MPL::PropertyEnum::thermal_conductivity));
            const double lambda_L = getScalar(
                liquid_phase.property(MPL::PropertyEnum::thermal_conductivity));
            const double lambda_S = getScalar(
                solid_phase.property(MPL::PropertyEnum::thermal_conductivity));

            const double lambda_eff = phi_G * lambda_G +
                                      phi_L * lambda_L +
                                      phi_S * lambda_S;

            // isotropic heat conductivity
            const auto conductivity_effective =
                lambda_eff * conductivity_tensor;

            /*
             * derived variables
             */
            const double k_over_mu_GR = k_rel_GR / mu_GR;
            const double k_over_mu_LR = k_rel_LR / mu_LR;

            const auto Sps = (alpha_B - phi) * beta_p_SR;
            const auto STs = (alpha_B - phi) * beta_T_SR;

            const double rho_cp_eff = phi_G * rho_GR * cp_G +
                    phi_L * rho_LR * cp_L + phi_S * rho_SR * cp_S;

            // darcy-velocities
            const auto w_LS =
                    (-permeability_tensor * k_over_mu_LR *
                    ( dNdx_p * (gas_phase_pressure - capillary_pressure) + rho_LR * b)).eval();

//            const auto dw_LS_dpGR =
//                    (-permeability_tensor * k_over_mu_LR *
//                            ( dNdx_p * p_LR + rho_LR * b)).eval();

            const auto dwGSdpGR = (permeability_tensor * dNdx_p).eval();
            const auto dwGSdpCap = (permeability_tensor * dNdx_p).eval();
            const auto dwGSdT = (permeability_tensor * dNdx_p).eval();
            const auto dwLSdpGR = (permeability_tensor * dNdx_p).eval();
            const auto dwLSdpCap = (permeability_tensor * dNdx_p).eval();
            const auto dwLSdT = (permeability_tensor * dNdx_p).eval();


            const auto w_GS =
                    (-permeability_tensor * k_over_mu_GR *
                    ( dNdx_p * gas_phase_pressure + rho_GR * b)).eval();

            // auxiliary operators:
            Laplace.noalias() =
                dNdx_p.transpose() * permeability_tensor * dNdx_p * w;

            auto const gravity_operator =
                (dNdx_p.transpose() * permeability_tensor * b * w).eval();

            // mass_operator = N_p_T * N_p * w;
            std::cout << " test 1\n";
            // output secondary variables
            ip_data.pressure_gas_linear = p_GR;
            ip_data.pressure_cap_linear = p_cap;
            ip_data.pressure_wet = p_LR;
            ip_data.density_gas = rho_GR;
            ip_data.density_liquid = rho_LR;
            ip_data.saturation = s_L;

            ip_data.velocity_liquid.noalias() = w_LS;
            ip_data.velocity_gas.noalias() = w_GS;


            std::cout << " test 1\n";
#ifdef DBG_OUTPUT
            std::cout << "==================================\n";
            std::cout << "            rho_SR : " << rho_SR << " \n";
            std::cout << "            rho_LR : " << rho_LR << " \n";
            std::cout << "            rho_GR : " << rho_GR << " \n";
            std::cout << "               rho : " << rho << " \n";
            std::cout << "==================================\n";
            std::cout << "               phi : " << phi << " \n";
            std::cout << "             phi_S : " << phi_S << " \n";
            std::cout << "==================================\n";
            std::cout << "             p_cap : " << p_GR << " \n";
            std::cout << "             p_cap : " << p_cap << " \n";
            std::cout << "               s_L : " << s_L << " \n";
            std::cout << "               s_G : " << s_G << " \n";
            std::cout << "             mu_LR : " << mu_LR << " \n";
            std::cout << "             mu_GR : " << mu_GR << " \n";
            std::cout << "==================================\n";
            std::cout << "    Gravity vector : \n";
            std::cout << "                 b : \n" << b << " \n";
            std::cout << "==================================\n";
            std::cout << "   volume strain e : " << e << " \n";
            std::cout << "==================================\n";
            std::cout << "                 C : " << "\n" << C << " \n\n";
            std::cout << "         sigma_eff : " << "\n" << sigma_eff << "\n\n";
            std::cout << "           alpha_B : " << alpha_B << " \n";
            std::cout << "      permeability : " << permeability << " \n";
            std::cout << "       perm_tensor : " << "\n" << permeability_tensor<< "\n\n";
            std::cout << "==================================\n";
            std::cout << "     drho_gr_dp_gr : " << drhoGRdpGR << " \n";
            std::cout << "     drho_lr_dp_lr : " << drhoLRdpLR << " \n";
            std::cout << "          drhoGRdT : " << drhoGRdT << " \n";
            std::cout << "          drhoLRdT : " << drhoLRdT << " \n";
            std::cout << "         beta_T_GR : " << beta_T_GR << " \n";
            std::cout << "         beta_T_LR : " << beta_T_LR << " \n";
            std::cout << "         beta_T_SR : " << beta_T_SR << " \n";
            std::cout << "==================================\n";
            std::cout << "            dsLdpc : " << dsLdpc << " \n";
            std::cout << "==================================\n";
            std::cout << "          k_rel_LR : " << k_rel_LR << " \n";
            std::cout << "          k_rel_GR : " << k_rel_GR << " \n";
            std::cout << "==================================\n";
            std::cout << "      k_over_mu_GR : " << k_over_mu_GR << " \n";
            std::cout << "      k_over_mu_LR : " << k_over_mu_LR << " \n";
            std::cout << "==================================\n";
            std::cout << "         beta_p_SR : " << beta_p_SR << " \n";
            std::cout << "         beta_p_PR : " << beta_p_PR << " \n";
            std::cout << "         beta_p_GR : " << beta_p_GR << " \n";
            std::cout << "         beta_p_LR : " << beta_p_LR << " \n";
            std::cout << "              cp_G : " << cp_G << " \n";
            std::cout << "              cp_L : " << cp_L << " \n";
            std::cout << "              cp_S : " << cp_S << " \n";
            std::cout << "            cp_eff : " << rho_cp_eff << " \n";
            std::cout << "==================================\n";
            std::cout << "        Velocities : \n";
            std::cout << "   ----------------------------------\n";
            std::cout << "              w_GS : \n" << w_GS << " \n";
            std::cout << "              w_LS : \n" << w_LS << " \n";
            std::cout << "==================================\n";
            std::cout << "           Laplace : \n";
            std::cout << "                 L :\n " << Laplace << " \n";
            std::cout << "==================================\n";

#endif
            Mgpg.noalias() += rho_GR * s_G * (phi * beta_p_GR + Sps) *
            		mass_operator;

            dMgpg_dpg.noalias() += s_G * drhoGRdpGR * (phi*beta_p_GR + Sps) *
            		mass_operator;
            dMgpg_dpc.noalias() -= dsLdpc * rho_GR * (phi*beta_p_GR + Sps) *
            		mass_operator;
            dMgpg_dT.noalias() += drhoGRdT  * (phi * beta_p_GR + Sps) *
            		mass_operator;

#ifdef DBG_OUTPUT
            std::cout << "   Mgpg:\n " << Mgpg << " \n";
            std::cout << "==================================\n";
#endif

            Mgpc.noalias() -=
            		rho_GR * (phi*dsLdpc + s_G * (s_L + p_cap*dsLdpc) * Sps) *
					mass_operator;

            dMgpc_dpg.noalias() -= drhoGRdpGR * (phi * dsLdpc + s_G * (s_L + p_cap*dsLdpc) * Sps) *
            		mass_operator;
            dMgpc_dpc.noalias() -= rho_GR * (phi * d2sLdpc2 + Sps *
            		(dsLdpc * (2 - 3 * s_L - p_cap * dsLdpc) + s_G * p_cap * d2sLdpc2)) *
            				mass_operator;
            dMgpc_dT.noalias()
 += drhoGRdT * (phi*dsLdpc + s_G * (s_L + p_cap*dsLdpc)* Sps) *
            		mass_operator;

#ifdef DBG_OUTPUT
            std::cout << "   Mgpc:\n " << Mgpc << " \n";
            std::cout << "==================================\n";
#endif



            MgT.noalias() -= N_p.transpose() * s_G * rho_GR *
            		(phi * beta_T_GR + beta_T_SR * (alpha_B - phi)) *
					N_p * w;

            dMgT_dpg.noalias() -= s_G * drhoGRdpGR * (phi * beta_T_GR * STs) *
            		mass_operator;
            dMgT_dpc.noalias() += dsLdpc * rho_GR * (phi * beta_T_GR * STs) *
            		mass_operator;
            dMgT_dT.noalias() -= s_G * drhoGRdT * (phi * beta_T_GR * STs) *
            		mass_operator;


#ifdef DBG_OUTPUT
            std::cout << "   MgT:\n " << MgT << " \n";
            std::cout << "==================================\n";
#endif

            Mgus.noalias() += s_G * rho_GR * alpha_B *
            		mass_operator2;
            dMgus_dpg.noalias() += s_G * drhoGRdpGR * alpha_B *
            		mass_operator2;
            dMgus_dpc.noalias() -= dsLdpc * rho_GR * alpha_B *
            		mass_operator2;
            dMgus_dpg.noalias() += s_G * drhoGRdT * alpha_B *
            		mass_operator2;


#ifdef DBG_OUTPUT
            std::cout << "   Mgus:\n " << Mgus << " \n";
            std::cout << "==================================\n";
#endif

            Mlpg.noalias() += N_p.transpose() * rho_LR * s_L *
                              (phi * beta_p_LR + Sps) *
                              N_p * w;
            dMlpg_dpg.noalias() += s_L * drhoLRdpLR * (phi * beta_p_LR + Sps) *
            		mass_operator;
            dMlpg_dpc.noalias() += (rho_LR * dsLdpc - s_L * drhoLRdpLR) * (phi * beta_p_LR + Sps) *
            		mass_operator;
            dMlpg_dT.noalias() += s_L * drhoLRdT * (phi * beta_p_LR + Sps) *
            		mass_operator;

#ifdef DBG_OUTPUT
            std::cout << "   Mlpg:\n " << Mlpg << " \n";
            std::cout << "==================================\n";
#endif

            Mlpc.noalias() += rho_LR * (dsLdpc * (phi - Sps*p_cap) - s_L*(phi*beta_p_LR + Sps)) *
            		mass_operator;
            dMlpc_dpg.noalias() += drhoLRdpLR  *
            		(phi * dsLdpc - phi_L * beta_p_LR - s_L * (s_L + p_cap*dsLdpc)* Sps) *
					mass_operator;
            dMlpc_dpc.noalias() += (rho_LR * (phi * d2sLdpc2 - dsLdpc *
            		( phi * beta_p_LR + Sps * ( 3 *s_L + p_cap * dsLdpc + s_L * p_cap* d2sLdpc2))) -
            		drhoLRdpLR * (dsLdpc * (phi - Sps*p_cap) - s_L*(phi*beta_p_LR + Sps))) *
					mass_operator;
            dMlpc_dT.noalias() += drhoLRdT*(dsLdpc * (phi - Sps*p_cap) - s_L*(phi*beta_p_LR + Sps)) *
            		mass_operator;

#ifdef DBG_OUTPUT
            std::cout << "     a_L: " << Sps << " \n";
            std::cout << "==================================\n";
            std::cout << "   Mlpc:\n " << Mlpc << " \n";
            std::cout << "==================================\n";
#endif

            MlT.noalias() -= s_L * rho_LR *
                             (phi * beta_T_LR + (alpha_B - phi) * beta_T_SR ) *
                             mass_operator;

            dMlT_dpg.noalias() -= s_L*drhoLRdpLR *
            		(phi * beta_T_LR +  (alpha_B - phi) * beta_T_SR) *
					mass_operator;
            dMlT_dpc.noalias() += (s_L * drhoLRdpLR - dsLdpc*rho_LR) *
            		(phi * beta_T_LR +  (alpha_B - phi) * beta_T_SR) *
					mass_operator;
            dMlT_dT.noalias() -= drhoLRdT * (phi * beta_T_LR +  (alpha_B - phi) * beta_T_SR) *
					mass_operator;

#ifdef DBG_OUTPUT
            std::cout << "   MlT:\n " << MlT << " \n";
            std::cout << "==================================\n";
#endif

            Mlus.noalias() += rho_LR * s_L * alpha_B *
            		mass_operator2;

            dMlus_dpg.noalias() += s_L * drhoLRdpLR * alpha_B *
            		mass_operator2;
            dMlus_dpc.noalias() += (dsLdpc*rho_LR - s_L * drhoLRdpLR) * alpha_B *
            		mass_operator2;
            dMlus_dT.noalias() += s_L * drhoLRdT * alpha_B *
            		mass_operator2;

#ifdef DBG_OUTPUT
            std::cout << "   Mlus:\n " << Mlus << " \n";
            std::cout << "==================================\n";
#endif

            Mepg.noalias() -= (phi_G*beta_T_GR + phi_L*beta_T_LR + phi_S*beta_T_SR) * T *
            		mass_operator;

            dMepg_dpg.noalias() += 0 * mass_operator;
            dMepg_dpc.noalias() += (phi * dsLdpc * (beta_T_GR - beta_T_LR) - phi_S * beta_T_SR) *
            		mass_operator;
            dMepg_dT.noalias() -= (phi * ( s_G * beta_T_GR + s_L * beta_T_LR) + phi_S * beta_T_SR) *
            		mass_operator;

#ifdef DBG_OUTPUT
            std::cout << "   Mepg:\n " << Mepg << " \n";
            std::cout << "=================================\n";
#endif
            Mepc.noalias() += (phi_L*beta_T_LR - (s_L + p_cap*dsLdpc)*phi_S*beta_T_SR) * T *
            		mass_operator;
            dMepc_dpg.noalias() += 0  *mass_operator;
            dMepc_dpc.noalias() += (phi * dsLdpc * beta_T_LR * T - ( 2 * dsLdpc + p_cap * d2sLdpc2)
            		* phi_S * beta_T_SR * T) * mass_operator;
            dMepc_dT.noalias() += (phi_L * beta_T_LR - (s_L + p_cap * dsLdpc) * phi_S * beta_T_SR) *
            		mass_operator;

#ifdef DBG_OUTPUT
            std::cout << "   Mepc:\n " << Mepc << " \n";
            std::cout << "=================================\n";
#endif

            MeT.noalias() += rho_cp_eff * mass_operator;
            dMeT_dpg.noalias() += (phi_G*drhoGRdpGR*cp_G + phi_L*drhoLRdpLR*cp_L + phi_S*drhoSRdpGR*cp_S) *
            		mass_operator;
            dMeT_dpc.noalias() += phi * (dsLdpc*(rho_LR*cp_L - rho_GR*cp_G) -
            		s_L * drhoLRdpLR*cp_L + phi_S*drhoSRdpCap*cp_S ) *
            		mass_operator;
            dMeT_dT.noalias() += (phi_G*drhoGRdT*cp_G + phi_L*drhoLRdT*cp_L + phi_S*drhoSRdT*cp_S) *
                        		mass_operator;

#ifdef DBG_OUTPUT
            std::cout << "   MeT:\n " << MeT << " \n";
            std::cout << "=================================\n";
#endif

            Lgpg.noalias() += rho_GR * k_over_mu_GR * Laplace;

            dLgpg_dpg.noalias() +=
            		k_over_mu_GR  * ( drhoGRdpGR - rho_GR/mu_GR*dmuGRdpGR)
            		* Laplace;
            dLgpg_dpc.noalias() +=
            		rho_GR/mu_GR*dkrelGdsL * dsLdpc * Laplace;
            dLgpg_dT.noalias() -= k_over_mu_GR * (drhoGRdT - rho_GR/mu_GR*dmuGRdT) *
            		Laplace;

#ifdef DBG_OUTPUT
            std::cout << "   Lgpg:\n " << Lgpg << " \n";
            std::cout << "==================================\n";
#endif
            Llpg.noalias() += rho_LR * k_over_mu_LR * Laplace;
            dLlpg_dpg.noalias() += k_over_mu_LR * (drhoLRdpLR - dmuLRdpLR*rho_LR/mu_LR) *
            		Laplace;
            dLlpg_dpc.noalias() += 1.0/mu_LR * ( (dkrelLdsL*dsLdpc - dmuLRdpLR*k_over_mu_LR) - drhoLRdpLR) *
            		Laplace;
            dLlpg_dT.noalias() = k_over_mu_LR * (drhoLRdT - dmuLRdT*rho_LR/mu_LR) *
            		Laplace;

#ifdef DBG_OUTPUT
            std::cout << "   Llpg:\n " << Llpg << " \n";
            std::cout << "==================================\n";
#endif

            Llpc.noalias() = - Llpg;
            dLlpc_dpg.noalias() = -dLlpg_dpg;
            dLlpc_dpc.noalias() = -dLlpg_dpc;
            dLlpc_dT.noalias() = -dLlpg_dT;

#ifdef DBG_OUTPUT
            std::cout << "   Llpc:\n " << Llpc << " \n";
            std::cout << "==================================\n";
#endif

            Aepg.noalias() -= N_p.transpose() *
                             (phi_G*beta_T_GR*w_GS.transpose() +
                              phi_L*beta_T_LR*w_LS.transpose()) * T *
                             dNdx_p * w;

            dAepgPg_dpg.noalias() -= (NpT*(beta_T_GR*w_GS.transpose() + beta_T_LR*w_LS.transpose())*T*dNdx_p+
            		NpT*beta_T_GR*T*gas_phase_pressure.transpose()*dNdx_p.transpose()*dwGSdpGR+
					NpT*beta_T_LR*T*gas_phase_pressure.transpose()*dNdx_p.transpose()*dwLSdpGR)*w;

            dAepgPg_dpc.noalias() -= (NpT*beta_T_GR*T*gas_phase_pressure.transpose()*dNdx_p.transpose()*dwGSdpCap+
					NpT*beta_T_LR*T*gas_phase_pressure.transpose()*dNdx_p.transpose()*dwLSdpCap)*w;

            dAepgPg_dT.noalias() -= (NpT*(beta_T_GR*w_GS.transpose() + beta_T_LR*w_LS.transpose())*
            		dNdx_p * gas_phase_pressure * Np +
					NpT*beta_T_GR*T*gas_phase_pressure.transpose()*dNdx_p.transpose()*dwGSdT+
					NpT*beta_T_LR*T*gas_phase_pressure.transpose()*dNdx_p.transpose()*dwLSdT)*w;

#ifdef DBG_OUTPUT
            std::cout << "   Aepg:\n " << Aepg << " \n";
            std::cout << "=================================\n";
#endif

            Aepc.noalias() =
            		N_p.transpose() * (phi_L*beta_T_LR*T*w_LS.transpose())
					* dNdx_p * w;

            dAepcPc_dpg += (NpT*beta_T_LR*T*capillary_pressure.transpose()*dNdx_p.transpose()*dwLSdpGR)*w;

            dAepcPc_dpc += (NpT*beta_T_LR*T*w_LS.transpose()*dNdx_p +
            		NpT*beta_T_LR*T*capillary_pressure.transpose()*dNdx_p.transpose()*dwLSdpGR)*w;

            dAepcPc_dT += (NpT*beta_T_LR*w_LS.transpose()*dNdx_p*capillary_pressure*Np +
                        		NpT*beta_T_LR*T*capillary_pressure.transpose()*dNdx_p.transpose()*dwLSdT)*w;

#ifdef DBG_OUTPUT
            std::cout << "   Aepc:\n " << Aepc << " \n";
            std::cout << "=================================\n";
#endif

            AeT.noalias() = N_p.transpose() *
                            (s_G*rho_GR*cp_G*w_GS.transpose() +
                             s_L*rho_LR*cp_L*w_LS.transpose()) *
                            dNdx_p * w;

            dAeTT_dpg += (NpT*(drhoGRdpGR*cp_G*w_GS.transpose()+drhoLRdpLR*cp_L*w_LS.transpose())*dNdx_p*temperature*Np +
            		NpT*rho_GR*cp_G*temperature.transpose()*dNdx_p.transpose()*dwGSdpGR +
					NpT*rho_LR*cp_L*temperature.transpose()*dNdx_p.transpose()*dwLSdpGR)*w;

            dAeTT_dpc += (NpT*drhoLRdpCap*cp_L*w_LS.transpose()*dNdx_p*temperature*Np +
            		NpT*rho_GR*cp_G*temperature.transpose()*dNdx_p.transpose()*dwGSdpCap +
					NpT*rho_LR*cp_L*temperature.transpose()*dNdx_p.transpose()*dwLSdpCap)*w;

            dAeTT_dT += (NpT*(rho_GR*cp_G*w_GS.transpose()+rho_LR*cp_L*w_LS.transpose())*dNdx_p +
            		NpT*(drhoGRdT*cp_G*w_GS.transpose()+drhoLRdT*cp_L*w_LS.transpose())*dNdx_p*temperature*Np +
					NpT*rho_GR*cp_G*temperature.transpose()*dNdx_p.transpose()*dwGSdT +
					NpT*rho_LR*cp_L*temperature.transpose()*dNdx_p.transpose()*dwLSdT)*w;

#ifdef DBG_OUTPUT
            std::cout << "   AeT:\n " << AeT << " \n";
            std::cout << "=================================\n";
#endif

            LeT.noalias() =
                dNdx_p.transpose() * conductivity_effective * dNdx_p * w;

            dLeT_dpc += (dNdx_p.transpose()*phi*dsLdpc*(lambda_L - lambda_G)*dNdx_p)*w;

#ifdef DBG_OUTPUT
            std::cout << "   LeT:\n " << LeT << " \n";
            std::cout << "=================================\n";
#endif

            Kupg.noalias() -= B.transpose() * alpha_B * identity2 * N_p * w;

#ifdef DBG_OUTPUT
            std::cout << "   Kupg:\n " << Kupg << " \n";
            std::cout << "==================================\n";
#endif
            Kupc.noalias() +=
                B.transpose() * alpha_B * identity2 * s_L * N_p * w;

            dKupc_dpc.noalias() -= (B.transpose()*alpha_B*identity2*dsLdpc*Np)*w;

#ifdef DBG_OUTPUT
            std::cout << "   Kupc:\n " << Kupc << " \n";
            std::cout << "==================================\n";
#endif

            fg.noalias() += rho_GR * rho_GR * k_over_mu_GR * gravity_operator;

            dfg_dpg += rho_GR * k_over_mu_GR * (2* drhoGRdpGR - rho_GR / mu_GR * dmuGRdpGR) *
            		gravity_operator;

            dfg_dpc += rho_GR * rho_GR / mu_GR * dkrelGdsL * dsLdpc *
            		gravity_operator;

            dfg_dT += rho_GR * k_over_mu_GR * (2*drhoGRdT + rho_GR/mu_GR*dmuGRdT) *
            		gravity_operator;

#ifdef DBG_OUTPUT
            std::cout << "   fg:\n " << fg << " \n";
            std::cout << "==================================\n";
#endif

            fl.noalias() += rho_LR * rho_LR * k_over_mu_LR * gravity_operator;

            dfl_dpg += rho_LR * k_over_mu_LR * (2 * drhoLRdpLR - rho_LR/mu_LR*dmuLRdpLR) *
            		gravity_operator;

            dfl_dpc += rho_LR / mu_LR * (rho_LR * (dkrelLdsL*dsLdpc + dmuLRdpLR*k_rel_LR/mu_LR) - 2*drhoLRdpLR*k_rel_LR) *
            		gravity_operator;

            dfl_dT += rho_LR * k_over_mu_LR * (2*drhoLRdT - rho_LR/mu_LR * dmuLRdT) *
            		gravity_operator;

#ifdef DBG_OUTPUT
            std::cout << "   fl:\n " << fl << " \n";
            std::cout << "==================================\n";
#endif

            fu.noalias() += N_u_op.transpose() * rho * b * w;

            dfu_dpg -= N_u_op.transpose() * drhodpGR * b * w;

            dfu_dpc -= N_u_op.transpose() * drhodpCap * b * w;

            dfu_dT -= (B.transpose() / beta_p_PR * beta_T_SR * identity2 -
            		N_u_op.transpose() * drhodT * b ) * w;

            dfu_du += B.transpose() * C * B * w;

#ifdef DBG_OUTPUT
            std::cout << "   fu:\n " << fu << " \n";
            std::cout << "==================================\n";
#endif

            // gas phase residuum
            rg += Mgpg * gas_phase_pressure_dot + Mgpc * capillary_pressure_dot +
            		MgT * temperature_dot + Mgus * displacement_dot + Lgpg * gas_phase_pressure
					- fg;

            // gas phase residua derivatives
            drg_dpg += (dMgpg_dpg * gas_phase_pressure_dot + dMgpc_dpg * capillary_pressure_dot +
            		dMgus_dpg * displacement_dot + dMgT_dpg * temperature_dot +
					dLgpg_dpg * gas_phase_pressure - dfg_dpg ) * Np +
							Lgpg + Mgpg / dt;

            drg_dpc += (dMgpg_dpc * gas_phase_pressure_dot + dMgpc_dpc * capillary_pressure_dot +
            		dMgus_dpc * displacement_dot + dMgT_dpc * temperature_dot +
					dLgpg_dpc * gas_phase_pressure - dfg_dpc ) * Np + Mgpc / dt;

            drg_dT += (dMgpg_dT * gas_phase_pressure_dot + dMgpc_dT * capillary_pressure_dot +
            		dMgus_dT * displacement_dot + dMgT_dT * temperature_dot +
					dLgpg_dT * gas_phase_pressure - dfg_dpc ) * Np + MgT / dt;

            drg_dus += Mgus / dt;

            // liquid phase residuum
            rl += Mlpg * gas_phase_pressure_dot + Mlpc * capillary_pressure_dot +
            		MlT * temperature_dot + Mlus * displacement_dot + Llpg * gas_phase_pressure +
					Llpc * capillary_pressure - fl;

            // liquid phase residua derivatives
            drl_dpg += (dMlpg_dpg * gas_phase_pressure_dot + dMlpc_dpg * capillary_pressure_dot +
            		dMlus_dpg * displacement_dot + dMlT_dpg * temperature_dot +
					dLlpg_dpg * gas_phase_pressure + dLlpc_dpg * capillary_pressure - dfl_dpg ) * Np +
							Llpg + Mlpg / dt;

            drl_dpc += (dMlpg_dpc * gas_phase_pressure_dot + dMlpc_dpc * capillary_pressure_dot +
            		dMlus_dpc * displacement_dot + dMlT_dpc * temperature_dot +
					dLlpg_dpc * gas_phase_pressure + dLlpc_dpc * capillary_pressure - dfl_dpc ) * Np +
							Llpc + Mlpc / dt;

            drl_dT += (dMlpg_dT * gas_phase_pressure_dot + dMlpc_dT * capillary_pressure_dot +
            		dMlus_dT * displacement_dot + dMlT_dT * temperature_dot +
					dLlpg_dT * gas_phase_pressure + dLlpc_dT * capillary_pressure - dfl_dT ) * Np +
							MlT / dt;

            drl_dus += Mlus / dt;

            // energy residuum
            re += Mepg * gas_phase_pressure_dot + Mepc * capillary_pressure_dot +
            		MeT * temperature_dot + Aepg * gas_phase_pressure + Aepc * capillary_pressure +
					(AeT + LeT) * temperature - fe;

            // energy residua derivatives
            dre_dpg += (dMepg_dpg * gas_phase_pressure_dot + dMepc_dpg * capillary_pressure_dot +
            		dMeT_dpg * temperature_dot - dfe_dpg) * Np + dAepgPg_dpg + dAepcPc_dpg +
            				dAeTT_dpg + Mepg /dt;

            dre_dpc += (dMepg_dpc * gas_phase_pressure_dot + dMepc_dpc* capillary_pressure_dot +
            		dMeT_dpc * temperature_dot + dLeT_dpc * temperature - dfe_dpc) * Np + dAepgPg_dpc + dAepcPc_dpc +
            				dAeTT_dpc + Mepc /dt;

            dre_dT += (dMepg_dT * gas_phase_pressure_dot + dMepc_dT* capillary_pressure_dot +
            		dMeT_dT * temperature_dot - dfe_dpc) * Np + dAepgPg_dT + dAepcPc_dT +
            				dAeTT_dT + LeT + Mepc /dt;


            std::cout << "   dre_dT :\n " << dre_dT << " \n";
            std::cout << "==================================\n";





            // displacement residuum
            ru += Kupg * gas_phase_pressure + Kupc * capillary_pressure - fu;

            // displacement residua derivatives
            dru_dpg += Kupg - dfu_dpg * Np;

            dru_dpc += ( dKupc_dpc * capillary_pressure - dfu_dpc )*  Np + Kupc;

            dru_dT -= dfu_dT * Np;

            dru_dus -= dfu_du;

        }

        std::cout << "#####################\n";
        std::cout << " Jacobian (analytical): \n";
        std::cout << J << "\n\n";
        std::cout << "=====================\n";
        std::cout << " residuum: \n\n";
        std::cout << r << "\n\n";
        std::cout << "=====================\n";


        OGS_FATAL("Intended stop.");


    }

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N_u = _secondary_data.N_u[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N_u.data(), N_u.size());
    }

private:
    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        static const int kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& sigma_eff = _ip_data[ip].sigma_eff;
            cache_mat.col(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_eff);
        }

        return cache;
    }

    virtual std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& eps = _ip_data[ip].eps;
            cache_mat.col(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(eps);
        }

        return cache;
    }

    std::vector<double> const& getIntPtSaturation(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        cache.resize(n_integration_points);
        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            cache[ip] = _ip_data[ip].saturation;
        }

        return cache;
    }

    std::vector<double> const& getIntPtWetPressure(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        cache.resize(n_integration_points);
        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            cache[ip] = _ip_data[ip].pressure_wet;
        }

        return cache;
    }

    std::vector<double> const& getIntPtDensityGas(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        cache.resize(n_integration_points);
        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            cache[ip] = _ip_data[ip].density_gas;
        }

        return cache;
    }

    std::vector<double> const& getIntPtDensityLiquid(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        cache.resize(n_integration_points);
        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            cache[ip] = _ip_data[ip].density_liquid;
        }

        return cache;
    }

    std::vector<double> const& getIntPtDarcyVelocityLiquid(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, DisplacementDim, n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            cache_mat.col(ip) = _ip_data[ip].velocity_liquid;
        }

        return cache;
    }

    std::vector<double> const& getIntPtPressureGas(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        cache.resize(n_integration_points);
        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            cache[ip] = _ip_data[ip].pressure_gas_linear;
        }

        return cache;
    }

    std::vector<double> const& getIntPtPressureCap(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        cache.resize(n_integration_points);
        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            cache[ip] = _ip_data[ip].pressure_cap_linear;
        }

        return cache;
    }

    std::vector<double> const& getIntPtDarcyVelocityGas(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, DisplacementDim, n_integration_points);
        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            cache_mat.col(ip) = _ip_data[ip].velocity_gas;
        }

        return cache;
    }

private:
    TH2MProcessData<DisplacementDim>& _process_data;

    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesTypeDisplacement,
                             ShapeMatricesTypePressure, DisplacementDim,
                             ShapeFunctionDisplacement::NPOINTS>;
    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    bool const _is_axially_symmetric;
    SecondaryData<
        typename ShapeMatricesTypeDisplacement::ShapeMatrices::ShapeType>
        _secondary_data;

    static const int Np_intPoints = ShapeFunctionPressure::NPOINTS;
    static const int gas_pressure_index = 0;
    static const int cap_pressure_index = 1 * Np_intPoints;
    static const int temperature_index = 2 * Np_intPoints;
    static const int displacement_index = 3 * Np_intPoints;

    static const int gas_pressure_size = Np_intPoints;
    static const int cap_pressure_size = Np_intPoints;
    static const int temperature_size = Np_intPoints;

    static const int Nu_intPoints =
        ShapeFunctionDisplacement::NPOINTS * DisplacementDim;
    static const int displacement_size = Nu_intPoints;
    static const int kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
};

}  // namespace TH2M
}  // namespace ProcessLib
