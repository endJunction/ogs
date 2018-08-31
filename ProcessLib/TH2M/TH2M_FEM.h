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
            // std::cout << " u:\n "
            //<< "???"
            //<< "\n";
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

            auto const s_L = 1 - 1.9722e-11 * pow (p_cap, 2.4279);

            if (s_L < 0.9)
            {
                OGS_FATAL ("s_L dropped below 0.9!");
            }
            auto const s_G = 1 - s_L;

            auto const s_e = (s_L - 0.2) / (1.0 - 0.2);



#ifdef DBG_OUTPUT
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

            ip_data.updateConstitutiveRelation(t, x_position, dt, displacement,
                    T, p_GR);

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
            std::cout << "==================================\n";
            std::cout << "   drho_gr_dp_gr : " << drhoGRdpGR << " \n";
            std::cout << "   drho_lr_dp_lr : " << drhoLRdpLR << " \n";
            std::cout << "==================================\n";
#endif
            auto const dsLdpc = MPL::getScalarDerivative(
                medium.property(MPL::PropertyEnum::saturation),
                variables, MPL::Variables::p_cap);
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

            auto const k_rel_LR = 1.0 - 2.207*pow((1.0 - s_L), 1.0121);

            auto const min_k_rel_GR = 0.000;

            auto const k_rel_GR = (1.0 - s_e) * (1 - s_e)
                    * (1.0 - pow(s_e, (5./3.))) + min_k_rel_GR;

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

            const double beta_p_SR = getScalar(
                solid_phase.property(MPL::PropertyEnum::compressibility));
            const double beta_p_PR =
                getScalar(medium.property(MPL::PropertyEnum::compressibility));
            const double beta_p_GR = 1. / rho_GR * drhoGRdpGR;
            const double beta_p_LR = 1. / rho_LR * drhoLRdpLR;

#ifdef DBG_OUTPUT
            std::cout << "   beta_p_SR: " << beta_p_SR << " \n";
            std::cout << "   beta_p_PR: " << beta_p_PR << " \n";
            std::cout << "   beta_p_GR: " << beta_p_GR << " \n";
            std::cout << "   beta_p_LR: " << beta_p_LR << " \n";
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

            const double beta_T_GR = drhoGRdT / rho_GR;

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
            MgT.noalias() += N_p.transpose() * s_G * rho_GR *
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

            Mlpc.noalias() +=
                N_p.transpose() * rho_LR *
                ((phi - s_L * (alpha_B - phi) * beta_p_SR * p_cap) * dsLdpc -
                 phi * s_L * beta_p_LR -
                 (alpha_B - phi) * beta_p_SR * s_L * s_L) *
                N_p * w;

#ifdef DBG_OUTPUT
            std::cout << "   Mlpc:\n " << Mlpc << " \n";
            std::cout << "==================================\n";
#endif

            MlT.noalias() += N_p.transpose() * s_L * rho_LR *
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

            const double cp_eff = phi_G * rho_GR * cp_G +
                                  phi_L * rho_LR * cp_L + phi_S * rho_SR * cp_S;

#ifdef DBG_OUTPUT
            std::cout << "   cp_G: " << cp_G << " \n";
            std::cout << "   cp_L: " << cp_L << " \n";
            std::cout << "   cp_S: " << cp_S << " \n";

            std::cout << " cp_eff: " << cp_eff << " \n";
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
                (phi_L*beta_T_LR * T + phi_S*s_L*beta_T_SR*T + p_cap*phi*dsLdpc) * N_p *
                w;

#ifdef DBG_OUTPUT
            std::cout << "   Mepc:\n " << Mepc << " \n";
            std::cout << "=================================\n";
#endif

            MeT.noalias() += N_p.transpose() * cp_eff * N_p * w;

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

            Aepg.noalias() = N_p.transpose() *
                             (s_G*beta_T_GR*w_GS.transpose() +
                              s_L*beta_T_LR*w_LS.transpose()) * T *
                             dNdx_p * w;

            Kepg.noalias() -= Aepg;

#ifdef DBG_OUTPUT
            std::cout << "   Aepg:\n " << Aepg << " \n";
            std::cout << "   Kepg:\n " << Kepg << " \n";
            std::cout << "=================================\n";
            // Aepc(0, 0) = 3.4;
#endif

            Aepc.noalias() =
                N_p.transpose() *
                ( p_GR*dsLdpc*w_GS.transpose() +
                        (s_L*beta_T_LR*T-p_LR*dsLdpc)*w_LS.transpose()) *
                        dNdx_p * w;

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

            const double lambda_eff = phi_G * rho_GR * lambda_G +
                                      phi_L * rho_LR * lambda_L +
                                      phi_S * rho_SR * lambda_S;

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

            Bu.noalias() -=
                (B.transpose() * sigma_eff - N_u_op.transpose() * rho * b) * w;

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

               OGS_FATAL ("##########################################");
#endif

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
