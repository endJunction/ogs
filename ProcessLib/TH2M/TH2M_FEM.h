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

        auto const material_id = _process_data.material_ids
                ? (*_process_data.material_ids)[_element.getID()]
                : 0;
        try
        {
            medium = _process_data.media.at(material_id).get();
            if (medium==nullptr)
            {
                OGS_FATAL("Medium for material ID %d was not created.",
                        material_id);
            }
        }
        catch(std::out_of_range)
        {
            OGS_FATAL("Requested material ID %d not found in MPL::media.",
                    material_id);
        }

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

        auto const element_id = _element.getID();
        x_position.setElementID(element_id);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        double const& dt = _process_data.dt;

        auto const& solid_phase = medium->phase(0);
        auto const& liquid_phase = medium->phase(1);
        auto const& gas_phase = medium->phase(2);

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

            static int const KelvinVectorSize =
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value;

            auto const& identity2 =
                MathLib::KelvinVector::Invariants<KelvinVectorSize>::identity2;

            auto const& mass_operator = N_p.transpose() * N_p * w;
            auto const& mass_operator2 = N_p.transpose() * identity2.transpose() * B * w;


            auto& eps = ip_data.eps;
            auto const& sigma_eff = ip_data.sigma_eff;


            // primary unknowns at integration point
            auto const p_GR = gas_phase_pressure.dot(N_p);
            auto const p_cap = capillary_pressure.dot(N_p);
            auto const T = temperature.dot(N_p);
            auto const u = N_u_op * displacement;

            typename ShapeMatricesTypeDisplacement::GlobalDimVectorType u_ip(
                    DisplacementDim);
            for (int i = 0; i < u_ip.size(); ++i)
            {
                NumLib::shapeFunctionInterpolate(
                        displacement.segment(i * ShapeFunctionDisplacement::NPOINTS,
                                ShapeFunctionDisplacement::NPOINTS),
                                N_u, u_ip.coeffRef(i));
            }
            eps.noalias() = B * displacement;
            typename ShapeMatricesTypeDisplacement::GlobalDimVectorType
            const grad_p_GR = dNdx_p * gas_phase_pressure;
            typename ShapeMatricesTypeDisplacement::GlobalDimVectorType
            const grad_p_Cap = dNdx_p * capillary_pressure;

            const double p_LR = p_GR - p_cap;

            // insert all primary variables into one object
            MPL::VariableArray variables;
            variables[MPL::Variables::capillary_pressure] = p_cap;
            variables[MPL::Variables::phase_pressure] = p_GR;
            variables[MPL::Variables::temperature] = T;
            // todo: displacement
            // insert all primary variables into one object for GPML

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


            double const permeability =
                MPL::getScalar(medium->property(MPL::permeability));

            GlobalDimMatrixType permeability_tensor =
                GlobalDimMatrixType::Zero(DisplacementDim, DisplacementDim);

            permeability_tensor.diagonal().setConstant(permeability);

            auto const alpha_B = MPL::getScalar(
                    medium->property(MPL::PropertyEnum::biot_coefficient),
                    variables);

            auto const phi = MPL::getScalar(
                medium->property(MPL::PropertyEnum::porosity), variables);
            auto const phi_S = 1. - phi;

            // Saturation
//          auto const s_L =
//              MPL::getScalar(medium->property(MPL::PropertyEnum::saturation),
//                               variables);
            auto const s_L_r = 0.2;
            auto const s_a = -1.9722e-11;
            auto const s_b = 2.4279;

//          auto const s_L = std::max(s_L_r,1. + s_a * std::pow (std::max(0.,p_cap), s_b));
            auto const s_L = p_cap*p_cap + 2*p_cap + 1;
            auto const dsLdpc = 2*p_cap + 2;
            auto const d2sLdpc2 = 2.0;

            auto const s_G = 1 - s_L;
            auto const s_e = (s_L - s_L_r) / (1.0 - s_L_r);

            const double phi_L = s_L * phi;
            const double phi_G = s_G * phi;

//          auto const dsLdpc = MPL::getScalarDerivative(
//              medium->property(MPL::PropertyEnum::saturation),
//              variables, MPL::Variables::p_cap);

//          auto const dsLdpc = s_a*s_b*std::pow(std::max(0.,p_cap), s_b - 1.0);
//            auto const d2sLdpc2 = s_a*(s_b*(s_b - 1)) * std::pow(std::max(0.,p_cap), s_b - 2.0);


//          auto const k_rel_LR = MPL::getPair(
//              mediumte.property(MPL::PropertyEnum::relative_permeability),
//              variables)[0];
//          auto const k_rel_GR = MPL::getPair(
//              medium->property(MPL::PropertyEnum::relative_permeability),
//              variables)[1];

//            auto const k_rel_LR = 1.0 - 2.207*std::pow((1.0 - s_L), 1.0121);
            auto const k_rel_LR = 1.0;
            auto const min_k_rel_GR = 0.0001;

            auto const k_rel_b = 5./3.;
//            auto const k_rel_GR = (1.0 - s_e) * (1 - s_e)
//                                                * (1.0 - std::pow(s_e, k_rel_b)) + min_k_rel_GR;


            auto const k_rel_GR = 2*s_L*s_L + 3*s_L + 3;

//            auto const dkrelGdsE = 1./s_e * (1 - s_e) *
//                  ((k_rel_b+2)*std::pow(s_e,k_rel_b+1)-k_rel_b*std::pow(s_e, k_rel_b) - 2*s_e);
//
//            auto const dsEdsL = 1. / (1. - s_L_r);
//
//            auto const dkrelGdsL = dkrelGdsE * dsEdsL;
            auto const dkrelGdsL = 4 * s_L + 3;
            auto const dkrelLdsL = 0.0;

            auto const rho_SR =
                    MPL::getScalar(solid_phase.property(MPL::PropertyEnum::density),
                            variables);
            auto const rho_LR = MPL::getScalar(
                    liquid_phase.property(MPL::PropertyEnum::density),
                    variables);
            auto const rho_GR_0 =
                    MPL::getScalar(gas_phase.property(MPL::PropertyEnum::density),
                            variables);

            const double beta_p_SR = getScalar(
                solid_phase.property(MPL::PropertyEnum::compressibility));
            const double beta_p_PR =
                getScalar(medium->property(MPL::PropertyEnum::compressibility));
            const double beta_p_GR = getScalar(
                gas_phase.property(MPL::PropertyEnum::compressibility));
            const double beta_p_LR = getScalar(
                liquid_phase.property(MPL::PropertyEnum::compressibility));

            // expansivities
            const double beta_T_GR = getScalar(
                gas_phase.property(MPL::PropertyEnum::thermal_expansivity));
            const double beta_T_LR = getScalar(
                liquid_phase.property(MPL::PropertyEnum::thermal_expansivity));
            const double beta_T_SR = getScalar(
                solid_phase.property(MPL::PropertyEnum::thermal_expansivity));


            const auto p_GR_0 = 100000;
            const auto T_GR_0 = 30;

            auto rho_GR = rho_GR_0 * (1 + beta_p_GR * (p_GR - p_GR_0) + beta_T_GR * (T - T_GR_0));

            auto const& b = _process_data.specific_body_force;

            variables[MPL::Variables::liquid_density] = rho_LR;
            variables[MPL::Variables::gas_density] = rho_GR;

//            auto const drhoGRdpGR = MPL::getScalarDerivative(
//                    gas_phase.property(MPL::density), variables, MPL::p_GR);
//            auto const drhoGRdpGR = rho_GR_0*beta_p_GR;
//            auto const drhoGRdpGR = rho_GR_0 * beta_p_GR;

//            auto const drhoLRdpLR =
//                                MPL::getScalarDerivative(liquid_phase.property(MPL::density),
//                                        variables, MPL::p_LR);
//            auto const drhoLRdpLR = 0.;
//
//            auto const drhoLRdpGR = drhoLRdpLR;
//            auto const drhoLRdpCap = - drhoLRdpLR;

//            auto const drhoGRdT = MPL::getScalarDerivative(
//                    gas_phase.property(MPL::density), variables, MPL::T);
//            auto const drhoGRdT = 0.;

//            auto const drhoLRdT = MPL::getScalarDerivative(
//                    liquid_phase.property(MPL::density), variables, MPL::T);
//            auto const drhoLRdT = 0.;
//
//            auto const drhoSRdpGR = 0.;
//            auto const drhoSRdpCap = 0.;
//            auto const drhoSRdT = 0.;

            // volume average density
            auto const rho =
                    phi_S * rho_SR + phi * (s_L * rho_LR + s_G * rho_GR);
//            auto const drhodpGR = phi_G*drhoGRdpGR + phi_L*drhoLRdpGR + phi_S*drhoSRdpGR;
//            auto const drhodpCap = phi*dsLdpc*(rho_LR-rho_GR) + phi_L*drhoLRdpCap + phi_S*drhoSRdpCap;
//            auto const drhodT = phi_G*drhoGRdT + phi_L*drhoLRdT + phi_S*drhoSRdT;
            // Viscosities
            auto const mu_LR = MPL::getScalar(
                liquid_phase.property(MPL::PropertyEnum::viscosity),
                variables);
            auto const mu_GR =
                MPL::getScalar(gas_phase.property(MPL::PropertyEnum::viscosity),
                               variables);

//            auto const dmuLRdpGR = 0.0;
//            auto const dmuLRdpCap = 0.0;
//            auto const dmuLRdT= 0.0;
//            auto const dmuGRdpGR = 0.0;
//            auto const dmuGRdT = 0.0;

            // heat capacities
            const double cp_G = getScalar(
                gas_phase.property(MPL::PropertyEnum::specific_heat_capacity));
            const double cp_L = getScalar(liquid_phase.property(
                MPL::PropertyEnum::specific_heat_capacity));
            const double cp_S = getScalar(solid_phase.property(
                MPL::PropertyEnum::specific_heat_capacity));

//            const double dcpGdpGR = 0.;
//            const double dcpLdpGR = 0.;
//            const double dcpSdpGR = 0.;
//            const double dcpLdpCap = 0.;
//            const double dcpSdpCap = 0.;
//            const double dcpGdT = 0.;
//            const double dcpLdT = 0.;
//            const double dcpSdT = 0.;

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

//            const double dLambdaGRdpGR = 0.;
//            const double dLambdaLRdpGR = 0.;
//            const double dLambdaLRdpCap = 0.;
//            const double dLambdaSRdpGR = 0.;
//            const double dLambdaSRdpCap = 0.;
//            const double dLambdaGRdT = 0.;
//            const double dLambdaLRdT = 0.;
//            const double dLambdaSRdT = 0.;

            const double lambda_eff = phi_G * lambda_G +
                                      phi_L * lambda_L +
                                      phi_S * lambda_S;

            // isotropic heat conductivity
            const auto conductivity_effective =
                lambda_eff * conductivity_tensor;

//            const auto dLambdadpGR = (phi_G * dLambdaGRdpGR +
//                    phi_L * dLambdaLRdpGR + phi_S * dLambdaSRdpGR) *
//                            conductivity_tensor;
//
//            const auto dLambdadpCap = (phi * dsLdpc * (lambda_L - lambda_G) +
//                    phi_L * dLambdaLRdpCap + phi_S * dLambdaSRdpCap) *
//                            conductivity_tensor;
//
//            const auto dLambdadT = (phi_G * dLambdaGRdT +
//                                phi_L * dLambdaLRdT + phi_S * dLambdaSRdT) *
//                                        conductivity_tensor;

            /*
             * derived variables
             */
            const double rho_cp_eff = phi_G * rho_GR * cp_G +
                    phi_L * rho_LR * cp_L + phi_S * rho_SR * cp_S;


//            const double drhocpGdpGR = phi_G * (drhoGRdpGR*cp_G + rho_GR*dcpGdpGR) +
//                    phi_L * (drhoLRdpGR * cp_L + rho_LR * dcpLdpGR) +
//                    phi_S * (drhoSRdpGR * cp_S + rho_SR * dcpSdpGR);
//
//            const double drhocpGdpCap = phi * dsLdpc * (rho_LR*cp_L + rho_GR*cp_G) +
//                    phi_L * (drhoLRdpCap * cp_L + rho_LR * dcpLdpCap) +
//                    phi_S * (drhoSRdpCap * cp_S + rho_SR * dcpSdpCap);
//
//            const double drhocpGdT = phi_G * (drhoGRdT*cp_G + rho_GR*dcpGdT) +
//                    phi_L * (drhoLRdT * cp_L + rho_LR * dcpLdT) +
//                    phi_S * (drhoSRdT * cp_S + rho_SR * dcpSdT);
//

            const double k_over_mu_GR = k_rel_GR / mu_GR;
            const double k_over_mu_LR = k_rel_LR / mu_LR;

            const auto Sps = (alpha_B - phi) * beta_p_SR;
            const auto STs = (alpha_B - phi) * beta_T_SR;

            // darcy-velocities
            const auto w_LS =
                    (-permeability_tensor * k_over_mu_LR *
                    ( dNdx_p * (gas_phase_pressure - capillary_pressure) - rho_LR * b)).eval();

            const auto w_GS =
                    (-permeability_tensor * k_over_mu_GR *
                    ( dNdx_p * gas_phase_pressure - rho_GR * b)).eval();

            // Elasticity matrix:
            auto C = ip_data.updateConstitutiveRelation(t, x_position, dt, displacement,
                    T, p_GR);

            // auxiliary operators:
            Laplace.noalias() =
                dNdx_p.transpose() * permeability_tensor * dNdx_p * w;

            auto const gravity_operator =
                (dNdx_p.transpose() * permeability_tensor * b * w).eval();

            // mass_operator = N_p_T * N_p * w;
            // output secondary variables
            ip_data.pressure_gas_linear = p_GR;
            ip_data.pressure_cap_linear = p_cap;
            ip_data.pressure_wet = p_LR;
            ip_data.rel_perm_gas = k_rel_GR;
            ip_data.rel_perm_liquid = k_rel_LR;
            ip_data.density_gas = rho_GR;
            ip_data.density_liquid = rho_LR;
            ip_data.saturation = s_L;

            ip_data.velocity_liquid.noalias() = w_LS;
            ip_data.velocity_gas.noalias() = w_GS;

#ifdef DBG_OUTPUT
            std::cout << "==================================\n";
            std::cout << "            rho_SR : " << rho_SR << " \n";
            std::cout << "            rho_LR : " << rho_LR << " \n";
            std::cout << "            rho_GR : " << rho_GR << " \n";
            std::cout << "               rho : " << rho << " \n";
            std::cout << "==================================\n";
            std::cout << "               phi : " << phi << " \n";
            std::cout << "             phi_G : " << phi_G << " \n";
            std::cout << "             phi_L : " << phi_L << " \n";
            std::cout << "             phi_S : " << phi_S << " \n";
            std::cout << "==================================\n";
            std::cout << "             p_GR  : " << p_GR << " \n";
            std::cout << "             p_cap : " << p_cap << " \n";
            std::cout << "----------------------------------\n";
            std::cout << "               s_L : " << s_L << " \n";
            std::cout << "               s_G : " << s_G << " \n";
            std::cout << "----------------------------------\n";
            std::cout << "             mu_LR : " << mu_LR << " \n";
            std::cout << "             mu_GR : " << mu_GR << " \n";
            std::cout << "==================================\n";
            std::cout << "    Gravity vector : \n";
            std::cout << "                 b : \n" << b << " \n";
            std::cout << "==================================\n";
            std::cout << "   volume strain e : " << /*e <<*/ " \n";
            std::cout << "==================================\n";
            std::cout << "                 C : " << "\n" << C << " \n\n";
            std::cout << "----------------------------------\n";
            std::cout << "         sigma_eff : " << "\n" << sigma_eff << "\n\n";
            std::cout << "----------------------------------\n";
            std::cout << "           alpha_B : " << alpha_B << " \n";
            std::cout << "----------------------------------\n";
            std::cout << "      permeability : " << permeability << " \n";
            std::cout << "----------------------------------\n";
            std::cout << "       perm_tensor : " << "\n" << permeability_tensor<< "\n\n";
            std::cout << "==================================\n";
            std::cout << "     drho_gr_dp_gr : "  << drhoGRdpGR << "\n";
            std::cout << "     drho_lr_dp_lr : "  << "? \n";
            std::cout << "          drhoGRdT : "  << "? \n";
            std::cout << "          drhoLRdT : "  << "? \n";
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
            std::cout << "         beta_p_PR : " << /*beta_p_PR << */ " \n";
            std::cout << "         beta_p_GR : " << beta_p_GR << " \n";
            std::cout << "         beta_p_LR : " << beta_p_LR << " \n";
            std::cout << "              cp_G : " << cp_G << " \n";
            std::cout << "              cp_L : " << cp_L << " \n";
            std::cout << "              cp_S : " << cp_S << " \n";
            std::cout << "            cp_eff : " << rho_cp_eff << " \n";
            std::cout << "==================================\n";
            std::cout << "          lambda_G : " << lambda_G << " \n";
            std::cout << "          lambda_L : " << lambda_L << " \n";
            std::cout << "          lambda_S : " << lambda_S << " \n";
            std::cout << "        lamdba_eff : " << lambda_eff << " \n";
            std::cout << "==================================\n";
            std::cout << "        Velocities : \n";
            std::cout << "   ----------------------------------\n";
            std::cout << "              w_GS : \n" << w_GS << " \n";
            std::cout << "              w_LS : \n" << w_LS << " \n";
            std::cout << "==================================\n";
            std::cout << "           Laplace : \n";
            std::cout << "                 L :\n " << Laplace << " \n";
            std::cout << "==================================\n";
            OGS_FATAL("CD");
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
                    ((phi + s_G * Sps * p_cap) * dsLdpc +
                            s_G * s_L * Sps) *
                            N_p * w;

#ifdef DBG_OUTPUT
            std::cout << "   Mgpc:\n " << Mgpc << " \n";
            std::cout << "==================================\n";
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
                              (phi * beta_p_LR + Sps) * N_p * w;

#ifdef DBG_OUTPUT
            std::cout << "   Mlpg:\n " << Mlpg << " \n";
            std::cout << "==================================\n";
#endif

            const auto c = phi * dsLdpc - phi_L*beta_p_LR - s_L *
                    (s_L + p_cap * dsLdpc) * Sps;

            Mlpc.noalias() +=
            		N_p.transpose() * rho_LR * c *  N_p * w;

#ifdef DBG_OUTPUT
     //       std::cout << "     a_L: " << a_L << " \n";
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

//            Aepg.noalias() = - N_p.transpose() *
//                             (phi_G*beta_T_GR*w_GS.transpose() +
//                              phi_L*beta_T_LR*w_LS.transpose()) * T *
//                             dNdx_p * w;
//
//            Kepg.noalias() += Aepg;

#ifdef DBG_OUTPUT
            std::cout << "   Aepg:\n " << Aepg << " \n";
            std::cout << "   Kepg:\n " << Kepg << " \n";
            std::cout << "=================================\n";
            // Aepc(0, 0) = 3.4;
#endif

//            Aepc.noalias() = N_p.transpose() * (phi_L*beta_T_LR*T*w_LS.transpose())
//					* dNdx_p * w;
//
//            Kepc.noalias() += Aepc;

#ifdef DBG_OUTPUT
            std::cout << "   Aepc:\n " << Aepc << " \n";
            std::cout << "   Kepc:\n " << Kepc << " \n";
            std::cout << "=================================\n";
#endif

//            AeT.noalias() = N_p.transpose() *
//                            (s_G*rho_GR*cp_G*w_GS.transpose() +
//                             s_L*rho_LR*cp_L*w_LS.transpose()) *
//                            dNdx_p * w;

#ifdef DBG_OUTPUT
            std::cout << "   AeT:\n " << AeT << " \n";
            std::cout << "=================================\n";
#endif

//            LeT.noalias() =
//                    dNdx_p.transpose() * conductivity_effective * dNdx_p * w;
//
//            KeT.noalias() += AeT + LeT;

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

            Bu.noalias() -= (B.transpose() * sigma_eff - N_u_op.transpose() *
                    rho * b) * w;

#ifdef DBG_OUTPUT
            std::cout << "   Bu:\n " << Bu << " \n";
            std::cout << "==================================\n";
#endif
        }


#ifdef DBG_OUTPUT
        std::cout << "#####################\n";
        std::cout << "Mgpg : \n" << Mgpg << "\n";
        std::cout << "Mgpc : \n" << Mgpc << "\n";
        std::cout << "MgT  : \n" << MgT << "\n";
        std::cout << "Mgus : \n" << Mgus << "\n";
        std::cout << "Lgpg : \n" << Kgpg << "\n";
        std::cout << "fg   : \n" << Bg << "\n";
        std::cout << "---------------------\n";
        std::cout << "Mlpg : \n" << Mlpg << "\n";
        std::cout << "Mlpc : \n" << Mlpc << "\n";
        std::cout << "MlT  : \n" << MlT << "\n";
        std::cout << "Mlus : \n" << Mlus << "\n";
        std::cout << "Llpg : \n" << Klpg << "\n";
        std::cout << "Llpc : \n" << Klpc << "\n";
        std::cout << "fl   : \n" << Bl << "\n";
        std::cout << "---------------------\n";
        std::cout << "Kupg : \n" << Kupg << "\n";
        std::cout << "Kupc : \n" << Kupc << "\n";
        std::cout << "fu   : \n" << Bu << "\n";
        std::cout << "---------------------\n";
        std::cout << "Mepg : \n" << MeT << "\n";
        std::cout << "Mepc : \n" << MeT << "\n";
        std::cout << "MeT  : \n" << MeT << "\n";
        std::cout << "Aepg : \n" << Aepg << "\n";
        std::cout << "Aepc : \n" << Aepc << "\n";
        std::cout << "Ket=(AeT+LeT)  : \n" << KeT << "\n";
#endif




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
                              std::vector<double>& local_M_data,
                              std::vector<double>& local_K_data,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override
    {
        const auto matrix_size = gas_pressure_size + cap_pressure_size +
                                 temperature_size + displacement_size;

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
//      auto dre_dus =
//              J.template block<temperature_size, displacement_size>(
//                      temperature_index, displacement_index);
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
        auto const element_id = _element.getID();
        x_position.setElementID(element_id);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

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
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                displacement_size> const>(local_x.data() + displacement_index,
                                          displacement_size);

        auto gas_phase_pressure_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                gas_pressure_size> const>(local_xdot.data() + gas_pressure_index,
                                          gas_pressure_size);
        auto capillary_pressure_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                cap_pressure_size> const>(local_xdot.data() + cap_pressure_index,
                                          cap_pressure_size);
        auto temperature_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                temperature_size> const>(local_xdot.data() + temperature_index,
                                          temperature_size);
        auto displacement_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                displacement_size> const>(local_xdot.data() + displacement_index,
                                          displacement_size);

        double const& dt = _process_data.dt;

        auto const& solid_phase = medium->phase(0);
        auto const& liquid_phase = medium->phase(1);
        auto const& gas_phase = medium->phase(2);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto& ip_data = _ip_data[ip];

            // shape matrices etc.
            auto const& w = ip_data.integration_weight;
            auto const& Np =  ip_data.N_p;
            auto const& NpT = Np.transpose().eval();
            auto const& gradNp = ip_data.dNdx_p;
            auto const& gradNpT = gradNp.transpose().eval();
            auto const& Nu = ip_data.N_u_op;
            auto const& NuT = Nu.transpose().eval();

            auto const x_coord =
                    interpolateXCoordinate<ShapeFunctionDisplacement,
                    ShapeMatricesTypeDisplacement>(_element,
                            ip_data.N_u);

            auto const B = LinearBMatrix::computeBMatrix<
                    DisplacementDim, ShapeFunctionDisplacement::NPOINTS,
                    typename BMatricesType::BMatrixType>(ip_data.dNdx_u,
                            ip_data.N_u, x_coord,_is_axially_symmetric);

            static int const KelvinVectorSize =
                    MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value;

            auto const& identity2 =
                    MathLib::KelvinVector::Invariants<KelvinVectorSize>::identity2;

            // primary unknowns
            auto const p_GR = gas_phase_pressure.dot(Np);
            auto const p_cap = capillary_pressure.dot(Np);
            auto const T = temperature.dot(Np);
            auto const u = Nu * displacement;

            auto const p_GR_dot = gas_phase_pressure_dot.dot(Np);
            auto const p_cap_dot = capillary_pressure_dot.dot(Np);
            auto const T_dot = temperature_dot.dot(Np);

            typename ShapeMatricesTypeDisplacement::GlobalDimVectorType
            const grad_p_GR = gradNp * gas_phase_pressure;
            typename ShapeMatricesTypeDisplacement::GlobalDimVectorType
            const grad_p_Cap = gradNp * capillary_pressure;

            double div_u_dot = identity2.transpose()*B*displacement_dot;

            const double p_LR = p_GR - p_cap;

            MPL::VariableArray variables;
            variables[MPL::Variables::capillary_pressure] = p_cap;
            variables[MPL::Variables::temperature] = T;

            auto const& b = _process_data.specific_body_force;

            // material parameters
            auto const alpha_B = MPL::getScalar(
                    medium->property(MPL::PropertyEnum::biot_coefficient),
                    variables);

            const double beta_p_SR = getScalar(
                solid_phase.property(MPL::PropertyEnum::compressibility),
                variables);
            const double beta_p_GR = getScalar(
                gas_phase.property(MPL::PropertyEnum::compressibility),
                variables);
            const double beta_p_LR = getScalar(
                liquid_phase.property(MPL::PropertyEnum::compressibility),
                variables);
            const double beta_p_S = getScalar(
                    medium->property(MPL::PropertyEnum::compressibility),
                    variables);

            const double beta_T_GR = getScalar(
                gas_phase.property(MPL::PropertyEnum::thermal_expansivity),
                variables);
            const double beta_T_LR = getScalar(
                liquid_phase.property(MPL::PropertyEnum::thermal_expansivity),
                variables);
            const double beta_T_SR = getScalar(
                solid_phase.property(MPL::PropertyEnum::thermal_expansivity),
                variables);

            auto const mu_LR = MPL::getScalar(
                liquid_phase.property(MPL::PropertyEnum::viscosity),
                variables);
            auto const dmuLRdpGR = MPL::getScalarDerivative(
                liquid_phase.property(MPL::PropertyEnum::viscosity),
                variables, MPL::Variables::phase_pressure);
            auto const dmuLRdpCap = MPL::getScalarDerivative(
                liquid_phase.property(MPL::PropertyEnum::viscosity),
                variables, MPL::Variables::capillary_pressure);
            auto const dmuLRdT = MPL::getScalarDerivative(
                liquid_phase.property(MPL::PropertyEnum::viscosity),
                variables, MPL::Variables::temperature);

            auto const mu_GR = MPL::getScalar(
                    gas_phase.property(MPL::PropertyEnum::viscosity),
                               variables);
            auto const dmuGRdpGR = MPL::getScalarDerivative(
                gas_phase.property(MPL::PropertyEnum::viscosity),
                variables, MPL::Variables::phase_pressure);
            auto const dmuGRdT = MPL::getScalarDerivative(
                gas_phase.property(MPL::PropertyEnum::viscosity),
                variables, MPL::Variables::temperature);

            auto const phi = MPL::getScalar(
                medium->property(MPL::PropertyEnum::porosity), variables);

            auto const sL = MPL::getScalar(medium->property(
                    MPL::PropertyEnum::saturation), variables);
            auto const sG = 1. - sL;
            auto const dsLdpc = MPL::getScalarDerivative(medium->property(
                    MPL::PropertyEnum::saturation), variables,
                    MPL::Variables::capillary_pressure);
            auto const d2sLdpc2 = MPL::getScalarDerivative(medium->property(
                    MPL::PropertyEnum::saturation), variables,
                    MPL::Variables::capillary_pressure,
                    MPL::Variables::capillary_pressure);

            auto const phi_S = 1. - phi;
            auto const phi_L = sL * phi;
            auto const phi_G = sG * phi;

            auto const rho_SR = MPL::getScalar(solid_phase.property(
                    MPL::PropertyEnum::density), variables);

            variables[MPL::Variables::phase_pressure] = p_GR;
            auto const rho_GR = MPL::getScalar(gas_phase.property(
                    MPL::PropertyEnum::density), variables);

            auto const drhoGRdpGR = MPL::getScalarDerivative(
                    gas_phase.property(MPL::PropertyEnum::density), variables,
                    MPL::Variables::phase_pressure);

            auto const drhoGRdT = MPL::getScalarDerivative(
                    gas_phase.property(MPL::PropertyEnum::density), variables,
                    MPL::Variables::temperature);

            variables[MPL::Variables::phase_pressure] = p_LR;
            auto const rho_LR = MPL::getScalar(liquid_phase.property(
                    MPL::PropertyEnum::density), variables);

            auto const drhoLRdpLR = MPL::getScalarDerivative(
                    liquid_phase.property(MPL::PropertyEnum::density),
                    variables, MPL::Variables::phase_pressure);

            auto const drhoLRdT = MPL::getScalarDerivative(
                    liquid_phase.property(MPL::PropertyEnum::density),
                    variables, MPL::Variables::temperature);

            auto const drhoLRdpGR = drhoLRdpLR;
            auto const drhoLRdpCap = - drhoLRdpLR;

            auto const drhoSRdpGR = 0.;
            auto const drhoSRdpCap = 0.;
            auto const drhoSRdT = 0.;

            double const permeability =
                MPL::getScalar(medium->property(MPL::permeability));

            GlobalDimMatrixType permeability_tensor =
                    GlobalDimMatrixType::Zero(DisplacementDim, DisplacementDim);

            permeability_tensor.diagonal().setConstant(permeability);

            auto const k_rel_LR = MPL::getPair(medium->property(
                    MPL::PropertyEnum::relative_permeability), variables)[0];
            auto const k_rel_GR = MPL::getPair(medium->property(
                    MPL::PropertyEnum::relative_permeability), variables)[1];

            auto const dkrelLdsL = MPL::getPairDerivative(
                                medium->property(MPL::PropertyEnum::relative_permeability),
                                variables, MPL::Variables::liquid_saturation)[0];
            auto const dkrelGdsL = MPL::getPairDerivative(
                                medium->property(MPL::PropertyEnum::relative_permeability),
                                variables, MPL::Variables::liquid_saturation)[1];

            const double cp_G = getScalar(gas_phase.property(
                    MPL::PropertyEnum::specific_heat_capacity), variables);
            const double dcpGdpGR = getScalarDerivative(gas_phase.property(
                    MPL::PropertyEnum::specific_heat_capacity),
                    variables, MPL::Variables::phase_pressure);
            const double dcpGdT = getScalarDerivative(gas_phase.property(
                    MPL::PropertyEnum::specific_heat_capacity),
                    variables, MPL::Variables::temperature);

            const double cp_L = getScalar(liquid_phase.property(
                    MPL::PropertyEnum::specific_heat_capacity), variables);
            const double dcpLdpGR = getScalarDerivative(liquid_phase.property(
                    MPL::PropertyEnum::specific_heat_capacity),
                    variables, MPL::Variables::phase_pressure);
            const double dcpLdpCap = getScalarDerivative(liquid_phase.property(
                    MPL::PropertyEnum::specific_heat_capacity),
                    variables, MPL::Variables::capillary_pressure);
            const double dcpLdT = getScalarDerivative(liquid_phase.property(
                    MPL::PropertyEnum::specific_heat_capacity),
                    variables, MPL::Variables::temperature);

            const double cp_S = getScalar(solid_phase.property(
                    MPL::PropertyEnum::specific_heat_capacity), variables);
            const double dcpSdpGR = getScalarDerivative(solid_phase.property(
                    MPL::PropertyEnum::specific_heat_capacity),
                    variables, MPL::Variables::phase_pressure);
            const double dcpSdpCap = getScalarDerivative(solid_phase.property(
                    MPL::PropertyEnum::specific_heat_capacity),
                    variables, MPL::Variables::capillary_pressure);
            const double dcpSdT = getScalarDerivative(solid_phase.property(
                    MPL::PropertyEnum::specific_heat_capacity),
                    variables, MPL::Variables::temperature);

            const double lambda_G = getScalar(gas_phase.property(
                    MPL::PropertyEnum::thermal_conductivity), variables);
            const double dlambdaGRdpGR = getScalarDerivative(
                    gas_phase.property(
                            MPL::PropertyEnum::thermal_conductivity),
                            variables, MPL::Variables::phase_pressure);
            const double dlambdaGRdT = getScalarDerivative(
                    gas_phase.property(
                            MPL::PropertyEnum::thermal_conductivity),
                            variables, MPL::Variables::temperature);

            const double lambda_L = getScalar(liquid_phase.property(
                    MPL::PropertyEnum::thermal_conductivity), variables);
            const double dlambdaLRdpGR = getScalarDerivative(
                    liquid_phase.property(
                            MPL::PropertyEnum::thermal_conductivity),
                            variables, MPL::Variables::phase_pressure);
            const double dlambdaLRdpCap = getScalarDerivative(
                    liquid_phase.property(
                            MPL::PropertyEnum::thermal_conductivity),
                            variables, MPL::Variables::capillary_pressure);
            const double dlambdaLRdT = getScalarDerivative(
                    liquid_phase.property(
                            MPL::PropertyEnum::thermal_conductivity),
                            variables, MPL::Variables::temperature);

            const double lambda_S = getScalar(solid_phase.property(
                    MPL::PropertyEnum::thermal_conductivity), variables);
            const double dlambdaSRdpGR = getScalarDerivative(
                    solid_phase.property(
                            MPL::PropertyEnum::thermal_conductivity),
                            variables, MPL::Variables::phase_pressure);
            const double dlambdaSRdpCap = getScalarDerivative(
                    solid_phase.property(
                            MPL::PropertyEnum::thermal_conductivity),
                            variables, MPL::Variables::capillary_pressure);
            const double dlambdaSRdT = getScalarDerivative(
                    solid_phase.property(
                            MPL::PropertyEnum::thermal_conductivity),
                            variables, MPL::Variables::temperature);

            const double lambda_eff = phi_G * lambda_G +
                    phi_L * lambda_L +
                    phi_S * lambda_S;

            Eigen::Matrix<double, DisplacementDim, DisplacementDim>
                conductivity_tensor;
            conductivity_tensor.setIdentity();

            const auto conductivity_effective =
                lambda_eff * conductivity_tensor;

            const auto dlambdadpGR = (phi_G * dlambdaGRdpGR +
                    phi_L * dlambdaLRdpGR + phi_S * dlambdaSRdpGR) *
                            conductivity_tensor;

            const auto dlambdadpCap = (phi * dsLdpc * (lambda_L - lambda_G) +
                    phi_L * dlambdaLRdpCap + phi_S * dlambdaSRdpCap) *
                            conductivity_tensor;

            const auto dlambdadT = (phi_G * dlambdaGRdT +
                                phi_L * dlambdaLRdT + phi_S * dlambdaSRdT) *
                                        conductivity_tensor;


            auto& eps = ip_data.eps;
            auto const& sigma_eff = ip_data.sigma_eff;

            eps.noalias() = B * displacement;
            auto C = ip_data.updateConstitutiveRelation(
                    t, x_position, dt, displacement,T, p_GR);


            auto const rho = phi_G * rho_GR + phi_L * rho_LR + phi_S * rho_SR;

            auto const drhodpGR = phi_G * drhoGRdpGR + phi_L * drhoLRdpGR +
                                phi_S * drhoSRdpGR;

            auto const drhodpCap = phi * dsLdpc * (rho_LR - rho_GR) + phi_L * drhoLRdpCap +
                                phi_S * drhoSRdpCap;
            auto const drhodT = phi_G * drhoGRdT + phi_L * drhoLRdT +
                                            phi_S * drhoSRdT;

            auto const rhocp = phi_G * rho_GR * cp_G + phi_L * rho_LR * cp_L +
                    phi_S * rho_SR * cp_S;

            auto const drhocpdpGR = phi_G *
                    (drhoGRdpGR * cp_G + rho_GR * dcpGdpGR) +
                    phi_L * (drhoLRdpGR * cp_L + rho_LR * dcpLdpGR) +
                    phi_S * (drhoSRdpGR * cp_S + rho_SR * dcpSdpGR);
            auto const drhocpdpCap = phi * dsLdpc *
                    (rho_LR * cp_L - rho_GR * cp_G) +
                    phi_L * (drhoLRdpCap * cp_L + rho_LR * dcpLdpCap) +
                    phi_S * (drhoSRdpCap * cp_S + rho_SR * dcpSdpCap);
            auto const drhocpdT = phi_G *
                    (drhoGRdT * cp_G + rho_GR * dcpGdT) +
                    phi_L * (drhoLRdT * cp_L + rho_LR * dcpLdT) +
                    phi_S * (drhoSRdT * cp_S + rho_SR * dcpSdT);

            auto const beta_T_R = phi_G * beta_T_GR + phi_L * beta_T_LR +
                    phi_S * beta_T_SR;

            auto const dbetaTRdCap = phi * dsLdpc * (beta_T_LR - beta_T_GR);

            auto const k_over_mu_GR = k_rel_GR / mu_GR;
            auto const k_over_mu_LR = k_rel_LR / mu_LR;

            // darcy-velocities
            const auto w_LS =
                    (-permeability_tensor * k_over_mu_LR *
                    ( gradNp * (gas_phase_pressure - capillary_pressure) - rho_LR * b)).eval();

            const auto w_GS =
                    (-permeability_tensor * k_over_mu_GR *
                    ( gradNp * gas_phase_pressure - rho_GR * b)).eval();

            const auto Sps = (alpha_B - phi) * beta_p_SR;
            const auto STs = (alpha_B - phi) * beta_T_SR;

            ip_data.rel_perm_gas = k_rel_GR;
            ip_data.rel_perm_liquid = k_rel_LR;

            /*
             *  Residuum and its derivatives
             */

            // Gas phase equation, gas pressure part
            const double cG1 = phi*beta_p_GR + Sps;
            const double cG2 = phi*beta_T_GR + STs;
            const double cG3 = phi*dsLdpc + sG*(sL + p_cap * dsLdpc) * Sps;

            const double dcG3dpc = phi * d2sLdpc2 + Sps *
                    (dsLdpc * (2 - 3 * sL - p_cap * dsLdpc) +
                            sG*p_cap*d2sLdpc2);

            rg += NpT * sG * rho_GR * cG1 * Np * w * gas_phase_pressure_dot; // G1
            rg -= NpT * sG * rho_GR * cG2 * Np * w * temperature_dot; // G2
            rg -= NpT * rho_GR * cG3 * Np * w * capillary_pressure_dot; // G3
            rg += NpT * sG * rho_GR * alpha_B *
                    identity2.transpose() * B * w * displacement_dot; // G4
            rg += gradNpT * rho_GR * k_over_mu_GR * permeability_tensor *
                    gradNp * w * gas_phase_pressure; // G5
            rg -= gradNpT * rho_GR * k_over_mu_GR * permeability_tensor *
                    rho_GR * b * w; // G6

            // residuum implemented as -r (?)
            rg *= -1.0;

            // Gas phase equation, gas pressure derivatives
            drg_dpg += NpT * sG * cG1 * (drhoGRdpGR * p_GR_dot + rho_GR/dt) *
                    Np * w; // G1
            drg_dpg -= NpT * sG * cG2 * drhoGRdpGR * T_dot * Np * w; // G2
            drg_dpg -= NpT * cG3 * drhoGRdpGR * p_cap_dot * Np * w; // G3
            drg_dpg += NpT * sG * drhoGRdpGR * alpha_B * div_u_dot *
                    Np * w; // G4
            drg_dpg += gradNpT * k_over_mu_GR * permeability_tensor *
                    (drhoGRdpGR - rho_GR/mu_GR * dmuGRdpGR) * grad_p_GR *
                    Np * w; // G5(1)
            drg_dpg += gradNpT * rho_GR * k_over_mu_GR * permeability_tensor *
                    gradNp * w; // G5(2)
            drg_dpg += gradNpT * rho_GR * k_over_mu_GR * permeability_tensor *
                    (rho_GR/mu_GR*dmuGRdpGR - 2 * drhoGRdpGR) * b *
                    Np * w; // G6

            // Gas phase equation, capillary pressure derivatives
            drg_dpc -= NpT * dsLdpc * rho_GR * cG1 * p_GR_dot * Np * w; // G1
            drg_dpc += NpT * dsLdpc * rho_GR * cG2 * T_dot * Np * w; // G2
            drg_dpc -= NpT * rho_GR * (cG3/dt + dcG3dpc * p_cap_dot) *
                    Np * w; // G3
            drg_dpc -= NpT * dsLdpc * rho_GR * alpha_B * div_u_dot *
                    Np * w; // G4
            drg_dpc += gradNpT * rho_GR * k_over_mu_GR * permeability_tensor *
                    grad_p_GR * Np * w; // G5
            drg_dpc -= gradNpT * rho_GR * rho_GR * dkrelGdsL * dsLdpc / mu_GR *
                    permeability_tensor * b * Np * w; // G6

            // Gas phase equation, displacement drivatives
            drg_dus +=  NpT * sG * rho_GR / dt * alpha_B *
                    identity2.transpose() * B * w; // G4

            // Gas phase equation, temperature drivatives
            drg_dT += NpT * sG * drhoGRdT * cG1 * p_GR_dot * Np * w; // G1
            drg_dT -= NpT * sG * cG2 * (drhoGRdT * T_dot + rho_GR/dt) *
                    Np * w; // G2
            drg_dT -= NpT * drhoGRdT * cG3 * p_cap_dot * Np * w; // G3
            drg_dT += NpT * sG * drhoGRdT * alpha_B * div_u_dot * Np * w; // G4
            drg_dT += gradNpT * k_over_mu_GR * permeability_tensor *
                    (drhoGRdT - rho_GR/mu_GR * dmuGRdT) * grad_p_GR *
                    Np * w; // G5
            drg_dT -= gradNpT * rho_GR * k_over_mu_GR * permeability_tensor *
                    (2*drhoGRdT - rho_GR/mu_GR * dmuGRdT) * b *
                    Np * w; // G6


            // liquid phase equation, gas pressure part
            const double cL1 = phi * beta_p_LR + Sps;
            const double cL2 = phi * beta_T_LR + STs;
            const double cL3 = phi * dsLdpc - phi_L*beta_p_LR - sL *
                    (sL + p_cap * dsLdpc) * Sps;

            rl += NpT * sL * rho_LR * cL1 * Np * w * gas_phase_pressure_dot; // L1
            rl -= NpT * sL * rho_LR * cL2 * Np * w * temperature_dot; // L2
            rl += NpT * rho_LR * cL3 * Np * w * capillary_pressure_dot; // L3
            rl += NpT * sL * rho_LR * alpha_B *
                    identity2.transpose() * B * w * displacement_dot; // L4
            rl += gradNpT * rho_LR * k_over_mu_LR * permeability_tensor *
                    gradNp * w * gas_phase_pressure; // L5
            rl -= gradNpT * rho_LR * k_over_mu_LR * permeability_tensor *
                    gradNp * w * capillary_pressure; // L6
            rl -= gradNpT * rho_LR * k_over_mu_LR * permeability_tensor *
                    rho_LR * b * w;

            drl_dpg += NpT * sL * cL1 * (rho_LR / dt + drhoLRdpGR * p_GR_dot) *
                    Np * w; // L1
            drl_dpg -= NpT * sL * drhoLRdpGR * cL2 * T_dot * Np * w; // L2
            drl_dpg += NpT * drhoLRdpGR * cL3 * p_cap_dot * Np * w; // L3
            drl_dpg += NpT * sL * drhoLRdpGR * alpha_B * div_u_dot *
                    Np * w; // L4
            drl_dpg += gradNpT * k_over_mu_LR * permeability_tensor  *
                    (drhoLRdpGR - rho_LR/mu_LR - dmuLRdpGR) * grad_p_GR *
                    Np * w; // L5a
            drl_dpg += gradNpT * rho_LR * k_over_mu_LR * permeability_tensor  *
                    gradNp * w; // L5b

            drl_dpg -= gradNpT * k_over_mu_LR * permeability_tensor *
                    (drhoLRdpGR - rho_LR/mu_LR * dmuLRdpGR) * grad_p_Cap *
                    Np * w; // L6
            drl_dpg -= gradNpT * rho_LR * k_over_mu_LR * permeability_tensor *
                    (2 * drhoLRdpGR - rho_LR/mu_LR * dmuLRdpGR) * b *
                    Np * w; // L7


            // liquid phase equation, capillary pressure part

            const double dcL3dpc = phi * (d2sLdpc2 - dsLdpc* beta_p_LR) -
                    Sps * (dsLdpc * (3*sL+p_cap*dsLdpc) + sL*p_cap*d2sLdpc2);

            drl_dpc += NpT * (dsLdpc * rho_LR + sL * drhoLRdpCap) * cL1 *
                    p_GR_dot * Np * w; // L1

            drl_dpc -= NpT * (dsLdpc * rho_LR + sL * drhoLRdpCap) * cL2 *
                    T_dot * Np * w; // L2

            drl_dpc += NpT * (( drhoLRdpCap * cL3 + rho_LR * dcL3dpc) *
                    p_cap_dot + rho_LR * cL3 / dt) * Np * w; // L3
            drl_dpc += NpT * sL * drhoLRdpCap * alpha_B * div_u_dot *
                    Np * w; // L4
            drl_dpc += gradNpT * permeability_tensor / mu_LR *
                    (drhoLRdpCap * k_rel_LR + rho_LR * dkrelLdsL * dsLdpc -
                            rho_LR * k_rel_LR / mu_LR * dmuLRdpCap) *
                            grad_p_GR * Np * w; // L5
            drl_dpc -= gradNpT * permeability_tensor / mu_LR *
                    (drhoLRdpCap * k_rel_LR + rho_LR * dkrelLdsL * dsLdpc -
                            rho_LR * k_rel_LR / mu_LR * dmuLRdpCap) *
                            grad_p_Cap * Np * w; // L6a
            drl_dpc -= gradNpT * rho_LR * k_over_mu_LR * permeability_tensor *
                    gradNp * w; // L6b
            drl_dpc -= gradNpT * rho_LR * permeability_tensor / mu_LR *
                    (2 * drhoLRdpCap * k_rel_LR + rho_LR * dkrelLdsL * dsLdpc -
                            rho_LR*k_rel_LR / mu_LR * dmuLRdpCap) * b *
                            Np * w;

            // liquid phase equation, displacement derivatives
            drl_dus += NpT * sL * rho_LR / dt * alpha_B *
                    identity2.transpose() * B * w; // L4

            // liquid phase equation, temperature derivatives
            drl_dT += NpT * sL * drhoLRdT * cL1 * p_GR_dot * Np * w; // L1
            drl_dT -= NpT * sL * cL2 * (rho_LR / dt + drhoLRdT * T_dot) *
                    Np * w; // L2
            drl_dT += NpT * drhoLRdT * cL3 * p_cap_dot * Np * w; // L3
            drl_dT += NpT * sL * drhoLRdT * alpha_B * div_u_dot * Np * w; // L4
            drl_dT += gradNpT * k_over_mu_LR * permeability_tensor *
                    (drhoLRdT - rho_LR/mu_LR*dmuLRdT) * grad_p_GR *
                    Np * w; // L5
            drl_dT -= gradNpT * k_over_mu_LR * permeability_tensor *
                    (drhoLRdT - rho_LR/mu_LR*dmuLRdT) * grad_p_Cap *
                    Np * w; // L6
            drl_dT -= gradNpT * rho_LR * k_over_mu_LR * permeability_tensor *
                    (2 * drhoLRdT - rho_LR/mu_LR*dmuLRdT) * b * Np * w; // L7

            // displacement equation, gas pressure derivatives
            ru += B.transpose() * sigma_eff * w;
            ru -= NuT * rho * b * w;
            ru -= B.transpose() * alpha_B * identity2 * Np * w *
                    gas_phase_pressure_dot; // U1
            ru += B.transpose() * alpha_B * identity2 * sL * Np * w *
                    capillary_pressure_dot; // U2

            dru_dpg -= B.transpose() * alpha_B * identity2 * Np * w; // U2
            dru_dpg -= NuT * drhodpGR * b * Np * w; // U4

            // displacement equation, gas pressure derivatives
            dru_dpc += B.transpose() * alpha_B * (dsLdpc*p_cap + sL) *
                    identity2 * Np * w; // U3
            dru_dpc -= NuT * drhodpCap * b * Np * w; // U4

            // displacement equation, displacement derivatives
            dru_dus += B.transpose() * C * B * w; // U1

            // displacement equation, temperature derivatives
//            dru_dT += B.transpose() / beta_p_S * beta_T_SR *identity2 *
//                    Np * w; // U1
            dru_dT -= NuT * drhodT * b * Np * w; // U4

            // energy equation, gas pressure derivatives
            const double cE2 = phi_L*beta_T_LR -
                    phi_S*beta_T_SR * (sL + p_cap*dsLdpc);
            const double dcE2dpc = phi * dsLdpc * beta_T_LR -
                    phi_S * beta_T_SR * (2*dsLdpc + p_cap * d2sLdpc2);

            re += NpT * rhocp * Np * w * temperature_dot; // E1
            re += NpT * cE2 * T * Np * w * capillary_pressure_dot; // E2
            re -= NpT * beta_T_R * T * Np * w * gas_phase_pressure_dot; // E3
            re += NpT * (rho_GR * cp_G * w_GS.transpose() +
                    rho_LR * cp_L * w_LS.transpose()) * gradNp *
                            w * temperature; // E4
            re -= NpT * (beta_T_GR * w_GS.transpose() +
                    beta_T_LR * w_LS.transpose()) * T * gradNpT *
                            w * gas_phase_pressure; // E5
            re += NpT * beta_T_LR * T * w_LS.transpose() * gradNpT *
                    w * capillary_pressure; // E6
            re += gradNpT * lambda_eff * gradNp * w * temperature; // E7

            dre_dpg += NpT * drhocpdpGR * T_dot * Np * w; // E1
            dre_dpg -= NpT * beta_T_R / dt * T * Np * w; // E3

            // energy equation, capillary pressure derivatives
            dre_dpc += NpT * drhocpdpCap * T_dot * Np * w; // E1
            dre_dpc += NpT * (cE2/dt + dcE2dpc * p_cap_dot) * T * Np * w; // E2
            dre_dpc -= NpT * dbetaTRdCap * T * p_GR_dot * Np * w; // E3

            // energy equation, temperature derivatives
            dre_dT += NpT * (rhocp/dt + drhocpdT * T_dot) * Np * w; // E1
            dre_dT += NpT * cE2 * p_cap_dot * Np * w; // E2
            dre_dT -= NpT * beta_T_R * p_GR_dot * Np * w; // E3

        }


#if 0

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

        // create M matrix
        auto M = MathLib::createZeroedMatrix<
            typename ShapeMatricesTypeDisplacement::template MatrixType<
                matrix_size, matrix_size>>(
            local_M_data, matrix_size, matrix_size);

        // create M matrix
        auto K = MathLib::createZeroedMatrix<
            typename ShapeMatricesTypeDisplacement::template MatrixType<
                matrix_size, matrix_size>>(
            local_K_data, matrix_size, matrix_size);

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


        // right hand side vectors
        PVector fg = PVector::Zero(Np_intPoints);
        PVector fl = PVector::Zero(Np_intPoints);
        PVector fe = PVector::Zero(Np_intPoints);
        UVector fu = UVector::Zero(Nu_intPoints);

        typename ShapeMatricesTypePressure::NodalMatrixType Laplace =
        		ShapeMatricesTypePressure::NodalMatrixType::Zero(Np_intPoints,
        				Np_intPoints);

        // M blocks:
        auto Mgpg =
                M.template block<gas_pressure_size, gas_pressure_size>(
                        gas_pressure_index, gas_pressure_index);
        auto Mgpc =
                M.template block<gas_pressure_size, cap_pressure_size>(
                        gas_pressure_index, cap_pressure_index);
        auto MgT =
                M.template block<gas_pressure_size, temperature_size>(
                        gas_pressure_index, temperature_index);
        auto Mgus =
                M.template block<gas_pressure_size, displacement_size>(
                        gas_pressure_index, displacement_index);

        auto Mlpg =
                M.template block<cap_pressure_size, gas_pressure_size>(
                        cap_pressure_index, gas_pressure_index);
        auto Mlpc =
                M.template block<cap_pressure_size, cap_pressure_size>(
                        cap_pressure_index, cap_pressure_index);
        auto MlT =
                M.template block<cap_pressure_size, temperature_size>(
                        cap_pressure_index, temperature_index);
        auto Mlus =
                M.template block<cap_pressure_size, displacement_size>(
                        cap_pressure_index, displacement_index);

        auto Mepg =
                M.template block<temperature_size, gas_pressure_size>(
                        temperature_index, gas_pressure_index);
        auto Mepc =
                M.template block<temperature_size, cap_pressure_size>(
                        temperature_index, cap_pressure_index);
        auto MeT =
                M.template block<temperature_size, temperature_size>(
                        temperature_index, temperature_index);

        // M blocks:
        auto Kgpg =
                K.template block<gas_pressure_size, gas_pressure_size>(
                        gas_pressure_index, gas_pressure_index);
        auto Klpg =
                K.template block<cap_pressure_size, gas_pressure_size>(
                        cap_pressure_index, gas_pressure_index);
        auto Klpc =
                K.template block<cap_pressure_size, cap_pressure_size>(
                        cap_pressure_index, cap_pressure_index);
        auto Kepg =
                K.template block<temperature_size, gas_pressure_size>(
                        temperature_index, gas_pressure_index);
        auto Kepc =
                K.template block<temperature_size, cap_pressure_size>(
                        temperature_index, cap_pressure_index);
        auto KeT =
                K.template block<temperature_size, temperature_size>(
                        temperature_index, temperature_index);
        auto Kupg =
                K.template block<displacement_size, gas_pressure_size>(
                        displacement_index, gas_pressure_index);
        auto Kupc =
                K.template block<displacement_size, cap_pressure_size>(
                        displacement_index, cap_pressure_index);


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
            auto const& NpT = Np.transpose().eval();
            auto const& gradNp = ip_data.dNdx_p;
            auto const& gradNpT = gradNp.transpose().eval();
            auto const& Nu = ip_data.N_u_op;
            auto const& NuT = Nu.transpose().eval();

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

            typename ShapeMatricesTypeDisplacement::GlobalDimVectorType
            const grad_p_GR = dNdx_p * gas_phase_pressure;
            typename ShapeMatricesTypeDisplacement::GlobalDimVectorType
            const grad_p_Cap = dNdx_p * capillary_pressure;

            double div_u_dot = identity2.transpose()*B*displacement_dot;

            const double p_LR = p_GR - p_cap;

            eps.noalias() = B * displacement;

//            auto const e = // volume strain
//                MathLib::KelvinVector::Invariants<KelvinVectorSize>::trace(eps);

#define nOUTPUT_IP
#define nDBG_OUTPUT


#ifdef OUTPUT_IP
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
            std::cout << "================================\n";
            std::cout << " mass_op :\n" <<  N_p.transpose() * N_p * w << "\n\n";
            std::cout << " mass    :\n" <<  mass_operator << "\n";

#endif

            // insert all primary variables into one object for GPML
            MPL::VariableArray variables;
            variables[MPL::Variables::p_cap] = p_cap;
            variables[MPL::Variables::p_GR] = p_GR;
            variables[MPL::Variables::p_LR] = p_LR;
            variables[MPL::Variables::T] = T;


            // get medium properties
            auto& medium = _process_data.medium;

            // get phase properties
            auto const& solid_phase = medium->phase(0);
            auto const& liquid_phase = medium->phase(1);
            auto const& gas_phase = medium->phase(2);

            double const permeability =
                MPL::getScalar(medium->property(MPL::permeability));

            GlobalDimMatrixType permeability_tensor =
                GlobalDimMatrixType::Zero(DisplacementDim, DisplacementDim);

            permeability_tensor.diagonal().setConstant(permeability);

            auto const alpha_B = MPL::getScalar(
                    medium->property(MPL::PropertyEnum::biot_coefficient),
                    variables);

            auto const phi = MPL::getScalar(
                medium->property(MPL::PropertyEnum::porosity), variables);
            auto const phi_S = 1. - phi;

            // Saturation
//          auto const s_L =
//              MPL::getScalar(medium->property(MPL::PropertyEnum::saturation),
//                               variables);
            auto const s_L_r = 0.2;
            auto const s_a = -1.9722e-11;
            auto const s_b = 2.4279;

//          auto const s_L = std::max(s_L_r,1. + s_a * std::pow (std::max(0.,p_cap), s_b));
            auto const s_L = p_cap*p_cap + 2*p_cap + 1;
            auto const dsLdpc = 2*p_cap + 2;
            auto const d2sLdpc2 = 2.0;

            auto const s_G = 1 - s_L;
            auto const s_e = (s_L - s_L_r) / (1.0 - s_L_r);

            const double phi_L = s_L * phi;
            const double phi_G = s_G * phi;

//          auto const dsLdpc = MPL::getScalarDerivative(
//              medium->property(MPL::PropertyEnum::saturation),
//              variables, MPL::Variables::p_cap);

//          auto const dsLdpc = s_a*s_b*std::pow(std::max(0.,p_cap), s_b - 1.0);
//            auto const d2sLdpc2 = s_a*(s_b*(s_b - 1)) * std::pow(std::max(0.,p_cap), s_b - 2.0);


//          auto const k_rel_LR = MPL::getPair(
//              mediumte.property(MPL::PropertyEnum::relative_permeability),
//              variables)[0];
//          auto const k_rel_GR = MPL::getPair(
//              medium->property(MPL::PropertyEnum::relative_permeability),
//              variables)[1];

//            auto const k_rel_LR = 1.0 - 2.207*std::pow((1.0 - s_L), 1.0121);
            auto const k_rel_LR = 1.0;
            auto const min_k_rel_GR = 0.0001;

            auto const k_rel_b = 5./3.;
//            auto const k_rel_GR = (1.0 - s_e) * (1 - s_e)
//                                                * (1.0 - std::pow(s_e, k_rel_b)) + min_k_rel_GR;


            auto const k_rel_GR = 2*s_L*s_L + 3*s_L + 3;

//            auto const dkrelGdsE = 1./s_e * (1 - s_e) *
//            		((k_rel_b+2)*std::pow(s_e,k_rel_b+1)-k_rel_b*std::pow(s_e, k_rel_b) - 2*s_e);
//
//            auto const dsEdsL = 1. / (1. - s_L_r);
//
//            auto const dkrelGdsL = dkrelGdsE * dsEdsL;
            auto const dkrelGdsL = 4 * s_L + 3;
            auto const dkrelLdsL = 0.0;

            auto const rho_SR =
                    MPL::getScalar(solid_phase.property(MPL::PropertyEnum::density),
                            variables);
            auto const rho_LR = MPL::getScalar(
                    liquid_phase.property(MPL::PropertyEnum::density),
                    variables);
            auto const rho_GR_0 =
                    MPL::getScalar(gas_phase.property(MPL::PropertyEnum::density),
                            variables);

            const double beta_p_SR = getScalar(
                solid_phase.property(MPL::PropertyEnum::compressibility));
            const double beta_p_PR =
                getScalar(medium->property(MPL::PropertyEnum::compressibility));
            const double beta_p_GR = getScalar(
                gas_phase.property(MPL::PropertyEnum::compressibility));
            const double beta_p_LR = getScalar(
                liquid_phase.property(MPL::PropertyEnum::compressibility));

            const auto p_GR_0 = 100000;
            auto rho_GR = rho_GR_0 * (1 + beta_p_GR * (p_GR - p_GR_0));

            auto const& b = _process_data.specific_body_force;

            variables[MPL::Variables::liquid_density] = rho_LR;
            variables[MPL::Variables::gas_density] = rho_GR;

//            auto const drhoGRdpGR = MPL::getScalarDerivative(
//                    gas_phase.property(MPL::density), variables, MPL::p_GR);
//            auto const drhoGRdpGR = rho_GR_0*beta_p_GR;
            auto const drhoGRdpGR = rho_GR_0 * beta_p_GR;

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

            auto const dmuLRdpGR = 0.0;
            auto const dmuLRdpCap = 0.0;
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

            const double dcpGdpGR = 0.;
            const double dcpLdpGR = 0.;
            const double dcpSdpGR = 0.;
            const double dcpLdpCap = 0.;
            const double dcpSdpCap = 0.;
            const double dcpGdT = 0.;
            const double dcpLdT = 0.;
            const double dcpSdT = 0.;


            // expansivities
            const double beta_T_GR = getScalar(
                gas_phase.property(MPL::PropertyEnum::thermal_expansivity));
            const double beta_T_LR = getScalar(
                liquid_phase.property(MPL::PropertyEnum::thermal_expansivity));
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

            const double dLambdaGRdpGR = 0.;
            const double dLambdaLRdpGR = 0.;
            const double dLambdaLRdpCap = 0.;
            const double dLambdaSRdpGR = 0.;
            const double dLambdaSRdpCap = 0.;
            const double dLambdaGRdT = 0.;
            const double dLambdaLRdT = 0.;
            const double dLambdaSRdT = 0.;

            const double lambda_eff = phi_G * lambda_G +
                                      phi_L * lambda_L +
                                      phi_S * lambda_S;

            // isotropic heat conductivity
            const auto conductivity_effective =
                lambda_eff * conductivity_tensor;

            const auto dLambdadpGR = (phi_G * dLambdaGRdpGR +
                    phi_L * dLambdaLRdpGR + phi_S * dLambdaSRdpGR) *
                            conductivity_tensor;

            const auto dLambdadpCap = (phi * dsLdpc * (lambda_L - lambda_G) +
                    phi_L * dLambdaLRdpCap + phi_S * dLambdaSRdpCap) *
                            conductivity_tensor;

            const auto dLambdadT = (phi_G * dLambdaGRdT +
                                phi_L * dLambdaLRdT + phi_S * dLambdaSRdT) *
                                        conductivity_tensor;

            /*
             * derived variables
             */
            const double rho_cp_eff = phi_G * rho_GR * cp_G +
                    phi_L * rho_LR * cp_L + phi_S * rho_SR * cp_S;

            const double drhocpGdpGR = phi_G * (drhoGRdpGR*cp_G + rho_GR*dcpGdpGR) +
                    phi_L * (drhoLRdpGR * cp_L + rho_LR * dcpLdpGR) +
                    phi_S * (drhoSRdpGR * cp_S + rho_SR * dcpSdpGR);

            const double drhocpGdpCap = phi * dsLdpc * (rho_LR*cp_L + rho_GR*cp_G) +
                    phi_L * (drhoLRdpCap * cp_L + rho_LR * dcpLdpCap) +
                    phi_S * (drhoSRdpCap * cp_S + rho_SR * dcpSdpCap);

            const double drhocpGdT = phi_G * (drhoGRdT*cp_G + rho_GR*dcpGdT) +
                    phi_L * (drhoLRdT * cp_L + rho_LR * dcpLdT) +
                    phi_S * (drhoSRdT * cp_S + rho_SR * dcpSdT);


            const double k_over_mu_GR = k_rel_GR / mu_GR;
            const double k_over_mu_LR = k_rel_LR / mu_LR;

            const auto Sps = (alpha_B - phi) * beta_p_SR;
            const auto STs = (alpha_B - phi) * beta_T_SR;

            // darcy-velocities
            const auto w_LS =
                    (-permeability_tensor * k_over_mu_LR *
                    ( dNdx_p * (gas_phase_pressure - capillary_pressure) - rho_LR * b)).eval();

            const auto w_GS =
                    (-permeability_tensor * k_over_mu_GR *
                    ( dNdx_p * gas_phase_pressure - rho_GR * b)).eval();

            // Elasticity matrix:
            auto C = ip_data.updateConstitutiveRelation(t, x_position, dt, displacement,
                    T, p_GR);

            // auxiliary operators:
            Laplace.noalias() =
                dNdx_p.transpose() * permeability_tensor * dNdx_p * w;

            auto const gravity_operator =
                (dNdx_p.transpose() * permeability_tensor * b * w).eval();

            // mass_operator = N_p_T * N_p * w;
            // output secondary variables
            ip_data.pressure_gas_linear = p_GR;
            ip_data.pressure_cap_linear = p_cap;
            ip_data.pressure_wet = p_LR;
            ip_data.rel_perm_gas = k_rel_GR;
            ip_data.rel_perm_liquid = k_rel_LR;
            ip_data.density_gas = rho_GR;
            ip_data.density_liquid = rho_LR;
            ip_data.saturation = s_L;

            ip_data.velocity_liquid.noalias() = w_LS;
            ip_data.velocity_gas.noalias() = w_GS;


#ifdef OUTPUT_IP
            std::cout << "==================================\n";
            std::cout << "            rho_SR : " << rho_SR << " \n";
            std::cout << "            rho_LR : " << rho_LR << " \n";
            std::cout << "            rho_GR : " << rho_GR << " \n";
            std::cout << "               rho : " << rho << " \n";
            std::cout << "==================================\n";
            std::cout << "               phi : " << phi << " \n";
            std::cout << "             phi_G : " << phi_G << " \n";
            std::cout << "             phi_L : " << phi_L << " \n";
            std::cout << "             phi_S : " << phi_S << " \n";
            std::cout << "==================================\n";
            std::cout << "             p_GR  : " << p_GR << " \n";
            std::cout << "             p_cap : " << p_cap << " \n";
            std::cout << "----------------------------------\n";
            std::cout << "               s_L : " << s_L << " \n";
            std::cout << "               s_G : " << s_G << " \n";
            std::cout << "----------------------------------\n";
            std::cout << "             mu_LR : " << mu_LR << " \n";
            std::cout << "             mu_GR : " << mu_GR << " \n";
            std::cout << "==================================\n";
            std::cout << "    Gravity vector : \n";
            std::cout << "                 b : \n" << b << " \n";
            std::cout << "==================================\n";
            std::cout << "   volume strain e : " << /*e <<*/ " \n";
            std::cout << "==================================\n";
            std::cout << "                 C : " << "\n" << C << " \n\n";
            std::cout << "----------------------------------\n";
            std::cout << "         sigma_eff : " << "\n" << sigma_eff << "\n\n";
            std::cout << "----------------------------------\n";
            std::cout << "           alpha_B : " << alpha_B << " \n";
            std::cout << "----------------------------------\n";
            std::cout << "      permeability : " << permeability << " \n";
            std::cout << "----------------------------------\n";
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
            std::cout << "          lambda_G : " << lambda_G << " \n";
            std::cout << "          lambda_L : " << lambda_L << " \n";
            std::cout << "          lambda_S : " << lambda_S << " \n";
            std::cout << "        lamdba_eff : " << lambda_eff << " \n";
            std::cout << "==================================\n";
            std::cout << "        Velocities : \n";
            std::cout << "   ----------------------------------\n";
            std::cout << "              w_GS : \n" << w_GS << " \n";
            std::cout << "              w_LS : \n" << w_LS << " \n";
            std::cout << "==================================\n";
            std::cout << "           Laplace : \n";
            std::cout << "                 L :\n " << Laplace << " \n";
            std::cout << "==================================\n";
            OGS_FATAL("AN");

#endif
/*
            Mgpg.noalias() += rho_GR * s_G * (phi * beta_p_GR + Sps) *
            		mass_operator;


#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   Mgpg:\n" << Mgpg << "\n";
            std::cout << "==================================\n";
#endif

            Mgpc.noalias() -=
            		rho_GR * (phi*dsLdpc + s_G * (s_L + p_cap*dsLdpc) * Sps) *
					mass_operator;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   Mgpc:\n" << Mgpc << "\n";
            std::cout << "==================================\n";
#endif



            MgT.noalias() -= N_p.transpose() * s_G * rho_GR *
            		(phi * beta_T_GR + beta_T_SR * (alpha_B - phi)) *
					N_p * w;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   MgT:\n" << MgT << "\n";
            std::cout << "==================================\n";
#endif

            Mgus.noalias() += s_G * rho_GR * alpha_B *
            		mass_operator2;


#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   Mgus:\n" << Mgus << "\n";
            std::cout << "==================================\n";
#endif

            Mlpg.noalias() += N_p.transpose() * rho_LR * s_L *
                              (phi * beta_p_LR + Sps) *
                              N_p * w;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   Mlpg:\n" << Mlpg << "\n";
            std::cout << "==================================\n";
#endif

            Mlpc.noalias() += rho_LR * (dsLdpc * (phi - Sps*p_cap) - s_L*(phi*beta_p_LR + Sps)) *
            		mass_operator;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   Mlpc:\n" << Mlpc << "\n";
            std::cout << "==================================\n";
#endif

            MlT.noalias() -= s_L * rho_LR *
                             (phi * beta_T_LR + (alpha_B - phi) * beta_T_SR ) *
                             mass_operator;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   MlT:\n" << MlT << "\n";
            std::cout << "==================================\n";
#endif

            Mlus.noalias() += rho_LR * s_L * alpha_B *
            		mass_operator2;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   Mlus:\n" << Mlus << "\n";
            std::cout << "==================================\n";
#endif

            Mepg.noalias() -= (phi_G*beta_T_GR + phi_L*beta_T_LR + phi_S*beta_T_SR) * T *
            		mass_operator;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   Mepg:\n" << Mepg << "\n";
            std::cout << "=================================\n";
#endif
            Mepc.noalias() += (phi_L*beta_T_LR - (s_L + p_cap*dsLdpc)*phi_S*beta_T_SR) * T *
            		mass_operator;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   Mepc:\n" << Mepc << "\n";
            std::cout << "=================================\n";
#endif

            MeT.noalias() += rho_cp_eff * mass_operator;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   MeT:\n" << MeT << "\n";
            std::cout << "=================================\n";
#endif

            Kgpg.noalias() += rho_GR * k_over_mu_GR * Laplace;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   Kgpg:\n" << Kgpg << "\n";
            std::cout << "==================================\n";
#endif
            Klpg.noalias() += rho_LR * k_over_mu_LR * Laplace;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   Klpg:\n" << Klpg << "\n";
            std::cout << "==================================\n";
#endif

            Klpc.noalias() = -Klpg;
//            dLlpc_dpg.noalias() = -dLlpg_dpg;
//            dLlpc_dpc.noalias() = -dLlpg_dpc;
//            dLlpc_dT.noalias() = -dLlpg_dT;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   Klpc:\n" << Klpc << "\n";
            std::cout << "==================================\n";
#endif

            Kepg.noalias() -= N_p.transpose() *
                             (phi_G*beta_T_GR*w_GS.transpose() +
                              phi_L*beta_T_LR*w_LS.transpose()) * T *
                             dNdx_p * w;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   Kepg:\n" << Kepg << "\n";
            std::cout << "=================================\n";
#endif

            Kepc.noalias() +=
            		N_p.transpose() * (phi_L*beta_T_LR*T*w_LS.transpose())
					* dNdx_p * w;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   KAepc:\n" << Kepc << "\n";
            std::cout << "=================================\n";
#endif

            KeT.noalias() += N_p.transpose() *
                            (s_G*rho_GR*cp_G*w_GS.transpose() +
                             s_L*rho_LR*cp_L*w_LS.transpose()) *
                            dNdx_p * w;

            KeT.noalias() +=
                dNdx_p.transpose() * conductivity_effective * dNdx_p * w;


#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   KeT:\n" << KeT << "\n";
            std::cout << "=================================\n";
#endif

            Kupg.noalias() -= B.transpose() * alpha_B * identity2 * N_p * w;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   Kupg:\n" << Kupg << "\n";
            std::cout << "==================================\n";
#endif
            Kupc.noalias() +=
                B.transpose() * alpha_B * identity2 * s_L * N_p * w;

//            dKupc_dpc.noalias() -= (B.transpose()*alpha_B*identity2*dsLdpc*Np)*w;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   Kupc:\n" << Kupc << "\n";
            std::cout << "==================================\n";
#endif

            fg.noalias() += rho_GR * rho_GR * k_over_mu_GR * gravity_operator;

//            dfg_dpg += rho_GR * k_over_mu_GR * (2* drhoGRdpGR - rho_GR / mu_GR * dmuGRdpGR) *
//            		gravity_operator;
//
//            dfg_dpc += rho_GR * rho_GR / mu_GR * dkrelGdsL * dsLdpc *
//            		gravity_operator;
//
//            dfg_dT += rho_GR * k_over_mu_GR * (2*drhoGRdT + rho_GR/mu_GR*dmuGRdT) *
//            		gravity_operator;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   fg:\n" << fg << "\n";
            std::cout << "==================================\n";
#endif

            fl.noalias() += rho_LR * rho_LR * k_over_mu_LR * gravity_operator;

//            dfl_dpg += rho_LR * k_over_mu_LR * (2 * drhoLRdpLR - rho_LR/mu_LR*dmuLRdpLR) *
//            		gravity_operator;
//
//            dfl_dpc += rho_LR / mu_LR * (rho_LR * (dkrelLdsL*dsLdpc + dmuLRdpLR*k_rel_LR/mu_LR) - 2*drhoLRdpLR*k_rel_LR) *
//            		gravity_operator;
//
//            dfl_dT += rho_LR * k_over_mu_LR * (2*drhoLRdT - rho_LR/mu_LR * dmuLRdT) *
//            		gravity_operator;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   fl:\n" << fl << "\n";
            std::cout << "==================================\n";
#endif

            fu.noalias() += (B.transpose() * sigma_eff - N_u_op.transpose() * rho * b )* w;

//            dfu_dpg -= N_u_op.transpose() * drhodpGR * b * w;
//
//            dfu_dpc -= N_u_op.transpose() * drhodpCap * b * w;
//
//            dfu_dT -= (B.transpose() / beta_p_PR * beta_T_SR * identity2 +
//            		N_u_op.transpose() * drhodT * b ) * w;
//
//            dfu_du += B.transpose() * C * B * w;

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "   fu:\n" << fu << "\n";
            std::cout << "==================================\n";
#endif
*/
            const auto F = 2* p_GR * p_GR + 3* p_GR + 1;
            Mgpg.noalias() += 1 * mass_operator;
            Mgpc.noalias() += F * mass_operator;


            // Gas phase equation, gas pressure part
            const double cG1 = phi*beta_p_GR + Sps;
            const double cG2 = phi*beta_T_GR + STs;
            const double cG3 = phi*dsLdpc + s_G*(s_L + p_cap * dsLdpc) * Sps;

//            drg_dpg += NpT * s_G * cG1 * (drhoGRdpGR * p_GR_dot + rho_GR/dt) * Np * w; // G1
//            drg_dpg -= NpT * s_G * cG2 * drhoGRdpGR * T_dot * Np * w; // G2
//            drg_dpg -= NpT * cG3 * drhoGRdpGR * p_cap_dot * Np * w; // G3



//            drg_dpg += NpT * s_G * drhoGRdpGR * alpha_B * div_u_dot * Np  * w; // G4
//            drg_dpg += gradNpT * k_over_mu_GR * permeability_tensor *
//            		(drhoGRdpGR - rho_GR/mu_GR * dmuGRdpGR) * grad_p_GR * Np * w; // G5(1)
//            drg_dpg += gradNpT * rho_GR * k_over_mu_GR * permeability_tensor * gradNp * w; // G5(2)
//            drg_dpg += gradNpT * rho_GR * k_over_mu_GR * permeability_tensor *
//            		(rho_GR/mu_GR*dmuGRdpGR - 2 * drhoGRdpGR) * b * Np * w; // G6

//            std::cout << "G1:\n" <<
//                    NpT * s_G * cG1 * (drhoGRdpGR * p_GR_dot + rho_GR/dt) * Np * w
//                    << "\n";
//            std::cout << "G2:\n" <<
//                    NpT * s_G * cG2 * drhoGRdpGR * T_dot * Np * w
//                    << "\n";
//            std::cout << "G3:\n" <<
//                    NpT * cG3 * drhoGRdpGR * p_cap_dot * Np * w
//                    << "\n";
//            std::cout << "G4:\n" <<
//                    NpT * s_G * drhoGRdpGR * alpha_B * div_u_dot * Np  * w
//                    << "\n";
//            std::cout << "G5(1):\n" <<
//                    gradNpT * k_over_mu_GR * permeability_tensor *
//                                        (drhoGRdpGR - rho_GR/mu_GR * dmuGRdpGR) * grad_p_GR * Np * w
//                    << "\n";
//            std::cout << "G5(2):\n" <<
//                    gradNpT * rho_GR * k_over_mu_GR * permeability_tensor * gradNp * w
//                    << "\n";
//            std::cout << "G6:\n" <<
//                    gradNpT * rho_GR * k_over_mu_GR * permeability_tensor *
//                                        (rho_GR/mu_GR*dmuGRdpGR - 2 * drhoGRdpGR) * b * Np * w
//                    << "\n";


            const double dcG3dpc = phi * d2sLdpc2 + Sps *
            		(dsLdpc * (2 - 3 * s_L - p_cap * dsLdpc) + s_G*p_cap*d2sLdpc2);
//
//            // Gas phase equation, capillary pressure part
//            drg_dpc -= NpT * dsLdpc * rho_GR * cG1 * p_GR_dot * Np * w; // G1
//            drg_dpc += NpT * dsLdpc * rho_GR * cG2 * T_dot * Np * w; // G2

//            drg_dpc -= NpT * rho_GR * (cG3/dt + dcG3dpc * p_cap_dot) * Np * w; // G3


//            drg_dpc -= NpT * dsLdpc * rho_GR * alpha_B * div_u_dot * Np * w; // G4
//            drg_dpc += gradNpT * rho_GR * k_over_mu_GR * permeability_tensor *
//            		grad_p_GR * Np * w; // G5
//            drg_dpc -= gradNpT * rho_GR * rho_GR * dkrelGdsL * dsLdpc / mu_GR *
//            		permeability_tensor * b * Np * w; // G6
//
//
//            // Gas phase equation, displacement part
//            drg_dus +=  NpT * s_G * rho_GR / dt * alpha_B *
//                    identity2.transpose() * B * w;
//
//            // liquid phase equation, gas pressure part
//            const double cL1 = phi * beta_p_LR * Sps;
//            const double cL2 = phi * beta_T_LR * STs;
//            const double cL3 = phi * dsLdpc - phi_L*beta_p_LR -
//                    s_L*(s_L + p_cap * dsLdpc) * Sps;
//
//            drl_dpg += NpT * s_L * cL1 * (rho_LR / dt + drhoLRdpGR * p_GR_dot) *
//                    Np * w; // L1
//            drl_dpg -= NpT * s_L * drhoLRdpGR * cL2 * T_dot * Np * w; // L2
//            drl_dpg += NpT * drhoLRdpGR * cL3 * p_cap_dot * Np * w; // L3
//            drl_dpg += NpT * s_L * drhoLRdpGR * alpha_B * div_u_dot * Np * w; // L4
//            drl_dpg += gradNpT * k_over_mu_LR * permeability_tensor  *
//                    (drhoLRdpGR - rho_LR/mu_LR - dmuLRdpGR) * grad_p_GR *
//                    Np * w; // L5
//            drl_dpg -= gradNpT * k_over_mu_LR * permeability_tensor *
//                    (drhoLRdpGR - rho_LR/mu_LR - dmuLRdpGR) * grad_p_Cap *
//                    Np * w; // L6
//            drl_dpg -= gradNpT * rho_LR * k_over_mu_LR * permeability_tensor *
//                    (2 * drhoLRdpGR - rho_LR/mu_LR * dmuLRdpGR) * b * Np * w; // L7
//
//
//            // liquid phase equation, capillary pressure part
//            const double dcL3dpc = phi * (d2sLdpc2 - dsLdpc) -
//                    Sps * (dsLdpc * (3*s_L+p_cap*dsLdpc) + s_L*p_cap * d2sLdpc2);
//
//            drl_dpc += NpT * (dsLdpc * rho_LR + s_L * drhoLRdpCap) * cL1 *
//                    p_GR_dot * Np * w; // L1
//            drl_dpc -= NpT * (dsLdpc * rho_LR + s_L * drhoLRdpCap) * cL2 *
//                    T_dot * Np * w; // L2
//            drl_dpc += NpT * (( drhoLRdpCap * cL3 + rho_LR * dcL3dpc) *
//                    p_cap_dot + rho_LR * cL3 / dt) * Np * w; // L3
//            drl_dpc += NpT * s_L * drhoLRdpCap * alpha_B * div_u_dot *
//                    Np * w; // L4
//            drl_dpc += gradNpT * permeability_tensor / mu_LR *
//                    (drhoLRdpCap * k_rel_LR + rho_LR * dkrelLdsL * dsLdpc -
//                            rho_LR * k_rel_LR / mu_LR * dmuLRdpCap) *
//                            grad_p_GR * Np * w; // L5
//            drl_dpc -= gradNpT * permeability_tensor / mu_LR *
//                    (drhoLRdpCap * k_rel_LR + rho_LR * dkrelLdsL * dsLdpc -
//                            rho_LR * k_rel_LR / mu_LR * dmuLRdpCap) *
//                            grad_p_Cap * Np * w; // L6a
//            drl_dpc -= gradNpT * rho_LR * k_over_mu_LR * permeability_tensor *
//                    gradNp * w; // L6b
//            drl_dpc -= gradNpT * rho_LR * permeability_tensor / mu_LR *
//                    (2 * drhoLRdpCap * k_rel_LR + rho_LR * dkrelLdsL * dsLdpc -
//                            rho_LR*k_rel_LR / mu_LR * dmuLRdpCap) * b * Np * w;
//
//            // liquid phase equation, displacement part
//            drl_dus += NpT * s_L * rho_LR / dt * alpha_B *
//                    identity2.transpose() * B * w; // L4
//
//            // displacement equation, gas pressure part
//            dru_dpg -= B.transpose() * alpha_B * identity2 *
//                    Np * w; // U2
//            dru_dpg -= NuT * drhodpGR * b * Np * w; // U4
//
//            // displacement equation, capillary pressure part
//            dru_dpc += B.transpose() * alpha_B * (dsLdpc*p_cap + s_L) *
//                    identity2 * Np * w; // U3
//            dru_dpc -= NuT * drhodpCap * b * Np * w; // U4
//
//            dru_dus += B.transpose() * C * B * w; // U1


#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "  rg  :\n" << rg << "\n";
            std::cout << "----------------------------------\n";
            std::cout << "  d_rg_dpg  :\n" << drg_dpg << "\n";
            std::cout << "  d_rg_dpc  :\n" << drg_dpc << "\n";
            std::cout << "  d_rg_dT   :\n" << drg_dT << "\n";
            std::cout << "  d_rg_dus  :\n" << drg_dus << "\n";
            std::cout << "==================================\n";
#endif



#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "  rl  :\n" << rl << "\n";
            std::cout << "----------------------------------\n";
            std::cout << "  drl_dpg  :\n" << drl_dpg << "\n";
            std::cout << "  drl_dpc  :\n" << drl_dpc << "\n";
            std::cout << "  drl_dT   :\n" << drl_dT << "\n";
            std::cout << "  drl_dus  :\n" << drl_dus << "\n";
            std::cout << "==================================\n";
#endif


#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "  re  :\n" << re << "\n";
            std::cout << "----------------------------------\n";
            std::cout << "  dre_dpg  :\n" << dre_dpg << "\n";
            std::cout << "  dre_dpc  :\n" << dre_dpc << "\n";
            std::cout << "  dre_dT   :\n" << dre_dT << "\n";
            std::cout << "==================================\n";
#endif

#ifdef OUTPUT_IP
            std::cout << "   IP: " << ip << "\n";
            std::cout << "  ru  :\n" << ru << "\n";
            std::cout << "----------------------------------\n";
            std::cout << "  dru_dpg  :\n" << dru_dpg << "\n";
            std::cout << "  dru_dpc  :\n" << dru_dpc << "\n";
            std::cout << "  dru_dT   :\n" << dru_dT << "\n";
            std::cout << "  dru_dus  :\n" << dru_dus << "\n";
            std::cout << "==================================\n";
#endif

        }

#ifdef mDBG_OUTPUT
        std::cout << "#####################\n";
        std::cout << "Mgpg : \n" << Mgpg << "\n";
        std::cout << "Mgpc : \n" << Mgpc << "\n";
        std::cout << "MgT  : \n" << MgT << "\n";
        std::cout << "Mgus : \n" << Mgus << "\n";
        std::cout << "Lgpg : \n" << Lgpg << "\n";
        std::cout << "fg   : \n" << fg << "\n";
        std::cout << "---------------------\n";
        std::cout << "Mlpg : \n" << Mlpg << "\n";
        std::cout << "Mlpc : \n" << Mlpc << "\n";
        std::cout << "MlT  : \n" << MlT << "\n";
        std::cout << "Mlus : \n" << Mlus << "\n";
        std::cout << "Llpg : \n" << Llpg << "\n";
        std::cout << "Llpc : \n" << Llpc << "\n";
        std::cout << "fl   : \n" << fl << "\n";
        std::cout << "---------------------\n";
        std::cout << "Kupg : \n" << Kupg << "\n";
        std::cout << "Kupc : \n" << Kupc << "\n";
        std::cout << "fu   : \n" << fu << "\n";
        std::cout << "---------------------\n";
        std::cout << "Mepg : \n" << MeT << "\n";
        std::cout << "Mepc : \n" << MeT << "\n";
        std::cout << "MeT  : \n" << MeT << "\n";
        std::cout << "Aepg : \n" << Aepg << "\n";
        std::cout << "Aepc : \n" << Aepc << "\n";
        std::cout << "Ket=(AeT+LeT)  : \n" <<  AeT+LeT << "\n";
        std::cout << "AeT  : \n" << AeT << "\n";
        std::cout << "LeT  : \n" << LeT << "\n";

#endif


#ifdef DBG_OUTPUT
        std::cout << "#####################\n";
        std::cout << " Jacobian (analytical): \n";
        std::cout << J << "\n\n";
        std::cout << "=====================\n";
//        std::cout << " residuum: \n\n";
//        std::cout << r << "\n\n";
//        std::cout << "=====================\n";
        OGS_FATAL("Intended stop.");
#endif

#endif



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

    std::vector<double> const& getIntPtRelPermGas(
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
            cache[ip] = _ip_data[ip].rel_perm_gas;
        }

        return cache;
    }


    std::vector<double> const& getIntPtRelPermLiquid(
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
            cache[ip] = _ip_data[ip].rel_perm_liquid;
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

    MPL::Medium const* medium = nullptr;

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
