/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>
#include <algorithm>

#include "MathLib/KelvinVector.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "LocalAssemblerInterface.h"
#include "TH2MProcessData.h"

#include <iostream>

namespace ProcessLib
{
namespace TH2M
{
template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatricesTypePressure, int DisplacementDim, int NPoints>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

    typename ShapeMatrixTypeDisplacement::template MatrixType<
        DisplacementDim, NPoints * DisplacementDim>
        N_u_op;
    typename BMatricesType::KelvinVectorType sigma_eff, sigma_eff_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;

    typename ShapeMatrixTypeDisplacement::NodalRowVectorType N_u;
    typename ShapeMatrixTypeDisplacement::GlobalDimNodalMatrixType dNdx_u;

    typename ShapeMatricesTypePressure::NodalRowVectorType N_p;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType dNdx_p;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;
    double integration_weight;

    void pushBackState()
    {
        eps_prev = eps;
        sigma_eff_prev = sigma_eff;
        material_state_variables->pushBackState();
    }

    template <typename DisplacementVectorType>
    typename BMatricesType::KelvinMatrixType updateConstitutiveRelation(
        double const t,
        SpatialPosition const& x_position,
        double const dt,
        DisplacementVectorType const& /*u*/)
    {
        auto&& solution = solid_material.integrateStress(
            t, x_position, dt, eps_prev, eps, sigma_eff_prev,
            *material_state_variables);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma_eff, material_state_variables, C) = std::move(*solution);

        return C;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

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
          _is_axially_symmetric(is_axially_symmetric),
          _saturation(std::vector<double>(_integration_method.getNumberOfPoints())),
          _pressure_wet(std::vector<double>(_integration_method.getNumberOfPoints()))
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
            _ip_data[ip].integration_weight =
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


        auto p_GR = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_x.data() + gas_pressure_index, gas_pressure_size);

        auto p_cap = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            cap_pressure_size> const>(local_x.data() + cap_pressure_index, cap_pressure_size);

        // create mass matrix
        auto local_M = MathLib::createZeroedMatrix<
            typename ShapeMatricesTypeDisplacement::template MatrixType<
            matrix_size, matrix_size>>(
            local_M_data, matrix_size, matrix_size);

        // create stiffness matrix:
        auto local_K = MathLib::createZeroedMatrix<
            typename ShapeMatricesTypeDisplacement::template MatrixType<
            matrix_size, matrix_size>>(
            local_K_data, matrix_size, matrix_size);

        //create rhs-vector:
        auto local_rhs = MathLib::createZeroedVector<
                    typename ShapeMatricesTypeDisplacement::template VectorType<
                    matrix_size>>(
                    local_rhs_data, matrix_size);




        // assign sub-matrices:
        // water component equation (aka non-wetting-pressure eq.)

        // mass-matrix, gas pressure
        auto Mpp =
                local_M.template block<gas_pressure_size, gas_pressure_size>(gas_pressure_index,
                                      gas_pressure_index);

        // mass-matrix, displacement
        auto Mpu =
                local_M.template block<gas_pressure_size, displacement_size>(gas_pressure_index,
                                        displacement_index);

        // stiffness-matrix, gas pressure
        auto Kpp =
                local_K.template block<gas_pressure_size, gas_pressure_size>(gas_pressure_index,
                            gas_pressure_index);

        // rhs-vector, gas pressure
        auto fp =
                local_rhs.template segment<gas_pressure_size>(gas_pressure_index);


        // displacement equation
        // stiffness-matrix, gas pressure
        auto Kup =
                local_K.template block<displacement_size, gas_pressure_size>(displacement_index,
                             gas_pressure_index);

        // rhs-vector, displacement
        auto fu =
                local_rhs.template segment<displacement_size>(displacement_index);



//         ********************************************************************
//        from H2:
//         ********************************************************************

        auto Mgp =
              local_M.template block<gas_pressure_size, gas_pressure_size>(
                  gas_pressure_index, gas_pressure_index);

          auto Mgpc = local_M.template block<gas_pressure_size, cap_pressure_size>(
              gas_pressure_index, cap_pressure_index);

          auto Mlpc = local_M.template block<cap_pressure_size, cap_pressure_size>(
              cap_pressure_index, cap_pressure_index);

//          NodalMatrixType laplace_operator =
//              NodalMatrixType::Zero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

          typename ShapeMatricesTypePressure::NodalMatrixType laplace_operator =
              ShapeMatricesTypePressure::NodalMatrixType::Zero(gas_pressure_size,
                                                               gas_pressure_size);


          auto Kgp =
              local_K.template block<gas_pressure_size, gas_pressure_size>(
                  gas_pressure_index, gas_pressure_index);

          auto Klp = local_K.template block<cap_pressure_size, gas_pressure_size>(
              cap_pressure_index, gas_pressure_index);

          auto Klpc = local_K.template block<cap_pressure_size, cap_pressure_size>(
              cap_pressure_index, cap_pressure_index);

          auto Bg = local_rhs.template segment<gas_pressure_size>(
              gas_pressure_index);

          auto Bl = local_rhs.template segment<cap_pressure_size>(cap_pressure_index);

//
//          auto gravity_operator =
//                  local_rhs.template segment<gas_pressure_size>(gas_pressure_index);


//         ********************************************************************



        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        double const& dt = _process_data.dt;

        auto u = Eigen::Map<typename ShapeMatricesTypeDisplacement::
                                                    template VectorType<displacement_size> const>(
                                local_x.data() + displacement_index, displacement_size);

//
//        const Eigen::MatrixXd& perm = _process_data.material->getPermeability(
//            material_id, t, pos, _element.getDimension());
//        assert(perm.rows() == _element.getDimension() || perm.rows() == 1);
//        GlobalDimMatrixType permeability = GlobalDimMatrixType::Zero(
//            _element.getDimension(), _element.getDimension());
//        if (perm.rows() == _element.getDimension())
//            permeability = perm;
//        else if (perm.rows() == 1)
//            permeability.diagonal().setConstant(perm(0, 0));

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;

            auto const& N_u_op = _ip_data[ip].N_u_op;

            auto const& N_u = _ip_data[ip].N_u;
            auto const& dNdx_u = _ip_data[ip].dNdx_u;

            auto const& N_p = _ip_data[ip].N_p;
            auto const& dNdx_p = _ip_data[ip].dNdx_p;

            auto const& mass_operator = N_p.transpose()*N_p*w;


            auto const x_coord =
                interpolateXCoordinate<ShapeFunctionDisplacement,
                                       ShapeMatricesTypeDisplacement>(_element,
                                                                      N_u);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunctionDisplacement::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx_u, N_u, x_coord,
                                                     _is_axially_symmetric);


            auto& eps = _ip_data[ip].eps;
            auto const& sigma_eff = _ip_data[ip].sigma_eff;

            double const S = _process_data.specific_storage(t, x_position)[0];
            double const K_over_mu =
                _process_data.intrinsic_permeability(t, x_position)[0] /
                _process_data.fluid_viscosity(t, x_position)[0];

            auto const permeability = _process_data.intrinsic_permeability(t, x_position)[0];

            using permTensorType = Eigen::Matrix<double, DisplacementDim, DisplacementDim>;

            permTensorType perm_tensor;

            perm_tensor.template block<2,
            2>(0,0).setIdentity();

            auto permeability_tensor = permeability * perm_tensor;


            auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
            auto const rho_sr = _process_data.solid_density(t, x_position)[0];
            auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
            // auto const porosity = _process_data.porosity(t, x_position)[0];
            double const porosity = 0.2975;
            auto const rho = rho_sr * (1 - porosity) + porosity * rho_fr;
            auto const& b = _process_data.specific_body_force;



            auto const& identity2 = MathLib::KelvinVector::Invariants<
                    MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value>::identity2;


            eps.noalias() = B * u;

            auto C = _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u);


//            // mass balance:
//
//            Mpp.noalias() += N_p.transpose() * S * N_p * w;
//            Mpu.noalias() += N_p.transpose() * alpha * identity2.transpose() * B * w;
//
//            Kpp.noalias() += dNdx_p.transpose() * K_over_mu * dNdx_p * w;
//
//            // Kpu = 0
//
//            fp.noalias() += dNdx_p.transpose() * rho_fr * K_over_mu * b * w;
//
//            // momentum balance:
//
//            //Muu = Mup = 0
//            Kup.noalias() -= B.transpose() * alpha * identity2 * N_p * w;
//            // Kuu = 0
//            fu.noalias() -= (B.transpose() * sigma_eff - N_u_op.transpose() * rho * b) * w;

            auto const pc_int_pt = p_cap.dot(N_p);
            auto const pn_int_pt = p_GR.dot(N_p);


            _pressure_wet[ip] = pn_int_pt - pc_int_pt;

            const double temperature = 293.15;
            const double molar_mass  = 0.02896;
            const double gas_constant = 8.3144621;

            double const rho_nonwet = molar_mass * pn_int_pt / gas_constant / temperature;

            double const rho_wet = 1000.;

            const double m_sw =-1.00623E-05;
            const double n_sw = 1.0;

            double const Sw = m_sw * pc_int_pt + n_sw;

            _saturation[ip] = Sw;
            double dSw_dpc = m_sw;

            double const drhononwet_dpn = molar_mass / gas_constant / temperature;

            const double _sw = std::min(std::max(0.0,Sw),0.8);

            double const Se = _sw/0.8;

            // gas
            double const k_rel_nonwet = std::max(0.0001,(1.-Se)*(1.-Se)*(1-std::pow(Se, 1.+2./3.)));

            double const mu_nonwet = 1.8e-5;
            double const lambda_nonwet = k_rel_nonwet / mu_nonwet;

            // liquid
            double const m_krel = 2.352941176;
            double const n_krel = -1.352941176;
            double const k_rel_wet = m_krel*Sw+n_krel;
            double const mu_wet = 1.0e-3;
            double const lambda_wet = k_rel_wet / mu_wet;



            // Assemble M matrix
            // nonwetting

            Mgp.noalias() +=
                    porosity * (1 - Sw) * drhononwet_dpn * N_p.transpose() * N_p * w;
            Mgpc.noalias() +=
                    -porosity * rho_nonwet * dSw_dpc * N_p.transpose() * N_p * w;

            Mlpc.noalias() +=
                    porosity * dSw_dpc * rho_wet * N_p.transpose() * N_p * w;

            laplace_operator.noalias() = dNdx_p.transpose() *
                    permeability_tensor * dNdx_p *
                    w;

            Kgp.noalias() += rho_nonwet * lambda_nonwet * laplace_operator;

            Klp.noalias() += rho_wet * lambda_wet * laplace_operator;
            Klpc.noalias() += -rho_wet * lambda_wet * laplace_operator;



// Das kompiliert, aber der erste Eintrag der Ergebnismatrix ist null!

//          auto const gravity_operator = dNdx_p.transpose() * permeability_tensor * b * w;
//
// Viel schlimmer noch, bei Ausgabe der einzelnen Komponenten mittels
// for r { for c { std::cout << gravity_operator(r,c);}} kommen andere
// Ergebnisse raus als bei std::cout << gravity_operator;



// Das geht nicht, weil ich angeblich Matrizen unterschiedlicher Größe
// verwenden würde (?)

//          Eigen::Matrix<double, 4, 1> gravity_operator;
//          gravity_operator.setZero();
//          gravity_operator=dNdx_p.transpose() * permeability_tensor * b * w;
//
//  btw.: die Multiplikation geht, der Fehler tritt erst bei der Berechnung
//  von Bg auf.
//
//            Bg.noalias() +=
//                    rho_nonwet * rho_nonwet * lambda_nonwet * gravity_operator;
//
//            Bl.noalias() += rho_wet * rho_wet * lambda_wet * gravity_operator;


//  Das hier geht immerhin.

            Bg.noalias() +=
                    rho_nonwet * rho_nonwet * lambda_nonwet *
                    dNdx_p.transpose() * permeability_tensor * b * w;

            Bl.noalias() += rho_wet * rho_wet * lambda_wet *
                    dNdx_p.transpose() * permeability_tensor * b * w;


//            std::cout << "------------------------\n";
//            std::cout << "Constitutive parameters:\n\n";
//            std::cout << "------------------------\n";
//
//            std::cout << "gas_pressure         : " <<  pn_int_pt << "\n";
//            std::cout << "capillary_pressure   : " <<  pc_int_pt << "\n";
//            std::cout << "liquid_pressure      : " << _pressure_wet[ip] << "\n";
//
//            std::cout << "temperature          : " << temperature << "\n";
//            std::cout << "density_gas          : " << rho_nonwet << "\n";
//            std::cout << "density_liquid       : " << rho_wet << "\n";
//
//            std::cout << "saturation_liquid    : " << Sw << "\n";
//
//            std::cout << "dSw/dpc              : " << dSw_dpc << "\n";
//            std::cout << "porosity             : " << porosity << "\n";
//
//            std::cout << "drho_gas/dp_gas      : " << drhononwet_dpn  << "\n";
//            std::cout << "rel_permeability_gas : " << k_rel_nonwet << "\n";
//            std::cout << "viscosity_gas        : " << mu_nonwet << "\n";
//
//            std::cout << "rel_permeability_liq : " << k_rel_wet << "\n";
//            std::cout << "viscosity_liqid      : " << mu_wet << "\n";
//
//            std::cout << "permeability_tensor  :\n";
//            std::cout << permeability_tensor << "\n";
//
//            std::cout << "specific_body_force  :\n";
//            std::cout << b << "\n";
//
//            std::cout << "\n";
//
//            std::cout << "-------------------------\n";
//            std::cout << "Element shape_functtions:\n\n";
//            std::cout << "-------------------------\n\n";
//
//            std::cout << "N_p                  :\n";
//            std::cout << _ip_data[ip].N_p << "\n";
//
//            std::cout << "dNdx_p               :\n";
//            std::cout << _ip_data[ip].dNdx_p << "\n";
//
//            std::cout << "dNdx_p.transpose()   :\n";
//            std::cout << _ip_data[ip].dNdx_p.transpose() << "\n";
//
//            std::cout << "integration_weight   :" << _ip_data[ip].integration_weight << "\n";
//            std::cout << "\n";
//
//
//            std::cout << "--------------------------------\n";
//            std::cout << "Auxiliary matrices and vectors :\n\n";
//            std::cout << "--------------------------------\n\n";
//
//            std::cout << "laplace_operator     :\n";
//            std::cout << laplace_operator << "\n";
//            std::cout << "\n";


//            std::cout << "-------------------------------\n";
//            std::cout << "Element mass matrix components:\n\n";
//            std::cout << "-------------------------------\n\n";
//
//            std::cout << "Mgpg:\n";
//            std::cout << Mgp << "\n";
//            std::cout << "Mgpc:\n";
//            std::cout << Mgpc << "\n";
//            std::cout << "Mlpc:\n";
//            std::cout << Mlpc << "\n";
//            std::cout << "\n";
//
//            std::cout << "------------------------------------\n";
//            std::cout << "Element stiffness matrix components:\n\n";
//            std::cout << "------------------------------------\n\n";
//
//            std::cout << "Kgpg:\n";
//            std::cout << Kgp << "\n";
//            std::cout << "Klpg:\n";
//            std::cout << Klp << "\n";
//            std::cout << "Klpc:\n";
//            std::cout << Klpc << "\n";
//            std::cout << "\n";
//
//            std::cout << "------------------------------\n";
//            std::cout << "Element rhs vector components:\n\n";
//            std::cout << "------------------------------\n\n";
//
//            std::cout << "fg:\n";
//            std::cout << Bg << "\n";
//            std::cout << "fl:\n";
//            std::cout << Bl << "\n";
//            std::cout << "\n";
//
//            std::cout << "------------------------------\n\n";


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
        assert(local_x.size() == gas_pressure_size + cap_pressure_size +
                temperature_size + displacement_size);

        auto p_GR =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                      gas_pressure_size);

        auto p_cap =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                cap_pressure_size> const>(local_x.data() + cap_pressure_index,
                                      cap_pressure_size);

        auto T =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                temperature_size> const>(local_x.data() + temperature_index,
                                      temperature_size);

        auto u = Eigen::Map<typename ShapeMatricesTypeDisplacement::
                                template VectorType<displacement_size> const>(
            local_x.data() + displacement_index, displacement_size);


        auto p_GR_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                gas_pressure_size> const>(local_xdot.data() + gas_pressure_index,
                                      gas_pressure_size);

        auto p_cap_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                cap_pressure_size> const>(local_x.data() + cap_pressure_index,
                                      cap_pressure_size);

        auto T_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                temperature_size> const>(local_x.data() + temperature_index,
                                      temperature_size);

        auto u_dot =
            Eigen::Map<typename ShapeMatricesTypeDisplacement::
                           template VectorType<displacement_size> const>(
                local_xdot.data() + displacement_index, displacement_size);


        const auto matrix_size = gas_pressure_size + cap_pressure_size +
                temperature_size + displacement_size;

        auto local_Jac = MathLib::createZeroedMatrix<
            typename ShapeMatricesTypeDisplacement::template MatrixType<
            matrix_size, matrix_size>>
            ( local_Jac_data, matrix_size, matrix_size);

        auto local_rhs = MathLib::createZeroedVector<
            typename ShapeMatricesTypeDisplacement::template VectorType<
            matrix_size>>(
            local_rhs_data, matrix_size);

        typename ShapeMatricesTypePressure::NodalMatrixType laplace_p =
            ShapeMatricesTypePressure::NodalMatrixType::Zero(gas_pressure_size,
                    gas_pressure_size);

        typename ShapeMatricesTypePressure::NodalMatrixType storage_p =
            ShapeMatricesTypePressure::NodalMatrixType::Zero(gas_pressure_size,
                    gas_pressure_size);

        typename ShapeMatricesTypeDisplacement::template MatrixType<
        Nu_intPoints, Np_intPoints>
            Kup = ShapeMatricesTypeDisplacement::template MatrixType<
                displacement_size, gas_pressure_size>::Zero(displacement_size,
                                                        gas_pressure_size);

        double const& dt = _process_data.dt;

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;

            auto const& N_u_op = _ip_data[ip].N_u_op;

            auto const& N_u = _ip_data[ip].N_u;
            auto const& dNdx_u = _ip_data[ip].dNdx_u;

            auto const& N_p = _ip_data[ip].N_p;
            auto const& dNdx_p = _ip_data[ip].dNdx_p;

            auto const x_coord =
                interpolateXCoordinate<ShapeFunctionDisplacement,
                                       ShapeMatricesTypeDisplacement>(_element,
                                                                      N_u);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunctionDisplacement::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx_u, N_u, x_coord,
                                                     _is_axially_symmetric);

            auto& eps = _ip_data[ip].eps;
            auto const& sigma_eff = _ip_data[ip].sigma_eff;

            double const S = _process_data.specific_storage(t, x_position)[0];
            double const K_over_mu =
                _process_data.intrinsic_permeability(t, x_position)[0] /
                _process_data.fluid_viscosity(t, x_position)[0];
            auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
            auto const rho_sr = _process_data.solid_density(t, x_position)[0];
            auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
            auto const porosity = _process_data.porosity(t, x_position)[0];
            auto const& b = _process_data.specific_body_force;
            auto const& identity2 = MathLib::KelvinVector::Invariants<
                    MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value>::identity2;

            //
            // displacement equation, displacement part
            //
            eps.noalias() = B * u;

            auto C =
                _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u);

            local_Jac
                .template block<displacement_size, displacement_size>(
                    displacement_index, displacement_index)
                .noalias() += B.transpose() * C * B * w;

            double const rho = rho_sr * (1 - porosity) + porosity * rho_fr;
            local_rhs.template segment<displacement_size>(displacement_index)
                .noalias() -=
                (B.transpose() * sigma_eff - N_u_op.transpose() * rho * b) * w;

            //
            // displacement equation, pressure part
            //
            Kup.noalias() += B.transpose() * alpha * identity2 * N_p * w;

            //
            // pressure equation, pressure part.
            //
            laplace_p.noalias() += dNdx_p.transpose() * K_over_mu * dNdx_p * w;

            storage_p.noalias() += N_p.transpose() * S * N_p * w;

            local_rhs.template segment<gas_pressure_size>(gas_pressure_index)
                .noalias() += dNdx_p.transpose() * rho_fr * K_over_mu * b * w;

            //
            // pressure equation, displacement part.
            //
            // Reusing Kup.transpose().
        }
        // displacement equation, pressure part
        local_Jac
            .template block<displacement_size, gas_pressure_size>(
                displacement_index, gas_pressure_index)
            .noalias() = -Kup;

        // pressure equation, pressure part.
        local_Jac
            .template block<gas_pressure_size, gas_pressure_size>(gas_pressure_index,
                    gas_pressure_index)
            .noalias() = laplace_p + storage_p / dt;

        // pressure equation, displacement part.
        local_Jac
            .template block<gas_pressure_size, displacement_size>(
                    gas_pressure_index, displacement_index)
            .noalias() = Kup.transpose() / dt;

        // pressure equation
        local_rhs.template segment<gas_pressure_size>(gas_pressure_index).noalias() -=
            laplace_p * p_GR + storage_p * p_GR_dot + Kup.transpose() * u_dot;

        // displacement equation
        local_rhs.template segment<displacement_size>(displacement_index)
            .noalias() += Kup * p_GR;



        local_Jac.template block<cap_pressure_size+temperature_size,
        cap_pressure_size+temperature_size>(cap_pressure_index, cap_pressure_index).setIdentity();

   //     local_rhs.template segment<cap_pressure_size+temperature_size>(cap_pressure_index).setOnes();

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
            auto const& sigma = _ip_data[ip].sigma;
            cache_mat.col(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma);
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
        std::vector<double>& /*cache*/) const override
    {
        assert(!_saturation.empty());
        return _saturation;
    }

    std::vector<double> const& getIntPtWetPressure(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_pressure_wet.empty());
        return _pressure_wet;
    }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        GlobalVector const& current_solution,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::vector<double>& cache) const override
    {
        auto const num_intpts = _ip_data.size();

        auto const indices = NumLib::getIndices(_element.getID(), dof_table);
        assert(!indices.empty());
        auto const local_x = current_solution.get(indices);

        cache.clear();
        auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, DisplacementDim, num_intpts);

        SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto p_GR =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                      gas_pressure_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            double const K_over_mu =
                _process_data.intrinsic_permeability(t, x_position)[0] /
                _process_data.fluid_viscosity(t, x_position)[0];

            auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
            auto const& b = _process_data.specific_body_force;

            // Compute the velocity
            auto const& dNdx_p = _ip_data[ip].dNdx_p;
            cache_matrix.col(ip).noalias() =
                -K_over_mu * dNdx_p * p_GR - K_over_mu * rho_fr * b;
        }

        return cache;
    }

private:
    std::vector<double> const& getIntPtSigma(std::vector<double>& cache,
                                             std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            if (component < 3)  // xx, yy, zz components
                cache.push_back(ip_data.sigma_eff[component]);
            else  // mixed xy, yz, xz components
                cache.push_back(ip_data.sigma_eff[component] / std::sqrt(2));
        }

        return cache;
    }

    std::vector<double> const& getIntPtEpsilon(
        std::vector<double>& cache, std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            cache.push_back(ip_data.eps[component]);
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


    // output vector for wetting phase saturation and pressure with
    // respect to each integration point
    std::vector<double> _saturation;
    std::vector<double> _pressure_wet;



    static const int Np_intPoints = ShapeFunctionPressure::NPOINTS;
    static const int gas_pressure_index = 0;
    static const int cap_pressure_index = 1*Np_intPoints;
    static const int temperature_index = 2*Np_intPoints;
    static const int displacement_index = 3*Np_intPoints;

    static const int gas_pressure_size = Np_intPoints;
    static const int cap_pressure_size = Np_intPoints;
    static const int temperature_size = Np_intPoints;

    static const int Nu_intPoints = ShapeFunctionDisplacement::NPOINTS * DisplacementDim;
    static const int displacement_size = Nu_intPoints;
    static const int kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
};

}  // namespace TH2M
}  // namespace ProcessLib
