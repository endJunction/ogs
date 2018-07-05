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

namespace MPL = MaterialPropertyLib;

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

   // typename BMatricesType::KelvinVectorType sigma_prev;

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
          _pressure_wet(std::vector<double>(_integration_method.getNumberOfPoints())),
          _density_gas(std::vector<double>(_integration_method.getNumberOfPoints())),
          _density_liquid(std::vector<double>(_integration_method.getNumberOfPoints())),
          _pressure_gas_linear(std::vector<double>(_integration_method.getNumberOfPoints())),
          _pressure_cap_linear(std::vector<double>(_integration_method.getNumberOfPoints())),
          _velocity_gas(std::vector<double>(DisplacementDim*_integration_method.getNumberOfPoints())),
          _velocity_liquid(std::vector<double>(DisplacementDim*_integration_method.getNumberOfPoints()))
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

        // primary variables
        auto gas_phase_pressure = Eigen::Map<typename
                ShapeMatricesTypePressure::template
                VectorType<gas_pressure_size> const>(local_x.data()
                        + gas_pressure_index, gas_pressure_size);

        auto capillary_pressure = Eigen::Map<typename
                ShapeMatricesTypePressure::template
                VectorType<cap_pressure_size> const>(local_x.data()
                        + cap_pressure_index, cap_pressure_size);

        auto temperature = Eigen::Map<typename
                ShapeMatricesTypePressure::template
                VectorType<temperature_size> const>(local_x.data()
                        + temperature_index, temperature_size);

        auto displacement = Eigen::Map<typename
                ShapeMatricesTypeDisplacement::template
                VectorType<displacement_size> const>(local_x.data()
                        + displacement_index, displacement_size);

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

//         ********************************************************************
//         ********************************************************************

          auto Mgpg = local_M.template block<gas_pressure_size, gas_pressure_size>(
                  gas_pressure_index, gas_pressure_index);

          auto Mgpc = local_M.template block<gas_pressure_size, cap_pressure_size>(
              gas_pressure_index, cap_pressure_index);

          auto MgT = local_M.template block<gas_pressure_size, temperature_size>(
              gas_pressure_index, temperature_index);

          auto Mgus = local_M.template block<gas_pressure_size, displacement_size>(
                        gas_pressure_index, displacement_index);

          auto Mlpg = local_M.template block<cap_pressure_size, gas_pressure_size>(
              cap_pressure_index, gas_pressure_index);

          auto Mlpc = local_M.template block<cap_pressure_size, cap_pressure_size>(
              cap_pressure_index, cap_pressure_index);

          auto MlT = local_M.template block<cap_pressure_size, temperature_size>(
              cap_pressure_index, temperature_index);

          auto Mlus = local_M.template block<cap_pressure_size, displacement_size>(
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

          auto Kgpg = local_K.template block<gas_pressure_size, gas_pressure_size>(
                          gas_pressure_index, gas_pressure_index);

          auto Klpg = local_K.template block<cap_pressure_size, gas_pressure_size>(
                  cap_pressure_index, gas_pressure_index);

          auto Klpc = local_K.template block<cap_pressure_size, cap_pressure_size>(
                  cap_pressure_index, cap_pressure_index);

//          auto KTT = local_K.template block<temperature_size, temperature_size>(
//                  temperature_index, temperature_index);
//
//          auto KTpg = local_K.template block<gas_pressure_size, temperature_size>(
//                                    gas_pressure_index, temperature_index);
//
//          auto KTpc = local_K.template block<cap_pressure_size, temperature_size>(
//                                    cap_pressure_index, temperature_index);


          auto Kupg =
                  local_K.template block<displacement_size, gas_pressure_size>(
                          displacement_index, gas_pressure_index);
          auto Kupc =
                  local_K.template block<displacement_size, cap_pressure_size>(
                          displacement_index, cap_pressure_index);



          auto Bg = local_rhs.template segment<gas_pressure_size>(
              gas_pressure_index);

          auto Bl = local_rhs.template segment<cap_pressure_size>(cap_pressure_index);

          auto Bu = local_rhs.template segment<displacement_size>(displacement_index);
//
//          auto gravity_operator =
//                  local_rhs.template segment<gas_pressure_size>(gas_pressure_index);


//         ********************************************************************



        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        double const& dt = _process_data.dt;


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


        _velocity_gas.clear();
        auto velocity_matrix_gas = MathLib::createZeroedMatrix<Eigen::Matrix<
                double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
                        _velocity_gas, DisplacementDim, _ip_data.size());

        _velocity_liquid.clear();
        auto velocity_matrix_liquid = MathLib::createZeroedMatrix<Eigen::Matrix<
                double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
                        _velocity_liquid, DisplacementDim, _ip_data.size());
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

//            double const S = _process_data.specific_storage(t, x_position)[0];
//            double const K_over_mu =
//                _process_data.intrinsic_permeability(t, x_position)[0] /
//                _process_data.fluid_viscosity(t, x_position)[0];



//            auto const permeability = _process_data.intrinsic_permeability(t, x_position)[0];


            // get porous medium properties from process data
            auto const& medium = _process_data.medium;

            // reset update status of all properties
            medium->resetPropertyUpdateStatus();


            auto const p_cap = capillary_pressure.dot(N_p);
            auto const p_GR = gas_phase_pressure.dot(N_p);
            auto const T = temperature.dot(N_p);
            auto const u = 0.;//displacement.dot(N_u);

            const double p_LR = p_GR - p_cap;

            std::cout << "==================================\n";
            std::cout << " Integration Point:\n";
            std::cout << " p_cap: " << p_cap << "\n";
            std::cout << " p_GR: " << p_GR << "\n";
            std::cout << " T: " << T << "\n";
            std::cout << " u: " << u << "\n";
            std::cout << "==================================\n";

            // insert all primary variables into one object
            MPL::VariableArray primaryVariables;
            primaryVariables[MPL::p_cap] = p_cap;
            primaryVariables[MPL::p_GR] = p_GR;
            primaryVariables[MPL::p_LR] = p_LR;
            primaryVariables[MPL::T] = 298.15;
            // todo: displacement

            // get fluid phase properties
            auto const& solid_phase = medium->phase(0);
            auto const& liquid_phase = medium->phase(1);
            auto const& gas_phase = medium->phase(2);

            // intrinsic permeability
            double const permeability =
                    MPL::getScalar(medium->property(MPL::permeability));

            std::cout << "   permeability: " << permeability << " \n";

            Eigen::Matrix<double, DisplacementDim, DisplacementDim> perm_tensor;
            perm_tensor.setIdentity();

            auto const permeability_tensor = permeability * perm_tensor;

            auto const alpha_B =
                    MPL::getScalar(medium->property(
                            MPL::PropertyEnum::biot_coefficient),
                            primaryVariables);

            std::cout << "   alpha_B: " << alpha_B << " \n";
            std::cout << "==================================\n";

            auto const rho_SR =
                    MPL::getScalar(solid_phase.property(
                            MPL::PropertyEnum::density),
                            primaryVariables);

            std::cout << "   rho_SR: " << rho_SR << " \n";

            auto const rho_LR =
                    MPL::getScalar(liquid_phase.property(
                            MPL::PropertyEnum::density),
                            primaryVariables);

            std::cout << "   rho_LR: " << rho_LR << " \n";

            auto const rho_GR =
                    MPL::getScalar(gas_phase.property(
                            MPL::PropertyEnum::density),
                            primaryVariables);

            std::cout << "   rho_GR: " << rho_GR << " \n";

            std::cout << "==================================\n";
            auto const phi =
                    MPL::getScalar(medium->property(
                            MPL::PropertyEnum::porosity),
                            primaryVariables);


            std::cout << "  phi : " << phi << " \n";

            auto const phi_S = 1. - phi;

            std::cout << "  phi_S : " << phi_S << " \n";

            std::cout << "==================================\n";

            auto const s_L =
                    MPL::getScalar(medium->property(
                            MPL::PropertyEnum::saturation),
                            primaryVariables);

            std::cout << "  s_L : " << s_L << " \n";

            auto const s_G = 1 - s_L;

            std::cout << "  s_G: " << s_G << " \n";

            auto const mu_LR = MPL::getScalar(liquid_phase.property(
                    MPL::PropertyEnum::viscosity),
                    primaryVariables);

            std::cout << "  mu_LR: " << mu_LR << " \n";

            auto const mu_GR = MPL::getScalar(gas_phase.property(
                    MPL::PropertyEnum::viscosity),
                    primaryVariables);

            std::cout << "  mu_GR: " << mu_GR << " \n";
            std::cout << "==================================\n";

            auto const rho = phi_S*rho_SR + phi*(s_L*rho_LR + s_G*rho_GR);

            std::cout << "  rho: " <<  rho<< " \n";

            auto const& b = _process_data.specific_body_force;

            static int const KelvinVectorSize =
                    MathLib::KelvinVector::KelvinVectorDimensions
                    <DisplacementDim>::value;

            auto const& identity2 = MathLib::KelvinVector::Invariants
                    <KelvinVectorSize>::identity2;

            eps.noalias() = B * displacement;

            auto const e = MathLib::KelvinVector::Invariants
                    <KelvinVectorSize>::trace(eps);

            std::cout << "   e: " << e << " \n";
            std::cout << "==================================\n";

            auto C = _ip_data[ip].updateConstitutiveRelation(
                    t, x_position, dt, displacement);

            auto const drhoGRdpGR = MPL::getScalarDerivative(
                    gas_phase.property(MPL::density),
                    primaryVariables, MPL::p_GR);

            auto const drhoLRdpLR = MPL::getScalarDerivative(
                    liquid_phase.property(MPL::density),
                    primaryVariables, MPL::p_LR);

            auto const molar_mass = getScalar(gas_phase.property(
                    MPL::PropertyEnum::molar_mass));

            std::cout << "   M_g : " << molar_mass << " \n";

            std::cout << "==================================\n";

            std::cout << "   drho_gr_dp_gr : " << drhoGRdpGR << " \n";
            std::cout << "   drho_lr_dp_lr : " << drhoLRdpLR << " \n";


            std::cout << "==================================\n";

            auto const dsLdpc =
                    MPL::getScalarDerivative(medium->property(
                            MPL::PropertyEnum::saturation), primaryVariables,
                            MPL::PrimaryVariables::p_cap);
            std::cout << "   dsLdpc : " << dsLdpc << " \n";
            std::cout << "==================================\n";


            auto const k_rel_LR = MPL::getPair(medium->property(
                    MPL::PropertyEnum::relative_permeability),
                    primaryVariables)[0];

            std::cout << "    k_rel_LR: " << k_rel_LR << " \n";

            auto const k_rel_GR = MPL::getPair(medium->property(
                    MPL::PropertyEnum::relative_permeability),
                    primaryVariables)[1];

            std::cout << "    k_rel_GR: " << k_rel_GR << " \n";

            std::cout << "==================================\n";

            // for output only
            _pressure_gas_linear[ip] = p_GR;
            _pressure_cap_linear[ip] = p_cap;
            _pressure_wet[ip] = p_LR;
            _density_gas[ip] = rho_GR;
            _density_liquid[ip] = rho_LR;
            _saturation[ip] = s_L;

            const double lambda_GR = k_rel_GR / mu_GR;
            const double lambda_LR = k_rel_LR / mu_LR;

            const double beta_p_SR = getScalar(solid_phase.property(
                    MPL::PropertyEnum::compressibility));
            std::cout << "    beta_p_SR: " << beta_p_SR << " \n";

            const double beta_p_PR = getScalar(medium->property(
                    MPL::PropertyEnum::compressibility));
            std::cout << "    beta_p_PR: " << beta_p_PR << " \n";

            const double beta_p_GR = 1./rho_GR*drhoGRdpGR;
            std::cout << "   beta_p_GR: " << beta_p_GR << " \n";

            const double beta_p_LR = 1./rho_LR*drhoLRdpLR;
            std::cout << "   beta_p_LR: " << beta_p_LR << " \n";

            std::cout << "==================================\n";

            velocity_matrix_liquid.col(ip).noalias() =
                    -permeability_tensor*lambda_LR * dNdx_p *
                    (gas_phase_pressure-capillary_pressure) +
                    permeability_tensor*lambda_LR * rho_LR * b;

            velocity_matrix_gas.col(ip).noalias() =
                    -permeability_tensor*lambda_GR * dNdx_p *
                    gas_phase_pressure +
                    permeability_tensor*lambda_GR * rho_GR * b;

            Laplace.noalias() = dNdx_p.transpose() * permeability_tensor * dNdx_p * w;

            Mgpg.noalias() += N_p.transpose()*
                    rho_GR*s_G*(phi*beta_p_GR+(alpha_B-phi)*beta_p_SR)*N_p*w;

            std::cout << "   Mgpg:\n " << Mgpg << " \n";
            std::cout << "==================================\n";

            Mgpc.noalias() -= N_p.transpose()*
                    rho_GR*( (phi + s_G*(alpha_B-phi)*beta_p_SR*p_cap)*dsLdpc
                    +s_G*s_L*(alpha_B-phi)*beta_p_SR)*N_p*w;

            std::cout << "   Mgpc:\n " << Mgpc << " \n";
            std::cout << "==================================\n";


//           Mgus.noalias() += N_p.transpose()*rho_GR*alpha_B*s_G*m.transpose()*B*w;
//
            Mlpg.noalias() += N_p.transpose()*
                    rho_LR*s_L*(phi*beta_p_LR+(alpha_B-phi)*beta_p_SR)*N_p*w;

            std::cout << "   Mlpg:\n " << Mlpg << " \n";
            std::cout << "==================================\n";

            Mlpc.noalias() += N_p.transpose()*
                    rho_LR*( (phi - s_L*(alpha_B-phi)*beta_p_SR*p_cap)*dsLdpc
                    -phi*s_L*beta_p_LR - (alpha_B-phi)*beta_p_SR*s_L*s_L)
                    *N_p*w;

            std::cout << "   Mlpc:\n " << Mlpc << " \n";
            std::cout << "==================================\n";

//            Mlus.noalias() += N_p.transpose()*rho_LR*s_L*alpha_B*m.transpose()*B*w;
//
            Kgpg.noalias() += rho_GR*lambda_GR*Laplace;

            std::cout << "   Kgpg:\n " << Kgpg << " \n";
            std::cout << "==================================\n";

            Klpg.noalias() += rho_LR*lambda_LR*Laplace;

            std::cout << "   Klpg:\n " << Klpg << " \n";
            std::cout << "==================================\n";


            Klpc.noalias() -= rho_LR*lambda_LR*Laplace;

            std::cout << "   Klpc:\n " << Klpc << " \n";
            std::cout << "==================================\n";


//            Kupg.noalias() += B.transpose()*alpha*m*N_p*w;
//            Kupc.noalias() -= B.transpose()*alpha*m*s_L*N_p*w;
//
//            // Gravity
//
//            auto const gravity_operator = (dNdx_p.transpose() *
//                                permeability_tensor * b * w).eval();
//
//            Bg.noalias() += rho_GR*rho_GR*lambda_GR*gravity_operator;
//            Bl.noalias() += rho_LR*rho_LR*lambda_LR*gravity_operator;
//            Bu.noalias() -= (B.transpose()*sigma_eff-N_u_op.transpose()*rho*b)*w;



           // OGS_FATAL("Intended halt. (numerical)");

        }


//        OGS_FATAL("Intended halt. (numerical)");
//        std::cout << "== Local M: ====\n";
//        std::cout << local_M << "\n";
//        std::cout << "================\n";
//        std::cout << "== Local K: ====\n";
//        std::cout << local_K << "\n";
//        std::cout << "================\n";
//        std::cout << "== Local f: ====\n";
//        std::cout << local_rhs << "\n";
//        std::cout << "================\n";
//
//        OGS_FATAL("Intended halt.");



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
                cap_pressure_size> const>(local_xdot.data() + cap_pressure_index,
                        cap_pressure_size);

        auto T_dot =
                Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                temperature_size> const>(local_xdot.data() + temperature_index,
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


        typename ShapeMatricesTypePressure::template MatrixType<
        Nu_intPoints, Np_intPoints> Kup = ShapeMatricesTypePressure::template MatrixType<
                displacement_size, gas_pressure_size>::Zero(displacement_size,
                                                        gas_pressure_size);



        typename ShapeMatricesTypePressure::template MatrixType<gas_pressure_size, gas_pressure_size>
        Mgpg = ShapeMatricesTypePressure::template  MatrixType
        <gas_pressure_size, gas_pressure_size>::Zero(gas_pressure_size, gas_pressure_size);

        typename ShapeMatricesTypePressure::template MatrixType<gas_pressure_size, gas_pressure_size>
        Mgpc = ShapeMatricesTypePressure::template  MatrixType
        <gas_pressure_size, gas_pressure_size>::Zero(gas_pressure_size, gas_pressure_size);

        typename ShapeMatricesTypePressure::template MatrixType<gas_pressure_size, gas_pressure_size>
        Mlpg = ShapeMatricesTypePressure::template  MatrixType
        <gas_pressure_size, gas_pressure_size>::Zero(gas_pressure_size, gas_pressure_size);

        typename ShapeMatricesTypePressure::template MatrixType<gas_pressure_size, gas_pressure_size>
        Mlpc = ShapeMatricesTypePressure::template  MatrixType
        <gas_pressure_size, gas_pressure_size>::Zero(gas_pressure_size, gas_pressure_size);

        typename ShapeMatricesTypePressure::template MatrixType<gas_pressure_size, gas_pressure_size>
        Kgpg = ShapeMatricesTypePressure::template  MatrixType
        <gas_pressure_size, gas_pressure_size>::Zero(gas_pressure_size, gas_pressure_size);

        typename ShapeMatricesTypePressure::template MatrixType<gas_pressure_size, gas_pressure_size>
        Klpg = ShapeMatricesTypePressure::template  MatrixType
        <gas_pressure_size, gas_pressure_size>::Zero(gas_pressure_size, gas_pressure_size);

        typename ShapeMatricesTypePressure::template MatrixType<gas_pressure_size, gas_pressure_size>
        Klpc = ShapeMatricesTypePressure::template  MatrixType
        <gas_pressure_size, gas_pressure_size>::Zero(gas_pressure_size, gas_pressure_size);


        typename ShapeMatricesTypePressure::template MatrixType<gas_pressure_size, displacement_size>
        Mgus = ShapeMatricesTypePressure::template  MatrixType
        <gas_pressure_size, displacement_size>::Zero(gas_pressure_size, displacement_size);


        typename ShapeMatricesTypePressure::template MatrixType<gas_pressure_size, displacement_size>
        Mlus = ShapeMatricesTypePressure::template  MatrixType
        <gas_pressure_size, displacement_size>::Zero(gas_pressure_size, displacement_size);


        typename ShapeMatricesTypePressure::template MatrixType<displacement_size, gas_pressure_size>
        Kupg = ShapeMatricesTypePressure::template  MatrixType
        <displacement_size, gas_pressure_size>::Zero(displacement_size, gas_pressure_size);

        typename ShapeMatricesTypePressure::template MatrixType<displacement_size, gas_pressure_size>
        Kupc = ShapeMatricesTypePressure::template  MatrixType
        <displacement_size, gas_pressure_size>::Zero(displacement_size, gas_pressure_size);

        typename ShapeMatricesTypePressure::template VectorType<gas_pressure_size>
        Bg = ShapeMatricesTypePressure::template  VectorType
        <gas_pressure_size>::Zero(gas_pressure_size);

        typename ShapeMatricesTypePressure::template VectorType<gas_pressure_size>
        Bl = ShapeMatricesTypePressure::template  VectorType
        <gas_pressure_size>::Zero(gas_pressure_size);


        typename ShapeMatricesTypePressure::template VectorType<displacement_size>
        Bu = ShapeMatricesTypePressure::template  VectorType
        <displacement_size>::Zero(displacement_size);


        double const& dt = _process_data.dt;

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        unsigned const n_integration_points =
                _integration_method.getNumberOfPoints();

        _velocity_liquid.clear();
        auto velocity_matrix_liquid = MathLib::createZeroedMatrix<Eigen::Matrix<
                double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
                        _velocity_liquid, DisplacementDim, _ip_data.size());

        _velocity_gas.clear();
        auto velocity_matrix_gas = MathLib::createZeroedMatrix<Eigen::Matrix<
                double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
                        _velocity_gas, DisplacementDim, _ip_data.size());


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

            auto const permeability = _process_data.intrinsic_permeability(t, x_position)[0];


//            std::cout << "analytical assembly of Jacobian:\n\n";
//            std::cout << "displacement_index : " << displacement_index << "\n";
//            std::cout << "gas_pressure_index : " << gas_pressure_index << "\n";
//            std::cout << "cap_pressure_index : " << cap_pressure_index << "\n";
//            std::cout << "\n";
//            std::cout << "displacement_size : " << displacement_size << "\n";
//            std::cout << "gas_pressure_size : " << gas_pressure_size << "\n";
//            std::cout << "cap_pressure_size : " << cap_pressure_size << "\n";
//            std::cout << "\n";


            const double temperature = 293.15;
            const double molar_mass  = 0.02896;
            const double gas_constant = 8.3144621;

            auto const pc_int_pt = p_cap.dot(N_p);
            auto const pn_int_pt = p_GR.dot(N_p);

            double const rho_GR = molar_mass * pn_int_pt / gas_constant / temperature;
            double const rho_LR = 1000.;
            double phi = porosity;

            const double p_LR = pn_int_pt - pc_int_pt;
            _pressure_wet[ip] = p_LR;
            _density_gas[ip] = rho_GR;


            const double m_sw =-1.00623E-05;
            const double n_sw = 1.0;

            double const Sw = m_sw * pc_int_pt + n_sw;

            const double sl = Sw;
            const double sg = 1.0 - sl;
            const double dsldpc = m_sw;


            double const drhononwet_dpn = molar_mass / gas_constant / temperature;

            auto const rho = rho_sr * (1 - porosity) + porosity *
                    (Sw * rho_LR + (1. - Sw)* rho_GR);

              const double beta_p_GR = 1./rho_GR*drhononwet_dpn;
              const double beta_p_SR = 0.000001;
              const double beta_p_LR = 0.;

              double alpha_B = alpha;

              const auto m = identity2;

              const double k_rel_GR = 1.0;
              const double k_rel_LR = 1.0;

              double const mu_GR = 1.8e-5;
              double const mu_LR = 1.0e-3;

              const double lambda_GR = k_rel_GR / mu_GR;
              const double lambda_LR = k_rel_LR / mu_LR;

              velocity_matrix_liquid.col(ip).noalias() =
                      -K_over_mu * dNdx_p * (p_GR-p_cap) - K_over_mu * rho_LR * b;

              velocity_matrix_gas.col(ip).noalias() =
                      -K_over_mu * dNdx_p * p_GR - K_over_mu * rho_GR * b;

              Eigen::Matrix<double, DisplacementDim, DisplacementDim> perm_tensor;
              perm_tensor.setIdentity();

              auto permeability_tensor = permeability * perm_tensor;

              auto const Laplace = (dNdx_p.transpose() * permeability_tensor * dNdx_p * w).eval();

              Mgpg.noalias() += N_p.transpose()*sg*rho_GR*(phi*beta_p_GR+(alpha_B-phi)*beta_p_SR)*N_p*w;
              Mgpc.noalias() -= N_p.transpose()*rho_GR*((alpha_B-phi)*beta_p_SR*sl*sg+(phi+(alpha_B-phi)*beta_p_SR*sg*pc_int_pt)*dsldpc)*N_p*w;
              Mgus.noalias() += N_p.transpose()*rho_GR*alpha_B*sg*m.transpose()*B*w;

              Mlpg.noalias() += N_p.transpose()*rho_LR*(phi*sl*beta_p_LR+(alpha_B-phi)*beta_p_SR*sl)*N_p*w;
              Mlpc.noalias() += N_p.transpose()*rho_LR*((phi-(alpha_B-phi)*beta_p_SR*sl*pc_int_pt)*dsldpc-(phi*sl*beta_p_LR+(alpha_B-phi)*beta_p_SR*sl*sl))*N_p*w;
              Mlus.noalias() += N_p.transpose()*rho_LR*sl*alpha_B*m.transpose()*B*w;

              Kgpg.noalias() += rho_GR*lambda_GR*Laplace;
              Klpg.noalias() += rho_LR*lambda_LR*Laplace;
              Klpc.noalias() -= rho_LR*lambda_LR*Laplace;

              Kupg.noalias() += B.transpose()*alpha*m*N_p*w;
              Kupc.noalias() -= B.transpose()*alpha*m*sl*N_p*w;

              auto const gravity_operator = (dNdx_p.transpose() *
                                  permeability_tensor * b * w).eval();

              Bg.noalias() += rho_GR*rho_GR*lambda_GR*gravity_operator;
              Bl.noalias() += rho_LR*rho_LR*lambda_LR*gravity_operator;
              Bu.noalias() -= (B.transpose()*sigma_eff-N_u_op.transpose()*rho*b)*w;

//              std::cout.precision(20);
//            std::cout   <<  "    Mgpg   :\n" <<  Mgpg   <<   "\n\n";
//            std::cout   <<  "    Mgpc   :\n" <<  Mgpc   <<   "\n\n";
//            std::cout   <<  "    Mgus   :\n" <<  Mgus   <<   "\n\n";
//
//            std::cout   <<  "    Mlpg   :\n" <<  Mlpg   <<   "\n\n";
//            std::cout   <<  "    Mlpc   :\n" <<  Mlpc   <<   "\n\n";
//            std::cout   <<  "    Mlus   :\n" <<  Mlus   <<   "\n\n";
//
//            std::cout   <<  "    Kgpg   :\n" <<  Kgpg   <<   "\n\n";
//            std::cout   <<  "    Klpg   :\n" <<  Klpg   <<   "\n\n";
//            std::cout   <<  "    Klpc   :\n" <<  Klpc   <<   "\n\n";
//            std::cout   <<  "    Kupg   :\n" <<  Kupg   <<   "\n\n";
//            std::cout   <<  "    Kupc   :\n" <<  Kupc   <<   "\n\n";
//
//            std::cout   <<  "    Bg     :\n" <<  Bg     <<   "\n\n";
//            std::cout   <<  "    Bl     :\n" <<  Bl     <<   "\n\n";
//            std::cout   <<  "    Bu     :\n" <<  Bu     <<   "\n\n";
//
//            std::cout   <<  "    Sigma  :\n" <<  sigma_eff  <<   "\n\n";
//
//            std::cout   <<  "   ________________________________\n\n\n";
//
//            OGS_FATAL("Intended halt (analytical).");


//             displacement equation, displacement part


            eps.noalias() = B * u;

            auto C =
                _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u);

            local_Jac
                .template block<displacement_size, displacement_size>(
                    displacement_index, displacement_index)
                .noalias() += B.transpose() * C * B * w;


//            local_rhs.template segment<displacement_size>(displacement_index)
//                .noalias() -=
//                (B.transpose() * sigma_eff - N_u_op.transpose() * rho * b) * w;

            //
            // displacement equation, pressure part
            //
//            Kup.noalias() += B.transpose() * alpha * identity2 * N_p * w;

            //
            // pressure equation, pressure part.
            //
//            laplace_p.noalias() += dNdx_p.transpose() * K_over_mu * dNdx_p * w;
//
//            storage_p.noalias() += N_p.transpose() * S * N_p * w;
//
//            local_rhs.template segment<gas_pressure_size>(gas_pressure_index)
//                .noalias() += dNdx_p.transpose() * rho_fr * K_over_mu * b * w;

            //
            // pressure equation, displacement part.
            //
            // Reusing Kup.transpose().
        }

//                    OGS_FATAL("Intended halt (analytical).");

        // gas phase equation, gas pressure part:
        local_Jac
            .template block<gas_pressure_size, gas_pressure_size>(
                gas_pressure_index, gas_pressure_index)
            .noalias() = Mgpg/dt+Kgpg;

        // gas phase equation, capillary pressure part:
        local_Jac
            .template block<gas_pressure_size, cap_pressure_size>(
                gas_pressure_index, cap_pressure_index)
            .noalias() = Mgpc/dt;

        // gas phase equation, displacement part:
        local_Jac
            .template block<gas_pressure_size, displacement_size>(
                gas_pressure_index, displacement_index)
            .noalias() = Mgus/dt;


        // liquid phase equation, gas pressure part:
        local_Jac
            .template block<cap_pressure_size, gas_pressure_size>(
                cap_pressure_index, gas_pressure_index)
            .noalias() = Mlpg/dt+Klpg;

        // liquid phase equation, capillary pressure part:
        local_Jac
            .template block<cap_pressure_size, cap_pressure_size>(
                cap_pressure_index, cap_pressure_index)
            .noalias() = Mlpc/dt+Klpc;

        // liquid phase equation, displacement part:
        local_Jac
            .template block<cap_pressure_size, displacement_size>(
                cap_pressure_index, displacement_index)
            .noalias() = Mlus/dt;


        // displacement equation, gas pressure part:
        local_Jac
            .template block<displacement_size, gas_pressure_size>(
                displacement_index, gas_pressure_index)
            .noalias() = Kupg;

        // displacement equation, capillary pressure part:
        local_Jac
            .template block<displacement_size, cap_pressure_size>(
                displacement_index, cap_pressure_index)
            .noalias() = Kupc;

        local_Jac.template block<temperature_size,
        temperature_size>(temperature_index, temperature_index).setIdentity();

        // residuum, gas phase equation

        local_rhs.template segment<gas_pressure_size>(gas_pressure_index)
                 .noalias() = Bg - Mgpg*p_GR_dot - Mgpc*p_cap_dot - Mgus*u_dot
                 - Kgpg*p_GR;

        // residuum, liquid phase equation
        local_rhs.template segment<cap_pressure_size>(cap_pressure_index)
                 .noalias() = Bl - Mlpg*p_GR_dot - Mlpc*p_cap_dot - Mlus*u_dot
                 - Klpg*p_GR - Klpc*p_cap;

        // rediduum, displacement equation
        local_rhs.template segment<displacement_size>(displacement_index)
                 .noalias() = Bu - Kupg*p_GR - Kupc*p_cap;

//        auto const rows = local_Jac.rows();
//        auto const cols = local_Jac.cols();
//        std::cout.precision(20);
//        std::cout << "Jacobi-Matrix:\n";
//
//        for (int r = 0; r < rows; r++)
//        {
//            for (int c = 0; c < cols; c++)
//            {
//                std::cout << local_Jac(r,c) << " ";
//            }
//            std::cout << "\n";
//        }
////
////        std::cout << "Residuum:\n";
////
////        for (int r = 0; r < rows; r++)
////        {
////            std::cout << local_rhs(r) << " ";
////        }
////
////        std::cout << "\n";
////
////        std::cout << "Analytical J.\n";
//        OGS_FATAL("Intended halt.");



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

    std::vector<double> const& getIntPtDensityGas(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_density_gas.empty());
        return _density_gas;
    }

    std::vector<double> const& getIntPtDensityLiquid(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_density_liquid.empty());
        return _density_liquid;
    }

    std::vector<double> const& getIntPtDarcyVelocityLiquid(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_velocity_liquid.empty());
        return _velocity_liquid;
    }

    std::vector<double> const& getIntPtPressureGas(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_pressure_gas_linear.empty());
        return _pressure_gas_linear;
    }

    std::vector<double> const& getIntPtPressureCap(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_pressure_cap_linear.empty());
        return _pressure_cap_linear;
    }

    std::vector<double> const& getIntPtDarcyVelocityGas(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_velocity_gas.empty());
        return _velocity_gas;
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
    std::vector<double> _density_gas;
    std::vector<double> _density_liquid;
    std::vector<double> _pressure_gas_linear;
    std::vector<double> _pressure_cap_linear;
    std::vector<double> _velocity_gas;
    std::vector<double> _velocity_liquid;



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
