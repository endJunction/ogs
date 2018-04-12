/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Implementation of ThermoplasticBDT plastic model.
 * see Parisio et al. 2018 "Early signs of volcanic eruption from the
 * brittle-ductile transition of rocks" for more details.
 *
 * Refer to "Single-surface benchmark of OpenGeoSys documentation
 * (https://docs.opengeosys.org/docs/benchmarks/small-deformations/mechanics-plasticity-single-surface)"
 * for more details for the tests.
 */

#pragma once

#ifndef NDEBUG
#include <ostream>
#endif

#include "BaseLib/Error.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/NewtonRaphson.h"
#include "ProcessLib/Parameter/Parameter.h"

#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{
namespace ThermoplasticBDT
{
/// material parameters in relation to Parisio et al. plastic model
struct MaterialPropertiesParameters
{
    using P = ProcessLib::Parameter<double>;

    MaterialPropertiesParameters(P const& G_, P const& K_, P const& fc_,
                                 P const& m_, P const& qp0_, P const& alpha_,
                                 P const& n_, P const& T0_)
        : G(G_),
          K(K_),
          fc(fc_),
          m(m_),
          qp0(qp0_),
          alpha(alpha_),
          n(n_),
          T0(T0_),
    {
    }
    // basic material parameters
    P const& G;  ///< shear modulus
    P const& K;  ///< bulk modulus

    P const& fc;
    P const& m;
    P const& qp0;
    P const& alpha;
    P const& n;
    P const& T0;
};

struct DamagePropertiesParameters
{
    using P = ProcessLib::Parameter<double>;
    P const& alpha_d;
    P const& beta_d;
    P const& h_d;
    double const m_d;
};

/// Evaluated MaterialPropertiesParameters container, see its documentation for
/// details.
struct MaterialProperties final
{
    MaterialProperties(double const t, ProcessLib::SpatialPosition const& x,
                       MaterialPropertiesParameters const& mp)
        : G(mp.G(t, x)[0]),
          K(mp.K(t, x)[0]),
          fc(mp.fc(t, x)[0]),
          m(mp.m(t, x)[0]),
          qp0(mp.qp0(t, x)[0]),
          alpha(mp.alpha(t, x)[0]),
          n(mp.n(t, x)[0]),
          T0(mp.T0(t, x)[0]),
    {
    }
    double const G;
    double const K;

    double const fc;
    double const m;
    double const qp0;
    double const alpha;
    double const n;
    double const T0;
};

/// Evaluated DamagePropertiesParameters container, see its documentation for
/// details.
struct DamageProperties
{
    DamageProperties(double const t,
                     ProcessLib::SpatialPosition const& x,
                     DamagePropertiesParameters const& dp)
        : alpha_d(dp.alpha_d(t, x)[0]),
          beta_d(dp.beta_d(t, x)[0]),
          h_d(dp.h_d(t, x)[0]),
          m_d(dp.m_d)
    {
    }
    double const alpha_d;
    double const beta_d;
    double const h_d;
    double const m_d;
};

template <typename KelvinVector>
struct PlasticStrain final
{
    PlasticStrain() : D(KelvinVector::Zero()) {}
    PlasticStrain(KelvinVector eps_p_D_, double const eps_p_V_,
                  double const eps_p_eff_)
        : D(std::move(eps_p_D_)), V(eps_p_V_), eff(eps_p_eff_){};

    KelvinVector D;  ///< deviatoric plastic strain
    double V = 0;    ///< volumetric strain
    double eff = 0;  ///< effective plastic strain

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

class Damage final
{
public:
    Damage() = default;
    Damage(double const kappa_d, double const value)
        : _kappa_d(kappa_d), _value(value)
    {
    }

    double kappa_d() const { return _kappa_d; }
    double value() const { return _value; }

private:
    double _kappa_d = 0;  ///< damage driving variable
    double _value = 0;    ///< isotropic damage variable
};

template <int DisplacementDim>
struct StateVariables
    : public MechanicsBase<DisplacementDim>::MaterialStateVariables
{
    StateVariables& operator=(StateVariables const&) = default;
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables& operator=(
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            state) noexcept override
    {
        assert(dynamic_cast<StateVariables const*>(&state) != nullptr);
        return operator=(static_cast<StateVariables const&>(state));
    }

    void setInitialConditions()
    {
        eps_p = eps_p_prev;
        damage = damage_prev;
    }

    void pushBackState() override
    {
        eps_p_prev = eps_p;
        damage_prev = damage;
    }

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

    PlasticStrain<KelvinVector> eps_p;  ///< plastic part of the state.
    Damage damage;                      ///< damage part of the state.

    // Initial values from previous timestep
    PlasticStrain<KelvinVector> eps_p_prev;  ///< \copydoc eps_p
    Damage damage_prev;                      ///< \copydoc damage

#ifndef NDEBUG
    friend std::ostream& operator<<(std::ostream& os,
                                    StateVariables<DisplacementDim> const& m)
    {
        os << "State:\n"
           << "eps_p_D: " << m.eps_p.D << "\n"
           << "eps_p_eff: " << m.eps_p.eff << "\n"
           << "kappa_d: " << m.damage.kappa_d() << "\n"
           << "damage: " << m.damage.value() << "\n"
           << "eps_p_D_prev: " << m.eps_p_prev.D << "\n"
           << "eps_p_eff_prev: " << m.eps_p_prev.eff << "\n"
           << "kappa_d_prev: " << m.damage_prev.kappa_d() << "\n"
           << "damage_prev: " << m.damage_prev.value() << "\n"
           << "lambda: " << m.lambda << "\n";
        return os;
    }
#endif  // NDEBUG

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <int DisplacementDim>
class SolidThermoPlasticBDT final : public MechanicsBase<DisplacementDim>
{
public:
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    static int const JacobianResidualSize =
        2 * KelvinVectorSize + 3;  // 2 is the number of components in the
                                   // jacobian/residual, not the space
                                   // dimension. And 3 is for additional
                                   // variables.
    using ResidualVectorType = Eigen::Matrix<double, JacobianResidualSize, 1>;
    using JacobianMatrix = Eigen::Matrix<double, JacobianResidualSize,
                                         JacobianResidualSize, Eigen::RowMajor>;

public:
    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() override
    {
        return std::make_unique<StateVariables<DisplacementDim>>();
    }

public:
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

public:
    explicit SolidThermoPlasticBDT(
        NumLib::NewtonRaphsonSolverParameters nonlinear_solver_parameters,
        MaterialPropertiesParameters material_properties,
        std::unique_ptr<DamagePropertiesParameters>&& damage_properties)
        : _nonlinear_solver_parameters(std::move(nonlinear_solver_parameters)),
          _mp(std::move(material_properties)),
          _damage_properties(std::move(damage_properties))
    {
    }

    double computeFreeEnergyDensity(
        double const /*t*/,
        ProcessLib::SpatialPosition const& /*x*/,
        double const /*dt*/,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override
    {
        assert(dynamic_cast<StateVariables<DisplacementDim> const*>(
                   &material_state_variables) != nullptr);

        auto const& eps_p = static_cast<StateVariables<DisplacementDim> const&>(
                                material_state_variables)
                                .eps_p;
        using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
        auto const& identity2 = Invariants::identity2;
        return (eps - eps_p.D - eps_p.V / 3 * identity2).dot(sigma) / 2;
    }

    boost::optional<std::tuple<KelvinVector,
                               std::unique_ptr<typename MechanicsBase<
                                   DisplacementDim>::MaterialStateVariables>,
                               KelvinMatrix>>
    integrateStress(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps_prev,
        KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) override;

    std::vector<typename MechanicsBase<DisplacementDim>::InternalVariable>
    getInternalVariables() const override;

    MaterialProperties evaluatedMaterialProperties(
        double const t, ProcessLib::SpatialPosition const& x) const
    {
        return MaterialProperties(t, x, _mp);
    }

    DamageProperties evaluatedDamageProperties(
        double const t, ProcessLib::SpatialPosition const& x) const
    {
        return DamageProperties(t, x, *_damage_properties);
    }

private:
    NumLib::NewtonRaphsonSolverParameters const _nonlinear_solver_parameters;

    MaterialPropertiesParameters _mp;
    std::unique_ptr<DamagePropertiesParameters> _damage_properties;
};

extern template class SolidThermoPlasticBDT<2>;
extern template class SolidThermoPlasticBDT<3>;
}  // namespace ThermoplasticBDT
}  // namespace Solids
}  // namespace MaterialLib
