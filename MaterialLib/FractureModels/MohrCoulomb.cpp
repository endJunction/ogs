/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MohrCoulomb.h"
#include "LogPenalty.h"

#include "BaseLib/Error.h"
#include "MathLib/MathTools.h"

namespace MaterialLib
{
namespace Fracture
{
namespace MohrCoulomb
{

namespace
{

struct MaterialPropertyValues
{
    double Kn = 0.0;
    double Ks = 0.0;
    double phi = 0.0; // friction angle
    double psi = 0.0; // dilation angle
    double c = 0.0;

    template <typename MaterialProperties>
    MaterialPropertyValues(
            MaterialProperties const& mp,
            double const t,
            ProcessLib::SpatialPosition const& x)
    {
        Kn = mp.normal_stiffness(t,x)[0];
        Ks = mp.shear_stiffness(t,x)[0];
        phi = MathLib::to_radians(mp.friction_angle(t,x)[0]);
        psi = MathLib::to_radians(mp.dilatancy_angle(t,x)[0]);
        c = mp.cohesion(t,x)[0];
    }
};

} // no namespace

template <int DisplacementDim>
void MohrCoulomb<DisplacementDim>::computeConstitutiveRelation(
    double const t,
    ProcessLib::SpatialPosition const& x,
    double const aperture0,
    Eigen::Ref<Eigen::VectorXd const>
        sigma0,
    Eigen::Ref<Eigen::VectorXd const>
        /*w_prev*/,
    Eigen::Ref<Eigen::VectorXd const>
        w,
    Eigen::Ref<Eigen::VectorXd const>
        /*sigma_prev*/,
    Eigen::Ref<Eigen::VectorXd>
        sigma,
    Eigen::Ref<Eigen::MatrixXd>
        Kep,
    typename FractureModelBase<DisplacementDim>::MaterialStateVariables&
        material_state_variables)
{
    assert(dynamic_cast<StateVariables<DisplacementDim> const*>(
               &material_state_variables) != nullptr);

    StateVariables<DisplacementDim>& state =
        static_cast<StateVariables<DisplacementDim>&>(material_state_variables);

    MaterialPropertyValues const mat(_mp, t, x);

    const int index_ns = DisplacementDim - 1;
    double const aperture = w[index_ns] + aperture0;

    Eigen::MatrixXd Ke;
    {  // Elastic tangent stiffness
        Ke = Eigen::MatrixXd::Zero(DisplacementDim, DisplacementDim);
        for (int i = 0; i < index_ns; i++)
            Ke(i, i) = mat.Ks;

        Ke(index_ns, index_ns) = mat.Kn;
    }

    // Total plastic aperture compression
    // NOTE: Initial condition sigma0 seems to be associated with an initial
    // condition of the w0 = 0. Therefore the initial state is not associated
    // with a plastic aperture change.
    {  // Exact elastic predictor
        sigma.noalias() = Ke * (w - state.w_p_prev);

        sigma.coeffRef(index_ns) *=
            logPenalty(aperture0, aperture, _penalty_aperture_cutoff);
    }

    sigma.noalias() += sigma0;

    // correction for an opening fracture
    if (_tension_cutoff && sigma[DisplacementDim - 1] >= 0)
    {
        Kep.setZero();
        sigma.setZero();
        material_state_variables.setTensileStress(true);
        return;

        // TODO (nagel); Update w_p for fracture opening and closing.
    }

    auto yieldFunction = [&mat](Eigen::VectorXd const& s) {
        double const sigma_n = s[DisplacementDim - 1];
        Eigen::VectorXd const sigma_s = s.head(DisplacementDim - 1);
        double const mag_tau = sigma_s.norm();  // magnitude
        return mag_tau + sigma_n * std::tan(mat.phi) - mat.c;
    };

    {  // Exit if still in elastic range by checking the shear yield function.
        double const Fs = yieldFunction(sigma);
        material_state_variables.setShearYieldFunctionValue(Fs);
        if (Fs < .0)
        {
            Kep = Ke;
            Kep(index_ns, index_ns) *= logPenaltyDerivative(
                aperture0, aperture, _penalty_aperture_cutoff);
            return;
        }
    }

    auto yieldFunction_derivative = [&mat](Eigen::VectorXd const& s) {
        Eigen::VectorXd dFs_dS(DisplacementDim);
        Eigen::VectorXd const sigma_s = s.head(DisplacementDim - 1);
        dFs_dS.head(DisplacementDim - 1).noalias() = sigma_s.normalized();
        dFs_dS[DisplacementDim - 1] = std::tan(mat.phi);
        return dFs_dS;
    };

    // plastic potential function: Qs = |tau| + Sn * tan da
    auto plasticPotential_derivative = [&mat](Eigen::VectorXd const& s) {
        Eigen::VectorXd dQs_dS(DisplacementDim);
        Eigen::VectorXd const sigma_s = s.head(DisplacementDim - 1);
        dQs_dS.head(DisplacementDim - 1).noalias() = sigma_s.normalized();
        dQs_dS[DisplacementDim - 1] = std::tan(mat.psi);
        return dQs_dS;
    };

    {  // Newton

        Eigen::FullPivLU<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>>
            linear_solver;
        using ResidualVectorType = Eigen::Matrix<double, 1, 1, Eigen::RowMajor>;
        using JacobianMatrix = Eigen::Matrix<double, 1, 1, Eigen::RowMajor>;

        JacobianMatrix jacobian;
        ResidualVectorType solution;
        solution << 0;

        auto const update_residual = [&](ResidualVectorType& residual) {
            residual[0] = yieldFunction(sigma);
        };

        auto const update_jacobian = [&](JacobianMatrix& jacobian) {
            jacobian(0, 0) = -yieldFunction_derivative(sigma).transpose() * Ke *
                             plasticPotential_derivative(sigma);
        };

        auto const update_solution = [&](ResidualVectorType const& increment) {
            solution += increment;
            /*DBUG("analytical = %g",
                 Fs / (mat.Ks + mat.Kn * std::tan(mat.psi) * std::tan(mat.phi)))
                 */
            state.w_p = state.w_p_prev +
                        solution[0] * plasticPotential_derivative(sigma);
            sigma.noalias() = sigma0 + (Ke * (w - state.w_p)) *
                                           logPenalty(aperture0, aperture,
                                                      _penalty_aperture_cutoff);
        };

        auto newton_solver =
            NumLib::NewtonRaphson<decltype(linear_solver), JacobianMatrix,
                                  decltype(update_jacobian), ResidualVectorType,
                                  decltype(update_residual),
                                  decltype(update_solution)>(
                linear_solver, update_jacobian, update_residual,
                update_solution, _nonlinear_solver_parameters);

        auto const success_iterations = newton_solver.solve(jacobian);

        if (!success_iterations)
            OGS_FATAL("MohrCoulomb nonlinear solver didn't converge.");

        // Solution containing lambda is not needed; w_p and sigma already
        // up to date.
    }

    {  // Update material state shear yield function value.
        double const Fs = yieldFunction(sigma);
        material_state_variables.setShearYieldFunctionValue(Fs);
    }

    Ke(index_ns, index_ns) *=
        logPenaltyDerivative(aperture0, aperture, _penalty_aperture_cutoff);
    Eigen::RowVectorXd const A = yieldFunction_derivative(sigma).transpose() *
                                 Ke /
                                 (yieldFunction_derivative(sigma).transpose() *
                                  Ke * plasticPotential_derivative(sigma));
    Kep = Ke - Ke * plasticPotential_derivative(sigma) * A;
}

template class MohrCoulomb<2>;
template class MohrCoulomb<3>;

}  // namespace MohrCoulomb
}  // namespace Fracture
}  // namespace MaterialLib
