/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigenvalues>
#include <boost/math/special_functions/pow.hpp>
#include <logog/include/logog.hpp>

#include "MaterialLib/SolidModels/Ehlers.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
/// Computes the damage internal material variable explicitly based on the
/// results obtained from the local stress return algorithm.
inline double calculateDamage(double const kappa_d, double const alpha_d,
                              double const beta_d)
{
    double const damage = (1 - beta_d) * (1 - std::exp(-kappa_d / alpha_d));

    if (damage < 0. || damage > 1.)
        ERR("Damage value %g outside of [0,1] interval.", damage);

    return damage;
}

/// \Returns a new kappa_d.
template <int DisplacementDim, typename KelvinVectorType>
double calculateDamageKappaD(
    double const eps_p_eff_diff,
    KelvinVectorType const& sigma,
    double const kappa_d_prev,
    MaterialLib::Solids::Ehlers::DamageProperties const& dp,
    MaterialLib::Solids::Ehlers::MaterialProperties const& mp)
{
    // Default case of the rate problem. Updated below if volumetric plastic
    // strain rate is positive (dilatancy).

    Eigen::Matrix<double, DisplacementDim, DisplacementDim> stress_mat =
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();

    for (int i = 0; i < DisplacementDim; ++i)
    {
        for (int j = 0; j < DisplacementDim; ++j)
        {
            if (i == j)
            {
                stress_mat(i, j) = sigma(i);
            }
            else
            {
                stress_mat(i, j) = sigma(i + j + 2);
            }
        }
    }

    Eigen::EigenSolver<decltype(stress_mat)> eigen_solver(stress_mat);
    auto const& principal_stress = eigen_solver.eigenvalues();
    // building kappa_d (damage driving variable)
    double prod_stress = 0.;
    for (int i = 0; i < DisplacementDim; ++i)
    {
        double const real_eigen_value = real(principal_stress(i, 0));
        prod_stress = prod_stress + boost::math::pow<2>(real_eigen_value);
    }

    // Brittleness decrease with confinement for the nonlinear flow rule.
    // ATTENTION: For linear flow rule -> constant brittleness.
    double const f_t =
        std::sqrt(3.0) * mp.kappa / (1 + std::sqrt(3.0) * mp.beta);
    double const r_s = std::sqrt(prod_stress) / f_t;

    double x_s = 0;
    if (r_s < 1)
    {
        x_s = 1;
    }
    else if (r_s >= 1 && r_s <= 2)
    {
        x_s = 1 + dp.h_d * (r_s - 1) * (r_s - 1);
    }
    else
    {
        x_s = 1 - 3 * dp.h_d + 4 * dp.h_d * std::sqrt(r_s - 1);
    }

    return kappa_d_prev + eps_p_eff_diff / x_s;
}
}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
