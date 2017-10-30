
/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

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

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
