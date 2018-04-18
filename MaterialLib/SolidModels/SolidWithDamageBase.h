/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/KelvinVector.h"

#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{

template <int DisplacementDim>
struct SolidWithDamageBase : public MechanicsBase<DisplacementDim>
{
    using KelvinVectorType =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

    virtual double getOvernonlocalGammaFactor(
        double const t,
        ProcessLib::SpatialPosition const& x_position) const = 0;

    /// Computes the damage internal material variable explicitly based on the
    /// results obtained from the local stress return algorithm.
    virtual double calculateDamage(
        double const t, ProcessLib::SpatialPosition const& x_position,
        double const kappa_d) const = 0;

    virtual double calculateDamageKappaD(
        double const t, ProcessLib::SpatialPosition const& x_position,
        double const eps_p_eff_diff, KelvinVectorType const& sigma,
        double const kappa_d_prev) const = 0;
};

}  // namespace Solids
}  // namespace MaterialLib
