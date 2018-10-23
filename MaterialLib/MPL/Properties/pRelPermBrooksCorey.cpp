/**
 * \author Norbert Grunwald
 * \date   02.07.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialLib/MPL/Properties/pRelPermBrooksCorey.h"
#include "MaterialLib/MPL/mpMedium.h"
#include "pUniversalConstants.h"

#include <iostream>
#include <algorithm>
#include <cmath>

namespace MaterialPropertyLib
{
BrooksCoreyRelPerm::BrooksCoreyRelPerm(Medium* m) : _medium(m){};
/// These constructors throw errors, since the property is not
/// implemented on phase or component scales.
BrooksCoreyRelPerm::BrooksCoreyRelPerm(Phase*) : _medium(0)
{
    notImplemented("BrooksCoreyRelPerm", "Phase");
}
BrooksCoreyRelPerm::BrooksCoreyRelPerm(Component*) : _medium(0)
{
    notImplemented("BrooksCoreyRelPerm", "Component");
}

/**
 */
PropertyDataType  BrooksCoreyRelPerm::value(VariableArray const& v)
{

    const double p_cap = getScalar(v[MaterialPropertyLib::capillary_pressure]);
    const double p_GR = getScalar(v[MaterialPropertyLib::phase_pressure]);

    const double s_L_res =
            getScalar(_medium->property(residual_liquid_saturation));
    const double s_L_max =
            1 - getScalar(_medium->property(residual_gas_saturation));
    const double s_L =
            getScalar(_medium->property(saturation), v);

    const double lambda =
            getScalar(_medium->property(brooks_corey_exponent));

    const double s_eff = (s_L - s_L_res)/(s_L_max - s_L_res);

    const double k_rel_LR = std::pow(s_eff, (2.+3.*lambda)/lambda);
    const double k_rel_GR = (1.-s_eff)*(1. - s_eff) * (1. - std::pow(s_eff, (2.+lambda)/lambda));

    const Pair value = {k_rel_LR, k_rel_GR};

    return value;
}

}  // MaterialPropertyLib
