/**
 * \author Norbert Grunwald
 * \date   27.06.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialLib/MPL/Properties/pSaturationBrooksCorey.h"
#include "MaterialLib/MPL/mpMedium.h"
#include "pUniversalConstants.h"

#include <algorithm>
#include <cmath>


namespace MaterialPropertyLib
{
BrooksCoreySaturation::BrooksCoreySaturation(Medium* m) : _medium(m){};
/// These constructors throw errors, since the property is not
/// implemented on phase or component scales.
BrooksCoreySaturation::BrooksCoreySaturation(Phase*) : _medium(0)
{
    notImplemented("BrooksCoreySaturation", "Phase");
}
BrooksCoreySaturation::BrooksCoreySaturation(Component*) : _medium(0)
{
    notImplemented("BrooksCoreySaturation", "Component");
}

/**
 */
PropertyDataType  BrooksCoreySaturation::value(VariableArray const& v)
{

    const double p_cap = getScalar(
            v[MaterialPropertyLib::Variables::capillary_pressure]);
    const double p = getScalar(
            v[MaterialPropertyLib::Variables::phase_pressure]);

    const double s_L_res =
            getScalar(_medium->property(
                    PropertyEnum::residual_liquid_saturation));
    const double s_L_max =
            1 - getScalar(_medium->property(
                    PropertyEnum::residual_gas_saturation));

    const double p_b =
            getScalar(_medium->property(PropertyEnum::entry_pressure));
    const double lambda =
            getScalar(_medium->property(PropertyEnum::brooks_corey_exponent));

    const double s_eff = std::pow(p_b/p_cap, lambda);
    return s_eff*(s_L_max - s_L_res) + s_L_res;
}

PropertyDataType BrooksCoreySaturation::dvalue(VariableArray const& v, Variables const pv)
{

    assert((pv == Variables::capillary_pressure) && "BrooksCoreySaturation::dvalue is implemented for "
            " derivatives with respect to capillary pressure only.");

    const double p_c = getScalar(v[pv]);
    const double s_L = getScalar(_medium->property(PropertyEnum::saturation));
    const double s_L_res =
            getScalar(_medium->property(PropertyEnum::residual_liquid_saturation));

    const double lambda =
            getScalar(_medium->property(PropertyEnum::brooks_corey_exponent));

    return -lambda / p_c * ( s_L - s_L_res);
}
}  // MaterialPropertyLib
