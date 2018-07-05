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

#include <MaterialLib/MPL/Properties/pBrooksCoreySaturation.h>
#include "../mpMedium.h"
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
    if (isUpdated())
        return _value;

    const double p_cap = getScalar(v[MaterialPropertyLib::PrimaryVariables::p_cap]);
    const double p_GR = getScalar(v[MaterialPropertyLib::PrimaryVariables::p_GR]);

    const double s_L_res =
            getScalar(_medium->property(PropertyEnum::residual_liquid_saturation));
    const double s_L_max =
            1 - getScalar(_medium->property(PropertyEnum::residual_gas_saturation));

    const double p_b =
            getScalar(_medium->property(PropertyEnum::entry_pressure));
    const double lambda =
            getScalar(_medium->property(PropertyEnum::brooks_corey_exponent));

    const double s_eff = std::pow(p_b/p_cap, lambda);

    _value= s_eff*(s_L_max - s_L_res) + s_L_res;

    double s_L = s_eff*(s_L_max - s_L_res) + s_L_res;

    isUpdated(true);

    return _value;
}

PropertyDataType BrooksCoreySaturation::dvalue(VariableArray const& v, PrimaryVariables const pv)
{

    assert((pv == PrimaryVariables::p_cap) && "BrooksCoreySaturation::dvalue is implemented for "
            " derivatives with respect to capillary pressure only.");

    const double p_c = getScalar(v[pv]);
    const double s_L = getScalar(_medium->property(PropertyEnum::saturation));
    const double s_L_res =
            getScalar(_medium->property(PropertyEnum::residual_liquid_saturation));

    const double lambda =
            getScalar(_medium->property(PropertyEnum::brooks_corey_exponent));

    _dvalue = -lambda / p_c * ( s_L - s_L_res);
    return _dvalue;

}
}  // MaterialPropertyLib
