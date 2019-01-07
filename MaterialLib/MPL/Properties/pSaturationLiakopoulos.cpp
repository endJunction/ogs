/**
 * \author Norbert Grunwald
 * \date   Oct 18, 2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <MaterialLib/MPL/Properties/pSaturationLiakopoulos.h>
#include "MaterialLib/MPL/mpMedium.h"
#include "pUniversalConstants.h"

#include <algorithm>
#include <cmath>

namespace MaterialPropertyLib
{
SaturationLiakopoulos::SaturationLiakopoulos(Medium* m) : _medium(m){};
/// These constructors throw errors, since the property is not
/// implemented on phase or component scales.
SaturationLiakopoulos::SaturationLiakopoulos(Phase*) : _medium(0)
{
    notImplemented("SaturationLiakopoulos", "Phase");
}
SaturationLiakopoulos::SaturationLiakopoulos(Component*) : _medium(0)
{
    notImplemented("SaturationLiakopoulos", "Component");
}

/**
 */
PropertyDataType SaturationLiakopoulos::value(VariableArray const& v)
{
    const double p_cap =
        getScalar(v[MaterialPropertyLib::Variables::capillary_pressure]);

    return std::max(s_L_r, 1. - s_a * std::pow(std::max(0., p_cap), s_b));
}

PropertyDataType SaturationLiakopoulos::dvalue(VariableArray const& v,
                                               Variables const pv)
{
    assert((pv == Variables::capillary_pressure) &&
           "SaturationLiakopoulos::dvalue is implemented for "
           " derivatives with respect to capillary pressure only.");

    const double p_cap = getScalar(v[pv]);
    const double s_L =
        getScalar(_medium->property(PropertyEnum::saturation), v);

    if (s_L <= s_L_r)
    {
        return 0.0;
    }
    else
    {
        return -s_a * s_b * std::pow(std::max(0., p_cap), s_b - 1.);
    }
}

PropertyDataType SaturationLiakopoulos::ddvalue(VariableArray const& v,
                                                Variables const pv1,
                                                Variables const pv2)
{
    assert((pv1 == Variables::capillary_pressure) &&
           (pv2 == Variables::capillary_pressure) &&
           "SaturationLiakopoulos::ddvalue is implemented for 2nd "
           "derivatives with respect to capillary pressure only.");

    const double p_cap = getScalar(v[pv1]);
    const double s_L =
        getScalar(_medium->property(PropertyEnum::saturation), v);

    if (s_L <= s_L_r)
    {
        return 0.0;
    }
    else
    {
        return -(s_b - 1.) * s_a * s_b *
               std::pow(std::max(0., p_cap), s_b - 2.);
    }
}
}  // namespace MaterialPropertyLib
