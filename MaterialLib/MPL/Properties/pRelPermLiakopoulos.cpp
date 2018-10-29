/**
 * \author Norbert Grunwald
 * \date   Oct 19, 2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <MaterialLib/MPL/Properties/pRelPermLiakopoulos.h>
#include "MaterialLib/MPL/mpMedium.h"
#include "pUniversalConstants.h"

#include <iostream>
#include <algorithm>
#include <cmath>

namespace MaterialPropertyLib
{
RelPermLiakopoulos::RelPermLiakopoulos(Medium* m) : _medium(m){};
/// These constructors throw errors, since the property is not
/// implemented on phase or component scales.
RelPermLiakopoulos::RelPermLiakopoulos(Phase*) : _medium(0)
{
    notImplemented("RelPermLiakopoulos", "Phase");
}
RelPermLiakopoulos::RelPermLiakopoulos(Component*) : _medium(0)
{
    notImplemented("RelPermLiakopoulos", "Component");
}

/**
 */
PropertyDataType  RelPermLiakopoulos::value(VariableArray const& v)
{

    const double s_L = getScalar(v[Variables::liquid_saturation]);

    const double s_eff = (s_L - s_L_r)/(1. - s_L_r);

    const double k_rel_LR = std::max(0., 1. - kA * std::pow((1. - s_L), kB));
    const double k_rel_GR = std::max(k_rel_G_min, (1.-s_eff)*(1. - s_eff) *
            (1. - std::pow(s_eff, (2.+lambda)/lambda)));

    const Pair value = {k_rel_LR, k_rel_GR};
    return value;
}

PropertyDataType  RelPermLiakopoulos::dvalue(VariableArray const& v,
        Variables const pv)
{

    assert((pv == Variables::liquid_saturation) &&
                "RelPermLiakopoulos::dvalue is implemented for "
                " derivatives with respect to liquid saturation only.");

    const double s_L = getScalar(v[Variables::liquid_saturation]);
    const double s_eff = (s_L - s_L_r)/(1. - s_L_r);

    const double k_rel_LR = std::max(0., 1. - kA * std::pow((1. - s_L), kB));
    const double k_rel_GR = std::max(k_rel_G_min, (1.-s_eff)*(1. - s_eff) *
            (1. - std::pow(s_eff, (2.+lambda)/lambda)));

    double dkrelLdsL;
    double dkrelGdsL;

    if (k_rel_LR <= 0.)
    {
        dkrelLdsL = 0.;
    }
    else
    {
        dkrelLdsL = kA * kB * std::pow((1. - s_L), kB - 1.);
    }
    if (k_rel_GR <= k_rel_G_min)
    {
        dkrelGdsL = 0.;
    }
    else
    {
        dkrelGdsL = 1. / ( lambda * (1. - s_L_r)) * (s_eff - 1.) *
                        ((2 + lambda - 2*s_eff - 3*lambda*s_eff) *
                        std::pow(s_eff, (2./lambda)) + 2*lambda) ;
    }

    const Pair value = {dkrelLdsL, dkrelGdsL};
    return value;
}

}  // MaterialPropertyLib



