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

#include "pBrooksCorey.h"
#include <algorithm>
#include <cmath>
#include "../../MathLib/poly34.h"
#include "../mpMedium.h"
#include "pUniversalConstants.h"

#include <iostream>


namespace MaterialPropertyLib
{
/// This constructor throws an error, since the property is not
/// implemented on the medium scale.
BrooksCorey::BrooksCorey(Medium* m) : _medium(m){};
BrooksCorey::BrooksCorey(Phase*) : _medium(0)
{
    notImplemented("BrooksCorey", "Phase");
}
BrooksCorey::BrooksCorey(Component*) : _medium(0)
{
    notImplemented("BrooksCorey", "Component");
}

/**
 */
PropertyDataType  BrooksCorey::value(VariableArray const& v)
{
    if (isUpdated())
        return _value;

    const double p_cap = getScalar(v[MaterialPropertyLib::p_cap]);
    const double p_GR = getScalar(v[MaterialPropertyLib::p_GR]);

    std::cout   <<  " In BrooksCorey: \n";
    std::cout   <<  "     p_cap " << p_cap << "\n";
    std::cout   <<  "     p_GR " << p_GR << "\n";

    const double s_L_res =
            getScalar(_medium->property(residual_liquid_saturation));
    const double s_L_max =
            1 - getScalar(_medium->property(residual_gas_saturation));

    const double p_b =
            getScalar(_medium->property(entry_pressure));
    const double lambda =
            getScalar(_medium->property(brooks_corey_exponent));

    const double s_eff = std::pow(p_b/p_cap, lambda);

    _value= s_eff*(s_L_max - s_L_res) + s_L_res;

    double s_L = s_eff*(s_L_max - s_L_res) + s_L_res;

    std::cout   <<  "     s_L_r:  " << s_L_res << "\n";
    std::cout   <<  "     s_L_max:  " << s_L_max << "\n";
    std::cout   <<  "     p_b:  " << p_b << "\n";
    std::cout   <<  "     lambda:  " << lambda << "\n";
    std::cout   <<  "     s_eff:  " << s_eff << "\n";
    std::cout   <<  "     s_L:  " << s_L << "\n";


    isUpdated(true);

    return _value;
}

}  // MaterialPropertyLib
