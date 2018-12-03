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
RelPermBrooksCorey::RelPermBrooksCorey(Medium* m,
        const double residual_liquid_saturation,
		const double residual_gas_saturation,
		const double k_rel_min_liquid,
		const double k_rel_min_gas,
        const double exponent)

: _medium(m),
  _residual_liquid_saturation(residual_liquid_saturation),
  _residual_gas_saturation(residual_gas_saturation),
  _k_rel_min_liquid(k_rel_min_liquid),
  _k_rel_min_gas(k_rel_min_gas),
  _exponent(exponent){};

/**
 */
PropertyDataType  RelPermBrooksCorey::value(VariableArray const& v)
{

    const double s_L = getScalar(v[MaterialPropertyLib::liquid_saturation]);

    const double s_L_res = _residual_liquid_saturation;
    const double s_L_max = 1. - _residual_gas_saturation;
    const double k_rel_min_LR = _k_rel_min_liquid;
    const double k_rel_min_GR = _k_rel_min_gas;

    const double lambda = _exponent;

    const double s_eff = (s_L - s_L_res)/(s_L_max - s_L_res);

    if (s_eff > 1.0)
    {
    	return Pair {{1.0, 0.0}};
    }
    if (s_eff < 0.0)
    {
    	return Pair {{0.0, 1.0}};
    }

    const double k_rel_LR = std::pow(s_eff, (2.+3.*lambda)/lambda);
    const double k_rel_GR = (1.-s_eff)*(1. - s_eff) * (1. - std::pow(s_eff, (2.+lambda)/lambda));

    const Pair kRel = {std::max(k_rel_LR, k_rel_min_LR),
    		std::max(k_rel_GR, k_rel_min_GR)};

    return kRel;
}

PropertyDataType  RelPermBrooksCorey::dvalue(VariableArray const& v,
        Variables const pv)
{
    assert((pv == Variables::liquid_saturation) &&
                "RelPermBrooksCorey::dvalue is implemented for "
                " derivatives with respect to liquid saturation only.");

    const double s_L = getScalar(v[MaterialPropertyLib::liquid_saturation]);

    const double s_L_res = _residual_liquid_saturation;
    const double s_L_max = 1. - _residual_gas_saturation;
    const double lambda =_exponent;

    const double s_eff = (s_L - s_L_res)/(s_L_max - s_L_res);
    const double dsedsL = 1. / (s_L_max - s_L_res);

    const double dk_rel_LRdse = (3*lambda + 2.)/lambda *
            std::pow(s_eff, 2./lambda + 2.);
    const double dk_rel_LRdsL = dk_rel_LRdse * dsedsL;

    const double _2L_L = (2.* lambda) / lambda;
    const double dk_rel_GRdse = -2.*(1-s_eff)*
            (1. - std::pow(s_eff, _2L_L)) - _2L_L * std::pow(s_eff, _2L_L-1.) *
            (1. - s_eff) * (1. - s_eff);
    const double dk_rel_GRdsL = dk_rel_GRdse * dsedsL;

    const Pair dkReldsL = {dk_rel_LRdsL, dk_rel_GRdsL};

    return dkReldsL;
}

}  // MaterialPropertyLib
