/**
 * \author Norbert Grunwald
 * \date   Oct 22, 2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "pSaturationFredlund.h"
#include "MaterialLib/MPL/mpMedium.h"
#include "MaterialLib/MPL/Properties/pUniversalConstants.h"
#include "MathLib/MathTools.h"

#include <iostream>

namespace MaterialPropertyLib
{
/// This constructor throws an error since it was used for a medium
/// property while it is a component property.
SaturationFredlund::SaturationFredlund(Medium* m,
                                               double const s_max,
                                               double const s_min,
                                               double const parameter_a,
                                               double const parameter_n,
                                               double const parameter_m,
                                               double const parameter_psi_r)
    : _medium(m),
      _s_max(s_max),
      _s_min(s_min),
      _parameter_a(parameter_a),
      _parameter_n(parameter_n),
      _parameter_m(parameter_m),
      _parameter_psi_r(parameter_psi_r)
{
}

/**
 */
PropertyDataType SaturationFredlund::value(VariableArray const& v)
{
    double const pc = std::max(0., getScalar(v[capillary_pressure]));

    const double s_max = _s_max;
    const double s_min = _s_min;
    const double a = _parameter_a;
    const double n = _parameter_n;
    const double m = _parameter_m;
    const double psi_r = _parameter_psi_r;

    const double A = eulersNumber + std::pow(pc/a, n);
    const double B = 1.0 / std::pow(std::log(A), m);
    const double C = 1. - std::log(1 + pc/psi_r) / std::log(1 + 1.e9/psi_r);

    const double s_eff = B*C;
    const double s_L = s_eff * (s_max - s_min) + s_min;
    return s_L;
}

PropertyDataType SaturationFredlund::dvalue(VariableArray const& v,
                                                Variables const pv)
{
    assert((pv == Variables::capillary_pressure) &&
           "SaturationFredlund::dvalue is implemented for derivatives with "
           "respect to capillary pressure only.");

    double const pc = std::max(0., getScalar(v[capillary_pressure]));

    const double s_max = _s_max;
    const double s_min = _s_min;
    const double a = _parameter_a;
    const double n = _parameter_n;
    const double m = _parameter_m;
    const double psi_r = _parameter_psi_r;

    const double A = eulersNumber + std::pow(pc/a, n);
    const double B = s_max * std::pow(std::log(A), -m);
    const double C = 1. - std::log(1 + pc/psi_r) / std::log(1 + 1.e9/psi_r);

    const double dA = n/pc * std::pow(pc/a, n);
    const double D = -m*s_max*dA/A;
    const double dB = D * std::pow(std::log(A),-m-1.);
    const double dC = -1./((pc+psi_r)*std::log((psi_r + 1.e9)/psi_r));

    const double ds_eff_dpc = dB*C + dC*B;
    const double ds_L_dpc = ds_eff_dpc * (s_max - s_min);
    return ds_L_dpc;
}

PropertyDataType SaturationFredlund::ddvalue(VariableArray const& v,
                                                Variables const pv1,
                                                Variables const pv2)
{
    assert((pv1 == Variables::capillary_pressure) &&
            (pv2 == Variables::capillary_pressure) &&
            "SaturationFredlund::ddvalue is implemented for 2nd derivatives"
            "with respect to capillary pressure only.");

    double const pc = std::max(0., getScalar(v[capillary_pressure]));

    const double s_max = _s_max;
    const double s_min = _s_min;
    const double a = _parameter_a;
    const double n = _parameter_n;
    const double m = _parameter_m;
    const double psi_r = _parameter_psi_r;

    const double A = eulersNumber + std::pow(pc/a, n);
    const double B = s_max * std::pow(std::log(A), -m);
    const double C = 1. - std::log(1 + pc/psi_r) / std::log(1 + 1.e9/psi_r);

    const double dA = n/pc * std::pow(pc/a, n);
    const double D = -m*s_max*dA/A;
    const double dB = D * std::pow(std::log(A),-m-1.);
    const double dC = -1./((pc+psi_r)*std::log((psi_r + 1.e9)/psi_r));
    const double ddA = (n-1)/(pc*pc)*n*std::pow(pc/a,n);
    const double dD = m*s_max*(dA*dA-A*ddA)/(A*A);
    const double g = -m-1.;
    const double ddB = std::pow(std::log(A),g) *
            (g*D*dA/(A*std::log(A))+dD);
    const double ddC = dC / (-(psi_r + pc));

    const double d2sldpc2 = ddB*C + 2*dC*dB + B*ddC;

    return d2sldpc2 * (s_max - s_min);
}

}  // namespace MaterialPropertyLib
