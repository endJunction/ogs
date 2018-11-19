/**
 * \author Norbert Grunwald
 * \date   16.11.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "pLinearPressure.h"
#include <iostream>
#include "../mpComponent.h"

namespace MaterialPropertyLib
{
LinearPressure::LinearPressure(Phase* phase,
        const double rho_0,
        const double p_0,
        const double beta_p) :
        _component(nullptr),
        _phase(phase),
        _rho_0(rho_0),
        _p_0(p_0),
        _beta_p(beta_p)
{};

LinearPressure::LinearPressure(Component* component,
        const double rho_0,
        const double p_0,
        const double beta_p) :
        _component(component),
        _phase(nullptr),
        _rho_0(rho_0),
        _p_0(p_0),
        _beta_p(beta_p)
{};

/**
 * This method computes a density by a linear correlation with temperature
 * by \f$\rho=\rho_0+\partial\rho/\partial T\left(T-T_\mathrm{ref}\right)\f$
 */
PropertyDataType LinearPressure::value(VariableArray const& v)
{
    const double pressure = getScalar(v[phase_pressure]);
    return _rho_0 * (1 + _beta_p * ( pressure - _p_0));
}

PropertyDataType LinearPressure::dvalue(VariableArray const& v,
        Variables const pv)
{
    assert((pv == Variables::phase_pressure) &&
           "LinearPressure::dvalue is implemented for derivatives with "
           "respect to phase pressure only.");

    return _rho_0 * _beta_p;
}


PropertyDataType LinearPressure::ddvalue(VariableArray const& v,
        Variables const pv1, Variables const pv2)
{
    assert((pv1 == Variables::phase_pressure) &&
           (pv2 == Variables::phase_pressure) &&
           "LinearPressure::ddvalue is implemented for derivatives with "
           "respect to phase pressure only.");
    return 0;
}


}  // namespace MaterialPropertyLib
