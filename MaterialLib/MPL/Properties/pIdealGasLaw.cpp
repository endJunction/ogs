/**
 * \author Norbert Grunwald
 * \date   18.09.2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "pIdealGasLaw.h"
#include "MaterialLib/MPL/mpComponent.h"
//#include "MaterialLib/MPL/mpEnums.h"

#include "../mpPhase.h"
#include "pUniversalConstants.h"

#include <boost/variant.hpp>
#include <iostream>

namespace MaterialPropertyLib
{
IdealGasLaw::IdealGasLaw(Medium* /*unused*/)
    : _phase(nullptr), _component(nullptr)
{
    notImplemented("IdealGasLaw", "Medium");
}

IdealGasLaw::IdealGasLaw(Phase* p) : _phase(p), _component(nullptr){};

IdealGasLaw::IdealGasLaw(Component* c) : _phase(nullptr), _component(c){};

/**
 */
PropertyDataType IdealGasLaw::value(VariableArray const& v)
{
    double molar_mass;
    if (_phase)
        molar_mass = getScalar(_phase->property(PropertyEnum::molar_mass));
    else
        molar_mass = getScalar(_component->property(PropertyEnum::molar_mass));

    const double R = gasConstant;
    const double p = getScalar(v[Variables::phase_pressure]);
    const double T = getScalar(v[Variables::temperature]);

    const double density = p*molar_mass / R / T;
    return density;
}

PropertyDataType IdealGasLaw::dvalue(VariableArray const& v, Variables const pv)
{

    const double density= getScalar(v[Variables::gas_density]);

    switch (pv)
    {
    case Variables::phase_pressure:
    {
        const double p = getScalar(v[pv]);
        return density / p;
    }
    case Variables::temperature:
    {
        const double T = getScalar(v[pv]);
        return -density / T;
    }
    default:
        return 0.;
    }
    return 0.;
}

}  // namespace MaterialPropertyLib
