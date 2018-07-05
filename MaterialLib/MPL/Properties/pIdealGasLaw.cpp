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
    if (isUpdated())
    {
        return _value;
    }

    double molar_mass;

    if (_phase)
        molar_mass = getScalar(_phase->property(PropertyEnum::molar_mass));
    else
        molar_mass = getScalar(_component->property(PropertyEnum::molar_mass));

    const double R = gasConstant;
    const double p = getScalar(v[PrimaryVariables::p_GR]);
    const double T = getScalar(v[PrimaryVariables::T]);

    _value = p*molar_mass / R / T;

    isUpdated(true);

    return _value;
}

PropertyDataType IdealGasLaw::dvalue(VariableArray const& v, PrimaryVariables const pv)
{

    double density;
    if (_phase)
        density = getScalar(_phase->property(PropertyEnum::density));
    else
        density = getScalar(_component->property(PropertyEnum::density));

    switch (pv)
    {
    case PrimaryVariables::p_GR:
    {
        const double p = getScalar(v[pv]);
        _dvalue = density / p;
        return _dvalue;
    }
    case PrimaryVariables::T:
    {
        const double T = getScalar(v[pv]);
        _dvalue = - density / T;
        return _dvalue;
    }
    default:
        return 0.;
    }
    return _dvalue;
}

}  // namespace MaterialPropertyLib
