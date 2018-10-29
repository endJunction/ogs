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
#include "pBilinearTemperaturePressure.h"
#include <iostream>
#include "../mpPhase.h"
#include "../mpComponent.h"


namespace MaterialPropertyLib
{
/// This constructor throws an error since it was used for a phase
/// property while it is a component property.
BilinearTemperaturePressure::BilinearTemperaturePressure(Phase* p,
        double const reference_density,
        double const reference_temperature,
        double const reference_pressure)
    : _phase(p),
      _component(nullptr),
      _reference_density(reference_density),
      _reference_temperature(reference_temperature),
      _reference_pressure(reference_pressure){};

/// This constructor copies the pointer of the component from the
/// arguments into the private attribute.
BilinearTemperaturePressure::BilinearTemperaturePressure(Component* c,
        double const reference_density,
        double const reference_temperature,
        double const reference_pressure)
    : _phase(nullptr),
      _component(c),
      _reference_density(reference_density),
      _reference_temperature(reference_temperature),
      _reference_pressure(reference_pressure){};

/**
 * This method computes a density by a linear correlation with temperature
 * by \f$\rho=\rho_0+\partial\rho/\partial T\left(T-T_\mathrm{ref}\right)\f$
 */
PropertyDataType BilinearTemperaturePressure::value(VariableArray const& v)
{
    const double p = getScalar(v[Variables::phase_pressure]);
    const double T = getScalar(v[Variables::temperature]);

    const double rho_0 = _reference_density;
    const double p_ref = _reference_pressure;
    const double T_ref = _reference_temperature;
    double beta_p;
    double beta_T;

    if (_phase)
    {
        beta_p = getScalar(_phase->property(compressibility));
        beta_T = getScalar(_phase->property(thermal_expansivity));
    }
    else
    {
        beta_p = getScalar(_component->property(compressibility));
        beta_T = getScalar(_component->property(thermal_expansivity));
    }

    const double density = rho_0 * (1 + beta_p * (p - p_ref) +
            beta_T * (T - T_ref));

    return density;
}

PropertyDataType BilinearTemperaturePressure::dvalue(VariableArray const& v,
        Variables const pv)
{

    const double p = getScalar(v[Variables::phase_pressure]);
    const double T = getScalar(v[Variables::temperature]);
    double derivative;

    const double rho_0 = _reference_density;
    double beta_p;
    double beta_T;

    if (_phase)
    {
        beta_p = getScalar(_phase->property(compressibility));
        beta_T = getScalar(_phase->property(thermal_expansivity));
    }
    else
    {
        beta_p = getScalar(_component->property(compressibility));
        beta_T = getScalar(_component->property(thermal_expansivity));
    }

    switch (pv)
    {
    case Variables::phase_pressure:
        derivative = rho_0 * beta_p;
        break;
    case Variables::temperature:
        derivative = rho_0 * beta_T;
        break;
    }
    return derivative;
}

}  // namespace MaterialPropertyLib



