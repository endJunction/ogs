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
#pragma once

#include "MaterialLib/MPL/mpProperty.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * This property (usually a component property) computes a linear
 * density function of temperature based on a reference density,
 * a slope, and a reference temperature.
 */
class BilinearTemperaturePressure final : public Property
{
private:
    Phase* _phase;
    Component* _component;
    double const _reference_density;
    double const _reference_temperature;
    double const _reference_pressure;

public:
    /// Constructor that passes a pointer to the phase.
    BilinearTemperaturePressure(Phase* p,
                                double const reference_density,
                                double const reference_temperature,
                                double const reference_pressure);
    /// Constructor that passes a pointer to the component.
    BilinearTemperaturePressure(Component* c,
                                double const reference_density,
                                double const reference_temperature,
                                double const reference_pressure);

    /// This method overrides the base class implementation and
    /// actually computes and sets the property _value.
    PropertyDataType value(VariableArray const& /*unused*/) override;
    PropertyDataType dvalue(VariableArray const& /*unused*/,
                            Variables const /*unused*/) override;
};

template <typename MaterialType>
std::unique_ptr<BilinearTemperaturePressure> createBilinearTemperaturePressure(
    BaseLib::ConfigTree const&, MaterialType)
{
    OGS_FATAL(
        "The BilinearTemperaturePressure property is implemented only on "
        "'phase' and 'component' scales.");
}

template <>
inline std::unique_ptr<BilinearTemperaturePressure>
createBilinearTemperaturePressure<Phase*>(BaseLib::ConfigTree const& config,
                                          Phase* phase)
{
    // check is reading the parameter, not peeking it...
    //! \ogs_file_param{prj__media__medium__properties__property__BilinearTemperaturePressure}
    // config.checkConfigParameter("type", "BilinearTemperaturePressure");
    DBUG("Create BilinearTemperaturePressure phase property");

    auto const reference_density =
        //! \ogs_file_param{prj__media__medium__properties__property__BilinearTemperaturePressure__reference_density}
        config.getConfigParameter<double>("reference_density");
    auto const reference_temperature =
        //! \ogs_file_param{prj__media__medium__properties__property__BilinearTemperaturePressure__reference_temperature}
        config.getConfigParameter<double>("reference_temperature");
    auto const reference_pressure =
        //! \ogs_file_param{prj__media__medium__properties__property__BilinearTemperaturePressure__reference_pressure}
        config.getConfigParameter<double>("reference_pressure");

    return std::make_unique<BilinearTemperaturePressure>(
        phase, reference_density, reference_temperature, reference_pressure);
}

template <>
inline std::unique_ptr<BilinearTemperaturePressure>
createBilinearTemperaturePressure<Component*>(BaseLib::ConfigTree const& config,
                                              Component* component)
{
    // check is reading the parameter, not peeking it...
    //! \ogs_file_param{prj__media__medium__properties__property__BilinearTemperaturePressure}
    // config.checkConfigParameter("type", "BilinearTemperaturePressure");
    DBUG("Create BilinearTemperaturePressure component property");

    auto const reference_density =
        //! \ogs_file_param{prj__media__medium__properties__property__BilinearTemperaturePressure__reference_density}
        config.getConfigParameter<double>("reference_density");
    auto const reference_temperature =
        //! \ogs_file_param{prj__media__medium__properties__property__BilinearTemperaturePressure__reference_temperature}
        config.getConfigParameter<double>("reference_temperature");
    auto const reference_pressure =
        //! \ogs_file_param{prj__media__medium__properties__property__BilinearTemperaturePressure__reference_pressure}
        config.getConfigParameter<double>("reference_pressure");

    return std::make_unique<BilinearTemperaturePressure>(
        component, reference_density, reference_temperature,
        reference_pressure);
}

}  // namespace MaterialPropertyLib
