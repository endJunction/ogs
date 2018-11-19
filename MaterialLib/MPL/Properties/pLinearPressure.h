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
class LinearPressure final : public Property
{
private:
    /// This property is (currently) implemented to obtain component
    /// properties only.
    Component* _component;
    Phase* _phase;
    const double _rho_0;
    const double _p_0;
    const double _beta_p;

public:

    /// Constructor that passes a pointer to the phase.
    LinearPressure(Phase* /*unused*/,
            const double /*rho_0*/,
            const double /*p_0*/,
            const double /*beta_p*/);
    /// Constructor that passes a pointer to the component.
    LinearPressure(Component* /*c*/,
            const double /*rho_0*/,
            const double /*p_0*/,
            const double /*beta_p*/);

    /// This method overrides the base class implementation and
    /// actually computes and sets the property _value.
    PropertyDataType value(VariableArray const& /*unused*/) override;
    PropertyDataType dvalue(VariableArray const& /*unused*/,
            Variables const) override;
    PropertyDataType ddvalue(VariableArray const& /*unused*/,
            Variables const, Variables const ) override;
};


template <typename MaterialType>
std::unique_ptr<LinearPressure> createLinearPressure(
    BaseLib::ConfigTree const&, MaterialType)
{
    OGS_FATAL(
        "The LinearPressure property is implemented on 'phase' and 'component'"
        " scale only.");
}

template <>
inline std::unique_ptr<LinearPressure> createLinearPressure<Phase*>(
    BaseLib::ConfigTree const& config, Phase* phase)
{
    // check is reading the parameter, not peeking it...
    //! \ogs_file_param{prj__media__medium__properties__property__LinearPressure}
    // config.checkConfigParameter("type", "LinearPressure");
    DBUG("Create LinearPressure phase property");

    auto const rho_0 =
        //! \ogs_file_param{prj__media__medium__properties__property__LinearPressure__rho_0}
        config.getConfigParameter<double>("rho_0");
    auto const p_0 =
        //! \ogs_file_param{prj__media__medium__properties__property__LinearPressure__p_0}
        config.getConfigParameter<double>("p_0");
    auto const beta_p =
        //! \ogs_file_param{prj__media__medium__properties__property__LinearPressure__beta_p}
        config.getConfigParameter<double>("beta_p");

    return std::make_unique<LinearPressure>(phase,rho_0,p_0,beta_p);
}

template <>
inline std::unique_ptr<LinearPressure> createLinearPressure<Component*>(
    BaseLib::ConfigTree const& config, Component* component)
{
    // check is reading the parameter, not peeking it...
    //! \ogs_file_param{prj__media__medium__properties__property__LinearPressure}
    // config.checkConfigParameter("type", "LinearPressure");
    DBUG("Create LinearPressure component property");

    auto const rho_0 =
        //! \ogs_file_param{prj__media__medium__properties__property__LinearPressure__rho_0}
        config.getConfigParameter<double>("rho_0");
    auto const p_0 =
        //! \ogs_file_param{prj__media__medium__properties__property__LinearPressure__p_0}
        config.getConfigParameter<double>("p_0");
    auto const beta_p =
        //! \ogs_file_param{prj__media__medium__properties__property__LinearPressure__beta_p}
        config.getConfigParameter<double>("beta_p");

    return std::make_unique<LinearPressure>(component,rho_0,p_0,beta_p);
}

}  // namespace MaterialPropertyLib
