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
#pragma once

#include "MaterialLib/MPL/mpProperty.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * \class AverageVolumeFraction
 * \brief A function averaging a property by volume fraction
 * \details This property is usually a medium property, it
 * computes the average of individual phase properties
 * weighted by volume fraction.
 */
class SaturationBrooksCorey final : public Property
{
private:
    /// A pointer to the phase object.
    Medium* _medium;
    const double _residual_liquid_saturation;
    const double _residual_gas_saturation;
    const double _exponent;
    const double _entry_pressure;



public:
    /// Constructor passing a pointer to the medium.
    SaturationBrooksCorey(Medium*, const double,
            const double, const double, const double);

    /// This method overrides the base class implementation and
    /// actually computes and sets the property _value.
    PropertyDataType value(VariableArray const&) override;
    PropertyDataType dvalue(VariableArray const&, Variables const ) override;
    PropertyDataType ddvalue(VariableArray const&, Variables const,
            Variables const) override;
};

template <typename MaterialType>
std::unique_ptr<SaturationBrooksCorey> createSaturationBrooksCorey(
    BaseLib::ConfigTree const&, MaterialType)
{
    OGS_FATAL(
        "The SaturationBrooksCorey property is implemented only on 'media' "
        "scale.");
}

template <>
inline std::unique_ptr<SaturationBrooksCorey> createSaturationBrooksCorey<Medium*>(
    BaseLib::ConfigTree const& config, Medium* medium)
{
    // check is reading the parameter, not peeking it...
    //! \ogs_file_param{prj__media__medium__properties__property__SaturationBrooksCorey}
    // config.checkConfigParameter("type", "SaturationBrooksCorey");
    DBUG("Create SaturationBrooksCorey medium property");

    auto const residual_liquid_saturation =
        //! \ogs_file_param{prj__media__medium__properties__property__SaturationBrooksCorey__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{prj__media__medium__properties__property__SaturationBrooksCorey__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const exponent =
        //! \ogs_file_param{prj__media__medium__properties__property__SaturationBrooksCorey__lambda}
        config.getConfigParameter<double>("lambda");
    auto const entry_pressure =
        //! \ogs_file_param{prj__media__medium__properties__property__SaturationBrooksCorey__entry_pressure}
        config.getConfigParameter<double>("entry_pressure");

    return std::make_unique<SaturationBrooksCorey>(medium,
            residual_liquid_saturation, residual_gas_saturation, exponent,
            entry_pressure);
}

}  // MaterialPropertyLib
