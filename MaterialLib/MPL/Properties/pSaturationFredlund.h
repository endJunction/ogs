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
#pragma once

#include "MaterialLib/MPL/mpProperty.h"


namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

class SaturationFredlund final : public Property
{

public:
    /// Constructor passing a pointer to the medium.
    SaturationFredlund(Medium* m, double const s_max,
            double const s_min,
            double const parameter_a,
            double const parameter_n,
            double const parameter_m,
            double const parameter_psi_r);

    /// This method overrides the base class implementation and
    /// actually computes and sets the property _value.
    PropertyDataType value(VariableArray const&) override;
    PropertyDataType dvalue(VariableArray const&,
            Variables const ) override;
    PropertyDataType ddvalue(VariableArray const&,
            Variables const,
            Variables const ) override;

private:
    Medium* _medium;
    double const _s_max;
    double const _s_min;
    double const _parameter_a;
    double const _parameter_n;
    double const _parameter_m;
    double const _parameter_psi_r;

};

template <typename MaterialType>
std::unique_ptr<SaturationFredlund> createSaturationFredlund(
    BaseLib::ConfigTree const&, MaterialType)
{
    OGS_FATAL(
        "The SaturationFredlund property is implemented only on 'media' "
        "scale.");
}

template <>
inline std::unique_ptr<SaturationFredlund> createSaturationFredlund<Medium*>(
    BaseLib::ConfigTree const& config, Medium* medium)
{
    // check is reading the parameter, not peeking it...
    //! \ogs_file_param{prj__media__medium__properties__property__SaturationFredlund}
    // config.checkConfigParameter("type", "SaturationFredlund");
    DBUG("Create SaturationFredlund medium property");

    auto const s_max =
        //! \ogs_file_param{prj__media__medium__properties__property__SaturationFredlund__s_max}
        config.getConfigParameter<double>("s_max");
    auto const s_min =
        //! \ogs_file_param{prj__media__medium__properties__property__SaturationFredlund__s_min}
        config.getConfigParameter<double>("s_min");
    auto const parameter_a =
        //! \ogs_file_param{prj__media__medium__properties__property__SaturationFredlund__parameter_a}
        config.getConfigParameter<double>("parameter_a");
    auto const parameter_n =
        //! \ogs_file_param{prj__media__medium__properties__property__SaturationFredlund__parameter_n}
        config.getConfigParameter<double>("parameter_n");
    auto const parameter_m =
        //! \ogs_file_param{prj__media__medium__properties__property__SaturationFredlund__parameter_m}
        config.getConfigParameter<double>("parameter_m");
    auto const parameter_psi_r =
        //! \ogs_file_param{prj__media__medium__properties__property__SaturationFredlund__parameter_psi_r}
        config.getConfigParameter<double>("parameter_psi_r");

    return std::make_unique<SaturationFredlund>(medium,
            s_max, s_min, parameter_a,parameter_n,parameter_m,parameter_psi_r );
}

}  // MaterialPropertyLib
