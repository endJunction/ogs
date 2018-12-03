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

class RelPermBrooksCorey final : public Property
{
private:
    /// A pointer to the phase object.
    Medium* _medium;
    const double _residual_liquid_saturation;
    const double _residual_gas_saturation;
    const double _k_rel_min_liquid;
    const double _k_rel_min_gas;
    const double _exponent;

public:
    /// Constructor passing a pointer to the medium.
    RelPermBrooksCorey(Medium*, const double, const double,
    		const double, const double, const double);

    /// This method overrides the base class implementation and
    /// actually computes and sets the property _value.
    PropertyDataType value(VariableArray const&) override;
    PropertyDataType dvalue(VariableArray const&, Variables const) override;

};

template <typename MaterialType>
std::unique_ptr<RelPermBrooksCorey> createRelPermBrooksCorey(
    BaseLib::ConfigTree const&, MaterialType)
{
    OGS_FATAL(
        "The RelPermBrooksCorey property is implemented only on 'media' "
        "scale.");
}

template <>
inline std::unique_ptr<RelPermBrooksCorey> createRelPermBrooksCorey<Medium*>(
    BaseLib::ConfigTree const& config, Medium* medium)
{
    // check is reading the parameter, not peeking it...
    //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey}
    // config.checkConfigParameter("type", "RelPermBrooksCorey");
    DBUG("Create RelPermBrooksCorey medium property");

    auto const residual_liquid_saturation =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const k_rel_min_liquid =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__k_rel_min_liquid}
        config.getConfigParameter<double>("k_rel_min_liquid");
    auto const k_rel_min_gas =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__k_rel_min_gas}
        config.getConfigParameter<double>("k_rel_min_gas");
    auto const exponent =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__lambda}
        config.getConfigParameter<double>("lambda");

    return std::make_unique<RelPermBrooksCorey>(medium,
            residual_liquid_saturation, residual_gas_saturation,
			k_rel_min_liquid, k_rel_min_gas, exponent);
}


}  // MaterialPropertyLib
