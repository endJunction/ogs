/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once
#include "RichardsFlowMaterialProperties.h"

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace RichardsFlow
{
struct RichardsFlowProcessData
{
    RichardsFlowProcessData(
        std::unique_ptr<RichardsFlowMaterialProperties>&& material_,
        Eigen::VectorXd const specific_body_force_,
        ParameterLib::Parameter<double> const& temperature_,
        bool const has_gravity_,
        bool const has_mass_lumping_)
        : material(std::move(material_)),
          specific_body_force(specific_body_force_),
          temperature(temperature_),
          has_gravity(has_gravity_),
          has_mass_lumping(has_mass_lumping_)
    {
    }

    RichardsFlowProcessData(RichardsFlowProcessData&&) = default;

    //! Copies are forbidden.
    RichardsFlowProcessData(RichardsFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(RichardsFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(RichardsFlowProcessData&&) = delete;

    std::unique_ptr<RichardsFlowMaterialProperties> material;
    Eigen::VectorXd const specific_body_force;
    ParameterLib::Parameter<double> const& temperature;
    double dt;
    bool const has_gravity;
    bool const has_mass_lumping;
};

}  // namespace RichardsFlow
}  // namespace ProcessLib
