/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter

#include "ThermoPlasticBDT.h"

namespace MaterialLib
{
namespace Solids
{
namespace ThermoPlasticBDT
{
inline NumLib::NewtonRaphsonSolverParameters
createNewtonRaphsonSolverParameters(BaseLib::ConfigTree const& config)
{
    DBUG("Create local nonlinear solver parameters.");
    auto const maximum_iterations =
        //! \ogs_file_param{material__solid__constitutive_relation__ThermoPlasticBDT__nonlinear_solver__maximum_iterations}
        config.getConfigParameter<int>("maximum_iterations");

    DBUG("\tmaximum_iterations: %d.", maximum_iterations);

    auto const error_tolerance =
        //! \ogs_file_param{material__solid__constitutive_relation__ThermoPlasticBDT__nonlinear_solver__error_tolerance}
        config.getConfigParameter<double>("error_tolerance");

    DBUG("\terror_tolerance: %g.", error_tolerance);

    return {maximum_iterations, error_tolerance};
}

inline std::unique_ptr<DamagePropertiesParameters> createDamageProperties(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param_special{material__solid__constitutive_relation__ThermoPlasticBDT__damage_properties__alpha_d}
    auto& alpha_d =
        ProcessLib::findParameter<double>(config, "alpha_d", parameters, 1);

    DBUG("Use \'%s\' as alpha_d.", alpha_d.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__ThermoPlasticBDT__damage_properties__beta_d}
    auto& beta_d =
        ProcessLib::findParameter<double>(config, "beta_d", parameters, 1);

    DBUG("Use \'%s\' as beta_d.", beta_d.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__ThermoPlasticBDT__damage_properties__h_d}
    auto& h_d = ProcessLib::findParameter<double>(config, "h_d", parameters, 1);

    DBUG("Use \'%s\' as h_d.", h_d.name.c_str());

    //! \ogs_file_param{material__solid__constitutive_relation__ThermoPlasticBDT__damage_properties__m_d}
    double const m_d = config.getConfigParameter<double>("m_d");
    DBUG("Over-non-local value m_d: %g", m_d);

    return std::make_unique<DamagePropertiesParameters>(
        DamagePropertiesParameters{alpha_d, beta_d, h_d, m_d});
}

template <int DisplacementDim>
std::unique_ptr<SolidThermoPlasticBDT<DisplacementDim>> createThermoPlasticBDT(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__solid__constitutive_relation__type}
    config.checkConfigParameter("type", "ThermoPlasticBDT");
    DBUG("Create ThermoPlasticBDT material");

    //! \ogs_file_param_special{material__solid__constitutive_relation__ThermoPlasticBDT__shear_modulus}
    auto& shear_modulus = ProcessLib::findParameter<double>(
        config, "shear_modulus", parameters, 1);

    DBUG("Use \'%s\' as shear modulus parameter.", shear_modulus.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__ThermoPlasticBDT__bulk_modulus}
    auto& bulk_modulus = ProcessLib::findParameter<double>(
        config, "bulk_modulus", parameters, 1);

    DBUG("Use \'%s\' as bulk modulus parameter.", bulk_modulus.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__ThermoPlasticBDT__fc}
    auto& fc = ProcessLib::findParameter<double>(config, "fc", parameters, 1);

    DBUG("Use \'%s\' as fc.", fc.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__ThermoPlasticBDT__m}
    auto& m = ProcessLib::findParameter<double>(config, "m", parameters, 1);

    DBUG("Use \'%s\' as m.", m.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__ThermoPlasticBDT__qp0}
    auto& qp0 = ProcessLib::findParameter<double>(config, "qp0", parameters, 1);

    DBUG("Use \'%s\' as qp0.", qp0.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__ThermoPlasticBDT__alpha}
    auto& alpha =
        ProcessLib::findParameter<double>(config, "alpha", parameters, 1);

    DBUG("Use \'%s\' as alpha.", alpha.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__ThermoPlasticBDT__n}
    auto& n = ProcessLib::findParameter<double>(config, "n", parameters, 1);

    DBUG("Use \'%s\' as n.", n.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__ThermoPlasticBDT__T0}
    auto& T0 = ProcessLib::findParameter<double>(config, "T0", parameters, 1);

    DBUG("Use \'%s\' as T0.", T0.name.c_str());

    MaterialPropertiesParameters mp{shear_modulus, bulk_modulus, fc, m,
                                    qp0,           alpha,        n,  T0};

    // Damage properties.
    std::unique_ptr<DamagePropertiesParameters>
        thermoplasticBDT_damage_properties;

    auto const& thermoplasticBDT_damage_config =
        //! \ogs_file_param{material__solid__constitutive_relation__ThermoPlasticBDT__damage_properties}
        config.getConfigSubtreeOptional("damage_properties");
    if (thermoplasticBDT_damage_config)
    {
        thermoplasticBDT_damage_properties =
            createDamageProperties(parameters, *thermoplasticBDT_damage_config);
    }

    auto const& nonlinear_solver_config =
        //! \ogs_file_param{material__solid__constitutive_relation__ThermoPlasticBDT__nonlinear_solver}
        config.getConfigSubtree("nonlinear_solver");
    auto const nonlinear_solver_parameters =
        createNewtonRaphsonSolverParameters(nonlinear_solver_config);

    return std::make_unique<SolidThermoPlasticBDT<DisplacementDim>>(
        nonlinear_solver_parameters,
        mp,
        std::move(thermoplasticBDT_damage_properties));
}

}  // namespace ThermoPlasticBDT
}  // namespace Solids
}  // namespace MaterialLib
