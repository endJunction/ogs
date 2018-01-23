/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateMohrCoulomb.h"

#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter
#include "MohrCoulomb.h"

namespace MaterialLib
{
namespace Fracture
{
// TODO (naumov); this is copy-paste from CreateEhlers => extract.
inline NumLib::NewtonRaphsonSolverParameters
createNewtonRaphsonSolverParameters(BaseLib::ConfigTree const& config)
{
    DBUG("Create local nonlinear solver parameters.");
    auto const maximum_iterations =
        //! \ogs_file_param{material__solid__constitutive_relation__Ehlers__nonlinear_solver__maximum_iterations}
        config.getConfigParameter<int>("maximum_iterations");

    DBUG("\tmaximum_iterations: %d.", maximum_iterations);

    auto const error_tolerance =
        //! \ogs_file_param{material__solid__constitutive_relation__Ehlers__nonlinear_solver__error_tolerance}
        config.getConfigParameter<double>("error_tolerance");

    DBUG("\terror_tolerance: %g.", error_tolerance);

    return {maximum_iterations, error_tolerance};
}

template <int DisplacementDim>
std::unique_ptr<FractureModelBase<DisplacementDim>>
createMohrCoulomb(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fracture_model__type}
    config.checkConfigParameter("type", "MohrCoulomb");
    DBUG("Create MohrCoulomb material");

    auto& Kn = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__fracture_model__MohrCoulomb__normal_stiffness}
        config, "normal_stiffness", parameters, 1);

    auto& Ks = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__fracture_model__MohrCoulomb__shear_stiffness}
        config, "shear_stiffness", parameters, 1);

    auto& friction_angle = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__fracture_model__MohrCoulomb__friction_angle}
        config, "friction_angle", parameters, 1);

    auto& dilatancy_angle = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__fracture_model__MohrCoulomb__dilatancy_angle}
        config, "dilatancy_angle", parameters, 1);

    auto& cohesion = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__fracture_model__MohrCoulomb__cohesion}
        config, "cohesion", parameters, 1);

    auto const penalty_aperture_cutoff =
        //! \ogs_file_param{material__fracture_model__MohrCoulomb__penalty_aperture_cutoff}
        config.getConfigParameter<double>("penalty_aperture_cutoff");

    auto const tension_cutoff =
        //! \ogs_file_param{material__fracture_model__MohrCoulomb__tension_cutoff}
        config.getConfigParameter<bool>("tension_cutoff");

    typename MohrCoulomb::MohrCoulomb<DisplacementDim>::MaterialProperties mp{
        Kn, Ks, friction_angle, dilatancy_angle, cohesion};

    auto const& nonlinear_solver_config =
        //! \ogs_file_param{material__solid__constitutive_relation__Ehlers__nonlinear_solver}
        config.getConfigSubtree("nonlinear_solver");
    auto const nonlinear_solver_parameters =
        createNewtonRaphsonSolverParameters(nonlinear_solver_config);

    return std::make_unique<MohrCoulomb::MohrCoulomb<DisplacementDim>>(
        nonlinear_solver_parameters, penalty_aperture_cutoff, tension_cutoff,
        mp);
}


template
std::unique_ptr<FractureModelBase<2>>
createMohrCoulomb(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

template
std::unique_ptr<FractureModelBase<3>>
createMohrCoulomb(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

}  // namespace Fracture
}  // namespace MaterialLib
