/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateNonequilibriumInitialState.h"
#include "NonequilibriumInitialState.h"

#include "ParameterLib/Utils.h"

namespace ProcessLib
{
std::unique_ptr<NonequilibriumInitialState> createNonequilibriumInitialState(
    boost::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::pair<std::string, int>> const& tag_names_and_components,
    std::vector<ProcessVariable> const& process_variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    MeshLib::Mesh const& mesh)
{
    if (!config)
    {
        return nullptr;
    }

    std::map<std::string, ParameterLib::Parameter<double> const*>
        name_parameter_map;
    for (auto [tag_name, n_components] : tag_names_and_components)
    {
        auto const parameter_name =
            //! \ogs_file_special
            config->getConfigParameterOptional<std::string>(tag_name);

        if (!parameter_name)
        {
            continue;
        }

        auto const& parameter = ParameterLib::findParameter<double>(
            *parameter_name, parameters, n_components, &mesh);

        if (name_parameter_map.find(tag_name) != cend(name_parameter_map))
        {
            OGS_FATAL("Tag name '%s' cannot be present multiple times.",
                      tag_name.c_str());
        }

        name_parameter_map[tag_name] = &parameter;
    }

    std::vector<std::reference_wrapper<ProcessVariable>> variables;
    const auto& equilibrate_process_variables_config =
        //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION__nonequilibrium_initial_state__equilibrate_process_variables}
        config->getConfigSubtree("equilibrate_process_variables");
    for (
        auto const& pv :
        //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION__nonequilibrium_initial_state__equilibrate_process_variables__process_variable}
        equilibrate_process_variables_config.getConfigSubtreeList(
            "process_variable"))
    {
        auto const name = pv.getValue<std::string>();
    }

    return std::make_unique<NonequilibriumInitialState>(
        NonequilibriumInitialState{std::move(name_parameter_map)});
}
}  // namespace ProcessLib
