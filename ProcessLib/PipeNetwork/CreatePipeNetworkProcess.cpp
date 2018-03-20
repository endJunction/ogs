/**
* \copyright
* Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#include "CreatePipeNetworkProcess.h"

#include "PipeNetworkProcess.h"
#include "PipeNetworkProcessData.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
    namespace PipeNetwork
    {
        std::unique_ptr<Process> createPipeNetworkProcess(
            MeshLib::Mesh& mesh,
            std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
            std::vector<ProcessVariable> const& variables,
            std::vector<std::unique_ptr<ParameterBase>> const& parameters,
            unsigned const integration_order,
            BaseLib::ConfigTree const& config,
            std::map<std::string, std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const& curves)
        {
            //! \ogs_file_param{prj__processes__process__type}
            config.checkConfigParameter("type", "PIPE_NETWORK");

            DBUG("Create PIPE_NETWORK Process.");

            //! \ogs_file_param{prj__processes__process__HEAT_CONDUCTION__process_variables}
            auto const pv_config = config.getConfigSubtree("process_variables");

            std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
                process_variables;
            auto per_process_variables = findProcessVariables(
                variables, pv_config,
                {//! \ogs_file_param_special{prj__processes__process__HEAT_CONDUCTION__process_variables__process_variable}
                    "process_variable" });
            process_variables.push_back(std::move(per_process_variables));

            // Process variable.

            // thermal conductivity parameter.
            auto& thermal_conductivity = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_CONDUCTION__thermal_conductivity}
                "thermal_conductivity", parameters, 1);

            DBUG("Use \'%s\' as thermal conductivity parameter.",
                thermal_conductivity.name.c_str());

            // heat capacity parameter.
            auto& heat_capacity = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_CONDUCTION__heat_capacity}
                "heat_capacity", parameters, 1);

            DBUG("Use \'%s\' as heat capacity parameter.", heat_capacity.name.c_str());

            // density parameter.
            auto& density = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_CONDUCTION__density}
                "density", parameters, 1);

            DBUG("Use \'%s\' as density parameter.", density.name.c_str());

            PipeNetworkProcessData process_data{ thermal_conductivity, heat_capacity,
                density };

            SecondaryVariableCollection secondary_variables;

            NumLib::NamedFunctionCaller named_function_caller(
                { "PipeNetwork_flow_velocity" }); // TODO be changed accordingly. 

            ProcessLib::createSecondaryVariables(config, secondary_variables,
                named_function_caller);

            return std::make_unique<PipeNetworkProcess>(
                mesh, std::move(jacobian_assembler), parameters, integration_order,
                std::move(process_variables), std::move(process_data),
                std::move(secondary_variables), std::move(named_function_caller));
        }

    }  // namespace PipeNetwork
}  // namespace ProcessLib
