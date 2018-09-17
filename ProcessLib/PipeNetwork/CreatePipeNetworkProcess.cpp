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
            std::map<std::string, std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const& /*curves*/)
        {
            //! \ogs_file_param{prj__processes__process__type}
            config.checkConfigParameter("type", "PIPE_NETWORK");

            DBUG("Create PIPE_NETWORK Process.");

            //! \ogs_file_param{prj__processes__process__PIPE_NETWORK__process_variables}
            auto const pv_config = config.getConfigSubtree("process_variables");

            std::vector<std::vector<std::reference_wrapper<ProcessVariable>>> process_variables;

            auto per_process_variables = findProcessVariables(variables, pv_config,
                {//! \ogs_file_param_special{prj__processes__process__PIPE_NETWORK__process_variables__flow_rate}
                "flow_rate",
                 //! \ogs_file_param_special{prj__processes__process__PIPE_NETWORK__process_variables__head}
                 "head" });
            process_variables.push_back(std::move(per_process_variables));

            // Process variable.

            // pipe diameter parameter.
            auto& pipe_diameter = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__PIPE_NETWORK__pipe_diameter}
                "pipe_diameter", parameters, 1);

            DBUG("Use \'%s\' as pipe diameter parameter.",
                pipe_diameter.name.c_str());

            // pipe roughness parameter.
            auto& pipe_roughness = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__PIPE_NETWORK__pipe_roughness}
                "pipe_roughness", parameters, 1);

            DBUG("Use \'%s\' as pipe roughness parameter.", pipe_roughness.name.c_str());

            // minorLoss coefficient parameter.
            auto& minorloss_coefficient = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__PIPE_NETWORK__minorLoss_coefficient}
                "minorloss_coefficient", parameters, 1);

            DBUG("Use \'%s\' as minorloss coefficient parameter.", minorloss_coefficient.name.c_str());

            PipeNetworkProcessData process_data{ pipe_diameter,
                pipe_roughness,
                minorloss_coefficient };

            SecondaryVariableCollection secondary_variables;

            NumLib::NamedFunctionCaller named_function_caller(
                { "PipeNetwork_flow_rate" }); // TODO be changed accordingly. 

            ProcessLib::createSecondaryVariables(config, secondary_variables,
                named_function_caller);

            return std::make_unique<PipeNetworkProcess>(
                mesh, std::move(jacobian_assembler), parameters, integration_order,
                std::move(process_variables), std::move(process_data),
                std::move(secondary_variables), std::move(named_function_caller));
        }

    }  // namespace PipeNetwork
}  // namespace ProcessLib
