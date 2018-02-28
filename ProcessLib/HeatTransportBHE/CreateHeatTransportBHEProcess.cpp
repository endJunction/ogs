/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#include "CreateHeatTransportBHEProcess.h"

#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "HeatTransportBHEProcess.h"
#include "HeatTransportBHEProcessData.h"
#include "BHE/BHEAbstract.h"
#include "BHE/BHE_1U.h"
#include "BHE/BHE_Net.h"
#include "BHE/CreateBHE1U.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/Density/createFluidDensityModel.h"
#include "MaterialLib/Fluid/Viscosity/createViscosityModel.h"
#include "MaterialLib/Fluid/SpecificHeatCapacity/CreateSpecificFluidHeatCapacityModel.h"
#include "MaterialLib/Fluid/ThermalConductivity/CreateFluidThermalConductivityModel.h"
#include "BaseLib/reorderVector.h"

namespace ProcessLib
{
    namespace HeatTransportBHE
    {
        std::unique_ptr<Process> createHeatTransportBHEProcess(
            MeshLib::Mesh& mesh,
            std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
            std::vector<ProcessVariable> const& variables,
            std::vector<std::unique_ptr<ParameterBase>> const& parameters,
            unsigned const integration_order,
            BaseLib::ConfigTree const& config,
            std::map<std::string, std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const& curves)
        {
            //! \ogs_file_param{prj__processes__process__type}
            config.checkConfigParameter("type", "HEAT_TRANSPORT_BHE");

            DBUG("Create HeatTransportBHE Process.");

            // Process variable.

            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__process_variables}
            auto const pv_config = config.getConfigSubtree("process_variables");

            std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
                process_variables;

            auto per_process_variables = findProcessVariables(
                variables, pv_config,
                {//! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__process_variables__process_variable}
                 "process_variable"});
            process_variables.push_back(std::move(per_process_variables));

            // solid phase thermal conductivity parameter.
            auto& thermal_conductivity_solid = findParameter<double>(
                config,
                //! \ogs_file_param_special prj__processes__process__HEAT_TRANSPORT_BHE__thermal_conductivity_solid}
                "thermal_conductivity_solid", parameters, 1);

            DBUG("Use \'%s\' as solid phase thermal conductivity parameter.",
                thermal_conductivity_solid.name.c_str());

            // solid phase thermal conductivity parameter.
            auto& thermal_conductivity_fluid = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__thermal_conductivity_fluid}
                "thermal_conductivity_fluid", parameters, 1);

            DBUG("Use \'%s\' as fluid phase thermal conductivity parameter.",
                thermal_conductivity_fluid.name.c_str());

            // gas phase thermal conductivity parameter.
            auto& thermal_conductivity_gas = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__thermal_conductivity_gas}
                "thermal_conductivity_gas", parameters, 1);

            DBUG("Use \'%s\' as gas phase thermal conductivity parameter.",
                thermal_conductivity_gas.name.c_str());

            // solid phase heat capacity parameter.
            auto& heat_capacity_solid = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__heat_capacity_solid}
                "heat_capacity_solid", parameters, 1);

            DBUG("Use \'%s\' as solid phase heat capacity parameter.", heat_capacity_solid.name.c_str());

            // fluid phase heat capacity parameter.
            auto& heat_capacity_fluid = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__heat_capacity_fluid}
                "heat_capacity_fluid", parameters, 1);

            DBUG("Use \'%s\' as fluid phase heat capacity parameter.", heat_capacity_fluid.name.c_str());

            // gas phase heat capacity parameter.
            auto& heat_capacity_gas = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__heat_capacity_gas}
                "heat_capacity_gas", parameters, 1);

            DBUG("Use \'%s\' as gas phase heat capacity parameter.", heat_capacity_gas.name.c_str());

            // solid phase density parameter.
            auto& density_solid = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__density_solid}
                "density_solid", parameters, 1);

            DBUG("Use \'%s\' as solid phase density parameter.", density_solid.name.c_str());

            // fluid phase density parameter.
            auto& density_fluid = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__density_fluid}
                "density_fluid", parameters, 1);

            DBUG("Use \'%s\' as fluid phase density parameter.", density_fluid.name.c_str());

            // gas phase density parameter.
            auto& density_gas = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__density_gas}
                "density_gas", parameters, 1);

            DBUG("Use \'%s\' as gas phase density parameter.", density_gas.name.c_str());

            // reading BHE parameters--------------------------------------------------------------
            std::vector<BHE::BHEAbstract*> vec_BHEs;
            BHE::BHE_Net BHE_network;

            // now read the BHE configurations
            auto const& bhe_configs =
                //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers}
                config.getConfigSubtree("borehole_heat_exchangers");

            for (
                auto const& bhe_conf :
                //! \ogs_file_param prj__processes__process___HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger}
                bhe_configs.getConfigSubtreeList("borehole_heat_exchanger"))
            {
                auto const bhe_id = bhe_conf.getConfigAttribute<int>("id");

                // read in the parameters
                using namespace BHE;
                const std::string bhe_type_str = bhe_conf.getConfigParameter<std::string>("bhe_type");

                // convert BHE type
                if (bhe_type_str == "BHE_TYPE_1U")
                {
                    // initialize the 1U type BHE
                    BHE::BHE_1U * m_bhe_1u = BHE::CreateBHE1U(config, bhe_conf, curves);

                    vec_BHEs.push_back(std::move(m_bhe_1u));
                    BHE_network.add_bhe_net_elem(m_bhe_1u);

                    // now adding a pipeline connecting the bottom of this BHE
                    BHE::BHE_Net_ELE_Pipe_Inner_1U * m_bhe_pipe_1u;
                    m_bhe_pipe_1u = new BHE::BHE_Net_ELE_Pipe_Inner_1U(m_bhe_1u->get_ele_name().append("_INNER_PIPE"), m_bhe_1u);
                    BHE_network.add_bhe_net_pipe(m_bhe_pipe_1u, m_bhe_1u->get_ele_name(), 0,
                        m_bhe_1u->get_ele_name(), 0);

                }
                else if (bhe_type_str == "BHE_TYPE_2U")
                {
                    // TODO
                }
                else if (bhe_type_str == "BHE_TYPE_CXC")
                {
                    // TODO
                }
                else if (bhe_type_str == "BHE_TYPE_CXA")
                {
                    // TODO
                }



            }
            // end of reading BHE parameters-------------------------------------------------------

            HeatTransportBHEProcessData process_data{ thermal_conductivity_solid,
                thermal_conductivity_fluid,
                thermal_conductivity_gas,
                heat_capacity_solid,
                heat_capacity_fluid,
                heat_capacity_gas,
                density_solid,
                density_fluid,
                density_gas };

            SecondaryVariableCollection secondary_variables;

            NumLib::NamedFunctionCaller named_function_caller(
                { "HeatConduction_temperature" });

            ProcessLib::createSecondaryVariables(config, secondary_variables,
                named_function_caller);

            return std::make_unique<HeatTransportBHEProcess>(
                mesh, std::move(jacobian_assembler), parameters, integration_order,
                std::move(process_variables), std::move(process_data),
                std::move(secondary_variables), std::move(named_function_caller));
        }

    }  // namespace HeatTransportBHE
}  // namespace ProcessLib
