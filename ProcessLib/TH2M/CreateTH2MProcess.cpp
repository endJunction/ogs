/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateTH2MProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"

#include "TH2MProcess.h"
#include "TH2MProcessData.h"

namespace ProcessLib
{
namespace TH2M
{
void validateTH2MProcessVariables(
    std::vector<std::reference_wrapper<ProcessVariable>> process_variables,
    int DisplacementDim)
{
    if (process_variables.size() != 4)
    {
        OGS_FATAL(
            "Wrong number of process variables found. TH2M requires four"
            " primary unknowns: gas pressure, capilalry pressure,"
            " temperature, and displacement!");
    }

    std::array<int, 4> const number_variable_components = {1, 1, 1,
                                                           DisplacementDim};
    const std::string process_variable_name[4] = {
        "gas pressure", "capillary pressure", "temperature", "displacement"};

    for (int i = indexGasPressure; i <= indexDisplacement; i++)
    {
        DBUG("Associate %s with process variable \'%s\'.",
             process_variable_name[i].c_str(),
             process_variables[i].get().getName().c_str());

        if (process_variables[i].get().getNumberOfComponents() !=
            number_variable_components[i])
        {
            OGS_FATAL(
                "Number of components of the process variable '%s' is "
                "different from the displacement dimension: got %d, expected "
                "%d",
                process_variables[i].get().getName().c_str(),
                process_variables[i].get().getNumberOfComponents(),
                number_variable_components[i]);
        }
    }
}

template <int DisplacementDim>
std::unique_ptr<Process> createTH2MProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    MaterialPropertyLib::Medium& medium)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "TH2M");
    DBUG("Create TH2MProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__TH2M__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__TH2M__process_variables__gas_pressure}
         "gas_pressure",
         //! \ogs_file_param_special{prj__processes__process__TH2M__process_variables__capillary_pressure}
         "capillary_pressure",
         //! \ogs_file_param_special{prj__processes__process__TH2M__process_variables__temperature}
         "temperature",
         //! \ogs_file_param_special{prj__processes__process__TH2M__process_variables__displacement}
         "displacement"});

    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;

    process_variables.push_back(std::move(per_process_variables));

    validateTH2MProcessVariables(process_variables[0], DisplacementDim);

    // Constitutive relation.
    // read type;
    auto const constitutive_relation_config =
        //! \ogs_file_param{prj__processes__process__TH2M__constitutive_relation}
        config.getConfigSubtree("constitutive_relation");

    auto const type =
        //! \ogs_file_param{prj__processes__process__TH2M__constitutive_relation__type}
        constitutive_relation_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material = nullptr;
    if (type == "LinearElasticIsotropic")
    {
        material =
            MaterialLib::Solids::createLinearElasticIsotropic<DisplacementDim>(
                parameters, constitutive_relation_config);
    }
    else
    {
        OGS_FATAL(
            "Cannot construct constitutive relation of given type \'%s\'.",
            type.c_str());
    }

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__TH2M__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != DisplacementDim)
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is %d, displacement "
                "dimension is %d",
                b.size(), DisplacementDim);

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    TH2MProcessData<DisplacementDim> process_data{std::move(material),
                                                  specific_body_force, medium};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"TH2M_gas_pressure", "TH2M_capillary_pressure", "TH2M_temperature",
         "TH2M_displacement"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<TH2MProcess<DisplacementDim>>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller));
}

template std::unique_ptr<Process> createTH2MProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    MaterialPropertyLib::Medium& medium);

template std::unique_ptr<Process> createTH2MProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    MaterialPropertyLib::Medium& medium);

}  // namespace TH2M
}  // namespace ProcessLib
