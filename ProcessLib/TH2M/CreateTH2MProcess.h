/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>
#include <map>

namespace BaseLib
{
class ConfigTree;
}
namespace MeshLib
{
class Mesh;
}
namespace ProcessLib
{
class AbstractJacobianAssembler;
struct ParameterBase;
class Process;
class ProcessVariable;
}  // namespace ProcessLib

namespace MaterialPropertyLib
{
class Medium;
}

namespace ProcessLib
{
namespace TH2M
{
template <int DisplacementDim>
std::unique_ptr<Process> createTH2MProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::unique_ptr<MaterialPropertyLib::Medium>> const& media);

extern template std::unique_ptr<Process> createTH2MProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::unique_ptr<MaterialPropertyLib::Medium>> const& media);

extern template std::unique_ptr<Process> createTH2MProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::unique_ptr<MaterialPropertyLib::Medium>> const& media);
}  // namespace TH2M
}  // namespace ProcessLib
