/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <utility>
#include <boost/optional.hpp>
#include <string>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}
namespace ParameterLib
{
struct ParameterBase;
}

namespace ProcessLib
{
struct NonequilibriumInitialState;
}

namespace ProcessLib
{
std::unique_ptr<NonequilibriumInitialState> createNonequilibriumInitialState(
    boost::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::pair<std::string, int>> const& tag_names_and_components,
    // std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    MeshLib::Mesh const& mesh);
}  // namespace ProcessLib
