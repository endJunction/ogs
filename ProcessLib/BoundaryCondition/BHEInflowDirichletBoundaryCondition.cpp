/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHEInflowDirichletBoundaryCondition.h"

#include <algorithm>
#include <logog/include/logog.hpp>
#include <vector>
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
void BHEInflowDirichletBoundaryCondition::getEssentialBCValues(
    const double /*t*/,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    SpatialPosition pos;

    bc_values.ids.clear();
    bc_values.values.clear();

    // convert mesh node ids to global index for the given component
    // TODO
    bc_values.ids.reserve(bc_values.ids.size() + _bc_values.ids.size());
    bc_values.values.reserve(bc_values.values.size() +
                             _bc_values.values.size());

    std::copy(_bc_values.ids.begin(), _bc_values.ids.end(),
              std::back_inserter(bc_values.ids));
    std::copy(_bc_values.values.begin(), _bc_values.values.end(),
              std::back_inserter(bc_values.values));
}

// update new values and corresponding indices.
void BHEInflowDirichletBoundaryCondition::preTimestep(const double /*t*/,
                                                      const GlobalVector& x)
{
    // for each BHE, the inflow temperature is dependent on
    // the ouflow temperature of the BHE.
    // Here the task is to get the outflow temperature and
    // save it locally
    double T_out(273.0), T_in_new(320.0);

    _bc_values.ids.clear();
    _bc_values.values.clear();

    auto const mesh_id = _bc_mesh.getID();
    auto const& nodes = _bc_mesh.getNodes();
    for (auto const* n : nodes)
    {
        std::size_t node_id = n->getID();
        MeshLib::Location l(mesh_id, MeshLib::MeshItemType::Node, node_id);
        // const auto g_idx =
        //      _dof_table_boundary.getGlobalIndex(l, _variable_id,
        //      _component_id);

        // read the T_out
        // TODO: notice that node_id + 1 is not always T_out!
        T_out = x[node_id + 1];

        // here calling the calculationg function in the
        // corresponding BHE class to get the updated
        // inflow temperature
        // T_in_new = // TODO
        // _bc_values.ids.emplace_back(g_idx);
        // _bc_values.values.emplace_back(T_in_new);
    }
}

std::unique_ptr<BHEInflowDirichletBoundaryCondition>
createBHEInflowDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    MeshLib::Mesh const& bc_mesh, int const variable_id,
    unsigned const integration_order, int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    DBUG(
        "Constructing BHEInflowDirichletBoundaryCondition from "
        "config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "BHEInflowDirichlet");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Dirichlet__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter %s", param_name.c_str());

    auto& param = findParameter<double>(param_name, parameters, 1);
    return std::make_unique<BHEInflowDirichletBoundaryCondition>(
        param, dof_table_bulk, bc_mesh, variable_id, integration_order,
        component_id);
}

}  // namespace ProcessLib
