/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHEBottomDirichletBoundaryCondition.h"

#include <algorithm>
#include <logog/include/logog.hpp>
#include <vector>
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
void BHEBottomDirichletBoundaryCondition::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    SpatialPosition pos;

    bc_values.ids.clear();
    bc_values.values.clear();

    bc_values.ids.resize(_bc_values.ids.size());
    bc_values.values.resize(_bc_values.values.size());

    const size_t n_nodes = _T_in_values.ids.size();
    double tmp_T_out(300.0); 
    for (size_t i=0; i<n_nodes; i++)
    {
        bc_values.ids[i] = _bc_values.ids[i];
        // here, the outflow temperature is always
        // the same as the inflow temperature
        // get the inflow temperature from here.
        tmp_T_out = x[_T_in_values.ids.at(i)];
        bc_values.values[i] = tmp_T_out;

    }
}

// update new values and corresponding indices.
void BHEBottomDirichletBoundaryCondition::preTimestep(
    const double t, const GlobalVector& x)
{
    // At the bottom of each BHE, the outflow temperature
    // is the same as the inflow temperature. 
    // Here the task is to get the inflow temperature and
    // save it locally
    auto const mesh_id = _bc_mesh.getID();
    auto const& nodes = _bc_mesh.getNodes();
    for (auto const* n : nodes)
    {
        std::size_t node_id = n->getID();
        auto g_idx = _T_in_values.ids.at(node_id);

        // read the T_in
        _T_in_values.values.at(node_id) = x[g_idx];
    }
	
}

std::unique_ptr<BHEBottomDirichletBoundaryCondition>
createBHEBottomDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    MeshLib::Mesh const& bc_mesh, int const variable_id,
    unsigned const integration_order, std::size_t const bulk_mesh_id,
    int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    DBUG(
        "Constructing BHEBottomDirichletBoundaryCondition from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter(
        "type", "BHEBottomDirichlet");


    // find the corresponding BHE id
    auto const bhe_id = config.getConfigParameter<std::size_t>("bhe_id");

    // TODO: try to locate the BHE instance here
    // if BHE instance cannot be found,
    // then complain about it.

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Dirichlet__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter %s", param_name.c_str());

    auto& param = findParameter<double>(param_name, parameters, 1);
    return std::make_unique<BHEBottomDirichletBoundaryCondition>(
        param, dof_table_bulk, bc_mesh, variable_id, integration_order,
        bulk_mesh_id, component_id);
}

}  // namespace ProcessLib
