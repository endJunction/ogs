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
    auto const n_nodes = _bc_values.ids.size();
    for (size_t i = 0; i < n_nodes; i++)
    {
        auto g_idx = _T_in_values.ids[i];

        // read the T_out
        _T_in_values.values.at(i) = x[g_idx];
    }
}

std::unique_ptr<BHEBottomDirichletBoundaryCondition>
createBHEBottomDirichletBoundaryCondition(
    GlobalIndexType global_idx_T_in_bottom, 
    GlobalIndexType global_idx_T_out_bottom, 
    MeshLib::Mesh const& bulk_mesh,
    std::vector<MeshLib::Node*> const& vec_outflow_bc_nodes,
    int const variable_id,
    unsigned const integration_order, std::size_t const bulk_mesh_id,
    int const component_id, unsigned const bhe_id)
{
    DBUG(
        "Constructing BHEBottomDirichletBoundaryCondition from config.");


    return std::make_unique<BHEBottomDirichletBoundaryCondition>(
        global_idx_T_in_bottom, global_idx_T_out_bottom, bulk_mesh,
        vec_outflow_bc_nodes, variable_id,
        integration_order, bulk_mesh_id,
        component_id, bhe_id);
}

}  // namespace ProcessLib
