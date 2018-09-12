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
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    SpatialPosition pos;

    // auto const& bulk_node_ids_map =
    //     *_bc_mesh.getProperties().getPropertyVector<std::size_t>(
    //         "bulk_node_ids");

    bc_values.ids.clear();
    bc_values.values.clear();

    bc_values.ids.resize(_bc_values.ids.size());
    bc_values.values.resize(_bc_values.values.size());

    const size_t n_nodes = _T_out_values.size();
    double tmp_T_in(320.0); 
    for (size_t i=0; i<n_nodes; i++)
    {
        bc_values.ids[i] = _bc_values.ids[i];
        // here call the corresponding BHE functions
        auto tmp_T_out = x[_T_out_indices[i]];
        tmp_T_in = _pt_bhe->get_Tin_by_Tout(tmp_T_out, t);
        bc_values.values[i] = tmp_T_in;
    }

}

// update new values and corresponding indices.
void BHEInflowDirichletBoundaryCondition::preTimestep(
    const double t, const GlobalVector& x)
{
    // for each BHE, the inflow temperature is dependent on
    // the ouflow temperature of the BHE. 
    // Here the task is to get the outflow temperature and
    // save it locally
    auto const n_nodes = _bc_values.ids.size();
    for (size_t i = 0; i < n_nodes; i++)
    {
        auto g_idx = _T_out_indices[i];

        // read the T_out
        _T_out_values[i] = x[_T_out_indices[i]];
    }
}

std::unique_ptr<BHEInflowDirichletBoundaryCondition>
createBHEInflowDirichletBoundaryCondition(
    GlobalIndexType global_idx_T_in_top, GlobalIndexType global_idx_T_out_top,
    MeshLib::Mesh const& bc_mesh,
    std::vector<MeshLib::Node*> const& vec_inflow_bc_nodes,
    int const variable_id, unsigned const integration_order,
    std::size_t const bulk_mesh_id, int const component_id,
    std::unique_ptr<ProcessLib::HeatTransportBHE::BHE::BHEAbstract> const&
        pt_bhe)
{
    DBUG(
        "Constructing BHEInflowDirichletBoundaryCondition.");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Dirichlet__parameter}

    return std::make_unique<BHEInflowDirichletBoundaryCondition>(
        global_idx_T_in_top, global_idx_T_out_top,
        bc_mesh, vec_inflow_bc_nodes,
        variable_id, integration_order, bulk_mesh_id, component_id, pt_bhe);
}

}  // namespace ProcessLib
