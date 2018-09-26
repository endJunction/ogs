/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BoundaryCondition.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHEAbstract.h"

namespace ProcessLib
{
class BHEBottomDirichletBoundaryCondition final : public BoundaryCondition
{
public:
    BHEBottomDirichletBoundaryCondition(
        GlobalIndexType global_idx_T_in_bottom,
        GlobalIndexType global_idx_T_out_bottom, 
        MeshLib::Mesh const& bulk_mesh,
        std::vector<MeshLib::Node*> const& vec_outflow_bc_nodes, 
        int const variable_id,
        unsigned const integration_order,
        std::size_t const bulk_mesh_id,
        int const component_id,
        unsigned const bhe_idx)        
        : _bulk_mesh_id(bulk_mesh_id),
          _bulk_mesh(bulk_mesh)
    {
        DBUG(
            "Found %d nodes for BHE bottom Dirichlet BCs for the variable %d "
            "and "
            "component %d",
            vec_outflow_bc_nodes.size(), variable_id, component_id);

        MeshLib::MeshSubset bc_mesh_subset{ _bulk_mesh, vec_outflow_bc_nodes };

        // create memory to store Tout values
        _T_in_values.ids.clear();
        _T_in_values.values.clear();

        SpatialPosition pos;

        _bc_values.ids.clear();
        _bc_values.values.clear();

        // convert mesh node ids to global index for the given component
        _bc_values.ids.reserve(bc_mesh_subset.getNumberOfNodes());
        _bc_values.values.reserve(bc_mesh_subset.getNumberOfNodes());
        for (auto const* const node : bc_mesh_subset.getNodes())
        {
            pos.setNodeID(node->getID());
            MeshLib::Location l(_bulk_mesh_id, MeshLib::MeshItemType::Node,
                                node->getID());
            // that might be slow, but only done once
            const auto g_T_out_idx = global_idx_T_out_bottom; 
            if (g_T_out_idx >= 0)
            {
                _bc_values.ids.emplace_back(g_T_out_idx);
                _bc_values.values.emplace_back(298.15);
            }

            const auto g_T_in_idx = global_idx_T_in_bottom; 

            if (g_T_in_idx >= 0)
            {
                _T_in_values.ids.emplace_back(g_T_in_idx);
                _T_in_values.values.emplace_back(298.15);
            }
        }
    }

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

    void preTimestep(const double t, const GlobalVector& x) override;

private:
    /// Local dof table, a subset of the global one restricted to the
    /// participating number of elements of the boundary condition.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table_boundary;

    /// Local dof table, a subset of the global one restricted to the
    /// participating number of elements of the boundary condition.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table_boundary_T_in;

    /// id of bulk mesh
    std::size_t const _bulk_mesh_id;

    /// the bulk mesh
    MeshLib::Mesh const& _bulk_mesh;

    /// Integration order for integration over the lower-dimensional elements
    // unsigned const _integration_order;

    NumLib::IndexValueVector<GlobalIndexType> _bc_values;

    NumLib::IndexValueVector<GlobalIndexType> _T_in_values;
};

std::unique_ptr<BHEBottomDirichletBoundaryCondition>
createBHEBottomDirichletBoundaryCondition(
    GlobalIndexType global_idx_T_in_bottom, 
    GlobalIndexType global_idx_T_out_bottom, 
    MeshLib::Mesh const& bulk_mesh,
    std::vector<MeshLib::Node*> const& vec_outflow_bc_nodes,
    int const variable_id,
    unsigned const integration_order, std::size_t const bulk_mesh_id,
    int const component_id, unsigned const bhe_id);

}  // namespace ProcessLib
