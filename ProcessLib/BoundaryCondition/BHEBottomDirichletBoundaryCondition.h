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
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        std::vector<MeshLib::Node*> const& bc_bottom_nodes,
        MeshLib::Mesh const& bulk_mesh,
        int const variable_id,
        unsigned const integration_order,
        std::size_t const bulk_mesh_id,
        int const component_id)
        : _variable_id(variable_id),
          _component_id(component_id),
          _bulk_mesh_id(bulk_mesh_id),
          _bulk_mesh(bulk_mesh),
          _integration_order(integration_order)
    {
        if (variable_id >=
                static_cast<int>(dof_table_bulk.getNumberOfVariables()) ||
            component_id >=
                dof_table_bulk.getNumberOfVariableComponents(variable_id))
        {
            OGS_FATAL(
                "Variable id or component id too high. Actual values: (%d, "
                "%d), "
                "maximum values: (%d, %d).",
                variable_id, component_id,
                dof_table_bulk.getNumberOfVariables(),
                dof_table_bulk.getNumberOfVariableComponents(variable_id));
        }

        DBUG(
            "Found %d nodes for BHE bottom Dirichlet BCs for the variable %d "
            "and "
            "component %d",
            bc_bottom_nodes.size(), variable_id, component_id);

        _bc_mesh = bc_mesh_subset.getMesh();

        // create memory to store Tout values
        _T_in_values.ids.clear();
        _T_in_values.values.clear();

        // Create local DOF table from intersected mesh subsets for the given
        // variable and component ids.
        _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
            variable_id, {component_id}, std::move(bc_mesh_subset)));

        // get the right T_in component dof_table
        int const component_id_T_in = 0;
        _dof_table_boundary_T_in.reset(
            dof_table_bulk.deriveBoundaryConstrainedMap(
                variable_id, {component_id_T_in}, std::move(bc_mesh_subset)));

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
            const auto g_T_out_idx = _dof_table_boundary->getGlobalIndex(
                l, _variable_id, _component_id);
            if (g_T_out_idx == NumLib::MeshComponentMap::nop)
                continue;
            // For the DDC approach (e.g. with PETSc option), the negative
            // index of g_idx means that the entry by that index is a ghost one,
            // which should be dropped. Especially for PETSc routines
            // MatZeroRows and MatZeroRowsColumns, which are called to apply the
            // Dirichlet BC, the negative index is not accepted like other
            // matrix or vector PETSc routines. Therefore, the following
            // if-condition is applied.
            if (g_T_out_idx >= 0)
            {
                _bc_values.ids.emplace_back(g_T_out_idx);
                _bc_values.values.emplace_back(298.15);
            }

            const auto g_T_in_idx = _dof_table_boundary_T_in->getGlobalIndex(
                l, _variable_id, component_id_T_in);
            if (g_T_in_idx == NumLib::MeshComponentMap::nop)
                continue;

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

    int const _variable_id;
    int const _component_id;

    /// Vector of (lower-dimensional) boundary elements on which the boundary
    /// condition is defined.
    MeshLib::Mesh& _bc_mesh;

    /// Integration order for integration over the lower-dimensional elements
    unsigned const _integration_order;

    std::size_t const _bulk_mesh_id;

    NumLib::IndexValueVector<GlobalIndexType> _bc_values;

    NumLib::IndexValueVector<GlobalIndexType> _T_in_values;

    HeatTransportBHE::BHE::BHEAbstract* _BHE_property;
};

std::unique_ptr<BHEBottomDirichletBoundaryCondition>
createBHEBottomDirichletBoundaryCondition(
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    std::vector<MeshLib::Node*> const& bc_inlet_nodes,
    MeshLib::Mesh const& bulk_mesh, int const variable_id,
    unsigned const integration_order, std::size_t const bulk_mesh_id,
    int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters);

}  // namespace ProcessLib
