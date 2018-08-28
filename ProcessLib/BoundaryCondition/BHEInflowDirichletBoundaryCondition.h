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
class BHEInflowDirichletBoundaryCondition final : public BoundaryCondition
{
public:
    BHEInflowDirichletBoundaryCondition(
        Parameter<double> const& parameter,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        MeshLib::Mesh const& bc_mesh,
        int const variable_id,
        unsigned const integration_order,
        std::size_t const bulk_mesh_id,
        int const component_id)
        : _parameter(parameter),
          _variable_id(variable_id),
          _component_id(component_id),
          _bulk_mesh_id(bulk_mesh_id),
          _bc_mesh(bc_mesh),
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

        std::vector<MeshLib::Node*> const& bc_nodes = _bc_mesh.getNodes();
        DBUG(
            "Found %d nodes for BHE Inflow Dirichlet BCs for the variable %d "
            "and "
            "component %d",
            bc_nodes.size(), variable_id, component_id);

        MeshLib::MeshSubset bc_mesh_subset{_bc_mesh, bc_nodes};

        // create memory to store Tout values
        _T_out_values.resize(bc_nodes.size());

        // Create local DOF table from intersected mesh subsets for the given
        // variable and component ids.
        _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
            variable_id, {component_id}, std::move(bc_mesh_subset)));

        SpatialPosition pos;

        auto const& bulk_node_ids_map =
            *_bc_mesh.getProperties().getPropertyVector<std::size_t>(
                "bulk_node_ids");

        _bc_values.ids.clear();
        _bc_values.values.clear();

        // convert mesh node ids to global index for the given component
        _bc_values.ids.reserve(_bc_mesh.getNumberOfNodes());
        _bc_values.values.reserve(_bc_mesh.getNumberOfNodes());
        for (auto const* const node : _bc_mesh.getNodes())
        {
            pos.setNodeID(node->getID());
            MeshLib::Location l(_bulk_mesh_id, MeshLib::MeshItemType::Node,
                                node->getID());
            // that might be slow, but only done once
            const auto g_idx = _dof_table_boundary->getGlobalIndex(
                l, _variable_id, _component_id);
            if (g_idx == NumLib::MeshComponentMap::nop)
                continue;
            // For the DDC approach (e.g. with PETSc option), the negative
            // index of g_idx means that the entry by that index is a ghost one,
            // which should be dropped. Especially for PETSc routines
            // MatZeroRows and MatZeroRowsColumns, which are called to apply the
            // Dirichlet BC, the negative index is not accepted like other
            // matrix or vector PETSc routines. Therefore, the following
            // if-condition is applied.
            if (g_idx >= 0)
            {
                _bc_values.ids.emplace_back(g_idx);
                _bc_values.values.emplace_back(
                    _parameter(0.0 /*using initial value*/, pos).front());
            }
        }
    }

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

    void preTimestep(const double t, const GlobalVector& x) override;

private:
    Parameter<double> const& _parameter;

    /// Local dof table, a subset of the global one restricted to the
    /// participating number of elements of the boundary condition.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table_boundary;

    int const _variable_id;
    int const _component_id;

    /// Vector of (lower-dimensional) boundary elements on which the boundary
    /// condition is defined.
    MeshLib::Mesh const& _bc_mesh;

    /// Integration order for integration over the lower-dimensional elements
    unsigned const _integration_order;

    std::size_t const _bulk_mesh_id;

    /// Stores the results of the outflow temperatures per boundary node.
    std::vector<double> _T_out_values;

    NumLib::IndexValueVector<GlobalIndexType> _bc_values;

    HeatTransportBHE::BHE::BHEAbstract* _BHE_property;
};

std::unique_ptr<BHEInflowDirichletBoundaryCondition>
createBHEInflowDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    MeshLib::Mesh const& bc_mesh, int const variable_id,
    unsigned const integration_order, std::size_t const bulk_mesh_id,
    int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters);

}  // namespace ProcessLib
