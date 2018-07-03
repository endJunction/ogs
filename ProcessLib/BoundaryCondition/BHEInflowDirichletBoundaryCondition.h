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
        int const component_id)
        : _parameter(parameter),
          _variable_id(variable_id),
          _component_id(component_id),
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
            "Found %d nodes for constraint Dirichlet BCs for the variable %d "
            "and "
            "component %d",
            bc_nodes.size(), variable_id, component_id);

        MeshLib::MeshSubset bc_mesh_subset{_bc_mesh, bc_nodes};

        // Create local DOF table from intersected mesh subsets for the given
        // variable and component ids.
        _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
            variable_id, {component_id}, std::move(bc_mesh_subset)));
    }

    void getEssentialBCValues(
        const double t,
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

    /// Stores the results of the outflow temperatures per boundary node.
    std::vector<double> _outflow_temperature;

    NumLib::IndexValueVector<GlobalIndexType> _bc_values;
};

std::unique_ptr<BHEInflowDirichletBoundaryCondition>
createBHEInflowDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    MeshLib::Mesh const& bc_mesh, int const variable_id,
    unsigned const integration_order, int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters);

}  // namespace ProcessLib
