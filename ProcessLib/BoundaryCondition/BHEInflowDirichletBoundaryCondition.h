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
        GlobalIndexType global_idx_T_in_top,
        GlobalIndexType global_idx_T_out_top,
        MeshLib::Mesh const& bc_mesh,
        std::vector<MeshLib::Node*> const& vec_inflow_bc_nodes,
        int const variable_id,
        unsigned const integration_order,
        std::size_t const bulk_mesh_id,
        int const component_id,
        std::unique_ptr<ProcessLib::HeatTransportBHE::BHE::BHEAbstract>& pt_bhe)
        : _variable_id(variable_id),
          _component_id(component_id),
          _bulk_mesh_id(bulk_mesh_id),
          _bc_mesh(bc_mesh),
          _integration_order(integration_order),
          _pt_bhe(pt_bhe)
    {

        DBUG(
            "Found %d nodes for BHE Inflow Dirichlet BCs for the variable %d "
            "and "
            "component %d",
            vec_inflow_bc_nodes.size(), variable_id, component_id);

        MeshLib::MeshSubset bc_mesh_subset{ _bc_mesh, vec_inflow_bc_nodes };

        // create memory to store Tout values
        _T_out_values.clear();
        _T_out_indices.clear();

        _bc_values.ids.clear();
        _bc_values.values.clear();

        // convert mesh node ids to global index for the given component
        _bc_values.ids.reserve(bc_mesh_subset.getNumberOfNodes());
        _bc_values.values.reserve(bc_mesh_subset.getNumberOfNodes());
        for (auto const* const node : bc_mesh_subset.getNodes())
        {
            // that might be slow, but only done once
            const auto g_idx_T_in = global_idx_T_in_top;
            const auto g_idx_T_out = global_idx_T_out_top;

            if (g_idx_T_in >= 0 && g_idx_T_out >= 0)
            {
                _T_out_indices.emplace_back(g_idx_T_out);
                _T_out_values.emplace_back(320.0 /*using initial value*/);
                _bc_values.ids.emplace_back(g_idx_T_in);
                _bc_values.values.emplace_back(320.0 /*using initial value*/);
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

    int const _variable_id;
    int const _component_id;

    MeshLib::Mesh const & _bc_mesh; 

    /// Integration order for integration over the lower-dimensional elements
    unsigned const _integration_order;

    std::size_t const _bulk_mesh_id;

    /// Stores the results of the outflow temperatures per boundary node.
    std::vector<double> _T_out_values;
    std::vector<GlobalIndexType> _T_out_indices;

    NumLib::IndexValueVector<GlobalIndexType> _bc_values;

    HeatTransportBHE::BHE::BHEAbstract* _BHE_property;

    std::unique_ptr<ProcessLib::HeatTransportBHE::BHE::BHEAbstract> const&
        _pt_bhe;
};

std::unique_ptr<BHEInflowDirichletBoundaryCondition>
createBHEInflowDirichletBoundaryCondition(
    GlobalIndexType global_idx_T_in_top, GlobalIndexType global_idx_T_out_top,
    MeshLib::Mesh const& bc_mesh,
    std::vector<MeshLib::Node*> const& vec_outflow_bc_nodes,
    int const variable_id, unsigned const integration_order,
    std::size_t const bulk_mesh_id, int const component_id,
    std::unique_ptr<ProcessLib::HeatTransportBHE::BHE::BHEAbstract>& pt_bhe);

}  // namespace ProcessLib
