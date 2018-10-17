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
#include "ProcessLib/HeatTransportBHE/BHE/BHETypes.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace ProcessLib
{
template <typename BHEType>
class BHEInflowDirichletBoundaryCondition final : public BoundaryCondition
{
public:
    BHEInflowDirichletBoundaryCondition(
        std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
        MeshLib::Mesh const& bc_mesh,
        std::vector<MeshLib::Node*> const& vec_inflow_bc_nodes,
        int const variable_id,
        int const component_id,
        BHEType& bhe);

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

    void preTimestep(const double t, const GlobalVector& x) override;

private:
    // TODO (haibing) re-organize as the bottom BC data structure
    MeshLib::Mesh const& _bc_mesh;

    /// Stores the results of the outflow temperatures per boundary node.
    std::vector<double> _T_out_values;
    std::vector<GlobalIndexType> _T_out_indices;

    NumLib::IndexValueVector<GlobalIndexType> _bc_values;

    BHEType& _bhe;
};

template <typename BHEType>
std::unique_ptr<BHEInflowDirichletBoundaryCondition<BHEType>>
createBHEInflowDirichletBoundaryCondition(
    std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
    MeshLib::Mesh const& bc_mesh,
    std::vector<MeshLib::Node*> const& vec_inflow_bc_nodes,
    int const variable_id, int const component_id, BHEType& bhe)
{
    DBUG("Constructing BHEInflowDirichletBoundaryCondition.");

    return std::make_unique<BHEInflowDirichletBoundaryCondition<BHEType>>(
        std::move(in_out_global_indices), bc_mesh, vec_inflow_bc_nodes,
        variable_id, component_id, bhe);
}
}  // namespace ProcessLib

#include "BHEInflowDirichletBoundaryCondition-impl.h"
