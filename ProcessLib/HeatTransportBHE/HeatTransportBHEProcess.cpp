/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HeatTransportBHEProcess.h"

#include <cassert>
#include <iostream>
#include "ProcessLib/HeatTransportBHE/BHE/MeshUtils.h"
#include "ProcessLib/HeatTransportBHE/LocalAssemblers/CreateLocalAssemblers.h"

#include "ProcessLib/HeatTransportBHE/LocalAssemblers/HeatTransportBHELocalAssemblerBHE.h"
#include "ProcessLib/HeatTransportBHE/LocalAssemblers/HeatTransportBHELocalAssemblerSoil.h"

#include "ProcessLib/BoundaryCondition/BHEBottomDirichletBoundaryCondition.h"
#include "ProcessLib/BoundaryCondition/BHEInflowDirichletBoundaryCondition.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
HeatTransportBHEProcess::HeatTransportBHEProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    HeatTransportBHEProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
    getBHEDataInMesh(mesh,
                     _vec_pure_soil_elements,
                     _vec_BHE_mat_IDs,
                     _vec_BHE_elements,
                     _vec_pure_soil_nodes,
                     _vec_BHE_nodes);

    if (_vec_BHE_mat_IDs.size() != _process_data._vec_BHE_property.size())
    {
        OGS_FATAL(
            "The number of the given BHE properties (%d) are not "
            "consistent"
            " with the number of BHE groups in a mesh (%d).",
            _process_data._vec_BHE_property.size(),
            _vec_BHE_mat_IDs.size());
    }

    // create a map from a material ID to a BHE ID
    auto max_BHE_mat_id =
        std::max_element(_vec_BHE_mat_IDs.begin(), _vec_BHE_mat_IDs.end());
    const std::size_t nBHEmatIDs = *max_BHE_mat_id + 1;
    _process_data._map_materialID_to_BHE_ID.resize(nBHEmatIDs);
    for (std::size_t i = 0; i < _vec_BHE_mat_IDs.size(); i++)
    {
        // by default, it is assumed that the soil compartment takes material ID
        // 0 and the BHE take the successive material group.
        _process_data._map_materialID_to_BHE_ID[_vec_BHE_mat_IDs[i]] = i;
    }

    // TODO check if materialIDs are there or not.
    MeshLib::PropertyVector<int> const* material_ids(
        mesh.getProperties().getPropertyVector<int>("MaterialIDs"));
    _process_data._mesh_prop_materialIDs = material_ids;
}

void HeatTransportBHEProcess::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());

    //------------------------------------------------------------
    // prepare mesh subsets to define DoFs
    //------------------------------------------------------------
    // first all the soil nodes
    _mesh_subset_pure_soil_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _vec_pure_soil_nodes);

    std::vector<int> vec_n_BHE_unknowns;
    vec_n_BHE_unknowns.reserve(_vec_BHE_nodes.size());
    // the BHE nodes need to be cherry-picked from the vector
    for (unsigned i = 0; i < _vec_BHE_nodes.size(); i++)
    {
        _mesh_subset_BHE_nodes.push_back(
            std::make_unique<MeshLib::MeshSubset const>(_mesh,
                                                        _vec_BHE_nodes[i]));
        int n_BHE_unknowns =
            _process_data._vec_BHE_property[i]->getNumUnknowns();
        vec_n_BHE_unknowns.emplace_back(n_BHE_unknowns);
    }

    // Collect the mesh subsets in a vector.
    // All the soil nodes has 1 temperature variable
    std::vector<MeshLib::MeshSubset> all_mesh_subsets{
        *_mesh_subset_pure_soil_nodes};

    // All the BHE nodes have additinal variables
    std::size_t count = 0;
    for (auto& ms : _mesh_subset_BHE_nodes)
    {
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        // Here the number of components equals to
                        // the number of unknowns on the BHE
                        vec_n_BHE_unknowns[count],
                        [&]() { return *ms; });
        count++;
    }

    std::vector<int> vec_n_components;
    // This is the soil temperature for first mesh subset.
    // 1, because for the soil part there is just one variable which is the soil
    // temperature.
    vec_n_components.push_back(1);
    // now the BHE subsets
    for (std::size_t i = 0; i < _vec_BHE_mat_IDs.size(); i++)
    {
        // Here the number of components equals to
        // the number of unknowns on the BHE
        vec_n_components.push_back(vec_n_BHE_unknowns[i]);
    }

    std::vector<std::vector<MeshLib::Element*> const*> vec_var_elements;
    // vec_var_elements.push_back(&_vec_pure_soil_elements);
    vec_var_elements.push_back(&(_mesh.getElements()));
    for (std::size_t i = 0; i < _vec_BHE_elements.size(); i++)
    {
        vec_var_elements.push_back(&_vec_BHE_elements[i]);
    }

    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets),
            vec_n_components,
            vec_var_elements,
            NumLib::ComponentOrder::BY_COMPONENT);

    // in case of debugging the dof table, activate the following line
    // std::cout << *_local_to_global_index_map << "\n";
}

void HeatTransportBHEProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    const int process_id = 0;

    ProcessLib::HeatTransportBHE::createLocalAssemblers<
        3, /*mesh.getDimension(),*/
        HeatTransportBHELocalAssemblerSoil, HeatTransportBHELocalAssemblerBHE>(
        mesh.getElements(), dof_table, _local_assemblers,
        _process_data._vec_ele_connected_BHE_IDs,
        _process_data._vec_BHE_property, mesh.isAxiallySymmetric(),
        integration_order, _process_data);

    // create BHE boundary conditions
    // for each BHE, one BC on the top
    // and one BC at the bottom.
    std::vector<std::unique_ptr<BoundaryCondition>> bc_collections =
        createBHEBoundaryConditionTopBottom();
    auto& current_process_BCs = _boundary_conditions[process_id];
    auto const bc_collections_size = bc_collections.size();
    for (unsigned i = 0; i < bc_collections_size; i++)
        current_process_BCs.addCreatedBC(std::move(bc_collections[i]));
}

void HeatTransportBHEProcess::assembleConcreteProcess(const double t,
                                                      GlobalVector const& x,
                                                      GlobalMatrix& M,
                                                      GlobalMatrix& K,
                                                      GlobalVector& b)
{
    DBUG("Assemble HeatTransportBHE process.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, x, M, K, b, _coupled_solutions);
}

void HeatTransportBHEProcess::assembleWithJacobianConcreteProcess(
    const double /*t*/, GlobalVector const& /*x*/, GlobalVector const& /*xdot*/,
    const double /*dxdot_dx*/, const double /*dx_dx*/, GlobalMatrix& /*M*/,
    GlobalMatrix& /*K*/, GlobalVector& /*b*/, GlobalMatrix& /*Jac*/)
{
    OGS_FATAL(
        "HeatTransportBHE: analytical Jacobian assembly is not implemented");
}

void HeatTransportBHEProcess::computeSecondaryVariableConcrete(
    const double t, GlobalVector const& x)
{
    DBUG("Compute heat flux for HeatTransportBHE process.");
    GlobalExecutor::executeMemberOnDereferenced(
        &HeatTransportBHELocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, *_local_to_global_index_map, t, x,
        _coupled_solutions);
}

std::vector<std::unique_ptr<BoundaryCondition>>
HeatTransportBHEProcess::createBHEBoundaryConditionTopBottom()
{
    /**
     * boundary conditions
     */
    std::vector<std::unique_ptr<BoundaryCondition>> bcs;

    int const n_BHEs = static_cast<int>(_process_data._vec_BHE_property.size());

    // bulk dof table
    NumLib::LocalToGlobalIndexMap dof_table_bulk = *_local_to_global_index_map;

    std::vector<MeshLib::Node*> bc_top_nodes;
    std::vector<MeshLib::Node*> bc_bottom_nodes;

    // for each BHE
    for (int bhe_i = 0; bhe_i < n_BHEs; bhe_i++)
    {
        // find the variable ID
        // the soil temperature is 0-th variable
        // the BHE temperature is therefore bhe_i + 1
        const int variable_id = bhe_i + 1;
        // find the node in mesh that are at the top
        auto const n_bhe_nodes = _vec_BHE_nodes[bhe_i].size();
        unsigned int idx_top = 0;
        unsigned int idx_bottom = _vec_BHE_nodes[bhe_i].size() - 1;
        double top_z_coord = _vec_BHE_nodes[bhe_i][idx_top]->getCoords()[2];
        double bottom_z_coord =
            _vec_BHE_nodes[bhe_i][idx_bottom]->getCoords()[2];
        for (unsigned bhe_node_i = 0; bhe_node_i < n_bhe_nodes; bhe_node_i++)
        {
            if (_vec_BHE_nodes[bhe_i][bhe_node_i]->getCoords()[2] >=
                top_z_coord)
                idx_top = bhe_node_i;

            if (_vec_BHE_nodes[bhe_i][bhe_node_i]->getCoords()[2] <=
                bottom_z_coord)
                idx_bottom = bhe_node_i;
        }
        bc_top_nodes.clear();
        bc_top_nodes.emplace_back(_vec_BHE_nodes[bhe_i][idx_top]);
        bc_bottom_nodes.clear();
        bc_bottom_nodes.emplace_back(_vec_BHE_nodes[bhe_i][idx_bottom]);

        MeshLib::MeshSubset bc_mesh_subset_top{_mesh, bc_top_nodes};
        MeshLib::Mesh const& bc_mesh_top = bc_mesh_subset_top.getMesh();

        MeshLib::MeshSubset bc_mesh_subset_bottom{_mesh, bc_bottom_nodes};
        MeshLib::Mesh const& bc_mesh_bottom = bc_mesh_subset_bottom.getMesh();

        auto get_global_bhe_bc_index = [&](std::size_t const mesh_id,
                                           std::size_t const node_id,
                                           int const component_id) {
            return dof_table_bulk.getGlobalIndex(
                {mesh_id, MeshLib::MeshItemType::Node, node_id}, variable_id,
                component_id);
        };

        // the create_BC function will be repeatedly used
        auto create_BC_top_inflow = [&](std::size_t const mesh_id,
                                        int const component_id_in,
                                        int const component_id_out) {
            auto const top_node_id = bc_top_nodes[0]->getID();
            auto const global_index_in =
                get_global_bhe_bc_index(mesh_id, top_node_id, component_id_in);
            auto const global_index_out =
                get_global_bhe_bc_index(mesh_id, top_node_id, component_id_out);
            return ProcessLib::createBHEInflowDirichletBoundaryCondition(
                global_index_in, global_index_out, _mesh, bc_top_nodes,
                variable_id, component_id_in,
                _process_data._vec_BHE_property[bhe_i]);
        };
        auto create_BC_bottom_outflow = [&](std::size_t const mesh_id,
                                            int const component_id_in,
                                            int const component_id_out) {
            auto const bottom_node_id = bc_bottom_nodes[0]->getID();
            auto const global_index_in = get_global_bhe_bc_index(
                mesh_id, bottom_node_id, component_id_in);
            auto const global_index_out = get_global_bhe_bc_index(
                mesh_id, bottom_node_id, component_id_out);
            return ProcessLib::createBHEBottomDirichletBoundaryCondition(
                global_index_in, global_index_out, _mesh, bc_top_nodes,
                variable_id, component_id_in);
        };

        auto create_top_bottom_bhe_BCs =
            [&](std::pair<int, int> const& in_out_component_id,
                auto bcs_output_iterator) {
                *bcs_output_iterator++ = create_BC_top_inflow(
                    bc_mesh_top.getID(), in_out_component_id.first,
                    in_out_component_id.second);
                *bcs_output_iterator++ = create_BC_bottom_outflow(
                    bc_mesh_bottom.getID(), in_out_component_id.first,
                    in_out_component_id.second);
            };
        for (auto const& in_out_component_id :
             _process_data._vec_BHE_property[bhe_i]
                 ->inflowOutflowBcComponentIds())
        {
            create_top_bottom_bhe_BCs(in_out_component_id,
                                      std::back_inserter(bcs));
        }
    }

    return bcs;
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
