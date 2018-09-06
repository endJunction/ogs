/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HeatTransportBHEProcess.h"

#include <cassert>

#include "ProcessLib/HeatTransportBHE/BHE/MeshUtils.h"
#include "ProcessLib/HeatTransportBHE/LocalAssemblers/CreateLocalAssemblers.h"

#include "ProcessLib/HeatTransportBHE/LocalAssemblers/HeatTransportBHELocalAssemblerBHE.h"
#include "ProcessLib/HeatTransportBHE/LocalAssemblers/HeatTransportBHELocalAssemblerSoil.h"

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
    _process_data._map_materialID_to_BHE_ID.resize(*max_BHE_mat_id + 1);
    for (unsigned i = 0; i < _vec_BHE_mat_IDs.size(); i++)
    {
        // by default, it is assumed that the soil compartment takes material ID
        // 0 and the BHE take the successive material group.
        _process_data._map_materialID_to_BHE_ID[_vec_BHE_mat_IDs[i]] = i;
    }

    /*
    // create a table of connected BHE IDs for each element
    _process_data._vec_ele_connected_BHE_IDs.resize(mesh.getNumberOfElements());
    for (unsigned i = 0; i < _vec_BHE_soil_elements.size(); i++)
    {
        for (auto e : _vec_BHE_soil_elements[i])
        {
            _process_data._vec_ele_connected_BHE_IDs[e->getID()].push_back(i);
        }
    }
    */

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

    std::vector<size_t> vec_n_BHE_unknowns;
    vec_n_BHE_unknowns.reserve(_vec_BHE_nodes.size());
    // the BHE nodes need to be cherry-picked from the vector
    for (unsigned i = 0; i < _vec_BHE_nodes.size(); i++)
    {
        _mesh_subset_BHE_nodes.push_back(
            std::make_unique<MeshLib::MeshSubset const>(_mesh,
                                                        _vec_BHE_nodes[i]));
        size_t n_BHE_unknowns =
            _process_data._vec_BHE_property[i]->get_n_unknowns();
        vec_n_BHE_unknowns.emplace_back(n_BHE_unknowns);
    }

    // Collect the mesh subsets in a vector.
    // All the soil nodes has 1 temperature variable
    std::vector<MeshLib::MeshSubset> all_mesh_subsets{
        *_mesh_subset_pure_soil_nodes};

    // All the BHE nodes have additinal variables
    int count = 0;
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
    // this is the soil temperature for first mesh subset
    // 1 because for the soil part ther is just one var which is the soile
    // temperatrure
    vec_n_components.push_back(1); 
    // now the BHE subsets
    for (unsigned i = 0; i < _vec_BHE_mat_IDs.size(); i++)
    {
        // Here the number of components equals to
        // the number of unknowns on the BHE
        vec_n_components.push_back(vec_n_BHE_unknowns[i]);
    }

    std::vector<std::vector<MeshLib::Element*> const*> vec_var_elements;
    // vec_var_elements.push_back(&_vec_pure_soil_elements);
    vec_var_elements.push_back(&(_mesh.getElements()));
    for (unsigned i = 0; i < _vec_BHE_elements.size(); i++)
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
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // this process can only run with 3-dimensional mesh
    ProcessLib::HeatTransportBHE::createLocalAssemblers<
        3, /*mesh.getDimension(),*/
        HeatTransportBHELocalAssemblerSoil,
        HeatTransportBHELocalAssemblerBHE>(
        mesh.getElements(), dof_table, _local_assemblers,
        _process_data._vec_ele_connected_BHE_IDs, 
        _process_data._vec_BHE_property, mesh.isAxiallySymmetric(),
        integration_order, _process_data);

    /*
    _secondary_variables.addSecondaryVariable(
        "heat_flux_x", makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &HeatTransportBHELocalAssemblerInterface::getIntPtHeatFluxX));

    if (mesh.getDimension() > 1)
    {
        _secondary_variables.addSecondaryVariable(
            "heat_flux_y", makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
                &HeatTransportBHELocalAssemblerInterface::getIntPtHeatFluxY));
    }
    if (mesh.getDimension() > 2)
    {
        _secondary_variables.addSecondaryVariable(
            "heat_flux_z", makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
                &HeatTransportBHELocalAssemblerInterface::getIntPtHeatFluxZ));
    }
    */
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
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian HeatTransportBHE process.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_table, t, x, xdot, dxdot_dx, dx_dx, M, K, b, Jac,
        _coupled_solutions);
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

}  // namespace HeatTransportBHE
}  // namespace ProcessLib
