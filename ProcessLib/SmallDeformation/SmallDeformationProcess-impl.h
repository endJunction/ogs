/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>

#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

#include "SmallDeformationFEM.h"

namespace ProcessLib
{
namespace SmallDeformation
{
template <int DisplacementDim>
SmallDeformationProcess<DisplacementDim>::SmallDeformationProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    SmallDeformationProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
    _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);
    _material_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "MaterialForces", MeshLib::MeshItemType::Node, DisplacementDim);
}

template <int DisplacementDim>
bool SmallDeformationProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, SmallDeformationLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables
    std::vector<MeshLib::MeshSubsets> all_mesh_subsets_single_component;
    all_mesh_subsets_single_component.emplace_back(
        _mesh_subset_all_nodes.get());
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);
    _nodal_forces->resize(DisplacementDim * mesh.getNumberOfNodes());

    _material_forces->resize(DisplacementDim * mesh.getNumberOfNodes());

    Base::_secondary_variables.addSecondaryVariable(
        "free_energy_density", 1,
        makeExtrapolator(getExtrapolator(), _local_assemblers,
                         &SmallDeformationLocalAssemblerInterface::
                             getIntPtFreeEnergyDensity));

    Base::_secondary_variables.addSecondaryVariable(
        "eps_p_V", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtEpsPV));
    Base::_secondary_variables.addSecondaryVariable(
        "eps_p_D_xx", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtEpsPDXX));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xx", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXX));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_yy", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaYY));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_zz", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xy", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXY));

    if (DisplacementDim == 3)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xz", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXZ));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_yz", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtSigmaYZ));
    }

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xx", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonXX));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_yy", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonYY));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_zz", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xy", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonXY));
    if (DisplacementDim == 3)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_yz", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonYZ));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xz", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonXZ));
    }

#ifdef PROTOBUF_FOUND
        Base::integration_point_writer = [this](
            MeshLib::PropertyVector<char>& output,
            MeshLib::PropertyVector<std::size_t>& offsets) {
            return writeIntegrationPointData(output, offsets);
        };
#endif  // PROTOBUF_FOUND
}

template <int DisplacementDim>
std::size_t SmallDeformationProcess<DisplacementDim>::writeIntegrationPointData(
    MeshLib::PropertyVector<char>& output,
    MeshLib::PropertyVector<std::size_t>& offsets)
{
    output.clear();
    offsets.clear();
    std::vector<char> local_data;
    std::size_t offset = 0;
    for (auto& la : _local_assemblers)
    {
        offsets.push_back(offset);
        std::size_t const local_offset =
            la->writeIntegrationPointData(local_data);
        std::copy_n(std::begin(local_data), local_offset,
                    std::back_inserter(output));
        offset += local_offset;
    }
    return offset;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, StaggeredCouplingTerm const& coupling_term)
{
    DBUG("Assemble SmallDeformationProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        *_local_to_global_index_map, t, x, M, K, b, coupling_term);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac,
        StaggeredCouplingTerm const& coupling_term)
{
    DBUG("AssembleWithJacobian SmallDeformationProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, *_local_to_global_index_map, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac, coupling_term);

    b.copyValues(*_nodal_forces);
    std::transform(_nodal_forces->begin(), _nodal_forces->end(),
                   _nodal_forces->begin(), [](double val) { return -val; });
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::
    setInitialConditionsConcreteProcess(double const t,
                                        GlobalVector const& x) override
{
    DBUG("SetInitialConditions SmallDeformationProcess.");

    if (!_mesh.getProperties().hasPropertyVector("integration_point_data"))
        return;
    if (!_mesh.getProperties().hasPropertyVector("integration_point_offsets"))
        OGS_FATAL(
            "integration_point_data field exists in the input but there is "
            "no integration_point_offsets cell data.");

    auto const& data = *_mesh.getProperties().template getPropertyVector<char>(
        "integration_point_data");
    assert(data.getMeshItemType() == MeshLib::MeshItemType::IntegrationPoint);

    auto const& offsets =
        *_mesh.getProperties().template getPropertyVector<std::size_t>(
            "integration_point_offsets");
    assert(offsets.getMeshItemType() == MeshLib::MeshItemType::Cell);

    std::vector<char> local_data;
    assert(_local_assemblers.size() == offsets.size());
    // Starting counting from one; the last cell is handled after the loop.
    std::size_t i = 0;
    for (; i < _local_assemblers.size() - 1; ++i)
    {
        std::size_t const size = offsets[i + 1] - offsets[i];
        local_data.resize(size);
        std::memcpy(local_data.data(), &data[offsets[i]], size);
        _local_assemblers[i]->readIntegrationPointData(local_data);
    }
    {  // last cell
        std::size_t const size = data.size() - offsets[i];
        local_data.resize(size);
        std::memcpy(local_data.data(), &data[offsets[i]], size);
        _local_assemblers[i]->readIntegrationPointData(local_data);
    }
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt)
{
    DBUG("PreTimestep SmallDeformationProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    GlobalExecutor::executeMemberOnDereferenced(
        &SmallDeformationLocalAssemblerInterface::preTimestep,
        _local_assemblers, *_local_to_global_index_map, x, t, dt);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::postTimestepConcreteProcess(
    GlobalVector const& x) override
{
    DBUG("PostTimestep SmallDeformationProcess.");

    ProcessLib::SmallDeformation::writeMaterialForces(
        *_material_forces, _local_assemblers, *_local_to_global_index_map, x);
    }

}  // namespace SmallDeformation
}  // namespace ProcessLib
