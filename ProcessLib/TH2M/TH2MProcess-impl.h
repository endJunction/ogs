/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <omp.h>
#include <cassert>

#include "CreateLocalAssemblers.h"
#include "MeshLib/Elements/Utils.h"
#include "ProcessLib/Process.h"

#include "TH2MProcessData.h"
#include "TH2M_FEM.h"

namespace ProcessLib
{
namespace TH2M
{
template <int DisplacementDim>
TH2MProcess<DisplacementDim>::TH2MProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    TH2MProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
}

template <int DisplacementDim>
bool TH2MProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::constructDofTable()
{
    DBUG("TH2M dof table construction..");
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, &_mesh.getNodes());
    // Create single component dof in the mesh's base nodes.
    _base_nodes = MeshLib::getBaseNodes(_mesh.getElements());
    _mesh_subset_base_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, &_base_nodes);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubsets> all_mesh_subsets_single_component;
    all_mesh_subsets_single_component.emplace_back(
        _mesh_subset_all_nodes.get());
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    if (_use_monolithic_scheme)
    {
        // Collect the mesh subsets in a vector.
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets;

        // For gas pressure, capillary pressure, and temperature:
        for (int i = 0; i <= indexTemperature; i++)
            all_mesh_subsets.emplace_back(_mesh_subset_base_nodes.get());

        // For displacement:
        const int monolithic_process_id = 0;
        std::generate_n(
                std::back_inserter(all_mesh_subsets),
                getProcessVariables(monolithic_process_id)[indexDisplacement].get().getNumberOfComponents(),
                [&]() { return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()}; });

        std::vector<int> const vec_n_components{1,1,1, DisplacementDim};
        _local_to_global_index_map =
                std::make_unique<NumLib::LocalToGlobalIndexMap>(
                        std::move(all_mesh_subsets), vec_n_components,
                        NumLib::ComponentOrder::BY_LOCATION);

        assert(_local_to_global_index_map);
    }
    else
    {
        OGS_FATAL("Staggered TH2M Process is not implemented!");
    }
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    const int mechanical_process_id = _use_monolithic_scheme ? 0 : 1;
   // const int deformation_variable_id = _use_monolithic_scheme ? 1 : 0;
    ProcessLib::TH2M::createLocalAssemblers<
        DisplacementDim, TH2MLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        // use displacement process variable to set shape function order
        getProcessVariables(mechanical_process_id)[indexDisplacement]
            .get()
            .getShapeFunctionOrder(),
        _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
        _process_data);

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_eff",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigma));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilon));

    Base::_secondary_variables.addSecondaryVariable(
            "saturation",
            makeExtrapolator(
                    1, getExtrapolator(), _local_assemblers,
                    &LocalAssemblerInterface::getIntPtSaturation));

    Base::_secondary_variables.addSecondaryVariable(
            "pressure_wet",
            makeExtrapolator(
                    1, getExtrapolator(), _local_assemblers,
                    &LocalAssemblerInterface::getIntPtWetPressure));

    Base::_secondary_variables.addSecondaryVariable(
            "density_gas",
            makeExtrapolator(
                    1, getExtrapolator(), _local_assemblers,
                    &LocalAssemblerInterface::getIntPtDensityGas));

    Base::_secondary_variables.addSecondaryVariable(
            "density_liquid",
            makeExtrapolator(
                    1, getExtrapolator(), _local_assemblers,
                    &LocalAssemblerInterface::getIntPtDensityLiquid));

    Base::_secondary_variables.addSecondaryVariable(
            "velocity_gas",
            makeExtrapolator(mesh.getDimension(), getExtrapolator(),
                    _local_assemblers,
                    &LocalAssemblerInterface::getIntPtDarcyVelocityGas));

    Base::_secondary_variables.addSecondaryVariable(
            "velocity_liquid",
            makeExtrapolator(mesh.getDimension(), getExtrapolator(),
                    _local_assemblers,
                    &LocalAssemblerInterface::getIntPtDarcyVelocityLiquid));
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble TH2MProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};

    // Call global assembler for each local assembly item.
    auto const n = _local_assemblers.size();

#pragma omp parallel for
    for (std::size_t i = 0; i < n; i++)
    {
        _global_assembler.assemble(i, *_local_assemblers[i], dof_table, t, x, M,
                                   K, b, _coupled_solutions);
    }
}


//template <int DisplacementDim>
//void TH2MProcess<DisplacementDim>::assembleWithJacobianConcreteProcess(
//    const double t, GlobalVector const& x, GlobalVector const& xdot,
//    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
//    GlobalVector& b, GlobalMatrix& Jac)
//{
//    DBUG("AssembleJacobian TH2MProcess.");
//    DBUG("x: %i", x.size());
//
//    // Call global assembler for each local assembly item.
//    GlobalExecutor::executeMemberDereferenced(
//        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
//        _local_assemblers, *_local_to_global_index_map, t, x, xdot, dxdot_dx,
//        dx_dx, M, K, b, Jac, _coupled_solutions);
//}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(const double t, GlobalVector const& x,
                                        GlobalVector const& xdot,
                                        const double dxdot_dx,
                                        const double dx_dx, GlobalMatrix& M,
                                        GlobalMatrix& K, GlobalVector& b,
                                        GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the monolithic scheme
    if (_use_monolithic_scheme)
    {
        DBUG(
            "Assemble the Jacobian of TH2M for the monolithic"
            " scheme.");
        dof_tables.emplace_back(*_local_to_global_index_map);
    }
    else
    {
        OGS_FATAL("Staggered TH2M Process not implemented.");
    }

    auto const n = _local_assemblers.size();

#pragma omp parallel for
    for (std::size_t i = 0; i < n; i++)
    {
        _global_assembler.assembleWithJacobian(
            i, *_local_assemblers[i], dof_tables, t, x, xdot, dxdot_dx, dx_dx,
            M, K, b, Jac, _coupled_solutions);
    }
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt,
    const int /*process_id*/)
{
    DBUG("PreTimestep TH2MProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        *_local_to_global_index_map, x, t, dt);
}

template <int DisplacementDim>
void TH2MProcess<DisplacementDim>::postTimestepConcreteProcess(
        GlobalVector const& x, const int process_id)
    {
        DBUG("PostTimestep HydroMechanicsProcess.");
        GlobalExecutor::executeMemberOnDereferenced(
            &LocalAssemblerInterface::postTimestep, _local_assemblers,
            getDOFTable(process_id), x);
    }


}  // namespace TH2M
}  // namespace ProcessLib
