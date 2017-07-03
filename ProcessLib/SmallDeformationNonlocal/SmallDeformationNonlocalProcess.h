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

#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformationCommon/Common.h"
#include "ProcessLib/SmallDeformationNonlocal/CreateLocalAssemblers.h"

#include "SmallDeformationNonlocalFEM.h"
#include "SmallDeformationNonlocalProcessData.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
template <int DisplacementDim>
class SmallDeformationNonlocalProcess final : public Process
{
    using Base = Process;

public:
    SmallDeformationNonlocalProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        SmallDeformationNonlocalProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller)
        : Process(mesh, std::move(jacobian_assembler), parameters,
                  integration_order, std::move(process_variables),
                  std::move(secondary_variables),
                  std::move(named_function_caller)),
          _process_data(std::move(process_data))
    {
        _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
            mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);
        _material_forces = MeshLib::getOrCreateMeshProperty<double>(
            mesh, "MaterialForces", MeshLib::MeshItemType::Node,
            DisplacementDim);
    }

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

private:
    using LocalAssemblerInterface =
        SmallDeformationNonlocalLocalAssemblerInterface;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override
    {
        ProcessLib::SmallDeformationNonlocal::createLocalAssemblers<
            DisplacementDim, LocalAssemblerData>(
            mesh.getDimension(), mesh.getElements(), dof_table,
            _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
            _process_data);

        // TODO move the two data members somewhere else.
        // for extrapolation of secondary variables
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets_single_component;
        all_mesh_subsets_single_component.emplace_back(
            _mesh_subset_all_nodes.get());
        _local_to_global_index_map_single_component.reset(
            new NumLib::LocalToGlobalIndexMap(
                std::move(all_mesh_subsets_single_component),
                // by location order is needed for output
                NumLib::ComponentOrder::BY_LOCATION));
        _nodal_forces->resize(DisplacementDim * mesh.getNumberOfNodes());
        _material_forces->resize(DisplacementDim * mesh.getNumberOfNodes());

        Base::_secondary_variables.addSecondaryVariable(
            "eps_p_V", 1,
            makeExtrapolator(getExtrapolator(), _local_assemblers,
                             &SmallDeformationNonlocalLocalAssemblerInterface::
                                 getIntPtEpsPV));
        Base::_secondary_variables.addSecondaryVariable(
            "eps_p_D_xx", 1,
            makeExtrapolator(getExtrapolator(), _local_assemblers,
                             &SmallDeformationNonlocalLocalAssemblerInterface::
                                 getIntPtEpsPDXX));

        Base::_secondary_variables.addSecondaryVariable(
            "damage", 1,
            makeExtrapolator(getExtrapolator(), _local_assemblers,
                             &SmallDeformationNonlocalLocalAssemblerInterface::
                                 getIntPtDamage));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xx", 1,
            makeExtrapolator(getExtrapolator(), _local_assemblers,
                             &SmallDeformationNonlocalLocalAssemblerInterface::
                                 getIntPtSigmaXX));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_yy", 1,
            makeExtrapolator(getExtrapolator(), _local_assemblers,
                             &SmallDeformationNonlocalLocalAssemblerInterface::
                                 getIntPtSigmaYY));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_zz", 1,
            makeExtrapolator(getExtrapolator(), _local_assemblers,
                             &SmallDeformationNonlocalLocalAssemblerInterface::
                                 getIntPtSigmaZZ));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xy", 1,
            makeExtrapolator(getExtrapolator(), _local_assemblers,
                             &SmallDeformationNonlocalLocalAssemblerInterface::
                                 getIntPtSigmaXY));

        if (DisplacementDim == 3)
        {
            Base::_secondary_variables.addSecondaryVariable(
                "sigma_xz", 1,
                makeExtrapolator(
                    getExtrapolator(), _local_assemblers,
                    &SmallDeformationNonlocalLocalAssemblerInterface::
                        getIntPtSigmaXZ));

            Base::_secondary_variables.addSecondaryVariable(
                "sigma_yz", 1,
                makeExtrapolator(
                    getExtrapolator(), _local_assemblers,
                    &SmallDeformationNonlocalLocalAssemblerInterface::
                        getIntPtSigmaYZ));
        }

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xx", 1,
            makeExtrapolator(getExtrapolator(), _local_assemblers,
                             &SmallDeformationNonlocalLocalAssemblerInterface::
                                 getIntPtEpsilonXX));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_yy", 1,
            makeExtrapolator(getExtrapolator(), _local_assemblers,
                             &SmallDeformationNonlocalLocalAssemblerInterface::
                                 getIntPtEpsilonYY));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_zz", 1,
            makeExtrapolator(getExtrapolator(), _local_assemblers,
                             &SmallDeformationNonlocalLocalAssemblerInterface::
                                 getIntPtEpsilonZZ));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xy", 1,
            makeExtrapolator(getExtrapolator(), _local_assemblers,
                             &SmallDeformationNonlocalLocalAssemblerInterface::
                                 getIntPtEpsilonXY));
        if (DisplacementDim == 3)
        {
            Base::_secondary_variables.addSecondaryVariable(
                "epsilon_yz", 1,
                makeExtrapolator(
                    getExtrapolator(), _local_assemblers,
                    &SmallDeformationNonlocalLocalAssemblerInterface::
                        getIntPtEpsilonYZ));

            Base::_secondary_variables.addSecondaryVariable(
                "epsilon_xz", 1,
                makeExtrapolator(
                    getExtrapolator(), _local_assemblers,
                    &SmallDeformationNonlocalLocalAssemblerInterface::
                        getIntPtEpsilonXZ));
        }

#ifdef PROTOBUF_FOUND
        Base::integration_point_writer =
            [this](MeshLib::PropertyVector<char>& output,
                   MeshLib::PropertyVector<std::size_t>& offsets) {
                return writeIntegrationPointData(output, offsets);
            };
#endif  // PROTOBUF_FOUND

        GlobalExecutor::executeMemberOnDereferenced(
            &SmallDeformationNonlocalLocalAssemblerInterface::nonlocal,
            _local_assemblers, _local_assemblers);
    }

    std::size_t writeIntegrationPointData(
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

    void assembleConcreteProcess(
        const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
        GlobalVector& b, StaggeredCouplingTerm const& coupling_term) override
    {
        DBUG("Assemble SmallDeformationNonlocalProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assemble,
            _local_assemblers, *_local_to_global_index_map, t, x, M, K, b,
            coupling_term);
    }

    void preAssembleConcreteProcess(const double t,
                                    GlobalVector const& x) override
    {
        DBUG("preAssemble SmallDeformationNonlocalProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::preAssemble,
            _local_assemblers, *_local_to_global_index_map, t, x);
    }

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac,
        StaggeredCouplingTerm const& coupling_term) override
    {
        DBUG("AssembleWithJacobian SmallDeformationNonlocalProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
            _local_assemblers, *_local_to_global_index_map, t, x, xdot,
            dxdot_dx, dx_dx, M, K, b, Jac, coupling_term);
    }

    void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                                    double const dt) override
    {
        DBUG("PreTimestep SmallDeformationNonlocalProcess.");

        _process_data.dt = dt;
        _process_data.t = t;

        GlobalExecutor::executeMemberOnDereferenced(
            &SmallDeformationNonlocalLocalAssemblerInterface::preTimestep,
            _local_assemblers, *_local_to_global_index_map, x, t, dt);
    }

    void postTimestepConcreteProcess(GlobalVector const& x) override
    {
        DBUG("PostTimestep SmallDeformationNonlocalProcess.");

        ProcessLib::SmallDeformation::writeMaterialForces(
            *_material_forces, _local_assemblers, *_local_to_global_index_map,
            x);
    }

private:
    SmallDeformationNonlocalProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;

    MeshLib::PropertyVector<double>* _nodal_forces = nullptr;
    MeshLib::PropertyVector<double>* _material_forces = nullptr;
};

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
