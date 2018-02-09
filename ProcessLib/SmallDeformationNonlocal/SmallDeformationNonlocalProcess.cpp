/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SmallDeformationNonlocalProcess.h"
#include <iostream>

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
template <int DisplacementDim>
SmallDeformationNonlocalProcess<DisplacementDim>::
    SmallDeformationNonlocalProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        SmallDeformationNonlocalProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
    _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);

    _integration_point_writer.emplace_back(
        std::make_unique<DamageIntegrationPointWriter>([this]() {
            // Result containing integration point data for each local
            // assembler.
            std::vector<std::vector<double>> result;
            result.resize(_local_assemblers.size());

            for (std::size_t i = 0; i < _local_assemblers.size(); ++i)
            {
                auto const& local_asm = *_local_assemblers[i];

                result[i] = local_asm.getDamage();
            }

            return result;
        }));
}

template <int DisplacementDim>
bool SmallDeformationNonlocalProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void SmallDeformationNonlocalProcess<DisplacementDim>::
    initializeConcreteProcess(NumLib::LocalToGlobalIndexMap const& dof_table,
                              MeshLib::Mesh const& mesh,
                              unsigned const integration_order)
{
    ProcessLib::SmallDeformationNonlocal::createLocalAssemblers<
        DisplacementDim, SmallDeformationNonlocalLocalAssembler>(
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

    Process::_secondary_variables.addSecondaryVariable(
        "sigma",
        makeExtrapolator(
            ProcessLib::KelvinVectorType<DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), _local_assemblers,
            &LocalAssemblerInterface::getIntPtSigma));
    Process::_secondary_variables.addSecondaryVariable(
        "eps_p_V",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsPV));
    Process::_secondary_variables.addSecondaryVariable(
        "eps_p_D_xx",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsPDXX));

    Process::_secondary_variables.addSecondaryVariable(
        "damage",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtDamage));
    // TODO remove the component-wise methods
    Process::_secondary_variables.addSecondaryVariable(
        "sigma_xx",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaXX));

    Process::_secondary_variables.addSecondaryVariable(
        "sigma_yy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaYY));

    Process::_secondary_variables.addSecondaryVariable(
        "sigma_zz",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaZZ));

    Process::_secondary_variables.addSecondaryVariable(
        "sigma_xy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaXY));

    if (DisplacementDim == 3)
    {
        Process::_secondary_variables.addSecondaryVariable(
            "sigma_xz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtSigmaXZ));

        Process::_secondary_variables.addSecondaryVariable(
            "sigma_yz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtSigmaYZ));
    }

    Process::_secondary_variables.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(
            ProcessLib::KelvinVectorType<DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), _local_assemblers,
            &LocalAssemblerInterface::getIntPtEpsilon));

    // TODO remove the component-wise methods
    Process::_secondary_variables.addSecondaryVariable(
        "epsilon_xx",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonXX));

    Process::_secondary_variables.addSecondaryVariable(
        "epsilon_yy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonYY));

    Process::_secondary_variables.addSecondaryVariable(
        "epsilon_zz",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonZZ));

    Process::_secondary_variables.addSecondaryVariable(
        "epsilon_xy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonXY));
    if (DisplacementDim == 3)
    {
        Process::_secondary_variables.addSecondaryVariable(
            "epsilon_yz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtEpsilonYZ));

        Process::_secondary_variables.addSecondaryVariable(
            "epsilon_xz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtEpsilonXZ));
    }

    // enable output of internal variables defined by material models
    auto const internal_variables =
        _process_data.material->getInternalVariables();
    for (auto const& internal_variable : internal_variables)
    {
        auto const& name = internal_variable.name;
        auto const& fct = internal_variable.getter;
        auto const num_components = internal_variable.num_components;
        DBUG("Registering internal variable %s.", name.c_str());

        auto getIntPtValues = BaseLib::easyBind(
            [fct, num_components](
                LocalAssemblerInterface const& loc_asm,
                const double /*t*/,
                GlobalVector const& /*current_solution*/,
                NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                std::vector<double>& cache) -> std::vector<double> const& {

                const unsigned num_int_pts =
                    loc_asm.getNumberOfIntegrationPoints();

                cache.clear();
                auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
                    double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
                    cache, num_components, num_int_pts);

                // TODO avoid the heap allocation (one per finite element)
                std::vector<double> cache_column(num_int_pts);

                for (unsigned i = 0; i < num_int_pts; ++i)
                {
                    auto const& state = loc_asm.getMaterialStateVariablesAt(i);

                    auto const& int_pt_values = fct(state, cache_column);
                    assert(int_pt_values.size() == num_components);
                    auto const int_pt_values_vec =
                        MathLib::toVector(int_pt_values);

                    cache_mat.col(i).noalias() = int_pt_values_vec;
                }

                return cache;
            });

        Process::_secondary_variables.addSecondaryVariable(
            name,
            makeExtrapolator(num_components, getExtrapolator(),
                             _local_assemblers, std::move(getIntPtValues)));
    }

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::nonlocal, _local_assemblers,
        _local_assemblers);


    // Set initial conditions for integration point data.
    for (auto const& ip_writer : _integration_point_writer)
    {
        auto const& name = ip_writer->name();
        if (!mesh.getProperties().existsPropertyVector<double>(name))
        {
            continue;
        }

        auto const& mesh_property =
            *mesh.getProperties().template getPropertyVector<double>(name);

        if (mesh_property.getMeshItemType() !=
            MeshLib::MeshItemType::IntegrationPoint)
        {
            continue;
        }

        auto const offsets_array_name = name + "_offsets";
        if (!mesh.getProperties().existsPropertyVector<std::size_t>(
                offsets_array_name))
        {
            OGS_FATAL(
                "Integration point data '%s' is present in the vtk field data "
                "but the corresponding '%s' array is not available.",
                name.c_str(), offsets_array_name.c_str());
        }

        auto const& mesh_property_offsets =
            *mesh.getProperties().template getPropertyVector<std::size_t>(
                offsets_array_name);

        if (mesh_property_offsets.getMeshItemType() !=
            MeshLib::MeshItemType::Cell)
        {
            OGS_FATAL(
                "Integration point data '%s' is present in the vtk field data "
                "but the corresponding '%s' array is not defined on cells.",
                name.c_str(), offsets_array_name.c_str());
        }

        // Now we have a properly named vtk's field data array and the
        // corresponding offsets array.
        for (std::size_t i = 0; i < _local_assemblers.size(); ++i)
        {
            auto& local_asm = _local_assemblers[i];
            auto const offset = mesh_property_offsets[i];

            // TODO (naumov) Check sizes / read size / etc.
            // OR reconstruct dimensions from size / component = ip_points
            local_asm->setIPDataInitialConditions(name, &mesh_property[offset]);
        }
    }
}

template <int DisplacementDim>
void SmallDeformationNonlocalProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble SmallDeformationNonlocalProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        *_local_to_global_index_map, t, x, M, K, b, _coupled_solutions);
}

template <int DisplacementDim>
void SmallDeformationNonlocalProcess<
    DisplacementDim>::preAssembleConcreteProcess(const double t,
                                                 GlobalVector const& x)
{
    DBUG("preAssemble SmallDeformationNonlocalProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::preAssemble,
        _local_assemblers, *_local_to_global_index_map, t, x);
}

template <int DisplacementDim>
void SmallDeformationNonlocalProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(const double t, GlobalVector const& x,
                                        GlobalVector const& xdot,
                                        const double dxdot_dx,
                                        const double dx_dx, GlobalMatrix& M,
                                        GlobalMatrix& K, GlobalVector& b,
                                        GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian SmallDeformationNonlocalProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, *_local_to_global_index_map, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac, _coupled_solutions);

    b.copyValues(*_nodal_forces);
    std::transform(_nodal_forces->begin(), _nodal_forces->end(),
                   _nodal_forces->begin(), [](double val) { return -val; });
}

template <int DisplacementDim>
void SmallDeformationNonlocalProcess<
    DisplacementDim>::preTimestepConcreteProcess(GlobalVector const& x,
                                                 double const t,
                                                 double const dt,
                                                 int const /*process_id*/)
{
    DBUG("PreTimestep SmallDeformationNonlocalProcess.");

    _process_data.dt = dt;
    _process_data.t = t;
    _process_data.injected_volume = _process_data.t;

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        *_local_to_global_index_map, x, t, dt);
}

template <int DisplacementDim>
void SmallDeformationNonlocalProcess<
    DisplacementDim>::postTimestepConcreteProcess(GlobalVector const& x)
{
    DBUG("PostTimestep SmallDeformationNonlocalProcess.");

    ProcessLib::SmallDeformation::writeMaterialForces(
        _material_forces, _local_assemblers, *_local_to_global_index_map, x);

    auto material_forces_property = MeshLib::getOrCreateMeshProperty<double>(
        _mesh, "MaterialForces", MeshLib::MeshItemType::Node, DisplacementDim);
    _material_forces->copyValues(*material_forces_property);
}

template <int DisplacementDim>
NumLib::IterationResult
SmallDeformationNonlocalProcess<DisplacementDim>::postIterationConcreteProcess(
    GlobalVector const& x)
{
    _process_data.crack_volume = 0.0;

    DBUG("PostNonLinearSolver crack volume computation.");

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::computeCrackIntegral, _local_assemblers,
        *_local_to_global_index_map, x, _process_data.crack_volume);

    INFO("Integral of crack: %g", _process_data.crack_volume);

    if (_process_data.propagating_crack)
    {
        _process_data.pressure_old = _process_data.pressure;
        _process_data.pressure = (_process_data.injected_volume - _process_data.crack_volume)*2e6;

        std::cout << "\n Pressure = " << _process_data.pressure << "\n";
        _process_data.pressure_error =
            _process_data.pressure == 0
                ? 0
                : std::abs(_process_data.pressure_old -
                           _process_data.pressure) /
                           _process_data.pressure;

        INFO("Internal pressure: %g and Pressure error: %.4e",
             _process_data.pressure, _process_data.pressure_error);

        // TODO (parisio) Is this correct?
        // Update displacement field
        //assert(_coupled_solutions == nullptr);
        //MathLib::LinAlg::scale(const_cast<GlobalVector&>(x),
        //                       _process_data.pressure);
    }

    // TODO (parisio) try this to enforce pressure convergence.
    /*
    if (_process_data.pressure_error > 1e-4)
    {
        return NumLib::IterationResult::REPEAT_ITERATION;
    }
    */

    return NumLib::IterationResult::SUCCESS;
}

template class SmallDeformationNonlocalProcess<2>;
template class SmallDeformationNonlocalProcess<3>;

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
