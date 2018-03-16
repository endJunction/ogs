/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VectorMatrixAssembler.h"

#include <cassert>
#include <functional>  // for std::reference_wrapper.

#include "NumLib/DOF/DOFTableUtil.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "LocalAssemblerInterface.h"

#include "CoupledSolutionsForStaggeredScheme.h"
#include "Process.h"

namespace ProcessLib
{
VectorMatrixAssembler::VectorMatrixAssembler(
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler)
    : _jacobian_assembler(std::move(jacobian_assembler))
{
}

void VectorMatrixAssembler::preAssemble(
    const std::size_t mesh_item_id, LocalAssemblerInterface& local_assembler,
    const NumLib::LocalToGlobalIndexMap& dof_table, const double t,
    const GlobalVector& x)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    local_assembler.preAssemble(t, local_x);
}

void VectorMatrixAssembler::assemble(
    const std::size_t mesh_item_id, LocalAssemblerInterface& local_assembler,
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
        dof_tables,
    const double t, const GlobalVector& x, MatrixCoordinateStorage& M,
    MatrixCoordinateStorage& K, VectorCoordinateStorage& b,
    CoupledSolutionsForStaggeredScheme const* const cpl_xs)
{
    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    for (std::size_t i = 0; i < dof_tables.size(); i++)
    {
        indices_of_processes.emplace_back(
            NumLib::getIndices(mesh_item_id, dof_tables[i].get()));
    }

    auto const& indices = (cpl_xs == nullptr)
                              ? indices_of_processes[0]
                              : indices_of_processes[cpl_xs->process_id];

    // TODO(naumov) Reserve guesses?
    std::vector<double> local_M_data;
    std::vector<double> local_K_data;
    std::vector<double> local_b_data;

    if (cpl_xs == nullptr)
    {
        auto const local_x = x.get(indices);
        local_assembler.assemble(t, local_x, local_M_data, local_K_data,
                                 local_b_data);
    }
    else
    {
        auto local_coupled_xs0 =
            getPreviousLocalSolutions(*cpl_xs, indices_of_processes);
        auto local_coupled_xs =
            getCurrentLocalSolutions(*cpl_xs, indices_of_processes);

        ProcessLib::LocalCoupledSolutions local_coupled_solutions(
            cpl_xs->dt, cpl_xs->process_id, std::move(local_coupled_xs0),
            std::move(local_coupled_xs));

        local_assembler.assembleForStaggeredScheme(t, local_M_data,
                                                   local_K_data, local_b_data,
                                                   local_coupled_solutions);
    }

    auto const r_c_indices =
        NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);

    M.add(r_c_indices, local_M_data);
    K.add(r_c_indices, local_K_data);
    assert(local_b_data.size() == indices.size());
    b.add(indices, local_b_data);
}

void VectorMatrixAssembler::assembleWithJacobian(
    std::size_t const mesh_item_id, LocalAssemblerInterface& local_assembler,
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
        dof_tables,
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac,
    CoupledSolutionsForStaggeredScheme const* const cpl_xs)
{
    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    for (std::size_t i = 0; i < dof_tables.size(); i++)
    {
        indices_of_processes.emplace_back(
            NumLib::getIndices(mesh_item_id, dof_tables[i].get()));
    }

    auto const& indices = (cpl_xs == nullptr)
                              ? indices_of_processes[0]
                              : indices_of_processes[cpl_xs->process_id];
    auto const local_xdot = xdot.get(indices);

    // TODO(naumov) Reserve guesses?
    std::vector<double> local_M_data;
    std::vector<double> local_K_data;
    std::vector<double> local_b_data;
    std::vector<double> local_Jac_data;

    if (cpl_xs == nullptr)
    {
        auto const local_x = x.get(indices);
        _jacobian_assembler->assembleWithJacobian(
            local_assembler, t, local_x, local_xdot, dxdot_dx, dx_dx,
            local_M_data, local_K_data, local_b_data, local_Jac_data);
    }
    else
    {
        auto local_coupled_xs0 =
            getPreviousLocalSolutions(*cpl_xs, indices_of_processes);
        auto local_coupled_xs =
            getCurrentLocalSolutions(*cpl_xs, indices_of_processes);

        ProcessLib::LocalCoupledSolutions local_coupled_solutions(
            cpl_xs->dt, cpl_xs->process_id, std::move(local_coupled_xs0),
            std::move(local_coupled_xs));

        _jacobian_assembler->assembleWithJacobianForStaggeredScheme(
            local_assembler, t, local_xdot, dxdot_dx, dx_dx, local_M_data,
            local_K_data, local_b_data, local_Jac_data,
            local_coupled_solutions);
    }

    auto const num_r_c = indices.size();
    auto const r_c_indices =
        NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);

    if (!local_M_data.empty())
    {
        auto const local_M = MathLib::toMatrix(local_M_data, num_r_c, num_r_c);
        M.add(r_c_indices, local_M);
    }
    if (!local_K_data.empty())
    {
        auto const local_K = MathLib::toMatrix(local_K_data, num_r_c, num_r_c);
        K.add(r_c_indices, local_K);
    }
    if (!local_b_data.empty())
    {
        assert(local_b_data.size() == num_r_c);
        b.add(indices, local_b_data);
    }
    if (!local_Jac_data.empty())
    {
        auto const local_Jac =
            MathLib::toMatrix(local_Jac_data, num_r_c, num_r_c);
        Jac.add(r_c_indices, local_Jac);
    }
    else
    {
        OGS_FATAL(
            "No Jacobian has been assembled! This might be due to programming "
            "errors in the local assembler of the current process.");
    }
}

}  // namespace ProcessLib
