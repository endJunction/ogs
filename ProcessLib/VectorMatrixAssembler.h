/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include "NumLib/NumericsConfig.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "AbstractJacobianAssembler.h"
#include "CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
struct VectorCoordinateStorage
{
    NumLib::LocalToGlobalIndexMap::RowColumnIndices::LineIndex indices;
    std::vector<double> entries;

    void add(NumLib::LocalToGlobalIndexMap::RowColumnIndices::LineIndex const&
                 local_indices,
             std::vector<double> const& local_vector)
    {
        if (local_vector.empty())
            return;
        indices.insert(indices.end(), local_indices.begin(),
                       local_indices.end());
        entries.insert(entries.end(), local_vector.cbegin(),
                       local_vector.cend());
    }

    void append(VectorCoordinateStorage const& other)
    {
        indices.insert(indices.end(), other.indices.cbegin(),
                       other.indices.cend());
        entries.insert(entries.end(), other.entries.cbegin(),
                       other.entries.cend());
    }
};

struct MatrixCoordinateStorage
{
    std::vector<GlobalIndexType> rows;
    std::vector<GlobalIndexType> columns;
    std::vector<double> entries;

    void add(NumLib::LocalToGlobalIndexMap::RowColumnIndices const& r_c_indices,
             std::vector<double> const& local_matrix)
    {
        if (local_matrix.empty())
            return;
        rows.insert(rows.end(), r_c_indices.rows.cbegin(),
                    r_c_indices.rows.cend());
        columns.insert(columns.end(), r_c_indices.columns.cbegin(),
                       r_c_indices.columns.cend());
        entries.insert(entries.end(), local_matrix.cbegin(),
                       local_matrix.cend());
    }

    void append(MatrixCoordinateStorage const& other)
    {
        rows.insert(rows.end(), other.rows.cbegin(), other.rows.cend());
        columns.insert(columns.end(), other.columns.cbegin(),
                       other.columns.cend());
        entries.insert(entries.end(), other.entries.cbegin(),
                       other.entries.cend());
    }
};

struct CoupledSolutionsForStaggeredScheme;

class LocalAssemblerInterface;

//! Utility class used to assemble global matrices and vectors.
//!
//! The methods of this class get the global matrices and vectors as input and
//! pass only local data on to the local assemblers.
class VectorMatrixAssembler final
{
public:
    explicit VectorMatrixAssembler(
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler);

    void preAssemble(const std::size_t mesh_item_id,
                     LocalAssemblerInterface& local_assembler,
                     const NumLib::LocalToGlobalIndexMap& dof_table,
                     const double t, const GlobalVector& x);

    //! Assembles\c M, \c K, and \c b.
    //! \remark Jacobian is not assembled here, see assembleWithJacobian().
    void assemble(std::size_t const mesh_item_id,
                  LocalAssemblerInterface& local_assembler,
                  std::vector<std::reference_wrapper<
                      NumLib::LocalToGlobalIndexMap>> const& dof_tables,
                  double const t, GlobalVector const& x, GlobalMatrix& M,
                  GlobalMatrix& K, GlobalVector& b,
                  CoupledSolutionsForStaggeredScheme const* const cpl_xs);

    //! Assembles \c M, \c K, \c b, and the Jacobian \c Jac of the residual.
    //! \note The Jacobian must be assembled.
    void assembleWithJacobian(
        std::size_t const mesh_item_id,
        LocalAssemblerInterface& local_assembler,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac,
        CoupledSolutionsForStaggeredScheme const* const cpl_xs);

private:
    //! Used to assemble the Jacobian.
    std::unique_ptr<AbstractJacobianAssembler> _jacobian_assembler;
};

}  // namespace ProcessLib
