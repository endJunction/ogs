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
#include "AbstractJacobianAssembler.h"
#include "CoupledSolutionsForStaggeredScheme.h"

namespace NumLib
{
class LocalToGlobalIndexMap;
}  // NumLib

namespace ProcessLib
{
struct TupleStorage
{
    NumLib::LocalToGlobalIndexMap::RowColumnIndices::LineIndex indices;
    std::vector<double> data;

    // vector add
    void add(NumLib::LocalToGlobalIndexMap::RowColumnIndices::LineIndex const&
                 local_indices,
             std::vector<double> const& local_vector)
    {
        if (local_vector.empty())
            return;
        indices.insert(indices.end(), local_indices.begin(),
                       local_indices.end());
        data.insert(data.end(), local_vector.begin(), local_vector.end());
    }

    void append(TupleStorage const& other)
    {
        indices.insert(indices.end(), other.indices.begin(),
                       other.indices.end());
        data.insert(data.end(), other.data.begin(), other.data.end());
    }
};

struct TripletStorage
{
    std::vector<Eigen::Triplet<double>> data;

    void add(NumLib::LocalToGlobalIndexMap::RowColumnIndices const& r_c_indices,
             std::vector<double> const& local_matrix)
    {
        if (local_matrix.empty())
            return;
        auto const n_rows = r_c_indices.rows.size();
        auto const n_cols = r_c_indices.columns.size();
        for (auto i = decltype(n_rows){0}; i < n_rows; i++)
        {
            auto const row = r_c_indices.rows[i];
            for (auto j = decltype(n_cols){0}; j < n_cols; j++)
            {
                auto const col = r_c_indices.columns[j];
                data.emplace_back(row, col, local_matrix[i * n_cols + j]);
            }
        }
    }

    void append(TripletStorage const& other)
    {
        data.insert(data.end(), other.data.begin(), other.data.end());
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
