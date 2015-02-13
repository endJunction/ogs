/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ROWCOLUMNINDICES_H_
#define ROWCOLUMNINDICES_H_

#include <vector>

namespace MathLib
{

template <typename IDX_TYPE>
class RowColumnIndices
{
public:
	typedef typename std::vector<IDX_TYPE> LineIndex;
	RowColumnIndices(LineIndex const& rows_, LineIndex const& columns_,
		const std::vector<IDX_TYPE>* non_ghost_local_ids_ = nullptr)
		: rows(rows_), columns(columns_), non_ghost_local_ids(non_ghost_local_ids_)
	{ }

	std::vector<IDX_TYPE> const* getNonGhostLocalIds() const
	{
		return non_ghost_local_ids;
	}

public:
	LineIndex const& rows;
	LineIndex const& columns;

private:
	std::vector<IDX_TYPE> const* non_ghost_local_ids;
};

} // MathLib

#endif  // ROWCOLUMNINDICES_H_
