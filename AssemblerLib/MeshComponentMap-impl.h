/**
 * \author Norihiro Watanabe
 * \author Wenqing Wang
 * \date   2013-04-16
 * \date   2014-11
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshComponentMap.h"

namespace AssemblerLib
{

using namespace detail;

template<RowDataType ROW_DATA_TYPE, typename T_DATA_TYPE>
std::vector<T_DATA_TYPE> MeshComponentMap
::getRowDataByLocation(const std::vector<Location> &l) const
{
    // Create vector of global indices sorted by location containing all
    // locations given in ls parameter.
    // or Create vector of ghost node flags.

    std::vector<T_DATA_TYPE> row_data;
    row_data.reserve(l.size());

    auto const &m = _dict.get<ByLocation>();
    for (auto lc = l.cbegin(); lc != l.cend(); ++lc)
    {
        auto const p = m.equal_range(Line(*lc));
        if(ROW_DATA_TYPE == RowDataType::ROW_GLOBAL_INDEX)
        {
            for (auto itr = p.first; itr != p.second; ++itr)
                row_data.push_back(itr->global_index);
        }
        else if(ROW_DATA_TYPE == RowDataType::ROW_GHOST_FLAG)
        {
            for (auto itr = p.first; itr != p.second; ++itr)
                row_data.push_back(itr->is_ghost);
        }
    }

    return row_data;
}

template<RowDataType ROW_DATA_TYPE, typename T_DATA_TYPE>
std::vector<T_DATA_TYPE> MeshComponentMap
::getRowDataByComponent(const std::vector<Location> &l) const
{
    // vector of (Component, global Index) pairs.
    typedef std::pair<std::size_t, std::size_t> CIPair;
    std::vector<CIPair> pairs;
    pairs.reserve(l.size());

    // Create a sub dictionary containing all lines with location from ls.
    auto const &m = _dict.get<ByLocation>();
    for (auto lc = l.cbegin(); lc != l.cend(); ++lc)
    {
        auto const p = m.equal_range(Line(*lc));
        if(ROW_DATA_TYPE == RowDataType::ROW_GLOBAL_INDEX)
        {
            for (auto itr = p.first; itr != p.second; ++itr)
                pairs.emplace_back(itr->comp_id, itr->global_index);
        }
        else if(ROW_DATA_TYPE == RowDataType::ROW_GHOST_FLAG)
        {
            for (auto itr = p.first; itr != p.second; ++itr)
                pairs.emplace_back(itr->comp_id, itr->is_ghost);
        }
    }

    auto CIPairLess = [](CIPair const& a, CIPair const& b)
    {
        return a.first < b.first;
    };

    // Create vector of global indices or ghost node flags from sub dictionary sorting by component.
    if (!std::is_sorted(pairs.begin(), pairs.end(), CIPairLess))
        std::stable_sort(pairs.begin(), pairs.end(), CIPairLess);

    std::vector<T_DATA_TYPE> row_data;
    row_data.reserve(pairs.size());
    for (auto p = pairs.cbegin(); p != pairs.cend(); ++p)
        row_data.push_back(p->second);

    return row_data;
}

}   // namespace AssemblerLib
