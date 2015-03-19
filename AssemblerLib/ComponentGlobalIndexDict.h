/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_COMPONENTGLOBALINDEXDICT_H_
#define ASSEMBLERLIB_COMPONENTGLOBALINDEXDICT_H_

#include <iostream>
#include <limits>

#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index_container.hpp>

#include "MeshLib/Location.h"

namespace AssemblerLib
{

/// \internal
namespace detail
{

struct Line
{
    MeshLib::Location location;

    // Physical component
    std::size_t comp_id;

    // Position in global matrix or vector
    std::size_t global_index;

    // Element entity status in a mesh partition, either ghost or real one in use
    bool is_ghost = false;

    Line(MeshLib::Location const& l, std::size_t c, std::size_t i,
        bool const ghost_flag = false)
    : location(l), comp_id(c), global_index(i), is_ghost(ghost_flag)
    {}

    Line(MeshLib::Location const& l, std::size_t c, bool const ghost_flag = false)
    : location(l), comp_id(c),
        global_index(std::numeric_limits<std::size_t>::max()), is_ghost(ghost_flag)
    {}

    explicit Line(MeshLib::Location const& l, bool const ghost_flag = false)
    : location(l),
        comp_id(std::numeric_limits<std::size_t>::max()),
        global_index(std::numeric_limits<std::size_t>::max()),
        is_ghost(ghost_flag)
    {}

    friend std::ostream& operator<<(std::ostream& os, Line const& l)
    {
        return os << l.location << ", " << l.comp_id << ", " << l.global_index
            << ", " << (l.is_ghost ? "ghost" : "");
    }
};

struct LineByLocationComparator
{
    bool operator()(Line const& a, Line const& b) const
    {
        return a.location < b.location;
    }
};

struct LineByLocationAndComponentAndGhostComparator
{
    bool operator()(Line const& a, Line const& b) const
    {
        if (a.location < b.location)
            return true;
        if (b.location < a.location)
            return false;
        // a.loc == b.loc

        if (a.comp_id < b.comp_id)
            return true;
        if (b.comp_id < a.comp_id)
            return false;
        // a.comp_id == b.comp_id

        return a.is_ghost < b.is_ghost;
    }
};

struct ByLocation {};
struct ByLocationAndComponent {};
struct ByComponent {};
struct ByGlobalIndex {};

typedef boost::multi_index::multi_index_container<
        Line,
        boost::multi_index::indexed_by
        <
            boost::multi_index::ordered_unique
            <
                boost::multi_index::tag<ByLocationAndComponent>,
                boost::multi_index::identity<Line>,
                LineByLocationAndComponentAndGhostComparator
            >,
            boost::multi_index::ordered_non_unique
            <
                boost::multi_index::tag<ByLocation>,
                boost::multi_index::identity<Line>,
                LineByLocationComparator
            >,
            boost::multi_index::ordered_non_unique
            <
                boost::multi_index::tag<ByComponent>,
                boost::multi_index::member<Line, std::size_t, &Line::comp_id>
            >,
            boost::multi_index::ordered_non_unique
            <
                boost::multi_index::tag<ByGlobalIndex>,
                boost::multi_index::member<Line, std::size_t, &Line::global_index>
            >
        >
    > ComponentGlobalIndexDict;

}    // namespace detail
}    // namespace AssemblerLib

#endif  // ASSEMBLERLIB_COMPONENTGLOBALINDEXDICT_H_
