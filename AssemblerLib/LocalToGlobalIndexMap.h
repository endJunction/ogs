/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_LOCALTOGLOBALINDEXMAP_H_
#define ASSEMBLERLIB_LOCALTOGLOBALINDEXMAP_H_

#ifndef NDEBUG
#include <iostream>
#endif  // NDEBUG

#include <vector>

#include "AssemblerLib/MeshComponentMap.h"
#include "MathLib/LinAlg/RowColumnIndices.h"
#include "MeshLib/MeshSubsets.h"

namespace AssemblerLib
{

/// Row and column indices in global linear algebra objects for each mesh item.
///
/// The row and column indices in global matrix and rhs vector for example,
/// required for addition of local contributions from every mesh item (like node
/// or cell) to global objects.
///
/// The number of rows should be equal to the number of mesh items and the
/// number of columns should be equal to the number of the components on that
/// mesh item.
class LocalToGlobalIndexMap
{
public:
    using IDX_TYPE = std::size_t;
    typedef MathLib::RowColumnIndices<IDX_TYPE> RowColumnIndices;
    typedef RowColumnIndices::LineIndex LineIndex;

public:
    /* \todo Extend the constructor for parallel meshes.
    LocalToGlobalIndexMap(
        std::vector<LineIndex> const& rows,
        std::vector<LineIndex> const & columns)
        : _rows(rows), _columns(columns)
    {
        assert(_rows.size() == _columns.size());
        assert(!_rows.empty());
    }
    */

    /// Creates a MeshComponentMap internally and stores the global indices for
    /// each mesh element of the given mesh_subsets.
    explicit LocalToGlobalIndexMap(
        std::vector<MeshLib::MeshSubsets*> const& mesh_subsets,
        AssemblerLib::ComponentOrder const order =
            AssemblerLib::ComponentOrder::BY_COMPONENT);

    /// Derive a LocalToGlobalIndexMap constrained to a set of mesh subsets and
    /// elements. A new mesh component map will be constructed using the passed
    /// mesh_subsets.
    ///
    /// \note The elements are not necessary those used in the mesh_subsets.
    LocalToGlobalIndexMap* deriveBoundaryConstrainedMap(
        std::vector<MeshLib::MeshSubsets*> const& mesh_subsets,
        std::vector<MeshLib::Element const*> const& elements,
        AssemblerLib::ComponentOrder const order =
            AssemblerLib::ComponentOrder::BY_COMPONENT) const;

    /// Returns total number of degrees of freedom.
    std::size_t dofSize() const;

    std::size_t size() const;

    RowColumnIndices operator[](std::size_t const mesh_item_id) const;

    LineIndex rowIndices(std::size_t const mesh_item_id) const;
    LineIndex columnIndices(std::size_t const mesh_item_id) const;

private:

    /// Private constructor used by internally created local-to-global index
    /// maps. The mesh_component_map is passed as argument instead of being
    /// created by the constructor.
    /// \attention The passed mesh_component_map is in undefined state after
    /// this construtor.
    explicit LocalToGlobalIndexMap(
        std::vector<MeshLib::MeshSubsets*> const& mesh_subsets,
        std::vector<MeshLib::Element const*> const& elements,
        AssemblerLib::MeshComponentMap&& mesh_component_map,
        AssemblerLib::ComponentOrder const order);

    template <typename ElementIterator>
    void
    findGlobalIndices(ElementIterator first, ElementIterator last,
        MeshLib::MeshSubset const* const mesh_subset,
        AssemblerLib::ComponentOrder const order)
    {
        std::size_t const mesh_id = mesh_subset->getMeshID();
        // For each element find the global indices for node/element
        // components.
        for (ElementIterator e = first; e != last; ++e)
        {
            std::vector<MeshLib::Location> vec_items;
            std::size_t const nnodes = (*e)->getNNodes();
            vec_items.reserve(nnodes);

            for (unsigned n = 0; n < nnodes; n++)
            {
                vec_items.emplace_back(
                    mesh_id,
                    MeshLib::MeshItemType::Node,
                    (*e)->getNode(n)->getID());
            }

            // Save a line of indices for the current element.
            switch (order)
            {
                case AssemblerLib::ComponentOrder::BY_LOCATION:
                    _rows.push_back(_mesh_component_map.getGlobalIndicesByLocation<IDX_TYPE>(vec_items));
                    break;
                case AssemblerLib::ComponentOrder::BY_COMPONENT:
                    _rows.push_back(_mesh_component_map.getGlobalIndicesByComponent<IDX_TYPE>(vec_items));
                    break;
            }
        }
    }

private:
    std::vector<MeshLib::MeshSubsets*> const& _mesh_subsets;
    AssemblerLib::MeshComponentMap _mesh_component_map;

    /// _rows contains for each element a vector of global indices to
    /// node/element process variables.
    std::vector<LineIndex> _rows;

    /// For non-parallel implementations the columns are equal to the rows.
    /// \todo This is to be overriden by any parallel implementation.
    std::vector<LineIndex> const& _columns = _rows;

#ifndef NDEBUG
    /// Prints first rows of the table, every line, and the mesh component map.
    friend std::ostream& operator<<(std::ostream& os, LocalToGlobalIndexMap const& map);
#endif  // NDEBUG

    /// Element wise local indices for the unknowns that are associated with non-ghost nodes.
    std::vector<std::vector<IDX_TYPE>> _element_non_ghost_local_ids;

    /// Set local indices for the unknowns that are associated with non-ghost nodes for each element.
    /// \param ghost_element Flag for ghost element.
    /// \param vec_items     Component info.
    template <ComponentOrder ORDER>
    void
    setLocalNonGhostIndices(const bool ghost_element,
        std::vector<MeshLib::Location> const& vec_items);

};

template <ComponentOrder ORDER>
void LocalToGlobalIndexMap::setLocalNonGhostIndices(const bool ghost_element,
                                                    const std::vector<MeshLib::Location>& vec_items)
{
    std::vector<IDX_TYPE> ele_local_ids;

    if(!ghost_element)
    {
        // Empty vector for non-ghost element.
        _element_non_ghost_local_ids.push_back(ele_local_ids);
        return;
    }

    std::vector<bool> ele_ghost_flags = _mesh_component_map.getGhostFlags<ORDER>(vec_items);

    for(std::size_t i = 0; i < ele_ghost_flags.size(); i++)
    {
        if(ele_ghost_flags[i])
            ele_local_ids.push_back(i);
    }

    _element_non_ghost_local_ids.push_back(ele_local_ids);
}

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_LOCALTOGLOBALINDEXMAP_H_
