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

#ifdef USE_PETSC
#include <petsc.h>
#include "MeshLib/NodePartitionedMesh.h"
#endif

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
#ifdef USE_PETSC
    using IDX_TYPE = PetscInt;
#else
    using IDX_TYPE = std::size_t;
#endif
    typedef MathLib::RowColumnIndices<IDX_TYPE> RowColumnIndices;
    typedef RowColumnIndices::LineIndex LineIndex;

public:
    /// Creates a MeshComponentMap internally and stores the global indices for
    /// each mesh element of the given mesh_subsets.
    explicit LocalToGlobalIndexMap(
        std::vector<MeshLib::MeshSubsets*> const& mesh_subsets,
        AssemblerLib::ComponentOrder const order =
            AssemblerLib::ComponentOrder::BY_COMPONENT,
            const bool is_linear_element = true);

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

    /// Returns total number of global degrees of freedom for DDC.
    std::size_t dofSizeGlobal() const
    {
        return _mesh_component_map.getNGlobalUnknowns();
    }

    std::size_t size() const;

    RowColumnIndices operator[](std::size_t const mesh_item_id) const;

    LineIndex rowIndices(std::size_t const mesh_item_id) const;
    LineIndex columnIndices(std::size_t const mesh_item_id) const;

#ifdef USE_PETSC
    const std::vector<IDX_TYPE>&
    getLocalNonGhostIndices(std::size_t const mesh_item_id) const
    {
        return _element_non_ghost_local_ids[mesh_item_id];
    }
#endif

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

    /// Vector contains for each element a vector of global row/or entry indices
    /// in the global stiffness matrix or vector
    std::vector<LineIndex> _rows;

    /// Vector alias to that contains for each element a vector of global column indices
    /// in the global stiffness matrix
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
