/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/MeshEnums.h"
#include "Element.h"
#include "EdgeReturn.h"
#include "CellRule.h"

namespace MeshLib
{

/**
 * This class represents a 3d tetrahedron element. The following sketch shows the node and edge numbering.
 * @anchor Tetrahedron4NodeAndEdgeNumbering
 * @code
 *          3
 *         /|\
 *        / | \
 *      3/  |  \5
 *      /   |4  \
 *     /    |    \
 *    0.....|.....2
 *     \    |  2 /
 *      \   |   /
 *      0\  |  /1
 *        \ | /
 *         \|/
 *          1
 *
 * @endcode
 */
class TetRule4 : public CellRule
{
public:
    /// Constant: The number of base nodes for this element
    static constexpr unsigned n_base_nodes = 4u;

    /// Constant: The number of all nodes for this element
    static constexpr unsigned n_all_nodes = 4u;

    /// Constant: The geometric type of the element
    static constexpr MeshElemType mesh_elem_type = MeshElemType::TETRAHEDRON;

    /// Constant: The FEM type of the element
    static constexpr CellType cell_type = CellType::TET4;

    /// Constant: The number of faces
    static constexpr unsigned n_faces = 4;

    /// Constant: The number of edges
    static constexpr unsigned n_edges = 6;

    /// Constant: The number of neighbors
    static constexpr unsigned n_neighbors = 4;

    /// Constant: Local node index table for faces
    static constexpr unsigned face_nodes[4][3] = {
        {0, 2, 1},  // Face 0
        {0, 1, 3},  // Face 1
        {1, 2, 3},  // Face 2
        {2, 0, 3}   // Face 3
    };

    /// Constant: Local node index table for edge
    static constexpr unsigned edge_nodes[6][2] = {
        {0, 1},  // Edge 0
        {1, 2},  // Edge 1
        {0, 2},  // Edge 2
        {0, 3},  // Edge 3
        {1, 3},  // Edge 4
        {2, 3}   // Edge 5
    };

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::LinearEdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

    /**
     * \copydoc MeshLib::Element::isPntInElement()
     * @param nodes the nodes of the element.
     */
    static bool isPntInElement(Node const* const* nodes,
                               MathLib::Point3d const& pnt, double eps);

    /**
     * Tests if the element is geometrically valid.
     */
    static ElementErrorCode validate(const Element* e);

    /// Returns the ID of a face given an array of nodes.
    static unsigned identifyFace(Node const* const*, Node* nodes[3]);

    /// Calculates the volume of the element
    static double computeVolume(Node const* const* _nodes);

}; /* class */

} /* namespace */
