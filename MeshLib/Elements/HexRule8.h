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
 * A 8-nodes Hexahedron Element.
 * @code
 *
 *  Hex:
 *                6
 *          7-----------6
 *         /:          /|
 *        / :         / |
 *      7/  :        /5 |
 *      / 11:       /   | 10
 *     /    : 4    /    |
 *    4-----------5     |
 *    |     :     | 2   |
 *    |     3.....|.....2
 *    |    .      |    /
 *  8 |   .       |9  /
 *    | 3.        |  / 1
 *    | .         | /
 *    |.          |/
 *    0-----------1
 *          0
 *
 * @endcode
 */
class HexRule8 : public CellRule
{
public:
    /// Constant: The number of base nodes for this element
    static constexpr unsigned n_base_nodes = 8u;

    /// Constant: The number of all nodes for this element
    static constexpr unsigned n_all_nodes = 8u;

    /// Constant: The geometric type of the element
    static constexpr MeshElemType mesh_elem_type = MeshElemType::HEXAHEDRON;

    /// Constant: The FEM type of the element
    static constexpr CellType cell_type = CellType::HEX8;

    /// Constant: The number of faces
    static constexpr unsigned n_faces = 6;

    /// Constant: The number of edges
    static constexpr unsigned n_edges = 12;

    /// Constant: The number of neighbors
    static constexpr unsigned n_neighbors = 6;

    /// Constant: Local node index table for faces
    static constexpr unsigned face_nodes[6][4] = {
        {0, 3, 2, 1},  // Face 0
        {0, 1, 5, 4},  // Face 1
        {1, 2, 6, 5},  // Face 2
        {2, 3, 7, 6},  // Face 3
        {3, 0, 4, 7},  // Face 4
        {4, 5, 6, 7}   // Face 5
    };

    /// Constant: Local node index table for edge
    static constexpr unsigned edge_nodes[12][2] = {
        {0, 1},  // Edge 0
        {1, 2},  // Edge 1
        {2, 3},  // Edge 2
        {0, 3},  // Edge 3
        {4, 5},  // Edge 4
        {5, 6},  // Edge 5
        {6, 7},  // Edge 6
        {4, 7},  // Edge 7
        {0, 4},  // Edge 8
        {1, 5},  // Edge 9
        {2, 6},  // Edge 10
        {3, 7}   // Edge 11
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

    /// Calculates the volume of a convex hexahedron by partitioning it into six tetrahedra.
    static double computeVolume(Node const* const* _nodes);

}; /* class */

} /* namespace */
