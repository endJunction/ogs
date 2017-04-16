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
#include "PrismRule6.h"

namespace MeshLib
{

/**
 * This class represents a 3d prism element with 15 nodes. The following sketch shows the node and edge numbering.
 * @anchor PrismNodeAndEdgeNumbering
 * @code
 *            5
 *           / \
 *          / : \
 *       11/  :  \10
 *        /   :14 \
 *       /    :    \
 *      3------9----4
 *      |     :     |
 *      |     2     |
 *      |    . .    |
 *    12|   .   .   |13
 *      | 8.     .7 |
 *      | .       . |
 *      |.         .|
 *      0-----------1
 *            6
 *
 * @endcode
 */
class PrismRule15 : public PrismRule6
{
public:
    /// Constant: The number of all nodes for this element
    static constexpr unsigned n_all_nodes = 15u;

    /// Constant: The FEM type of the element
    static constexpr CellType cell_type = CellType::PRISM15;

    /// Constant: Local node index table for faces
    static constexpr unsigned face_nodes[5][8] = {
        {0, 2, 1, 8, 7, 6, 99, 99},    // Face 0
        {0, 1, 4, 3, 6, 10, 12, 9},    // Face 1
        {1, 2, 5, 4, 7, 11, 13, 10},   // Face 2
        {2, 0, 3, 5, 8, 9, 14, 11},    // Face 3
        {3, 4, 5, 12, 13, 14, 99, 99}  // Face 4
    };

    /// Constant: Local node index table for edge
    static constexpr unsigned edge_nodes[9][3] = {
        {0, 1, 6},   // Edge 0
        {1, 2, 7},   // Edge 1
        {0, 2, 8},   // Edge 2
        {0, 3, 9},   // Edge 3
        {1, 4, 10},  // Edge 4
        {2, 5, 11},  // Edge 5
        {3, 4, 12},  // Edge 6
        {4, 5, 13},  // Edge 7
        {3, 5, 14}   // Edge 8
    };

    /// Constant: Table for the number of nodes for each face
    static constexpr unsigned n_face_nodes[5] = {6, 8, 8, 8, 6};

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::QuadraticEdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

}; /* class */

} /* namespace */
