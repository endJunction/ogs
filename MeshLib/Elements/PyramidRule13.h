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
#include "PyramidRule5.h"

namespace MeshLib
{

/**
 * This class represents a 3d pyramid element. The following sketch shows the node and edge numbering.
 * @anchor PyramidNodeAndEdgeNumbering
 * @code
 *
 *               4
 *             //|\
 *            // | \
 *         11//  |  \10
 *          //   |9  \
 *         //    |    \
 *        3/.... |.....2
 *       ./      |  7 /
 *      ./12     |   /
 *    8./        |  /6
 *    ./         | /
 *   ./          |/
 *  0------------1
 *        5
 * @endcode

 */
class PyramidRule13 : public PyramidRule5
{
public:
    /// Constant: The number of all nodes for this element
    static constexpr unsigned n_all_nodes = 13u;

    /// Constant: The FEM type of the element
    static constexpr CellType cell_type = CellType::PYRAMID13;

    /// Constant: Local node index table for faces
    static constexpr unsigned face_nodes[5][8] = {
        {0, 1, 4, 5, 10, 9, 99, 99},   // Face 0
        {1, 2, 4, 6, 11, 10, 99, 99},  // Face 1
        {2, 3, 4, 7, 12, 11, 99, 99},  // Face 2
        {3, 0, 4, 8, 9, 12, 99, 99},   // Face 3
        {0, 3, 2, 1, 8, 7, 6, 5}       // Face 4
    };

    /// Constant: Local node index table for edge
    static constexpr unsigned edge_nodes[8][3] = {
        {0, 1, 5},   // Edge 0
        {1, 2, 6},   // Edge 1
        {2, 3, 7},   // Edge 2
        {0, 3, 8},   // Edge 3
        {0, 4, 9},   // Edge 4
        {1, 4, 10},  // Edge 5
        {2, 4, 11},  // Edge 6
        {3, 4, 12}   // Edge 7
    };

    /// Constant: Table for the number of nodes for each face
    static constexpr unsigned n_face_nodes[5] = {6, 6, 6, 6, 8};

    /// Returns the i-th edge of the element.
    using EdgeReturn = MeshLib::QuadraticEdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

}; /* class */

} /* namespace */
