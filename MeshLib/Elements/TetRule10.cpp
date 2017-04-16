/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TetRule10.h"

#include <array>

#include <logog/include/logog.hpp>

#include "MeshLib/Node.h"
#include "Tri.h"

namespace MeshLib
{
constexpr unsigned TetRule10::face_nodes[4][6];
constexpr unsigned TetRule10::edge_nodes[6][3];

const Element* TetRule10::getFace(const Element* e, unsigned i)
{
    if (i<n_faces)
    {
        std::array<Node*,6> nodes;
        for (unsigned j=0; j<6; j++)
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));
        return new Tri6(nodes);
    }
    ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
    return nullptr;
}

} // end namespace MeshLib
