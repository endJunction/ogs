/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HexRule20.h"

#include <array>

#include <logog/include/logog.hpp>

#include "MeshLib/Node.h"
#include "Quad.h"
#include "Line.h"

namespace MeshLib {
constexpr unsigned HexRule20::face_nodes[6][8];
constexpr unsigned HexRule20::edge_nodes[12][3];

const Element* HexRule20::getFace(const Element* e, unsigned i)
{
    if (i < n_faces)
    {
        std::array<Node*, 8> nodes;
        for (unsigned j=0; j<8; j++)
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));
        return new Quad8(nodes);
    }
    ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
    return nullptr;
}

} // end namespace MeshLib
