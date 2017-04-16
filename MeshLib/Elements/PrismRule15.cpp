/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PrismRule15.h"

#include <logog/include/logog.hpp>

#include "MeshLib/Node.h"
#include "Quad.h"
#include "Tri.h"

namespace MeshLib {

constexpr unsigned PrismRule15::face_nodes[5][8];
constexpr unsigned PrismRule15::edge_nodes[9][3];
constexpr unsigned PrismRule15::n_face_nodes[5];

const Element* PrismRule15::getFace(const Element* e, unsigned i)
{
    if (i < n_faces)
    {
        unsigned nFaceNodes(PrismRule15::n_face_nodes[i]);
        auto** nodes = new Node*[nFaceNodes];
        for (unsigned j=0; j<nFaceNodes; j++)
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));

        if (i == 0 || i == 4)
            return new Tri6(nodes);

        return new Quad8(nodes);
    }
    ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
    return nullptr;
}

} // end namespace MeshLib
