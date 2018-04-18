/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshUtils.h"

#include "BaseLib/makeVectorUnique.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{

namespace HeatTransportBHE
{
void getBHEDataInMesh(
    MeshLib::Mesh const& mesh,
    std::vector<MeshLib::Element*>& vec_soil_elements,
    std::vector<int>& vec_BHE_mat_IDs,
    std::vector<std::vector<MeshLib::Element*>>& vec_BHE_elements,
    std::vector<std::vector<MeshLib::Node*>>& vec_BHE_nodes,
    std::vector<std::vector<MeshLib::Element*>>& vec_BHE_soil_elements)
{
    // IsCrackTip isCrackTip(mesh);

    // get vectors of matrix elements and fracture elements
    vec_soil_elements.reserve(mesh.getNumberOfElements());
    std::vector<MeshLib::Element*> all_BHE_elements;
    for (MeshLib::Element* e : mesh.getElements())
    {
        if (e->getDimension() == mesh.getDimension())
            vec_soil_elements.push_back(e);
        else if (e->getDimension() == (unsigned int)1)
            all_BHE_elements.push_back(e);
    }
    DBUG("-> found total %d soil elements and %d BHE elements",
         vec_soil_elements.size(), all_BHE_elements.size());

    // get BHE material IDs
    auto opt_material_ids(
        mesh.getProperties().getPropertyVector<int>("MaterialIDs"));
    for (MeshLib::Element* e : all_BHE_elements)
        vec_BHE_mat_IDs.push_back((*opt_material_ids)[e->getID()]);
    BaseLib::makeVectorUnique(vec_BHE_mat_IDs);
    DBUG("-> found %d BHE material groups", vec_BHE_mat_IDs.size());

    // create a vector of BHE elements for each group
    vec_BHE_elements.resize(vec_BHE_mat_IDs.size());
    for (unsigned bhe_id = 0; bhe_id < vec_BHE_mat_IDs.size(); bhe_id++)
    {
        const auto bhe_mat_id = vec_BHE_mat_IDs[bhe_id];
        std::vector<MeshLib::Element*>& vec_elements = vec_BHE_elements[bhe_id];
        std::copy_if(all_BHE_elements.begin(), all_BHE_elements.end(),
                     std::back_inserter(vec_elements),
                     [&](MeshLib::Element* e) {
                         return (*opt_material_ids)[e->getID()] == bhe_mat_id;
                     });
        DBUG("-> found %d elements on the BHE_%d", vec_elements.size(), bhe_id);
    }

    // get a vector of BHE nodes
    vec_BHE_nodes.resize(vec_BHE_mat_IDs.size());
    for (unsigned bhe_id = 0; bhe_id < vec_BHE_mat_IDs.size(); bhe_id++)
    {
        std::vector<MeshLib::Node*>& vec_nodes = vec_BHE_nodes[bhe_id];
        for (MeshLib::Element* e : vec_BHE_elements[bhe_id])
        {
            for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
            {
                vec_nodes.push_back(const_cast<MeshLib::Node*>(e->getNode(i)));
            }
        }
        BaseLib::makeVectorUnique(
            vec_nodes, [](MeshLib::Node* node1, MeshLib::Node* node2) {
                return node1->getID() < node2->getID();
            });
        DBUG("-> found %d nodes on the BHE_%d", vec_nodes.size(), bhe_id);
    }

    // get a vector of soil elements connected with BHEs
    vec_BHE_soil_elements.resize(vec_BHE_mat_IDs.size());
    for (unsigned bhe_id = 0; bhe_id < vec_BHE_mat_IDs.size(); bhe_id++)
    {
        std::vector<MeshLib::Element*>& vec_elements =
            vec_BHE_soil_elements[bhe_id];
        for (MeshLib::Node* node : vec_BHE_nodes[bhe_id])
        {
            for (unsigned i = 0; i < node->getNumberOfElements(); i++)
            {
                // make sure it is not the 1D BHE element
                if (node->getElement(i)->getDimension() != 1)
                    vec_elements.push_back(
                        const_cast<MeshLib::Element*>(node->getElement(i)));
            }
        }
        // remove any redundant elements
        BaseLib::makeVectorUnique(
            vec_elements, [](MeshLib::Element* elem1, MeshLib::Element* elem2) {
                return elem1->getID() < elem2->getID();
            });

        DBUG("-> found %d soil elements connected with BHE_%d",
             vec_elements.size(), bhe_id);
    }  // end of for
    }

}  // end of namespace HeatTransportBHE

}  // namespace ProcessLib


