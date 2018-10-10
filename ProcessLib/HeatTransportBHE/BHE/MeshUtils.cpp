/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshUtils.h"

#include "BaseLib/Algorithm.h"
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
    std::vector<MeshLib::Node*>& vec_soil_nodes,
    std::vector<std::vector<MeshLib::Node*>>& vec_BHE_nodes)
{
    // partition all the mesh elements, and copy them into
    // two seperate vectors, one with only matrix elements
    // and the other only BHE elements
    vec_soil_elements.reserve(mesh.getNumberOfElements());
    std::vector<MeshLib::Element*> all_BHE_elements;
    auto& const all_mesh_elements = mesh.getElements();
    std::partition_copy(
        std::begin(all_mesh_elements),
        std::end(all_mesh_elements),
        std::back_inserter(all_BHE_elements),
        std::back_inserter(vec_soil_elements),
        [](MeshLib::Element* e) { return e->getDimension() == 1; });

    // get BHE material IDs
    // TODO (haibing) Use MeshLib::materialIDs(), Explicitly abort on error.
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

    // get all the nodes into the pure soil nodes vector
    vec_soil_nodes.reserve(mesh.getNumberOfNodes());
    for (MeshLib::Node* n : mesh.getNodes())
    {
        // All elements are counted as a soil element
        vec_soil_nodes.push_back(n);
    }

    // finalLy counting two types of elements
    // They are (i) soil, and (ii) BHE type of elements
    DBUG("-> found total %d soil elements and %d BHE elements",
         vec_soil_elements.size(),
         all_BHE_elements.size());

    // TODO (haibing) return a struct with proper names of the result, don't
    // change the input arguments.
}
}  // end of namespace HeatTransportBHE
}  // namespace ProcessLib
