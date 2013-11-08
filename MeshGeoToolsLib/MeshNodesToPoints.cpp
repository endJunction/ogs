/**
 * @file MeshNodesToPoints.cpp
 * @date Nov 7, 2013
 * @brief
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "MeshNodesToPoints.h"

// GeoLib
#include "Point.h"

// MeshLib
#include "Node.h"

namespace MeshGeoTools
{

MeshNodesToPoints::MeshNodesToPoints(MeshLib::Mesh const& mesh,
		std::vector<std::size_t> const& mesh_node_ids,
		GeoLib::GEOObjects &geo_obj,
		std::string &geo_name)
{
	std::vector<MeshLib::Node*> const& nodes(mesh.getNodes());
	std::vector<GeoLib::Point*> *points(new std::vector<GeoLib::Point*>);
	points->reserve(mesh_node_ids.size());
	for (auto it(mesh_node_ids.cbegin()); it != mesh_node_ids.cend(); it++) {
		points->push_back(new GeoLib::Point(nodes[*it]->getCoords()));
	}

	geo_obj.addPointVec(points, geo_name, nullptr, 0.0);
}

} // end namespace MeshGeoTools
