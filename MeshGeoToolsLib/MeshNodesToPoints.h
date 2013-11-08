/**
 * @file MeshNodesToPoints.h
 * @date Nov 7, 2013
 * @brief
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */
#ifndef MESHNODESTOPOINTS_H_
#define MESHNODESTOPOINTS_H_

#include <vector>
#include <string>

// GeoLib
#include "GEOObjects.h"

// MeshLib
#include "Mesh.h"

namespace MeshGeoTools
{

/**
 * The class extracts mesh nodes to GeoLib::Point objects. This
 * can be used, for example, for visualization or debugging purposes.
 */
class MeshNodesToPoints
{
public:
	/**
	 * Extracts mesh nodes to GeoLib::Point objects. The points are
	 * stored in a geometry within an instance of the class
	 * GeoLib::GEOObjects.
	 * @param mesh The mesh containing the mesh nodes.
	 * @param mesh_node_ids The ids of the mesh nodes that should be
	 * extracted.
	 * @param geo_obj GeoLib::GEOObjects instance managing the created points
	 * @param geo_name the name of the new created geometry within the
	 * GeoLib::GEOObjects instance
	 */
	MeshNodesToPoints(MeshLib::Mesh const& mesh,
			std::vector<std::size_t> const& mesh_node_ids,
			GeoLib::GEOObjects &geo_obj,
			std::string &geo_name);
};

} // end namespace MeshGeoTools

#endif /* MESHNODESTOPOINTS_H_ */
