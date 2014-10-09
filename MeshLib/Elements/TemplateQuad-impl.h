/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the Quad class-
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <array>

#include "MeshLib/Node.h"

#include "MathLib/MathTools.h"
#include "GeoLib/AnalyticalGeometry.h"

namespace MeshLib
{

template <unsigned NNODES, CellType CELLQUADTYPE>
const unsigned TemplateQuad<NNODES, CELLQUADTYPE>::n_all_nodes;

template <unsigned NNODES, CellType CELLQUADTYPE>
const unsigned TemplateQuad<NNODES, CELLQUADTYPE>::n_base_nodes;

template <unsigned NNODES, CellType CELLQUADTYPE>
const unsigned TemplateQuad<NNODES, CELLQUADTYPE>::_edge_nodes[4][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{2, 3}, // Edge 2
	{0, 3}  // Edge 3
};

template <unsigned NNODES, CellType CELLQUADTYPE>
TemplateQuad<NNODES,CELLQUADTYPE>::TemplateQuad(Node* nodes[NNODES], unsigned value, std::size_t id)
	: Face(value, id)
{
	_nodes = nodes;

	_neighbors = new Element*[4];
	std::fill(_neighbors, _neighbors + 4, nullptr);

	this->_area = this->computeVolume();
}

template<unsigned NNODES, CellType CELLQUADTYPE>
TemplateQuad<NNODES,CELLQUADTYPE>::TemplateQuad(std::array<Node*, NNODES> const& nodes,
                                                unsigned value, std::size_t id)
	: Face(value, id)
{
	_nodes = new Node*[NNODES];
	std::copy(nodes.begin(), nodes.end(), _nodes);

	_neighbors = new Element*[4];
	std::fill(_neighbors, _neighbors + 4, nullptr);

	this->_area = this->computeVolume();
}

template <unsigned NNODES, CellType CELLQUADTYPE>
TemplateQuad<NNODES,CELLQUADTYPE>::TemplateQuad(const TemplateQuad<NNODES,CELLQUADTYPE> &quad)
	: Face(quad.getValue(), quad.getID())
{
	_nodes = new Node*[NNODES];
	for (unsigned i=0; i<NNODES; i++) {
		_nodes[i] = quad._nodes[i];
	}

	_neighbors = new Element*[4];
	for (unsigned i=0; i<4; i++) {
		_neighbors[i] = quad._neighbors[i];
	}

	_area = quad.getArea();
}

template <unsigned NNODES, CellType CELLQUADTYPE>
TemplateQuad<NNODES,CELLQUADTYPE>::~TemplateQuad()
{
}

template <unsigned NNODES, CellType CELLQUADTYPE>
double TemplateQuad<NNODES,CELLQUADTYPE>::computeVolume()
{
	return GeoLib::calcTriangleArea(*_nodes[0], *_nodes[1], *_nodes[2])
         + GeoLib::calcTriangleArea(*_nodes[2], *_nodes[3], *_nodes[0]);
}

template <unsigned NNODES, CellType CELLQUADTYPE>
bool TemplateQuad<NNODES,CELLQUADTYPE>::isEdge(unsigned idx1, unsigned idx2) const
{
	for (unsigned i(0); i<4; i++)
	{
		if (_edge_nodes[i][0]==idx1 && _edge_nodes[i][1]==idx2) return true;
		if (_edge_nodes[i][1]==idx1 && _edge_nodes[i][0]==idx2) return true;
	}
	return false;
}

template <unsigned NNODES, CellType CELLQUADTYPE>
bool TemplateQuad<NNODES,CELLQUADTYPE>::isPntInElement(GeoLib::Point const& pnt, double eps) const
{
	return (GeoLib::isPointInTriangle(pnt, *_nodes[0], *_nodes[1], *_nodes[2], eps) ||
	        GeoLib::isPointInTriangle(pnt, *_nodes[0], *_nodes[2], *_nodes[3], eps));
}

template <unsigned NNODES, CellType CELLQUADTYPE>
Element* TemplateQuad<NNODES,CELLQUADTYPE>::clone() const
{
	return new TemplateQuad(*this);
}

template <unsigned NNODES, CellType CELLQUADTYPE>
unsigned TemplateQuad<NNODES,CELLQUADTYPE>::identifyFace(Node* nodes[3]) const
{
	for (unsigned i=0; i<4; i++)
	{
		unsigned flag(0);
		for (unsigned j=0; j<2; j++)
			for (unsigned k=0; k<2; k++)
				if (_nodes[_edge_nodes[i][j]] == nodes[k])
					flag++;
		if (flag==2)
			return i;
	}
	return std::numeric_limits<unsigned>::max();
}

template <unsigned NNODES, CellType CELLQUADTYPE>
ElementErrorCode TemplateQuad<NNODES,CELLQUADTYPE>::validate() const
{
	ElementErrorCode error_code;
	error_code[ElementErrorFlag::ZeroVolume]  = this->hasZeroVolume();
	error_code[ElementErrorFlag::NonCoplanar] = (!GeoLib::isCoplanar(*_nodes[0], *_nodes[1], *_nodes[2], *_nodes[3]));
	// for collapsed quads (i.e. reduced to a line) this test might result "false" as all four points are actually located on a line.
	if (!error_code[ElementErrorFlag::ZeroVolume]) 
		error_code[ElementErrorFlag::NonConvex]   = (!(GeoLib::dividedByPlane(*_nodes[0], *_nodes[2], *_nodes[1], *_nodes[3]) &&
			                                           GeoLib::dividedByPlane(*_nodes[1], *_nodes[3], *_nodes[0], *_nodes[2])));
	error_code[ElementErrorFlag::NodeOrder]  = !this->testElementNodeOrder();
	return error_code;
}

}

