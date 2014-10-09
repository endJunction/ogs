/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the TemplatePyramid class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// Thirdparty
#include "logog/include/logog.hpp"

#include "GeoLib/AnalyticalGeometry.h"
#include "MeshLib/Node.h"

#include "Tri.h"
#include "Tet.h"
#include "Quad.h"

namespace MeshLib {

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
const unsigned TemplatePyramid<NNODES, CELLPYRAMIDTYPE>::n_all_nodes;

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
const unsigned TemplatePyramid<NNODES, CELLPYRAMIDTYPE>::n_base_nodes;

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
const unsigned TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::_face_nodes[5][4] =
{
	{0, 1, 4, 99}, // Face 0
	{1, 2, 4, 99}, // Face 1
	{2, 3, 4, 99}, // Face 2
	{3, 0, 4, 99}, // Face 3
	{0, 3, 2,  1}  // Face 4
};

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
const unsigned TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::_edge_nodes[8][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{2, 3}, // Edge 2
	{0, 3}, // Edge 3
	{0, 4}, // Edge 4
	{1, 4}, // Edge 5
	{2, 4}, // Edge 6
	{3, 4}  // Edge 7
};

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
const unsigned TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::_n_face_nodes[5] = { 3, 3, 3, 3, 4 };

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::TemplatePyramid(Node* nodes[NNODES], unsigned value, std::size_t id)
	: Cell(value, id)
{
	_nodes = nodes;

	_neighbors = new Element*[5];
	std::fill(_neighbors, _neighbors + 5, nullptr);

	this->_volume = this->computeVolume();
}

template<unsigned NNODES, CellType CELLPYRAMIDTYPE>
TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::TemplatePyramid(std::array<Node*, NNODES> const& nodes,
                                                         unsigned value, std::size_t id)
	: Cell(value, id)
{
	_nodes = new Node*[NNODES];
	std::copy(nodes.begin(), nodes.end(), _nodes);

	_neighbors = new Element*[5];
	std::fill(_neighbors, _neighbors + 5, nullptr);

	this->_volume = this->computeVolume();
}

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::TemplatePyramid(const TemplatePyramid<NNODES,CELLPYRAMIDTYPE> &pyramid)
	: Cell(pyramid.getValue(), pyramid.getID())
{
	_nodes = new Node*[NNODES];
	for (unsigned i=0; i<NNODES; i++) {
		_nodes[i] = pyramid._nodes[i];
	}

	_neighbors = new Element*[5];
	for (unsigned i=0; i<5; i++) {
		_neighbors[i] = pyramid._neighbors[i];
	}

	_volume = pyramid.getVolume();
}

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::~TemplatePyramid()
{
}

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
double TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::computeVolume()
{
	return GeoLib::calcTetrahedronVolume(_nodes[0]->getCoords(), _nodes[1]->getCoords(), _nodes[2]->getCoords(), _nodes[4]->getCoords())
		 + GeoLib::calcTetrahedronVolume(_nodes[2]->getCoords(), _nodes[3]->getCoords(), _nodes[0]->getCoords(), _nodes[4]->getCoords());
}

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
const Element* TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::getFace(unsigned i) const
{
	if (i<this->getNFaces())
	{
		unsigned nFaceNodes (this->getNFaceNodes(i));
		Node** nodes = new Node*[nFaceNodes];
		for (unsigned j=0; j<nFaceNodes; j++)
			nodes[j] = _nodes[_face_nodes[i][j]];

		if (i<4)
			return new Tri(nodes);
		else
			return new Quad(nodes);
	}
	ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
	return NULL;
}

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
unsigned TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::getNFaceNodes(unsigned i) const
{
	if (i<5)
		return _n_face_nodes[i];
	ERR("Error in MeshLib::Element::getNFaceNodes() - Index %d does not exist.", i);
	return 0;
}

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
bool TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::isEdge(unsigned idx1, unsigned idx2) const
{
	for (unsigned i(0); i<8; i++)
	{
		if (_edge_nodes[i][0]==idx1 && _edge_nodes[i][1]==idx2) return true;
		if (_edge_nodes[i][1]==idx1 && _edge_nodes[i][0]==idx2) return true;
	}
	return false;
}

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
bool TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::isPntInElement(GeoLib::Point const& pnt, double eps) const
{
    return (GeoLib::isPointInTetrahedron(pnt, *_nodes[0], *_nodes[1], *_nodes[2], *_nodes[4], eps) ||
            GeoLib::isPointInTetrahedron(pnt, *_nodes[0], *_nodes[2], *_nodes[3], *_nodes[4], eps));
}

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
Element* TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::clone() const
{
	return new TemplatePyramid<NNODES,CELLPYRAMIDTYPE>(*this);
}

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
unsigned TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::identifyFace(Node* nodes[3]) const
{
	for (unsigned i=0; i<5; i++)
	{
		unsigned flag(0);
		for (unsigned j=0; j<4; j++)
			for (unsigned k=0; k<3; k++)
				if (_face_nodes[i][j] != 99 && _nodes[_face_nodes[i][j]] == nodes[k])
					flag++;
		if (flag==3)
			return i;
	}
	return std::numeric_limits<unsigned>::max();
}

template <unsigned NNODES, CellType CELLPYRAMIDTYPE>
ElementErrorCode TemplatePyramid<NNODES,CELLPYRAMIDTYPE>::validate() const
{
	ElementErrorCode error_code;
	error_code[ElementErrorFlag::ZeroVolume] = this->hasZeroVolume();

	const MeshLib::Quad* base (dynamic_cast<const MeshLib::Quad*>(this->getFace(4)));
	if (base)
	{
		error_code |= base->validate();
		error_code[ElementErrorFlag::NodeOrder] = !this->testElementNodeOrder();
	}
	else
		error_code.set(ElementErrorFlag::NodeOrder);
	delete base;

	return error_code;
}

} // end namespace MeshLib
