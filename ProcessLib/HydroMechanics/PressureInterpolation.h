/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/Fem/CoordinatesMapping/NaturalNodeCoordinates.h"

namespace ProcessLib
{
namespace HydroMechanics
{
// For each higher order node evaluate the shape matrices for the lower
// order element (the base nodes)
template <typename ShapeMatricesTypePressure, typename ShapeFunctionPressure,
          typename MeshElement, int GlobalDim>
void interpolatePressureOnLowerOrderElement(
    typename ShapeMatricesTypePressure::NodalVectorType const& p,
    MeshLib::Element const& element,
    bool const is_axially_symmetric,
    MeshLib::PropertyVector<double>& nodal_pressure)
{
    using FemType = NumLib::TemplateIsoparametric<ShapeFunctionPressure,
                                                  ShapeMatricesTypePressure>;

    FemType fe(*static_cast<const typename ShapeFunctionPressure::MeshElement*>(
        &element));
    int const number_base_nodes = element.getNumberOfBaseNodes();
    int const number_all_nodes = element.getNumberOfNodes();

    for (int n = 0; n < number_base_nodes; ++n)
    {
        std::size_t const global_index = element.getNodeIndex(n);
        nodal_pressure[global_index] = p[n];
    }

    for (int n = number_base_nodes; n < number_all_nodes; ++n)
    {
        // Evaluated at higher order nodes' coordinates.
        typename ShapeMatricesTypePressure::ShapeMatrices shape_matrices_p{
            ShapeFunctionPressure::DIM, GlobalDim,
            ShapeFunctionPressure::NPOINTS};

        fe.computeShapeFunctions(
            NumLib::NaturalCoordinates<MeshElement>::coordinates[n].data(),
            shape_matrices_p, GlobalDim, is_axially_symmetric);

        auto const& N_p = shape_matrices_p.N;

        std::size_t const global_index = element.getNodeIndex(n);
        nodal_pressure[global_index] = N_p * p;
    }
};
}  // namespace HydroMechanics
}  // namespace ProcessLib
