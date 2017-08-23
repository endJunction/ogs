/**
 * \file
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>

namespace NumLib
{
/// Create an element's points coordinates matrix.
/// \returns the NPoints-by-GlobalDim matrix containing the point coordinates of
/// the given element.
template <int NPoints, int GlobalDim>
Eigen::Matrix<double, NPoints, GlobalDim> createElementPointsCoordinatesMatrix(
    MeshLib::Element const& element)
{
    Eigen::Matrix<double, NPoints, GlobalDim> coordinates;
    auto const n_nodes = element.getNumberOfNodes();
    assert(NPoints == n_nodes);
    std::cout << "----\n";
    for (unsigned n = 0; n < n_nodes; ++n)
    {
        auto const& node = *element.getNode(n);
        for (int d = 0; d < GlobalDim; ++d)
            coordinates(n, d) = node[d];
        for (int d = 0; d < 3; ++d)
            std::cout << node[d] << " ";
    }
    std::cout << "\n";

    std::cout << GlobalDim << "-dim; for element " << element.getID() << "; nodes:\n"
              << coordinates << "\n";

    return coordinates;
}

/// Computes the normal vector of an element for the given dN/dr and the
/// coordinates matrix.
template <int GlobalDim>
struct ComputeNormalVector
{
    template <typename GlobalDimVector,
              typename CoordinatesMatrix,
              typename DNDRMatrix>
    static GlobalDimVector calculate(DNDRMatrix const& /*dNdr*/,
                                     CoordinatesMatrix const& /*coordinates*/)
    {
        // TODO: Remove unnecessary instantiations.
        OGS_FATAL("Invalid for the given GlobalDim.");
    }
};

template <>
struct ComputeNormalVector<3>
{
    template <typename GlobalDimVector,
              typename CoordinatesMatrix,
              typename DNDRMatrix>
    static GlobalDimVector calculate(DNDRMatrix const& dNdr,
                                     CoordinatesMatrix const& coordinates)
    {
        auto const tangents = (dNdr * coordinates).eval();

        GlobalDimVector const normal =
            tangents.row(0).template cross(tangents.row(1));
        return normal / normal.norm();
    }
};

template <>
struct ComputeNormalVector<2>
{
    template <typename GlobalDimVector,
              typename CoordinatesMatrix,
              typename DNDRMatrix>
    static GlobalDimVector calculate(DNDRMatrix const& dNdr,
                                     CoordinatesMatrix const& coordinates)
    {
        auto const tangents = (dNdr * coordinates).eval();

        GlobalDimVector normal(2);
        normal << -tangents(1), tangents(0);
        return normal / normal.norm();
    }
};

}  // NumLib
