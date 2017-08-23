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
#include <Eigen/Dense>

namespace NumLib
{
/// Create an element's points coordinates matrix.
/// \returns the NPoints-by-GlobalDim matrix containing the point coordinates of
/// the given element.
template <typename GlobalDimNodalMatrixType>
GlobalDimNodalMatrixType createElementPointsCoordinatesMatrix(
    MeshLib::Element const& element)
{
    GlobalDimNodalMatrixType coordinates;
    auto const n_nodes = element.getNumberOfNodes();
    assert(GlobalDimNodalMatrixType::RowsAtCompileTime == n_nodes);
    std::cout << "----\n";
    for (unsigned n = 0; n < n_nodes; ++n)
    {
        auto const& node = *element.getNode(n);
        for (int d = 0; d < GlobalDimNodalMatrixType::ColsAtCompileTime; ++d)
            coordinates(d, n) = node[d];
        for (int d = 0; d < 3; ++d)
            std::cout << node[d] << " ";
    }
    std::cout << "\n";

    std::cout << GlobalDimNodalMatrixType::ColsAtCompileTime
              << "-dim; for element " << element.getID() << "; nodes:\n"
              << coordinates << "\n";

    return coordinates;
}

/// Computes a unit-length normal vector of an element for the given dN/dr and
/// the coordinates matrix. For each of the three dimensions there are
/// specializations \see ComputeNormalVector<1>, ComputeNormalVector<2>, and
/// ComputeNormalVector<3>.
template <int GlobalDim>
struct ComputeNormalVector;

/// \copydoc ComputeNormalVector
///
/// The 3d-case computes the tangents and uses cross-product to compute the
/// resulting normal.
template <>
struct ComputeNormalVector<3>
{
    template <typename GlobalDimVector,
              typename CoordinatesMatrix,
              typename DNDRMatrix>
    static Eigen::Matrix<double,
                         DNDRMatrix::ColsAtCompileTime,
                         CoordinatesMatrix::RowsAtCompileTime>
    calculate(DNDRMatrix const& dNdr, CoordinatesMatrix const& coordinates)
    {
        // The following cross product call requires vectors of fixed size 3.
        Eigen::Matrix<double, 33, 3> const tangents =
            dNdr.transpose() * coordinates;

        GlobalDimVector const normal =
            tangents.row(0).template cross(tangents.row(1));
        return normal / normal.norm();
    }
};

/// \copydoc ComputeNormalVector
///
/// The 2d-case computes the tangent (along single local coordinate) and rotates
/// it 90 degrees to obtain a normal vector.
template <>
struct ComputeNormalVector<2>
{
    template <typename GlobalDimVector,
              typename CoordinatesMatrix,
              typename DNDRMatrix>
    static GlobalDimVector calculate(DNDRMatrix const& dNdr,
                                     CoordinatesMatrix const& coordinates)
    {
        auto const tangents = (dNdr.transpose() * coordinates).eval();

        GlobalDimVector normal(2);
        normal << -tangents(1), tangents(0);
        return normal / normal.norm();
    }
};

/// \copydoc ComputeNormalVector
///
/// In 1d (i.e. for the end points of lines) a 1d unity vector is always
/// returned.
template <>
struct ComputeNormalVector<1>
{
    template <typename GlobalDimVector,
              typename CoordinatesMatrix,
              typename DNDRMatrix>
    static GlobalDimVector calculate(DNDRMatrix const& /*dNdr*/,
                                     CoordinatesMatrix const& /*coordinates*/)
    {
        GlobalDimVector normal(1);
        normal << 1.0;
        return normal;
    }
};
}  // NumLib
