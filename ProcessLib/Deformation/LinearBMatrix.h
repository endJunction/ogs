/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cmath>

namespace ProcessLib
{
namespace LinearBMatrix
{
namespace detail
{
template <int NPOINTS, typename DNDX_Type, typename BMatrixType>
void fillBMatrix2DCartesianPart(DNDX_Type const& dNdx, BMatrixType& B)
{
    for (int i = 0; i < NPOINTS; ++i)
    {
        B(1, NPOINTS + i) = dNdx(1, i);
        B(3, i) = dNdx(1, i) / std::sqrt(2);
        B(3, NPOINTS + i) = dNdx(0, i) / std::sqrt(2);
        B(0, i) = dNdx(0, i);
    }
}
}  // detail

/// Fills a B-matrix based on given shape function dN/dx values.
template <int DisplacementDim,
          int NPOINTS,
          typename BMatrixType,
          typename N_Type,
          typename DNDX_Type>
BMatrixType computeBMatrix(DNDX_Type const& dNdx,
                           N_Type const& N,
                           const double radius,
                           const bool is_axially_symmetric)
{
    static_assert(0 < DisplacementDim && DisplacementDim <= 3,
                  "LinearBMatrix::computeBMatrix: DisplacementDim must be in "
                  "range [1,3].");

    BMatrixType B =
        BMatrixType::Zero(KelvinVectorDimensions<DisplacementDim>::value,
                          NPOINTS * DisplacementDim);

    switch (DisplacementDim)
    {
        case 3:
            for (int i = 0; i < NPOINTS; ++i)
            {
                B(2, 2 * NPOINTS + i) = dNdx(2, i);
                B(4, NPOINTS + i) = dNdx(2, i) / std::sqrt(2);
                B(4, 2 * NPOINTS + i) = dNdx(1, i) / std::sqrt(2);
                B(5, i) = dNdx(2, i) / std::sqrt(2);
                B(5, 2 * NPOINTS + i) = dNdx(0, i) / std::sqrt(2);
            }
            detail::fillBMatrix2DCartesianPart<NPOINTS>(dNdx, B);
            break;
        case 2:
            detail::fillBMatrix2DCartesianPart<NPOINTS>(dNdx, B);
            if (is_axially_symmetric)
            {
                for (int i = 0; i < NPOINTS; ++i)
                {
                    B(2, i) = N[i] / radius;
                }
            }
            break;
        default:
            break;
    }

    return B;
}

/// Fills a G-matrix based on given shape function dN/dx values.
template <int DisplacementDim,
          int NPOINTS,
          typename N_Type,
          typename DNDX_Type,
          typename BMatrixType>
void computeGMatrix(DNDX_Type const& dNdx,
                    BMatrixType& g_matrix,
                    const bool is_axially_symmetric,
                    N_Type const& N,
                    const double radius)
{
    static_assert(0 < DisplacementDim && DisplacementDim <= 3,
                  "LinearGMatrix::computeGMatrix: DisplacementDim must be in "
                  "range [1,3].");

    g_matrix.setZero();

    switch (DisplacementDim)
    {
        case 3:
            // The gradient coordinates are organized in the following order:
            // (1,1), (1,2), (1,3)
            // (2,1), (2,2), (2,3)
            // (3,1), (3,2), (3,3)
            for (int d = 0; d < DisplacementDim; ++d){
                for (int i = 0; i < NPOINTS; ++i)
                {
                    g_matrix(d + 0 * DisplacementDim, i + 0 * NPOINTS) =
                        dNdx(d, i);
                    g_matrix(d + 1 * DisplacementDim, i + 1 * NPOINTS) =
                        dNdx(d, i);
                    g_matrix(d + 2 * DisplacementDim, i + 2 * NPOINTS) =
                        dNdx(d, i);
                }
            }
            break;
        case 2:
            // The gradient coordinates are organized in the following order:
            // (1,1), (1,2)
            // (2,1), (2,2)
            // (3,3)
            for (int d = 0; d < DisplacementDim; ++d){
                for (int i = 0; i < NPOINTS; ++i)
                {
                    g_matrix(d, i) = dNdx(d, i);
                    g_matrix(d + DisplacementDim, i + NPOINTS) = dNdx(d, i);
                }
            }
            if (is_axially_symmetric)
            {
                for (int i = 0; i < NPOINTS; ++i)
                {
                    g_matrix(4, i) = N[i] / radius;
                }
            }
            break;
        default:
            break;
    }
}

}  // namespace LinearBMatrix
}  // namespace ProcessLib
