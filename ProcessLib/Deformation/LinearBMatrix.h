/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Deformation/BMatrixPolicy.h"

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

    BMatrixType B = BMatrixType::Zero(
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value,
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

// Overload for 2D case.
template <int DisplacementDim,
          int NPOINTS,
          typename N_Type,
          typename DNDX_Type,
          typename std::enable_if_t<DisplacementDim == 2>* = nullptr>
MathLib::KelvinVector::KelvinVectorType<2> computeStrain(
    const Eigen::Ref<Eigen::Matrix<double, NPOINTS * 2, 1> const>& u,
    DNDX_Type const& dNdx,
    N_Type const& N,
    const double radius,
    const bool is_axially_symmetric)
{
    MathLib::KelvinVector::KelvinVectorType<2> strain;
    strain[0] = dNdx.row(0).dot(u.template segment<NPOINTS>(0 * NPOINTS));
    strain[1] = dNdx.row(1).dot(u.template segment<NPOINTS>(1 * NPOINTS));
    strain[2] = is_axially_symmetric
                    ? N.dot(u.template segment<NPOINTS>(0 * NPOINTS)) / radius
                    : 0;
    strain[3] = dNdx.row(1).dot(u.template segment<NPOINTS>(0 * NPOINTS)) /
                    std::sqrt(2.) +
                dNdx.row(0).dot(u.template segment<NPOINTS>(1 * NPOINTS)) /
                    std::sqrt(2.);

    return strain;
}

// Overload for 3D case.
template <int DisplacementDim,
          int NPOINTS,
          typename N_Type,
          typename DNDX_Type,
          typename std::enable_if_t<DisplacementDim == 3>* = nullptr>
MathLib::KelvinVector::KelvinVectorType<3> computeStrain(
    const Eigen::Ref<Eigen::Matrix<double, NPOINTS * 3, 1> const>& u,
    DNDX_Type const& dNdx,
    N_Type const& N,
    const double radius,
    const bool is_axially_symmetric)
{
    MathLib::KelvinVector::KelvinVectorType<3> strain;

    return strain;
}

// Overload for 2D case.
template <int DisplacementDim,
          int NPOINTS,
          typename N_Type,
          typename DNDX_Type,
          typename std::enable_if_t<DisplacementDim == 2>* = nullptr>
Eigen::Matrix<double, DisplacementDim * NPOINTS, 1> const computeInternalForces(
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& sigma,
    DNDX_Type const& dNdx,
    N_Type const& N,
    const double radius,
    const bool is_axially_symmetric)
{
    Eigen::Matrix<double, DisplacementDim * NPOINTS, 1> forces;
    forces.template segment<NPOINTS>(0 * NPOINTS) =
        (dNdx.row(0) * sigma[0] + dNdx.row(1) * (sigma[3] / std::sqrt(2.)))
            .transpose();
    forces.template segment<NPOINTS>(1 * NPOINTS) =
        (dNdx.row(1) * sigma[1] + dNdx.row(0) * (sigma[3] / std::sqrt(2.)))
            .transpose();

    return forces;
}

// Overload for 3D case.
template <int DisplacementDim,
          int NPOINTS,
          typename N_Type,
          typename DNDX_Type,
          typename std::enable_if_t<DisplacementDim == 3>* = nullptr>
Eigen::Matrix<double, DisplacementDim * NPOINTS, 1> const computeInternalForces(
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& sigma,
    DNDX_Type const& dNdx,
    N_Type const& N,
    const double radius,
    const bool is_axially_symmetric)
{
    Eigen::Matrix<double, DisplacementDim * NPOINTS, 1> forces;

    return forces;
}

// Overload for 2D case.
template <int DisplacementDim,
          int NPOINTS,
          typename N_Type,
          typename DNDX_Type,
          typename std::enable_if_t<DisplacementDim == 2>* = nullptr>
Eigen::Matrix<double, DisplacementDim * NPOINTS, 1> const computeBCBProduct(
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> const& C,
    DNDX_Type const& dNdx,
    N_Type const& N,
    const double radius,
    const bool is_axially_symmetric)
{
    // Intermediate matrix holding the B^T * C product.
    Eigen::Matrix<double, DisplacementDim * NPOINTS, 4> bc;
    for (int i = 0; i < 4; ++i)
    {
        // TODO if is_axially_symmetric
        bc.template block<NPOINTS, 1>(0 * NPOINTS, i) =
            dNdx.row(0) * C(0, i) + dNdx.row(1) * (C(3, i) / std::sqrt(2.));

        bc.template block<NPOINTS, 1>(1 * NPOINTS, i) =
            dNdx.row(0) * (C(3, i) / std::sqrt(2.)) + dNdx.row(1) * C(1, i);
    }

    Eigen::Matrix<double, DisplacementDim * NPOINTS, DisplacementDim * NPOINTS> bcb;
    forces.template segment<NPOINTS>(0 * NPOINTS) =
        (dNdx.row(0) * sigma[0] + dNdx.row(1) * (sigma[3] / std::sqrt(2.)))
            .transpose();
    forces.template segment<NPOINTS>(1 * NPOINTS) =
        (dNdx.row(1) * sigma[1] + dNdx.row(0) * (sigma[3] / std::sqrt(2.)))
            .transpose();

    return bcb;
}

// Overload for 3D case.
template <int DisplacementDim,
          int NPOINTS,
          typename N_Type,
          typename DNDX_Type,
          typename std::enable_if_t<DisplacementDim == 3>* = nullptr>
Eigen::Matrix<double, DisplacementDim * NPOINTS, 1> const computeInternalForces(
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& sigma,
    DNDX_Type const& dNdx,
    N_Type const& N,
    const double radius,
    const bool is_axially_symmetric)
{
    Eigen::Matrix<double, DisplacementDim * NPOINTS, 1> forces;

    return forces;
}
}  // namespace LinearBMatrix
}  // namespace ProcessLib
