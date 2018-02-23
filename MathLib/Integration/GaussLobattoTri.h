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

#include <array>

#include "mathlib_export.h"

namespace MathLib
{
/// Gauss-Lobatto quadrature on triangles
///
/// \tparam ORDER   integration order.
template <unsigned ORDER>
struct GaussLobattoTri
{
    static MATHLIB_EXPORT const unsigned Order = ORDER;
    static MATHLIB_EXPORT const unsigned NPoints = ORDER;
    static MATHLIB_EXPORT const std::array<std::array<double, 2>, NPoints> X;
    static MATHLIB_EXPORT const double W[NPoints];
};

template <>
struct GaussLobattoTri<2>
{
    static MATHLIB_EXPORT const unsigned Order = 2;
    static MATHLIB_EXPORT const unsigned NPoints = 3;
    static MATHLIB_EXPORT const std::array<std::array<double, 2>, NPoints> X;
    static MATHLIB_EXPORT const double W[NPoints];
};

template <>
struct GaussLobattoTri<3>
{
    static MATHLIB_EXPORT const unsigned Order = 3;
    static MATHLIB_EXPORT const unsigned NPoints = 4;
    static MATHLIB_EXPORT const std::array<std::array<double, 2>, NPoints> X;
    static MATHLIB_EXPORT const double W[NPoints];
};

#ifndef _MSC_VER  // The following explicit instantatiation declaration does not
                  // compile on that particular compiler but is necessary.
template <>
const std::array<std::array<double, 2>, GaussLobattoTri<1>::NPoints>
    GaussLobattoTri<1>::X;
template <>
double const GaussLobattoTri<1>::W[1];
#endif
}  // namespace MathLib
