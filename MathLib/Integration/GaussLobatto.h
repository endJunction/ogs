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

#include "mathlib_export.h"

#include "WeightedSum.h"

namespace MathLib
{
/// Gauss-Lobatto quadrature method
///
/// \tparam ORDER   integration order.
template <unsigned ORDER>
struct GaussLobatto
{
    static MATHLIB_EXPORT const unsigned Order = ORDER;
    static MATHLIB_EXPORT const double X[Order + 2];
    static MATHLIB_EXPORT const double W[Order + 2];
};

#ifndef _MSC_VER  // The following explicit instantatiation declaration does not
                  // compile on that particular compiler.
template <>
double const GaussLobatto<1>::X[3];
template <>
double const GaussLobatto<1>::W[3];
template <>
double const GaussLobatto<2>::X[4];
template <>
double const GaussLobatto<2>::W[4];
template <>
double const GaussLobatto<3>::X[5];
template <>
double const GaussLobatto<3>::W[5];
template <>
double const GaussLobatto<4>::X[6];
template <>
double const GaussLobatto<4>::W[6];
template <>
double const GaussLobatto<5>::X[7];
template <>
double const GaussLobatto<5>::W[7];
template <>
double const GaussLobatto<6>::X[8];
template <>
double const GaussLobatto<6>::W[8];
#endif

}  // namespace MathLib
