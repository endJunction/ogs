/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
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

/// Gauss-Legendre quadrature method
///
/// \tparam ORDER   integration order.
template <unsigned ORDER>
struct GaussLegendre { };

template <>
struct GaussLegendre<1>
{
    static constexpr unsigned Order = 1;
    static constexpr double X[1] = {0.};
    static constexpr double W[1] = {2.};
};

template<>
struct GaussLegendre<2>
{
    static constexpr unsigned Order = 2;
    static constexpr double X[2] = {0.577350269189626, -0.577350269189626};
    static constexpr double W[2] = {1., 1.};
};

template<>
struct GaussLegendre<3>
{
    static constexpr unsigned Order = 3;
    static constexpr double X[3] = {0.774596669241483, 0., -0.774596669241483};
    static constexpr double W[3] = {5./9, 8./9, 5./9};
};

template<>
struct GaussLegendre<4>
{
    static constexpr unsigned Order = 4;
    static constexpr double X[4] =
        {-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053};
    static constexpr double W[4] =
        { 0.347854845137454,  0.652145154862546, 0.652145154862546, 0.347854845137454};
};

}  // namespace MathLib
