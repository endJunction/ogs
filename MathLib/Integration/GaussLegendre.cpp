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

#include "GaussLegendre.h"

namespace MathLib
{
constexpr unsigned GaussLegendre<1>::Order;
constexpr double GaussLegendre<1>::X[1];
constexpr double GaussLegendre<1>::W[1];

constexpr unsigned GaussLegendre<2>::Order;
constexpr double GaussLegendre<2>::X[2];
constexpr double GaussLegendre<2>::W[2];

constexpr unsigned GaussLegendre<3>::Order;
constexpr double GaussLegendre<3>::X[3];
constexpr double GaussLegendre<3>::W[3];

constexpr unsigned GaussLegendre<4>::Order;
constexpr double GaussLegendre<4>::X[4];
constexpr double GaussLegendre<4>::W[4];


}   // namespace MathLib
