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

#include "GaussLobattoTri.h"

namespace MathLib
{
template <>
const std::array<std::array<double, 2>, GaussLobattoTri<1>::NPoints>
    GaussLobattoTri<1>::X = {{{{1. / 3., 1. / 3.}}}};
template <>
double const GaussLobattoTri<1>::W[1] = {1.0};

const std::array<std::array<double, 2>, GaussLobattoTri<2>::NPoints>
    GaussLobattoTri<2>::X = {
        {{{1. / 6., 1. / 6.}}, {{2. / 3., 1. / 6.}}, {{1. / 6., 2. / 3.}}}};
double const GaussLobattoTri<2>::W[3] = {1. / 3., 1. / 3., 1. / 3.};

const std::array<std::array<double, 2>, GaussLobattoTri<3>::NPoints>
    GaussLobattoTri<3>::X = {{{{1. / 3., 1. / 3.}},
                              {{1. / 5., 3. / 5.}},
                              {{1. / 5., 1. / 5.}},
                              {{3. / 5., 1. / 5.}}}};
double const GaussLobattoTri<3>::W[4] = {-27. / 48., 25. / 48., 25. / 48.,
                                         25. / 48.};

}  // namespace MathLib
