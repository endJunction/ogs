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

#include "GaussLobatto.h"

namespace MathLib
{
template <>
double const GaussLobatto<1>::X[3] = {-1., 0., 1.};
template <>
double const GaussLobatto<1>::W[3] = {1. / 3, 4. / 3, 1. / 3};

template <>
double const GaussLobatto<2>::X[4] = {-1., -1. / 5 * std::sqrt(5.),
                                      1. / 5 * std::sqrt(5.), 1.};

template <>
double const GaussLobatto<2>::W[4] = {1. / 6, 5. / 6, 5. / 6, 1. / 6};

template <>
double const GaussLobatto<3>::X[5] = {-1., -1. / 7 * std::sqrt(21.), 0.,
                                      1. / 7 * std::sqrt(21.), 1};
template <>
double const GaussLobatto<3>::W[5] = {0.1, 49. / 90, 32. / 45, 49. / 90, 0.1};

template <>
double const GaussLobatto<4>::X[6] = {
    -1.,
    -std::sqrt(1. / 21 * (7 + 2 * std::sqrt(7))),
    -std::sqrt(1. / 21 * (7 - 2 * std::sqrt(7))),
    std::sqrt(1. / 21 * (7 - 2 * std::sqrt(7))),
    std::sqrt(1. / 21 * (7 + 2 * std::sqrt(7))),
    1.};
template <>
double const GaussLobatto<4>::W[6] = {1. / 15,
                                      1. / 30 * (14 - std::sqrt(7)),
                                      1. / 30 * (14 + std::sqrt(7)),
                                      1. / 30 * (14 + std::sqrt(7)),
                                      1. / 30 * (14 - std::sqrt(7)),
                                      1. / 15};

}  // namespace MathLib
