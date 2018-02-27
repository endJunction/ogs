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

template <>
double const GaussLobatto<5>::X[7] = {
    -1, -0.83022389627856692987203221, -0.46884879347071421380377188,
    0,  0.4688487934707142138037719,   0.83022389627856692987203221,
    1};
template <>
double const GaussLobatto<5>::W[7] = {
    0.047619047619047619047619048, 0.2768260473615659480107004,
    0.43174538120986262341787102,  0.48761904761904761904761905,
    0.43174538120986262341787102,  0.27682604736156594801070041,
    0.047619047619047619047619048};

template <>
double const GaussLobatto<6>::X[8] = {-1,
                                      -0.87174014850960661533744576,
                                      -0.59170018143314230214451073,
                                      -0.20929921790247886876865726,
                                      0.20929921790247886876865726,
                                      0.59170018143314230214451073,
                                      0.87174014850960661533744576,
                                      1};

template <>
double const GaussLobatto<6>::W[8] = {
    0.035714285714285714285714286, 0.2107042271435060393829921,
    0.3411226924835043647642407,   0.41245879465870388156705297,
    0.41245879465870388156705297,  0.3411226924835043647642407,
    0.2107042271435060393829921,   0.035714285714285714285714286};

}  // namespace MathLib
