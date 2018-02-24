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

#include <cassert>

namespace NumLib
{
template <>
inline std::array<unsigned, 1>
IntegrationGaussLobattoRegular<1>::getPositionIndices(unsigned order,
                                                      unsigned igp)
{
    return {{igp}};
}

template <>
inline std::array<unsigned, 2>
IntegrationGaussLobattoRegular<2>::getPositionIndices(unsigned order,
                                                      unsigned igp)
{
    auto const n_points_per_dim = order + 2;
    assert(igp < n_points_per_dim * n_points_per_dim);
    std::array<unsigned, 2> result;
    result[0] = igp / n_points_per_dim;
    result[1] = igp % n_points_per_dim;
    return result;
}

template <>
inline std::array<unsigned, 3>
IntegrationGaussLobattoRegular<3>::getPositionIndices(unsigned order,
                                                      unsigned igp)
{
    auto const n_points_per_dim = order + 2;
    assert(igp < n_points_per_dim * n_points_per_dim * n_points_per_dim);
    unsigned const gp_r = igp / (n_points_per_dim * n_points_per_dim);
    unsigned const gp_s = igp % (n_points_per_dim * n_points_per_dim);
    std::array<unsigned, 3> result;
    result[0] = gp_r;
    result[1] = gp_s / n_points_per_dim;
    result[2] = gp_s % n_points_per_dim;
    return result;
}

template <unsigned N_DIM>
inline MathLib::TemplateWeightedPoint<double, double, N_DIM>
IntegrationGaussLobattoRegular<N_DIM>::getWeightedPoint(unsigned order,
                                                        unsigned igp)
{
    assert(igp < std::pow(order + 2, N_DIM));
    std::array<unsigned, N_DIM> const pos = getPositionIndices(order, igp);

    switch (order)
    {
        case 1:
            return getWeightedPoint<MathLib::GaussLobatto<1>>(pos);
        case 2:
            return getWeightedPoint<MathLib::GaussLobatto<2>>(pos);
        case 3:
            return getWeightedPoint<MathLib::GaussLobatto<3>>(pos);
        case 4:
            return getWeightedPoint<MathLib::GaussLobatto<4>>(pos);
    }

    return MathLib::TemplateWeightedPoint<double, double, N_DIM>(
        std::array<double, N_DIM>(), 0);
}

template <unsigned N_DIM>
template <typename Method>
inline MathLib::TemplateWeightedPoint<double, double, N_DIM>
IntegrationGaussLobattoRegular<N_DIM>::getWeightedPoint(
    std::array<unsigned, N_DIM> const& pos)
{
    std::array<double, N_DIM> coords;
    double weight = 1;
    for (unsigned d = 0; d < N_DIM; d++)
    {
        coords[d] = Method::X[pos[d]];
        weight *= Method::W[pos[d]];
    }

    return MathLib::TemplateWeightedPoint<double, double, N_DIM>(coords,
                                                                 weight);
}
}  // namespace NumLib
