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

#include "boost/math/constants/constants.hpp"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE  // namespace of borehole heat exchanger
{
inline double pipeFlowVelocity(double const& flow_rate,
                               double const& pipe_radius)
{
    constexpr double pi = boost::math::constants::pi<double>();
    return flow_rate / (pi * pipe_radius * pipe_radius);
}

// Flow velocity in circular annulus defined by r_inner and r_outer.
inline double annulusFlowVelocity(double const flow_rate, double const r_outer,
                                  double const r_inner)
{
    constexpr double pi = boost::math::constants::pi<double>();
    return flow_rate / (pi * (r_outer * r_outer - r_inner * r_inner));
}

inline double prandtlNumber(double const& viscosity,
                            double const& heat_capacity,
                            double const& heat_conductivity)
{
    return viscosity * heat_capacity / heat_conductivity;
}

inline double reynoldsNumber(double const velocity_norm,
                             double const pipe_diameter,
                             double const viscosity,
                             double const density)
{
    return velocity_norm * pipe_diameter / (viscosity / density);
}

inline double nusseltNumber(double const reynolds_number,
                            double const prandtl_number,
                            double const pipe_diameter,
                            double const pipe_length)
{
    if (reynolds_number < 2300.0)
    {
        return 4.364;
    }
    if (reynolds_number < 10000.0)
    {
        double const gamma = (reynolds_number - 2300) / (10000 - 2300);

        return (1.0 - gamma) * 4.364 +
               gamma *
                   ((0.0308 / 8.0 * 1.0e4 * prandtl_number) /
                    (1.0 + 12.7 * std::sqrt(0.0308 / 8.0) *
                               (std::pow(prandtl_number, 2.0 / 3.0) - 1.0)) *
                    (1.0 + std::pow(pipe_diameter / pipe_length, 2.0 / 3.0)));
    }

    double const xi = std::pow(1.8 * std::log10(reynolds_number) - 1.5, -2.0);
    return (xi / 8.0 * reynolds_number * prandtl_number) /
           (1.0 + 12.7 * std::sqrt(xi / 8.0) *
                      (std::pow(prandtl_number, 2.0 / 3.0) - 1.0)) *
           (1.0 + std::pow(pipe_diameter / pipe_length, 2.0 / 3.0));
}

}  // end of namespace BHE
}  // end of namespace HeatTransportBHE
}  // end of namespace ProcessLib
