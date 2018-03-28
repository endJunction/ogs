/**
* \copyright
* Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

namespace ProcessLib
{
    template <typename T>
    struct Parameter;

    namespace PipeNetwork
    {

        struct PipeNetworkProcessData
        {
            PipeNetworkProcessData(Parameter<double> const& thermal_conductivity_,
                Parameter<double> const& heat_capacity_,
                Parameter<double> const& density_)
                : thermal_conductivity(thermal_conductivity_),
                heat_capacity(heat_capacity_),
                density(density_)
            {
            }

            PipeNetworkProcessData(PipeNetworkProcessData&& other)
                : thermal_conductivity(other.thermal_conductivity),
                heat_capacity(other.heat_capacity),
                density(other.density)
            {
            }

            //! Copies are forbidden.
            PipeNetworkProcessData(PipeNetworkProcessData const&) = delete;

            //! Assignments are not needed.
            void operator=(PipeNetworkProcessData const&) = delete;

            //! Assignments are not needed.
            void operator=(PipeNetworkProcessData&&) = delete;

            // TODO
            // these 3 parameters need to be changed accordingly. 
            Parameter<double> const& thermal_conductivity;
            Parameter<double> const& heat_capacity;
            Parameter<double> const& density;
            
        };

    }  // namespace PipeNetwork
}  // namespace ProcessLib
