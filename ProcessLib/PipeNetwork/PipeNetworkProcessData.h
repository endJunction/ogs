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
            PipeNetworkProcessData(
                Parameter<double> const& pipe_diameter_,
                Parameter<double> const& pipe_roughness_,
                Parameter<double> const& minorloss_coefficient_)
                : pipe_diameter(pipe_diameter_),
                  pipe_roughness(pipe_roughness_),
                  minorloss_coefficient(minorloss_coefficient_)
            {
            }

            PipeNetworkProcessData(PipeNetworkProcessData&& other)
                : pipe_diameter(other.pipe_diameter),
                  pipe_roughness(other.pipe_roughness),
                  minorloss_coefficient(other.minorloss_coefficient)
            {
            }

            //! Copies are forbidden.
            PipeNetworkProcessData(PipeNetworkProcessData const&) = delete;

            //! Assignments are not needed.
            void operator=(PipeNetworkProcessData const&) = delete;

            //! Assignments are not needed.
            void operator=(PipeNetworkProcessData&&) = delete;

            Parameter<double> const& pipe_diameter;
            Parameter<double> const& pipe_roughness;
            Parameter<double> const& minorloss_coefficient;
        };

    }  // namespace PipeNetwork
}  // namespace ProcessLib
