/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <string>
#include <map>

#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
struct NonequilibriumInitialState
{
    std::map<std::string, ParameterLib::Parameter<double> const*> const
        parameter_map;
};
}  // namespace ProcessLib
