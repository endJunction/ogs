/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Parameter.h"

namespace MaterialPropertyLib
{
Parameter::Parameter(ParameterLib::Parameter<double> const& parameter)
    : _parameter(parameter)
{
    _dvalue = boost::apply_visitor(
        [](auto const& value) -> PropertyDataType { return decltype(value){}; },
        v);
};
}  // namespace MaterialPropertyLib
