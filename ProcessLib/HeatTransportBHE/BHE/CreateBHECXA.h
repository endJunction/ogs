/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include "ProcessLib/Process.h"
#include "BHE_CXA.h"

namespace BHE  // namespace of borehole heat exchanger
{
    BHE::BHE_CXA *
        CreateBHECXA(BaseLib::ConfigTree const& config,
                    BaseLib::ConfigTree const& bhe_conf, 
                    std::map<std::string,
                    std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const& curves);


}  // end of namespace
