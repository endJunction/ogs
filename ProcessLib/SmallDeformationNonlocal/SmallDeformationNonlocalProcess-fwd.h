/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "SmallDeformationNonlocalProcess.h"

extern template class ProcessLib::SmallDeformationNonlocal::
    SmallDeformationNonlocalProcess<2>;
extern template class ProcessLib::SmallDeformationNonlocal::
    SmallDeformationNonlocalProcess<3>;
