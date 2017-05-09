/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SmallDeformationNonlocalProcess-fwd.h"
#include "SmallDeformationNonlocalProcess.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{

template class SmallDeformationNonlocalProcess<2>;
template class SmallDeformationNonlocalProcess<3>;

}   // namespace SmallDeformationNonlocal
}   // namespace ProcessLib
