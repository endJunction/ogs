/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <utility>
#include <vector>

#include "BaseLib/Error.h"

#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "ProcessLib/LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
class HeatTransportBHELocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    HeatTransportBHELocalAssemblerInterface() : _dofIndex_to_localIndex{} {}
    HeatTransportBHELocalAssemblerInterface(
        std::size_t n_local_size, std::vector<unsigned> dofIndex_to_localIndex)
        : _dofIndex_to_localIndex(std::move(dofIndex_to_localIndex))
    {
    }

private:

    std::vector<unsigned> const _dofIndex_to_localIndex;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
