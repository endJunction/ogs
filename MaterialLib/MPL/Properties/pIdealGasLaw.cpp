/**
 * \author Norbert Grunwald
 * \date   18.09.2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "pIdealGasLaw.h"

namespace MaterialPropertyLib
{
IdealGasLaw::IdealGasLaw(Medium*) : _phase(nullptr), _component(nullptr)
{
    notImplemented("IdealGasLaw", "Medium");
}

IdealGasLaw::IdealGasLaw(Phase* p) : _phase(p), _component(nullptr){};

IdealGasLaw::IdealGasLaw(Component* c) : _phase(nullptr), _component(c){};

/**
 */
PropertyDataType IdealGasLaw::value(VariableArray const&)
{
    if (isUpdated())
        return _value;

    /// \todo: implementation of IdealGasLaw.
    return _value;
}
}  // MaterialPropertyLib
