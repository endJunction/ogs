/**
 * \author Norbert Grunwald
 * \date   27.06.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "MaterialLib/MPL/mpProperty.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

class BrooksCoreyRelPerm final : public Property
{
private:
    /// A pointer to the phase object.
    Medium* _medium;


public:
    /// Constructor passing a pointer to the medium.
    BrooksCoreyRelPerm(Medium*);
    /// Constructor passing a pointer to a phase.
    BrooksCoreyRelPerm(Phase*);
    /// Constructor passing a pointer to a component.
    BrooksCoreyRelPerm(Component*);

    /// This method overrides the base class implementation and
    /// actually computes and sets the property _value.
    PropertyDataType value(VariableArray const&) override;
};


}  // MaterialPropertyLib
