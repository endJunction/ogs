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
/**
 * \class AverageVolumeFraction
 * \brief A function averaging a property by volume fraction
 * \details This property is usually a medium property, it
 * computes the average of individual phase properties
 * weighted by volume fraction.
 */
class BrooksCoreySaturation final : public Property
{
private:
    /// A pointer to the phase object.
    Medium* _medium;


public:
    /// Constructor passing a pointer to the medium.
    BrooksCoreySaturation(Medium*);
    /// Constructor passing a pointer to a phase.
    BrooksCoreySaturation(Phase*);
    /// Constructor passing a pointer to a component.
    BrooksCoreySaturation(Component*);

    /// This method overrides the base class implementation and
    /// actually computes and sets the property _value.
    PropertyDataType value(VariableArray const&) override;
    PropertyDataType dvalue(VariableArray const&, Variables const ) override;
};


}  // MaterialPropertyLib
