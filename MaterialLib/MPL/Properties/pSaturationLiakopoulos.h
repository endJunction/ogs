/**
 * \author Norbert Grunwald
 * \date   Oct 18, 2018
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

class SaturationLiakopoulos final : public Property
{
private:
    /// A pointer to the phase object.
    Medium* _medium;
    const double s_L_r = 0.2;
    const double s_a = 1.9722e-11;
    const double s_b = 2.4279;

public:
    /// Constructor passing a pointer to the medium.
    SaturationLiakopoulos(Medium*);
    /// Constructor passing a pointer to a phase.
    SaturationLiakopoulos(Phase*);
    /// Constructor passing a pointer to a component.
    SaturationLiakopoulos(Component*);

    /// This method overrides the base class implementation and
    /// actually computes and sets the property _value.
    PropertyDataType value(VariableArray const&) override;
    PropertyDataType dvalue(VariableArray const&, Variables const) override;
    PropertyDataType ddvalue(VariableArray const&, Variables const,
                             Variables const) override;
};

}  // namespace MaterialPropertyLib
