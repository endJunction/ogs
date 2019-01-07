/**
 * \author Norbert Grunwald
 * \date   Oct 19, 2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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

class RelPermLiakopoulos final : public Property
{
private:
    /// A pointer to the phase object.
    Medium* _medium;

    const double s_L_r = 0.2;
    const double lambda = 3.;
    const double kA = 2.207;
    const double kB = 1.0121;
    const double k_rel_G_min = 1.0e-4;

public:
    /// Constructor passing a pointer to the medium.
    RelPermLiakopoulos(Medium*);
    /// Constructor passing a pointer to a phase.
    RelPermLiakopoulos(Phase*);
    /// Constructor passing a pointer to a component.
    RelPermLiakopoulos(Component*);

    /// This method overrides the base class implementation and
    /// actually computes and sets the property _value.
    PropertyDataType value(VariableArray const&) override;
    PropertyDataType dvalue(VariableArray const&, Variables const) override;
};

}  // namespace MaterialPropertyLib
