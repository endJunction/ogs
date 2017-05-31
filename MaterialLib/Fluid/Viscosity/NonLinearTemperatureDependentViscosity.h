/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <string>
#include "BaseLib/ConfigTree.h"
#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/// Non-linear temperature dependent viscosity model.
class NonLinearTemperatureDependentViscosity final : public FluidProperty
{
public:
    /// Get model name.
    std::string getName() const override
    {
        return "Non-linear temperature dependent viscosity";
    }

    /// Get viscosity value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariableType.
    double getValue(const ArrayType& var_vals) const override
    {
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];

        return 1 + b * T + c * T * T + d * T * T * T + e * T * T * T * T +
               f * T * T * T * T * T + g * T * T * T * T * T * T;
    }

    /// Get the partial temperature derivative of viscosity; Other derivatives
    /// are 0;
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariableType.
    /// \param var       Variable type.
    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var) const override
    {
        if (var != PropertyVariableType::T)
            return 0;

        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];

        return b + c * T + d * T * T + e * T * T * T + f * T * T * T * T +
               g * T * T * T * T * T;
    }

private:
    static constexpr double a = 0.00175879595266029;
    static constexpr double b = -0.0000516976926689464 / a;
    static constexpr double c = 8.60121220346358E-07 / a;
    static constexpr double d = -8.17972869539407E-09 / a;
    static constexpr double e = 4.33281560663312E-11 / a;
    static constexpr double f = -1.18247503562765E-13 / a;
    static constexpr double g = 1.29200607741534E-16 / a;
};

}  // namespace Fluid
}  // namespace MaterialLib
