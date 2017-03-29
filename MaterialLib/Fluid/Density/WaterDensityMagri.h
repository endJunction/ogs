/**
 * \file
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <string>

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/// Nonlinear temperature dependent density model for water as described in the
/// FEFLOW White Papers Vol. III by Fabien Magri (2009).
class WaterDensityMagri final : public FluidProperty
{
public:

    /// Get model name.
    std::string getName() const override
    {
        return "Magri nonlinear temperature dependent density";
    }

    /// Get density value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariableType.
    double getValue(const ArrayType& var_vals) const override
    {
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        const double p = var_vals[static_cast<int>(PropertyVariableType::p)];
        return _a + _b * p + _c * p * p + (_d + _e * p + _f * p * p) * T +
               (_g + _h * p + _i * p * p) * T * T +
               (_j + _k * p + _l * p * p) * T * T * T +
               (_m + _n * p + _o * p * p) * T * T * T * T +
               (_p + _q * p + _r * p * p) * T * T * T * T * T +
               (_s + _t * p + _u * p * p) * T * T * T * T * T * T;
    }

    /// Get the partial differential of the density with respect to temperature.
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                  is given in enum class PropertyVariableType.
    /// \param var       Variable type.
    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var) const override
    {
        if (var != PropertyVariableType::T)
            return 0.0;

        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        const double p = var_vals[static_cast<int>(PropertyVariableType::p)];
        if (var == PropertyVariableType::T)
            return (_d + _e * p + _f * p * p) + (_g + _h * p + _i * p * p) * T +
                   (_j + _k * p + _l * p * p) * T * T +
                   (_m + _n * p + _o * p * p) * T * T * T +
                   (_p + _q * p + _r * p * p) * T * T * T * T +
                   (_s + _t * p + _u * p * p) * T * T * T * T * T;

        return 0;
    }

private:
    const double _a = 9.99792877961606e+02;
    const double _b = 5.07605113140940e-04;
    const double _c = -5.28425478164183e-10;
    const double _d = 5.13864847162196e-02;
    const double _e = -3.61991396354483e-06;
    const double _f = 7.97204102509724e-12;
    const double _g = -7.53557031774437e-03;
    const double _h = 6.32712093275576e-08;
    const double _i = -1.66203631393248e-13;
    const double _j = 4.60380647957350e-05;
    const double _k = -5.61299059722121e-10;
    const double _l = 1.80924436489400e-15;
    const double _m = -2.26651454175013e-07;
    const double _n = 3.36874416675978e-12;
    const double _o = -1.30352149261326e-17;
    const double _p = 6.14889851856743e-10;
    const double _q = -1.06165223196756e-14;
    const double _r = 4.75014903737416e-20;
    const double _s = -7.39221950969522e-13;
    const double _t = 1.42790422913922e-17;
    const double _u = -7.13130230531541e-23;
};

}  // end namespace
}  // end namespace
