/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct FlowAndTemperature
{
    double const flow_rate;
    double const temperature;
};

template <typename ControlBase>
struct FlowAndTemperatureControl : ControlBase
{
    template <typename... Args>
    FlowAndTemperatureControl(Args&&... args)
        : _control_base{std::forward<Args>(args)...}
    {
    }

    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        return _control_base.flowAndTemperature(T_out, time);
    }

    ControlBase _control_base;
};

struct TemperatureCurveConstantFlow
{
    FlowAndTemperature flowAndTemperature(double const /*T_out*/,
                                          double const time) const
    {
        return {flow_rate, temperature_curve->getValue(time)};
    }
    double flow_rate;
    MathLib::PiecewiseLinearInterpolation const* temperature_curve;
};

struct FixedPowerConstantFlow
{
    FlowAndTemperature flowAndTemperature(double const T_out,
                                          double const /*time*/) const
    {
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    double flow_rate;
    double power;  // Value is expected to be in Watt.
    double heat_capacity;
    double density;
};

struct FixedPowerFlowCurve
{
    FlowAndTemperature flowAndTemperature(double const T_out,
                                          double const time) const
    {
        double const flow_rate = flow_curve->getValue(time);
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    MathLib::PiecewiseLinearInterpolation const* flow_curve;

    double power;   // Value is expected to be in Watt.
    double heat_capacity;
    double density;
};

// TODO (haibing) Providing test cases or unit tests is required to enable the
// following sections
#if 0
inline double computeCopCurvePowerCoolingHeating(
    double const building_power, double const T_out,
    MathLib::PiecewiseLinearInterpolation* const heating_cop_curve,
    MathLib::PiecewiseLinearInterpolation* const cooling_cop_curve)
{
    if (building_power <= 0.0)
    {
        // get COP value based on T_out in the curve
        double const COP = heating_cop_curve->getValue(T_out);

        // now calculate how much power needed from BHE
        return building_power * (COP - 1.0) / COP;
        // also how much power from electricity
        // power_elect_tmp = building_power_tmp - power_tmp;
        // print the amount of power needed
        // std::cout << "COP: " << COP << ", Q_bhe: " << power_tmp
        // << ", Q_elect: " << power_elect_tmp << std::endl;
    }
    else
    {
        // get COP value based on T_out in the curve
        double const COP = cooling_cop_curve->getValue(T_out);

        // now calculate how much power needed from BHE
        return building_power * (COP + 1.0) / COP;
        // also how much power from electricity
        // power_elect_tmp = -building_power_tmp + power_tmp;
        // print the amount of power needed
        // std::cout << "COP: " << COP << ", Q_bhe: " << power_tmp
        // << ", Q_elect: " << power_elect_tmp << std::endl;
    }
}

struct FixedTemperatureDifferenceConstantFlow
{
    FlowAndTemperature operator()(double const T_out,
                                  double const /*time*/) const override
    {
        return {flow_rate, T_out + temperature_difference};
    }
    double flow_rate;
    double temperature_difference;
};

struct PowerCurveWithFixedTemperatureDifference
{
    FlowAndTemperature operator()(double const T_out,
                                  double const time) const override
    {
        double const power = power_curve->getValue(time);
        return powerBasedFlowAndTemperature(power);
    }
    MathLib::PiecewiseLinearInterpolation * const power_curve;
    double const power_threshold;
    double const heat_capacity;
    double const density;
    double const minimum_flow_rate;
    double const temperature_difference;

protected:
    FlowAndTemperature powerBasedFlowAndTemperature(double const power) const
    {

        if (std::fabs(power) < power_threshold)
        {
            return {minimum_flow_rate, T_out};
        }
        double const temperature_increment =
            (power < 0 ? -1 : 1) * temperature_difference;
        // calculate the corresponding flow rate needed using the defined
        // delta_T value
        double const flow_rate =
            power / temperature_increment / heat_capacity / density;
        // calculate the new T_in
        return {flow_rate, T_out + temperature_increment};
    }
};

/// Coefficient of performance scaled power curve with fixed temperature
/// difference. The scaling is heating/cooling dependent and is given as a
/// curve in both cases.
struct CopScaledPowerCurveWithFixedTemperatureDifference
    : PowerCurveWithFixedTemperatureDifference
{
    FlowAndTemperature operator()(double const T_out,
                                  double const time) const override
    {
        double const power = power_curve->getValue(time);
        double const scaled_power = computeCopCurvePowerCoolingHeating(
            power, T_out, heating_cop_curve, cooling_cop_curve);

        return powerBasedFlowAndTemperature(scaled_power);
    }
    MathLib::PiecewiseLinearInterpolation* const heating_cop_curve;
    MathLib::PiecewiseLinearInterpolation* const cooling_cop_curve;
};

struct PowerCurveWithFixedFlowRate
{
    FlowAndTemperature operator()(double const T_out,
                                  double const time) const override
    {
        double const power = power_curve->getValue(current_time);
        return powerBasedFlowAndTemperature(power);
    }
    MathLib::PiecewiseLinearInterpolation * const power_curve;
    double const power_threshold;
    double const heat_capacity;
    double const density;
    double const minimum_flow_rate;
    double const flow_rate;

protected:
    FlowAndTemperature powerBasedFlowAndTemperature(double const power,
                                                    double const T_out) const
    {
        if (std::fabs(power) < power_threshold)
        {
            return {minimum_flow_rate, T_out};
        }
        // calculate the dT value based on fixed flow rate
        double const dT = power / flow_rate / heat_capacity / density;
        // calcuate the new T_in
        return {flow_rate, T_out + dT};
    }
};

struct CopScaledPowerCurveWithFixedFlowRate
{
    FlowAndTemperature operator()(double const T_out,
                                  double const time) const override
    {
        double const power = power_curve->getValue(time);
        double const scaled_power = computeCopCurvePowerCoolingHeating(
            power, T_out, heating_cop_curve, cooling_cop_curve);
        return powerBasedFlowAndTemperature(scaled_power, T_out);
    }
    MathLib::PiecewiseLinearInterpolation* const heating_cop_curve;
    MathLib::PiecewiseLinearInterpolation* const cooling_cop_curve;
};
#endif
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
