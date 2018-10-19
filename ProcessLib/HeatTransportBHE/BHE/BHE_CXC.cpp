/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHE_CXC.h"

using namespace ProcessLib::HeatTransportBHE::BHE;

constexpr std::pair<int, int> BHE_CXC::inflow_outflow_bc_component_ids[];

/**
 * calculate thermal resistance
 */
void ProcessLib::HeatTransportBHE::BHE::BHE_CXC::initialize()
{
    flow_properties_in = calculateThermoMechanicalFlowPropertiesPipe(
        pipe_param, borehole_geometry.length, refrigerant_param, Q_r);

    flow_properties_out = calculateThermoMechanicalFlowPropertiesAnnulus(
        pipe_param, borehole_geometry.length, refrigerant_param, Q_r);

    calcThermalResistances();
}

void BHE_CXC::calcThermalResistances()
{
    double const& Nu_in = flow_properties_in.nusselt_number;
    double const& Nu_out = flow_properties_out.nusselt_number;
    double const& lambda_r = refrigerant_param.lambda_r;

    double const _R_ff = calculateThermalResistanceFf(Nu_in, Nu_out, lambda_r);
    double const _R_fog = calculateThermalResistanceFog(Nu_out, lambda_r);
    double const _R_gs = calculateThermalResistanceGroutSoil(Nu_out, lambda_r);

    calcHeatTransferCoefficients(_R_fog, _R_ff, _R_gs);
}

double BHE_CXC::calculateThermalResistanceGrout(double const Nu_out,
                                                double const lambda_r) const
{
    if (extern_Ra_Rb.use_extern_Ra_Rb)
    {
        // thermal resistance due to advective flow of refrigerant in the pipes
        // Eq. 58, 59, and 60 in Diersch_2011_CG
        double const R_adv_b_o1 = thermalResistanceMagicalOuverture(
            Nu_out, lambda_r,
            pipe_param.r_outer - (pipe_param.r_inner + pipe_param.b_in),
            pipe_param.r_outer);
        // thermal resistance due to thermal conductivity of the pip wall
        // material Eq. 66 in Diersch_2011_CG
        double const R_con_o1 = thermalResistanceMagicalMur(
            pipe_param.r_outer + pipe_param.b_out, pipe_param.r_outer,
            pipe_param.lambda_p_o);

        return extern_Ra_Rb.ext_Rb - R_adv_b_o1 - R_con_o1;
    }
    // Eq. 69
    return thermalResistanceMagicalMur(borehole_geometry.diameter / 2,
                                       pipe_param.r_outer + pipe_param.b_out,
                                       grout_param.lambda_g);
}

// Eq. 67
double BHE_CXC::calculateThermalResistanceFf(double const Nu_in,
                                             double const Nu_out,
                                             double const lambda_r) const
{
    if (extern_Ra_Rb.use_extern_Ra_Rb)
    {
        return extern_Ra_Rb.ext_Ra;
    }
    if (extern_def_thermal_resistances.if_use_defined_therm_resis)
    {
        return extern_def_thermal_resistances
            .ext_Rgg1;  // Attention! Here ext_Rgg1 is treated as Rff for
                        // coaxial type
    }
    // thermal resistance due to advective flow of refrigerant in the pipes
    // Eq. 58, 59, and 60 in Diersch_2011_CG
    double const R_adv_i1 =
        thermalResistanceMagicalIntroduction(Nu_in, lambda_r);
    double const R_adv_a_o1 = thermalResistanceMagicalOuverture(
        Nu_out, lambda_r,
        pipe_param.r_outer - (pipe_param.r_inner + pipe_param.b_in),
        pipe_param.r_inner + pipe_param.b_in);
    // thermal resistance due to thermal conductivity of the pip wall material
    // Eq. 66 in Diersch_2011_CG
    double const R_con_i1 =
        thermalResistanceMagicalMur(pipe_param.r_inner + pipe_param.b_in,
                                    pipe_param.r_inner, pipe_param.lambda_p_i);
    // Eq. 56
    return R_adv_i1 + R_adv_a_o1 + R_con_i1;
}

// Eq. 57
double BHE_CXC::calculateThermalResistanceFog(double const Nu_out,
                                              double const lambda_r) const
{
    if (extern_def_thermal_resistances.if_use_defined_therm_resis)
        return extern_def_thermal_resistances.ext_Rfog;

    // thermal resistance due to advective flow of refrigerant in the pipes
    // Eq. 58, 59, and 60 in Diersch_2011_CG
    double const R_adv_b_o1 = thermalResistanceMagicalOuverture(
        Nu_out, lambda_r,
        pipe_param.r_outer - (pipe_param.r_inner + pipe_param.b_in),
        pipe_param.r_outer);
    // thermal resistance due to thermal conductivity of the pip wall
    // material Eq. 66 in Diersch_2011_CG
    double const R_con_o1 =
        thermalResistanceMagicalMur(pipe_param.r_outer + pipe_param.b_out,
                                    pipe_param.r_outer, pipe_param.lambda_p_o);
    // thermal resistance due to the grout transition
    // Eq. 68
    double const chi = chimereDimensionlessFactor(
        borehole_geometry.diameter,
        2.0 * (pipe_param.r_outer + pipe_param.b_out));

    // thermal resistances of the grout
    double const R_g = calculateThermalResistanceGrout(Nu_out, lambda_r);
    double const R_con_b = chi * R_g;
    return R_adv_b_o1 + R_con_o1 + R_con_b;
}

// thermal resistance due to grout-soil exchange
double BHE_CXC::calculateThermalResistanceGroutSoil(double const Nu_out,
                                                    double const lambda_r) const
{
    if (extern_def_thermal_resistances.if_use_defined_therm_resis)
    {
        return extern_def_thermal_resistances.ext_Rgs;
    }
    // thermal resistance due to the grout transition
    // Eq. 68
    double const chi = chimereDimensionlessFactor(
        borehole_geometry.diameter,
        2.0 * (pipe_param.r_outer + pipe_param.b_out));

    // thermal resistances of the grout
    double const R_g = calculateThermalResistanceGrout(Nu_out, lambda_r);
    return (1 - chi) * R_g;
}

/**
 * calculate heat transfer coefficient
 */
void BHE_CXC::calcHeatTransferCoefficients(double const R_fog,
                                           double const R_ff, double const R_gs)
{
    boundary_heat_exchange_coefficients[0] = 1.0 / R_fog;
    boundary_heat_exchange_coefficients[1] = 1.0 / R_ff;

    if (!std::isfinite(R_gs))
    {
        OGS_FATAL(
            "Error!!! Grout Thermal Resistance is an infinite number! The "
            "simulation will be stopped! \n");
    }
    boundary_heat_exchange_coefficients[2] = 1.0 / R_gs;
}

double BHE_CXC::getTinByTout(double const T_out, double const current_time)
{
    double const& rho_r = refrigerant_param.rho_r;
    double const& heat_cap_r = refrigerant_param.heat_cap_r;
    double T_in(0.0);
    double power_tmp(0.0);
    double Q_r_tmp(0.0);

    switch (this->boundary_type)
    {
        case BHE_BOUNDARY_TYPE::POWER_IN_WATT_BOUNDARY:
            T_in = power_in_watt_val / Q_r / heat_cap_r / rho_r + T_out;
            break;
        case BHE_BOUNDARY_TYPE::FIXED_TEMP_DIFF_BOUNDARY:
            T_in = T_out + delta_T_val;
            break;
        case BHE_BOUNDARY_TYPE::POWER_IN_WATT_CURVE_FIXED_DT_BOUNDARY:
            // get the power value in the curve
            // power_tmp = GetCurveValue(power_in_watt_curve_idx, 0,
            // current_time, &flag_valid);
            power_tmp = power_in_watt_curve->getValue(current_time);

            // if power value exceeds threshold, calculate new values
            if (std::fabs(power_tmp) > threshold)
            {
                // calculate the corresponding flow rate needed
                // using the defined delta_T value
                Q_r_tmp = power_tmp / delta_T_val / heat_cap_r / rho_r;
                // update all values dependent on the flow rate
                Q_r = Q_r_tmp;
                initialize();
                // calculate the new T_in
                T_in = T_out + delta_T_val;
            }
            else
            {
                Q_r_tmp = 1.0e-06;  // this has to be a small value to avoid
                                    // division by zero
                // update all values dependent on the flow rate
                Q_r = Q_r_tmp;
                initialize();
                // calculate the new T_in
                T_in = T_out;
            }
            break;
        case BHE_BOUNDARY_TYPE::POWER_IN_WATT_CURVE_FIXED_FLOW_RATE_BOUNDARY:
            // get the power value in the curve
            // power_tmp = GetCurveValue(power_in_watt_curve_idx, 0,
            // current_time, &flag_valid);
            power_tmp = power_in_watt_curve->getValue(current_time);

            // calculate the dT value based on fixed flow rate
            delta_T_val = power_tmp / Q_r / heat_cap_r / rho_r;
            // calcuate the new T_in
            T_in = T_out + delta_T_val;
            break;
        default:
            T_in = T_out;
            break;
    }

    return T_in;
}
