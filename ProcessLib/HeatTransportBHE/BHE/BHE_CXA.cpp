/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHE_CXA.h"
#include "Physics.h"

using namespace ProcessLib::HeatTransportBHE::BHE;

constexpr std::pair<int, int> BHE_CXA::inflow_outflow_bc_component_ids[];

namespace {
std::pair<double, double> calcReynoldsNumber(
    double const u_out, double const u_in, PipeParameters const& pipe_param,
    RefrigerantParameters const& refrigerant_param)
{
    double const& r_outer = pipe_param.r_outer;
    double const& r_inner = pipe_param.r_inner;
    double const& b_in = pipe_param.b_in;
    double const& mu_r = refrigerant_param.mu_r;
    double const& rho_r = refrigerant_param.rho_r;

    double const d_o = 2.0 * r_inner;
    double const d_i = 2.0 * (r_outer - (r_inner + b_in));

    double const Re_o1 = reynoldsNumber(u_out, d_o, mu_r, rho_r);
    double const Re_i1 = reynoldsNumber(u_in, d_i, mu_r, rho_r);
    return {Re_o1, Re_i1};
}
}  // namespace

/**
 * calculate thermal resistance
 */
void ProcessLib::HeatTransportBHE::BHE::BHE_CXA::initialize()
{
    double const u_in = annulusFlowVelocity(
        Q_r, pipe_param.r_outer, pipe_param.r_inner + pipe_param.b_in);
    double const u_out = pipeFlowVelocity(Q_r, pipe_param.r_inner);
    _u(0) = u_in;
    _u(1) = u_out;

    auto const Re =
        calcReynoldsNumber(u_out, u_in, pipe_param, refrigerant_param);
    double const Pr = prandtlNumber(refrigerant_param.mu_r,
                       refrigerant_param.heat_cap_r,
                       refrigerant_param.lambda_r);

    double const Nu_out = nusseltNumber(Re.first, Pr, 2 * pipe_param.r_inner,
                                        borehole_geometry.length);
    _Nu(1) = Nu_out;
    double const diameter_ratio =
        (pipe_param.r_inner + pipe_param.b_in) / pipe_param.r_outer;
    double const pipe_aspect_ratio =
        2.0 * (pipe_param.r_outer - (pipe_param.r_inner + pipe_param.b_in)) /
        borehole_geometry.length;
    double const Nu_in =
        nusseltNumberAnnulus(Re.second, Pr, diameter_ratio, pipe_aspect_ratio);
    _Nu(0) = Nu_in;

    calcThermalResistances();
    calcHeatTransferCoefficients();
}

void BHE_CXA::calcThermalResistances()
{
    double Nu_in, Nu_out;
    double d_o1, d_i1, d_h;
    double chi;
    double _R_con_i1, _R_con_o1;
    double const& D = borehole_geometry.diameter;
    double const& r_outer = pipe_param.r_outer;
    double const& r_inner = pipe_param.r_inner;
    double const& b_in = pipe_param.b_in;
    double const& b_out = pipe_param.b_out;
    double const& lambda_r = refrigerant_param.lambda_r;
    double const& lambda_g = grout_param.lambda_g;
    double const& lambda_p_i = pipe_param.lambda_p_i;
    double const& lambda_p_o = pipe_param.lambda_p_o;

    Nu_in = _Nu(0);
    Nu_out = _Nu(1);
    d_i1 = 2.0 * r_outer;
    d_o1 = 2.0 * (r_inner + b_in);
    d_h = 2.0 * (r_outer - (r_inner + b_in));

    constexpr double PI = boost::math::constants::pi<double>();
    // thermal resistance due to advective flow of refrigerant in the pipes
    // Eq. 58, 59, and 60 in Diersch_2011_CG
    _R_adv_o1 = 1.0 / (Nu_out * lambda_r * PI);
    _R_adv_a_i1 = 1.0 / (Nu_in * lambda_r * PI) * (d_h / d_o1);
    _R_adv_b_i1 = 1.0 / (Nu_in * lambda_r * PI) * (d_h / d_i1);

    // thermal resistance due to thermal conductivity of the pip wall material
    // Eq. 66 in Diersch_2011_CG
    _R_con_i1 = std::log((r_outer + b_out) / r_outer) / (2.0 * PI * lambda_p_o);
    _R_con_o1 = std::log((r_inner + b_out) / r_inner) / (2.0 * PI * lambda_p_i);

    // thermal resistance due to the grout transition
    d_i1 = 2.0 * (r_outer + b_out);
    // Eq. 68
    chi = std::log(std::sqrt(D * D + d_i1 * d_i1) / std::sqrt(2) / d_i1) /
          std::log(D / d_i1);
    if (extern_Ra_Rb.use_extern_Ra_Rb)
    {
        _R_g = extern_Ra_Rb.ext_Rb - _R_adv_b_i1 - _R_con_i1;
    }
    else
    {
        // Eq. 69
        _R_g = std::log(D / d_i1) / (2.0 * PI * lambda_g);
    }
    // Eq. 67
    _R_con_b = chi * _R_g;
    if (extern_Ra_Rb.use_extern_Ra_Rb)
    {
        _R_ff = extern_Ra_Rb.ext_Ra;
    }
    else if (extern_def_thermal_resistances.if_use_defined_therm_resis)
    {
        _R_ff = extern_def_thermal_resistances
                    .ext_Rgg1;  // Attention! Here ext_Rgg1 is treated as Rff
                                // for coaxial type
    }
    else
    {
        // Eq. 56
        _R_ff = _R_adv_o1 + _R_adv_a_i1 + _R_con_o1;
    }

    // Eq. 57
    if (extern_def_thermal_resistances.if_use_defined_therm_resis)
        _R_fig = extern_def_thermal_resistances.ext_Rfig;
    else
        _R_fig = _R_adv_b_i1 + _R_con_i1 + _R_con_b;

    // thermal resistance due to grout-soil exchange
    if (extern_def_thermal_resistances.if_use_defined_therm_resis)
        _R_gs = extern_def_thermal_resistances.ext_Rgs;
    else
        _R_gs = (1 - chi) * _R_g;

    if (!std::isfinite(_R_gs))
    {
        OGS_FATAL(
            "Error!!! Grout Thermal Resistance is an infinite number! The "
            "simulation will be stopped! \n");
    }
}

/**
 * calculate heat transfer coefficient
 */
void BHE_CXA::calcHeatTransferCoefficients()
{
    _PHI_fig = 1.0 / _R_fig;
    _PHI_ff = 1.0 / _R_ff;
    _PHI_gs = 1.0 / _R_gs;
}

double BHE_CXA::getBoundaryHeatExchangeCoeff(std::size_t idx_unknown) const
{
    // Here we calculates the boundary heat exchange coefficients
    // in the governing equations of BHE.
    // These governing equations can be found in
    // 1) Diersch (2013) FEFLOW book on page 958, M.3, or
    // 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 90-97.

    double exchange_coeff(0);
    double const PHI = _PHI_fig;
    switch (idx_unknown)
    {
        case 0:
            exchange_coeff = PHI;
            break;
        case 1:
            // PHI_ff
            exchange_coeff = _PHI_ff;
            break;
        case 2:
            // PHI_gs
            exchange_coeff = _PHI_gs;
            break;
        default:
            OGS_FATAL(
                "Error !!! The index passed to "
                "getBoundaryHeatExchangeCoeff for BHE is not correct. ");
            break;
    }
    return exchange_coeff;
}

double BHE_CXA::getTinByTout(double T_out, double current_time = -1.0)
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
            if (fabs(power_tmp) > threshold)
            {
                // calculate the corresponding flow rate needed
                // using the defined delta_T value
                Q_r_tmp = power_tmp / delta_T_val / heat_cap_r / rho_r;
                // update all values dependent on the flow rate
                updateFlowRate(Q_r_tmp);
                // calculate the new T_in
                T_in = T_out + delta_T_val;
            }
            else
            {
                Q_r_tmp = 1.0e-06;  // this has to be a small value to avoid
                                    // division by zero
                // update all values dependent on the flow rate
                updateFlowRate(Q_r_tmp);
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
