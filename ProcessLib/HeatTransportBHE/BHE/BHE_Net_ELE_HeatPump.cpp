/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHE_Net_ELE_HeatPump.h"
#include "ProcessLib/Utils/ProcessUtils.h"

using namespace ProcessLib::HeatTransportBHE::BHE;

BHE_Net_ELE_HeatPump::BHE_Net_ELE_HeatPump(std::string& name)
    : BHE_Net_ELE_Abstract(name, BHE_NET_ELE::BHE_NET_HEATPUMP, 1, 1)
{
    _heat_pump_BC_type = HEAT_PUMP_BOUND_POWER_FIXED_DT;
}

double BHE_Net_ELE_HeatPump::set_BC(double T_in, double /*current_time*/)
{
    double T_out = 0.0;
    double power_hp = 0.0;
    double power_bhe = 0.0;
    // int flag_valid = false;
    double COP = 0.0;
    double power_el = 0.0;
    double delta_T = 0.0;

    switch (_heat_pump_BC_type)
    {
        case HEAT_PUMP_BOUND_POWER_FIXED_DT:
            break;
        case HEAT_PUMP_BOUND_POWER_FIXED_FLOWRATE:
            double rho_cp_u = _fluid_density * _fluid_heat_capacity * _flowrate;

            power_hp = _power_val;
            COP = _cop_curve->getValue(T_in);
            power_bhe = power_hp * (COP - 1.0) / COP;
            // also how much power from electricity
            power_el = power_hp - power_bhe;
            if (fabs(power_hp) < 0.1)
            {
                T_out = T_in;
            }
            else
            {
                delta_T = -power_bhe / rho_cp_u;
                T_out = T_in - delta_T;
            }

            DBUG("Heat pump: \'%s\', T_in: %.2d, T_out: %.2d. \n",
                 get_ele_name().c_str(), T_in, T_out);
            DBUG("COP: %.2d, Q_bhe: %.2d, Q_electricity: %.2d. \n", COP,
                 power_bhe, power_el);
            break;
    }

    return T_out;
}

double BHE_Net_ELE_HeatPump::get_RHS_value()
{
    double rt_RHS_val = 0.0;

    // depending on the boundary condition,
    // calculate the RHS value
    switch (_heat_pump_BC_type)
    {
        case HEAT_PUMP_BOUND_POWER_FIXED_DT:
            rt_RHS_val = _delta_T_val;
            break;
        case HEAT_PUMP_BOUND_POWER_FIXED_FLOWRATE:
            // TODO
            break;
    }
    return rt_RHS_val;
}
