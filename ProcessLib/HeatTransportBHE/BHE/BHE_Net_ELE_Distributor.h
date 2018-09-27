/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BHE_Net_ELE_Abstract.h"
namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE  // namespace of borehole heat exchanger
{
class BHE_Net_ELE_Distributor : public BHE_Net_ELE_Abstract
{
public:
    BHE_Net_ELE_Distributor(std::string& name, Eigen::VectorXd& vec_Inlet_Ratio,
                            Eigen::VectorXd& vec_Outlet_Ratio);

    double set_BC(double /*T_in*/, double /*current_time*/) { return 0; }

    double get_flowrate() { return 0; }

private:
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
