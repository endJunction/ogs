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
class BHE_Net_ELE_Pipe : public BHE_Net_ELE_Abstract
{
public:
    /**
     * constructor
     */
    BHE_Net_ELE_Pipe(std::string& name,
                     BHE_NET_ELE::type type = BHE_NET_ELE::BHE_NET_PIPE);

    double get_RHS_value();

    /*
    double set_BC(double T_in, double current_time)
    {
        return 0;
    }

    double get_flowrate()
    {
        return 0;
    }
    */
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
