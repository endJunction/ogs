/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHE_Net_ELE_Pipe_Inner_2U.h"

using namespace ProcessLib::HeatTransportBHE::BHE;

BHE_Net_ELE_Pipe_Inner_2U::BHE_Net_ELE_Pipe_Inner_2U(
    std::string& name, ProcessLib::HeatTransportBHE::BHE::BHEAbstract* m_BHE)
    : BHE_Net_ELE_Pipe(name, BHE_NET_ELE::BHE_NET_PIPE_INNER_2U)
{
    // configure the penalty factor
    this->set_penalty_factor(1.0e6);
}

std::size_t ProcessLib::HeatTransportBHE::BHE::BHE_Net_ELE_Pipe_Inner_2U::
    get_global_idx_in()
{
    return _global_idx_in;
}

std::size_t ProcessLib::HeatTransportBHE::BHE::BHE_Net_ELE_Pipe_Inner_2U::
    get_global_idx_out()
{
    return _global_idx_out;
}
