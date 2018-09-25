/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include "ProcessLib/Process.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "BHE_2U.h"
namespace ProcessLib
{
	namespace HeatTransportBHE
	{
		namespace BHE  // namespace of borehole heat exchanger
		{
			BHE::BHE_2U *
				CreateBHE2U(BaseLib::ConfigTree const& config,
					BaseLib::ConfigTree const& bhe_conf,
					std::map<std::string,
					std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const& curves,
                    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const&
                        bhe_refrigerant_density,
                    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const&
                        bhe_refrigerant_viscosity,
                    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const&
                        bhe_refrigerant_heat_capacity,
                    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const&
                        bhe_regrigerant_heat_conductivity);

		}  // end of namespace
	}
}
