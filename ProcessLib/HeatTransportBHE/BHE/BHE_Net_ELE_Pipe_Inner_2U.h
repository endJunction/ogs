/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include "BHE_Net_ELE_Abstract.h"
#include "BHEAbstract.h"
#include "BHE_Net_ELE_Pipe.h"

namespace ProcessLib
{
	namespace HeatTransportBHE
	{
		namespace BHE  // namespace of borehole heat exchanger
		{
			class BHE_Net_ELE_Pipe_Inner_2U : public BHE_Net_ELE_Pipe {

			public:
				/**
				* constructor
				*/
				BHE_Net_ELE_Pipe_Inner_2U(std::string & name, BHE::BHEAbstract * m_BHE);

			protected:
				/**
				  * obtain the global index at the pipeline inlet
				  */
				std::size_t get_global_idx_in();

				/**
				  * obtain the global index at the pipeline outlet
				  */
				std::size_t get_global_idx_out();

			private:
				/**
				  * the global index at the pipeline inlet
				  */
				std::size_t _global_idx_in;

				/**
				  * the global index at the pipeline outlet
				  */
				std::size_t _global_idx_out;
			};

		}
	}
}
