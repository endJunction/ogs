/**
* \copyright
* Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include <Eigen/Eigen>

#include "ProcessLib/HeatTransportBHE/BHE/BHEAbstract.h"

namespace ProcessLib
{
    namespace HeatTransportBHE
    {
        using namespace BHE;

        struct IntegrationPointDataBHE final
        {
            explicit IntegrationPointDataBHE(BHEAbstract& bhe_instance)
                : _bhe_instance(bhe_instance)
            {
                // depending on the type of BHE
                // switch (_bhe_instance.)
                // _BHE_element_prop.lambda_g = _bhe_instance.
            }

            // typename HMatricesType::HMatrixType _h_matrices;
            // typename HMatricesType::ForceVectorType _sigma, _sigma_prev;
            // typename HMatricesType::ForceVectorType _w, _w_prev;
            // double _aperture = 0.0;
            // double _aperture_prev = 0.0;
            // double _aperture0 = 0.0;

            BHEAbstract& _bhe_instance;

            /*
            std::unique_ptr<typename MaterialLib::Fracture::FractureModelBase<
                DisplacementDim>::MaterialStateVariables>
                _material_state_variables;
            */

            // Eigen::MatrixXd _C;
            double integration_weight;

            // product of refrigerant density and heat capacity
            double rho_r_cp_r;

            // hydrothermal dispersion cofficient of refrigerant
            double Lambda;

            // product of grout density and heat capacity
            double rho_g_cp_g;

            // grout thermla conductivity
            double lambda_g;

            // vector of refrigerant flow velocity in the downward pipe
            Eigen::Vector3d vec_flow_velocity_in_1;

            // vector of refrigerant flow velocity in the upward pipe
            Eigen::Vector3d vec_flow_velocity_out_1;

            // vector of refrigerant flow velocity in the downward pipe
            Eigen::Vector3d vec_flow_velocity_in_2;

            // vector of refrigerant flow velocity in the upward pipe
            Eigen::Vector3d vec_flow_velocity_out_2;

            void pushBackState()
            {
                // _w_prev = _w;
                // _sigma_prev = _sigma;
                // _aperture_prev = _aperture;
            }

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        };

    }  // namespace HeatTransportBHE
}  // namespace ProcessLib
