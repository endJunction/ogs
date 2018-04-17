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
    explicit IntegrationPointDataBHE(BHE::BHEAbstract& bhe_instance)
        : _bhe_instance(bhe_instance)
    {
    }

    // typename HMatricesType::HMatrixType _h_matrices;
    // typename HMatricesType::ForceVectorType _sigma, _sigma_prev;
    // typename HMatricesType::ForceVectorType _w, _w_prev;
    // double _aperture = 0.0;
    // double _aperture_prev = 0.0;
    // double _aperture0 = 0.0;

    BHE::BHEAbstract& _bhe_instance;
    /*
    std::unique_ptr<typename MaterialLib::Fracture::FractureModelBase<
        DisplacementDim>::MaterialStateVariables>
        _material_state_variables;
    */

    // Eigen::MatrixXd _C;
    double integration_weight;

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
