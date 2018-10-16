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

template <typename ShapeMatrixType, typename BHEType>
struct IntegrationPointDataBHE final
{
    explicit IntegrationPointDataBHE(BHEType const& bhe) : _bhe(bhe)
    {
        // depending on the type of BHE
        const int unknown_size = _bhe.getNumUnknowns();
        // initialization
        for (int i = 0; i < unknown_size; i++)
        {
            Eigen::VectorXd vec_advection(3);
            vec_advection.setZero();
            _vec_Advection_vectors.push_back(vec_advection);
        }

        // set parameter values
        for (int j = 0; j < unknown_size; j++)
        {
            // advection vector
            _bhe.getAdvectionVector(j, _vec_Advection_vectors[j]);
        }
    }

    BHEAbstract const& _bhe;

    typename ShapeMatrixType::NodalRowVectorType N;
    typename ShapeMatrixType::GlobalDimNodalMatrixType dNdx;
    double integration_weight;

    // Advection vectors
    std::vector<Eigen::VectorXd> _vec_Advection_vectors;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
