/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHEAbstract.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
template <typename ShapeMatricesType, typename BMatricesType, int GlobalDim>
struct IntegrationPointDataSoilNearBHE final
{
    explicit IntegrationPointDataSoilNearBHE(BHEAbstract& bhe_instance)
        : _bhe_instance(bhe_instance)
    {
        // depending on the type of BHE
        const int unknown_size = _bhe_instance.get_n_unknowns();
        /*
        // initialization
        _vec_mass_coefficients.resize(unknown_size);
        for (int i = 0; i < unknown_size; i++)
        {
            Eigen::MatrixXd mat_laplace(3, 3);
            mat_laplace.setZero();
            _vec_mat_Laplace.push_back(mat_laplace);
            Eigen::VectorXd vec_advection(3);
            vec_advection.setZero();
            _vec_Advection_vectors.push_back(vec_advection);
        }

        // set parameter values
        for (int j = 0; j < unknown_size; j++)
        {
            // mass matrix coefficients
            _vec_mass_coefficients[j] = _bhe_instance.get_mass_coeff(j);
            // laplace matrix
            _bhe_instance.get_laplace_matrix(j, _vec_mat_Laplace[j]);
            // advection vector
            _bhe_instance.get_advection_vector(
                j, _vec_Advection_vectors[j]);
        }
        */
    }

    double integration_weight;

    typename ShapeMatricesType::NodalRowVectorType N;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx;

    void pushBackState() {}

    BHEAbstract& _bhe_instance;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace HeatTransportBHE
}  // namespace ProcessLib
