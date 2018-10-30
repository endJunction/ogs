/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "TH2MProcessData.h"

namespace ProcessLib
{
namespace TH2M
{
namespace MPL = MaterialPropertyLib;

template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatricesTypePressure, int DisplacementDim, int NPoints>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

    typename ShapeMatrixTypeDisplacement::template MatrixType<
        DisplacementDim, NPoints * DisplacementDim>
        N_u_op;
    typename BMatricesType::KelvinVectorType sigma_eff, sigma_eff_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;

    // typename BMatricesType::KelvinVectorType sigma_prev;

    typename ShapeMatrixTypeDisplacement::NodalRowVectorType N_u;
    typename ShapeMatrixTypeDisplacement::GlobalDimNodalMatrixType dNdx_u;

    typename ShapeMatricesTypePressure::NodalRowVectorType N_p;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType dNdx_p;

    using GlobalDimVectorType =
        typename ShapeMatricesTypePressure::GlobalDimVectorType;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;
    double integration_weight;

    // output vector for wetting phase saturation and pressure with
    // respect to each integration point
    double saturation;
    double pressure_wet;
    double rel_perm_gas;
    double rel_perm_liquid;
    double density_gas;
    double density_liquid;
    double pressure_gas_linear;
    double pressure_cap_linear;
    GlobalDimVectorType velocity_gas{DisplacementDim};
    GlobalDimVectorType velocity_liquid{DisplacementDim};
    void pushBackState()
    {
        eps_prev = eps;
        sigma_eff_prev = sigma_eff;
        material_state_variables->pushBackState();
    }

    template <typename DisplacementVectorType>
    typename BMatricesType::KelvinMatrixType updateConstitutiveRelation(
        double const t,
        SpatialPosition const& x_position,
        double const dt,
        DisplacementVectorType const& /*u*/,
        double const T,
        double const /*p*/)
    {
        auto&& solution = solid_material.integrateStress(
            t, x_position, dt, eps_prev, eps, sigma_eff_prev,
            *material_state_variables, T);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma_eff, material_state_variables, C) = std::move(*solution);

        return C;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace TH2M
}  // namespace ProcessLib
