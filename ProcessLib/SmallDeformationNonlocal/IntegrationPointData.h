/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "IntegrationPointDataNonlocalInterface.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
template <typename BMatricesType, typename ShapeMatricesType,
          int DisplacementDim>
struct IntegrationPointData final : public IntegrationPointDataNonlocalInterface
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
        if (auto const msv = dynamic_cast<typename MaterialLib::Solids::Ehlers::
                                              StateVariables<DisplacementDim>*>(
                material_state_variables.get()))
        {
            eps_p_V = &msv->eps_p.V;
            eps_p_D_xx = &(msv->eps_p.D[0]);
        }

        /* TODO_MATERIAL_FORCES
        material_force.setZero(DisplacementDim);
        */
    }

#if defined(_MSC_VER) && _MSC_VER < 1900
    // The default generated move-ctor is correctly generated for other
    // compilers.
    explicit IntegrationPointData(IntegrationPointData&& other)
        : b_matrices(std::move(other.b_matrices)),
          sigma(std::move(other.sigma)),
          sigma_prev(std::move(other.sigma_prev)),
          eps(std::move(other.eps)),
          eps_prev(std::move(other.eps_prev)),
          solid_material(other.solid_material),
          material_state_variables(std::move(other.material_state_variables)),
          C(std::move(other.C)),
          integration_weight(std::move(other.integration_weight)),
    {
    }
#endif  // _MSC_VER

    typename BMatricesType::BMatrixType b_matrices;
    typename BMatricesType::KelvinVectorType sigma, sigma_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;
    double free_energy_density = 0;
    double damage = 0;       ///< isotropic damage
    double damage_prev = 0;  ///< \copydoc damage
    double kappa_d_prev = 0;  ///< \copydoc kappa_d
    // double nonlocal_kappa_d = 0;
    // double nonlocal_kappa_d_prev = 0;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    typename BMatricesType::KelvinMatrixType C;
    typename ShapeMatricesType::NodalRowVectorType N;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx;

    double const* eps_p_V;
    double const* eps_p_D_xx;

    void pushBackState()
    {
        eps_prev = eps;
        sigma_prev = sigma;
        damage_prev = damage;
        kappa_d_prev = kappa_d;
        // nonlocal_kappa_d_prev = nonlocal_kappa_d;
        material_state_variables->pushBackState();
    }

    // Unused double getLocalRateKappaD() const { return kappa_d - kappa_d_prev; }
};

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
