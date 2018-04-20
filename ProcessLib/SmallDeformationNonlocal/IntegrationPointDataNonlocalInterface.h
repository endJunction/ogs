/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
struct IntegrationPointDataNonlocalInterface;

struct IntegrationPointDataNonlocalInterface
{
    virtual ~IntegrationPointDataNonlocalInterface() = default;

    std::vector<IntegrationPointDataNonlocalInterface*> ip_l_pointer;
    std::vector<double> alpha_kl_times_w_l;

    double kappa_d = 0;      ///< damage driving variable.
    double integration_weight;
    /* TODO_MATERIAL_FORCES
    double nonlocal_internal_length;
    typename ShapeMatricesType::GlobalDimVectorType material_force;
    */
    bool active_self = false;
    bool activated = false;
};
}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
