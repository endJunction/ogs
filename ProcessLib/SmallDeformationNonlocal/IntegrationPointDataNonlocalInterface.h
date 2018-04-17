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
struct IntegrationPointDataNonlocalInterface
{
    virtual ~IntegrationPointDataNonlocalInterface() = default;

    std::vector<std::tuple<
        // element's local assembler
        IntegrationPointDataNonlocalInterface* const,
        double,  // squared distance to current integration point
        double   // alpha_kl * w_l
        >>
        non_local_assemblers;

    double kappa_d = 0;      ///< damage driving variable.
    double integration_weight;
    double nonlocal_internal_length;
    /* TODO_MATERIAL_FORCES
    typename ShapeMatricesType::GlobalDimVectorType material_force;
    */
    Eigen::Vector3d coordinates;
    bool active_self = false;
    bool activated = false;
};
}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
