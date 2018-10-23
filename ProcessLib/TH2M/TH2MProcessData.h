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

#include <Eigen/Dense>

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}  // namespace MaterialLib
namespace MaterialPropertyLib
{
class Medium;
}
namespace ProcessLib
{
namespace TH2M
{
template <int DisplacementDim>
struct TH2MProcessData
{
    TH2MProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            material_,
        Eigen::Matrix<double, DisplacementDim, 1>
            specific_body_force_,
            std::map<int, std::unique_ptr<MaterialPropertyLib::Medium>> const& media_,
            MeshLib::PropertyVector<int> const* material_ids_)
        : material{std::move(material_)},
          specific_body_force(std::move(specific_body_force_)),
          media(media_),
          material_ids(material_ids_)
    {
    }

    /// The constitutive relation for the mechanical part.
    /// \note Linear elasticity is the only supported one in the moment.
    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    std::map<int, std::unique_ptr<MaterialPropertyLib::Medium>> const& media;

    MeshLib::PropertyVector<int> const* const material_ids;

    double dt = 0.0;
    double t = 0.0;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace TH2M
}  // namespace ProcessLib
