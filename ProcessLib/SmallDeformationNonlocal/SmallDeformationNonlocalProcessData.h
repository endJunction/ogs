/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
template <int DisplacementDim>
struct SmallDeformationNonlocalProcessData
{
    SmallDeformationNonlocalProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            material,
        double const internal_length_)
        : material{std::move(material)},
          internal_length_squared(internal_length_ * internal_length_)
    {
    }

    SmallDeformationNonlocalProcessData(
        SmallDeformationNonlocalProcessData&& other)
        : material{std::move(other.material)},
          dt{other.dt},
          t{other.t},
          internal_length_squared{other.internal_length_squared}
    {
    }

    //! Copies are forbidden.
    SmallDeformationNonlocalProcessData(
        SmallDeformationNonlocalProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationNonlocalProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationNonlocalProcessData&&) = delete;

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    double dt = 0;
    double t = 0;
    double const internal_length_squared;

    double injected_volume = 0.0;
    double crack_volume = 0.0;
    bool propagating_crack = true;

    double pressure = 0.0;
    double pressure_old = 0.0;
    double pressure_error = 0.0;
};

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
