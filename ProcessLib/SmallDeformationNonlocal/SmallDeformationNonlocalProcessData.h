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
        : material{std::move(material)}, internal_length(internal_length_)
    {
    }

    SmallDeformationNonlocalProcessData(
        SmallDeformationNonlocalProcessData&& other)
        : material{std::move(other.material)},
          dt{other.dt},
          t{other.t},
          internal_length{other.internal_length}
    {
    }

    //! Copies are forbidden.
    SmallDeformationNonlocalProcessData(SmallDeformationNonlocalProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationNonlocalProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationNonlocalProcessData&&) = delete;

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    double dt = 0;
    double t = 0;
    double const internal_length;
};

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
