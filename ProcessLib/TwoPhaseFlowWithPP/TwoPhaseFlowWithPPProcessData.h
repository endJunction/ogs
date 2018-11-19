/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once
#include "TwoPhaseFlowWithPPMaterialProperties.h"
namespace MaterialPropertyLib
{
class Medium;
}
namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace TwoPhaseFlowWithPP
{
struct TwoPhaseFlowWithPPProcessData
{
    TwoPhaseFlowWithPPProcessData(
        Eigen::VectorXd const specific_body_force_,
        bool const has_gravity_,
        bool const has_mass_lumping_,
        Parameter<double> const& temperature_,
        std::unique_ptr<TwoPhaseFlowWithPPMaterialProperties>&& material_,
        std::map<int, std::unique_ptr<MaterialPropertyLib::Medium>> const& media_,
        MeshLib::PropertyVector<int> const* material_ids_)
        : specific_body_force(specific_body_force_),
          has_gravity(has_gravity_),
          has_mass_lumping(has_mass_lumping_),
          temperature(temperature_),
          material(std::move(material_)),
          media(media_),
          material_ids(material_ids_)

    {
    }

    TwoPhaseFlowWithPPProcessData(TwoPhaseFlowWithPPProcessData&& other)
        : specific_body_force(other.specific_body_force),
          has_gravity(other.has_gravity),
          has_mass_lumping(other.has_mass_lumping),
          temperature(other.temperature),
          material(std::move(other.material)),
          media(other.media),
          material_ids(other.material_ids)
    {
    }

    //! Copies are forbidden.
    TwoPhaseFlowWithPPProcessData(TwoPhaseFlowWithPPProcessData const&) =
        delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseFlowWithPPProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseFlowWithPPProcessData&&) = delete;

    //! Specific body forces applied to solid and fluid.
    //! It is usually used to apply gravitational forces.
    //! A vector of displacement dimension's length.
    Eigen::VectorXd const specific_body_force;

    bool const has_gravity;

    //! Enables lumping of the mass matrix.
    bool const has_mass_lumping;
    Parameter<double> const& temperature;
    std::unique_ptr<TwoPhaseFlowWithPPMaterialProperties> material;
    std::map<int, std::unique_ptr<MaterialPropertyLib::Medium>> const& media;
    MeshLib::PropertyVector<int> const* const material_ids;
};

}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib
