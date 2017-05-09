/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#ifdef PROTOBUF_FOUND
#include "SerializationLib/integration_point.pb.h"
#endif  // PROTOBUF_FOUND

namespace ProcessLib
{
template <typename LocalAssembler>
void readSmallDeformationIntegrationPointData(
    std::vector<char> const& data, LocalAssembler& local_assembler)
{
#ifdef PROTOBUF_FOUND
    OGS::ElementData element_data;
    if (!element_data.ParseFromArray(data.data(), data.size()))
        OGS_FATAL("Parsing ElementData protobuf failed.");

    // check element number
    if (local_assembler._element.getID() != element_data.element_id())
        OGS_FATAL("Reading input failed somewhat. Mesh item id does not match");

    // check number of integration points
    unsigned const n_integration_points =
        local_assembler._integration_method.getNumberOfPoints();
    if (n_integration_points != element_data.n_integration_points())
        OGS_FATAL(
            "Reading input failed somewhat. The value of "
            "n_integration_points does not match");

    // The actual data load
    if (!element_data.has_small_deformation())
        OGS_FATAL(
            "Reading input data failed: Expected SmallDeformation message "
            "is not set.");
    auto const small_deformation_data = element_data.small_deformation();

    // sigma
    assert(n_integration_points == small_deformation_data.sigma_size());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto sigma = small_deformation_data.sigma(ip);
        if (LocalAssembler::DisplacementDim != sigma.dimension())
            OGS_FATAL("Dimension of a Kelvin vector do not match.");
        assert(local_assembler._ip_data[ip].sigma.size() ==
               sigma.value_size());

        for (int i = 0; i < local_assembler._ip_data[ip].sigma.size(); ++i)
            local_assembler._ip_data[ip].sigma[i] = sigma.value(i);
    }

    // epsilon
    assert(n_integration_points == small_deformation_data.eps_size());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto eps = small_deformation_data.eps(ip);
        if (LocalAssembler::DisplacementDim != eps.dimension())
            OGS_FATAL("Dimension of a Kelvin vector do not match.");
        assert(local_assembler._ip_data[ip].eps.size() == eps.value_size());

        for (int i = 0; i < local_assembler._ip_data[ip].eps.size(); ++i)
            local_assembler._ip_data[ip].eps[i] = eps.value(i);
    }
#else   // PROTOBUF_FOUND
    (void)data;  // Unused argument
#endif  // PROTOBUF_FOUND
}

template <typename LocalAssembler>
std::size_t writeSmallDeformationIntegrationPointData(
    std::vector<char>& data, LocalAssembler const& local_assembler)
{
#ifdef PROTOBUF_FOUND
    unsigned const n_integration_points =
        local_assembler._integration_method.getNumberOfPoints();

    OGS::ElementData element_data;
    element_data.set_element_id(local_assembler._element.getID());
    element_data.set_n_integration_points(n_integration_points);

    auto small_deformation_data = element_data.mutable_small_deformation();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto sigma = small_deformation_data->add_sigma();
        sigma->set_dimension(LocalAssembler::DisplacementDim);
        for (int i = 0; i < local_assembler._ip_data[ip].sigma.size(); ++i)
            sigma->add_value(local_assembler._ip_data[ip].sigma[i]);
    }

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto eps = small_deformation_data->add_eps();
        eps->set_dimension(LocalAssembler::DisplacementDim);
        for (int i = 0; i < local_assembler._ip_data[ip].eps.size(); ++i)
            eps->add_value(local_assembler._ip_data[ip].eps[i]);
    }

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& ip_data = local_assembler._ip_data[ip];
        auto material_state = small_deformation_data->add_material_state();
        material_state->CopyFrom(ip_data.solid_material.writeMaterialState(
            *ip_data.material_state_variables));
    }

    data.resize(element_data.ByteSize());
    element_data.SerializeToArray(data.data(), element_data.ByteSize());

    return element_data.ByteSize();
#else   // PROTOBUF_FOUND
    (void)data;  // Unused argument
    return 0;    // Dummy value needed for compilation. Code is not executed
                 // because the integration_point_writer is not created in
                 // absence of protobuffer.
#endif  // PROTOBUF_FOUND
};

}  // namespace ProcessLib
