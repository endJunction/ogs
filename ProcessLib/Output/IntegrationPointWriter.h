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
#include <memory>
#include <vector>

#pragma once

namespace MeshLib
{
class Mesh;
}

namespace ProcessLib
{
struct IntegrationPointWriter
{
    virtual ~IntegrationPointWriter() = default;

    virtual int numberOfComponents() const = 0;
    virtual int integrationOrder() const = 0;
    virtual std::string name() const = 0;
    virtual std::vector<std::vector<double>> values() const = 0;
};

/// Add integration point data the the mesh's properties.
///
/// Integration point data stored as field data (contrary to point or cell
/// data). The data is supplemented with information in JSON format.
void addIntegrationPointWriter(
    MeshLib::Mesh& mesh,
    std::vector<std::unique_ptr<IntegrationPointWriter>> const&
        integration_point_writer);

struct IntegrationPointMetaData
{
    std::string const name;
    int const n_components;
    int const integration_order;
};

}  // namespace ProcessLib
