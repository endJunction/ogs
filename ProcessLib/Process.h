/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PROCESS_H_
#define PROCESS_LIB_PROCESS_H_

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{

class Process
{
public:
    Process(MeshLib::Mesh const& mesh)
        : _mesh(mesh)
    { }

    virtual ~Process() = default;

    virtual void initialize() = 0;
    virtual void solve() = 0;

    /// Postprocessing after solve().
    /// The file_name is indicating the name of possible output file.
    virtual void post(std::string const& file_name) = 0;

protected:
    MeshLib::Mesh const& _mesh;
};

}   // namespace ProcessLib

//
// Include all known processes here.
//
#include "GroundwaterFlowProcess.h"

#endif  // PROCESS_LIB_PROCESS_H_
