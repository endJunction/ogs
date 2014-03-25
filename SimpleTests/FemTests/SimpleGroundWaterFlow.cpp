/**
 * \date   2013-10-08
 * \brief  Implementation of simple ground water flow process
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cstdlib>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// ThirdParty/tclap
#include "ThirdParty/tclap/CmdLine.h"

// BaseLib
#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/FileTools.h"
#include "Configure.h"

// FileIO
#include "XmlIO/Boost/BoostXmlCndInterface.h"
#include "XmlIO/Boost/BoostXmlGmlInterface.h"
#include "readMeshFromFile.h"

// GeoLib
#include "GEOObjects.h"

// MeshLib
#include "Mesh.h"

// OGS
#include "ProjectData.h"

int main(int argc, char *argv[])
{
	// logog
	LOGOG_INITIALIZE();
	BaseLib::LogogSimpleFormatter *custom_format(new BaseLib::LogogSimpleFormatter);
	logog::Cout *logog_cout(new logog::Cout);
	logog_cout->SetFormatter(*custom_format);

	// tclap
	TCLAP::CmdLine cmd(
			"Simple ground water flow test, reading mesh (only 2d quad elements), geometry and bc and simulate ground water flow",
			' ', "0.1");

	TCLAP::ValueArg<std::string> mesh_arg("m", "mesh", "file name of the mesh", true, "", "string");
	cmd.add(mesh_arg);

	TCLAP::ValueArg<std::string> geometry_arg("g", "geometry", "file name of the geometry", true,
			"", "string");
	cmd.add(geometry_arg);

	TCLAP::ValueArg<std::string> bc_arg("", "boundary_condition",
			"file name of the boundary condition", true, "", "string");
	cmd.add(bc_arg);

	cmd.parse(argc, argv);

	ProjectData project_data;

	// *** read geometry
	{
		FileIO::BoostXmlGmlInterface geo_io(project_data);
		geo_io.readFile(geometry_arg.getValue());
	}
	std::string unique_name;

	// *** read mesh
	std::string mesh_name(mesh_arg.getValue());
	project_data.addMesh(FileIO::readMeshFromFile(mesh_name));
	mesh_name = BaseLib::extractBaseNameWithoutExtension(mesh_name);

	// *** read boundary conditions
	FileIO::BoostXmlCndInterface xml_io(project_data);
	xml_io.readFile(bc_arg.getValue());

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
