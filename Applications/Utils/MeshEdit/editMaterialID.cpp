/**
 * @copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "BaseLib/LogogSimpleFormatter.h"

// FileIO
#include "FileIO/Legacy/MeshIO.h"
#include "FileIO/readMeshFromFile.h"
#include "FileIO/writeMeshToFile.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshEditing/ElementValueModification.h"

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Edit material IDs of mesh elements.", ' ', "0.1");
	TCLAP::SwitchArg replaceArg("r", "replace", "replace material IDs", false);
	TCLAP::SwitchArg condenseArg("c", "condense", "condense material IDs", false);
	TCLAP::SwitchArg specifyArg("s", "specify", "specify material IDs by element types (-e)", false);
	std::vector<TCLAP::Arg*> vec_xors;
	vec_xors.push_back(&replaceArg);
	vec_xors.push_back(&condenseArg);
	vec_xors.push_back(&specifyArg);
	cmd.xorAdd(vec_xors);
	TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
	                                     "the name of the file containing the input mesh", true,
	                                     "", "file name");
	cmd.add(mesh_in);
	TCLAP::ValueArg<std::string> mesh_out("o", "mesh-output-file",
	                                      "the name of the file the mesh will be written to", true,
	                                      "", "file name");
	cmd.add(mesh_out);
	TCLAP::MultiArg<unsigned> matIDArg("m", "current-material-id",
	                                      "current material id to be replaced", false, "number");
	cmd.add(matIDArg);
	TCLAP::ValueArg<unsigned> newIDArg("n", "new-material-id",
	                                      "new material id", false, 0, "number");
	cmd.add(newIDArg);
	std::vector<std::string> eleList(MeshLib::getMeshElemTypeStringsShort());
	TCLAP::ValuesConstraint<std::string> allowedVals(eleList);
	TCLAP::ValueArg<std::string> eleTypeArg("e", "element-type",
	                                      "element type", false, "", &allowedVals);
	cmd.add(eleTypeArg);

	cmd.parse(argc, argv);

	if (!replaceArg.isSet() && !condenseArg.isSet() && !specifyArg.isSet()) {
		INFO("Please select editing mode: -r or -c or -s");
		return 0;
	} else if (replaceArg.isSet() && condenseArg.isSet()) {
		INFO("Please select only one editing mode: -r or -c or -s");
		return 0;
	} else if (replaceArg.isSet()) {
		if (!matIDArg.isSet() || !newIDArg.isSet()) {
			INFO("current and new material IDs must be provided for replacement");
			return 0;
		}
	} else if (specifyArg.isSet()) {
		if (!eleTypeArg.isSet() || !newIDArg.isSet()) {
			INFO("element type and new material IDs must be provided to specify elements");
			return 0;
		}
	}

	MeshLib::Mesh* mesh (FileIO::readMeshFromFile(mesh_in.getValue()));
	INFO("Mesh read: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());

	if (condenseArg.isSet()) {
		INFO("Condensing material ID...");
		MeshLib::ElementValueModification::condense(*mesh);
	} else if (replaceArg.isSet()) {
		INFO("Replacing material ID...");
		const auto vecOldID = matIDArg.getValue();
		const unsigned newID = newIDArg.getValue();
		for (auto oldID : vecOldID) {
			INFO("%d -> %d", oldID, newID);
			MeshLib::ElementValueModification::replace(*mesh, oldID, newID, true);
		}
	} else if (specifyArg.isSet()) {
		INFO("Specifying material ID...");
		const std::string eleTypeName(eleTypeArg.getValue());
		const MeshLib::MeshElemType eleType = MeshLib::String2MeshElemType(eleTypeName);
		const unsigned newID = newIDArg.getValue();
		unsigned cnt = MeshLib::ElementValueModification::setByElementType(*mesh, eleType, newID);
		INFO("updated %d elements", cnt);
	}

	// write into a file
	FileIO::writeMeshToFile(*mesh, mesh_out.getValue());

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}

