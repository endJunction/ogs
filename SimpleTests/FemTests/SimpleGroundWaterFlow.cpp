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
#include <functional>
#include <unordered_map>
#include <typeinfo>
#include <typeindex>

#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>
#endif

// AssemblerLib
#include "AssemblerLib/GlobalSetup.h"
#include "AssemblerLib/SimpleAssembler.h"
#include "AssemblerLib/VectorMatrixAssembler.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// ThirdParty/tclap
#include "ThirdParty/tclap/CmdLine.h"

// BaseLib
#include "BaseLib/metaprogramming.h"
#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/FileTools.h"
#include "Configure.h"

// FileIO
#include "XmlIO/Boost/BoostXmlCndInterface.h"
#include "XmlIO/Boost/BoostXmlGmlInterface.h"
#include "XmlIO/Boost/BoostVtuInterface.h"
#include "readMeshFromFile.h"

// GeoLib
#include "GEOObjects.h"
#include "GeoObject.h"

// MathLib
#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#include "MathLib/TemplateWeightedPoint.h"

// MeshGeoToolsLib
#include "MeshNodeSearcher.h"
#include "MeshNodesToPoints.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSubsets.h"
#include "MeshLib/Elements/Quad.h"

// NumLib
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

// OGS
#include "BoundaryCondition.h"
#include "ProjectData.h"

const std::size_t NumLib::ShapeQuad4::DIM;
const std::size_t NumLib::ShapeQuad4::NPOINTS;

// The following defines types depending on Lis or DirectSolver choice.
// Right includes are also made.

//#undef LIS
#define LIS

#ifdef LIS
	#include "AssemblerLib/SerialExecutor.h"
	#include "AssemblerLib/SerialLisVectorMatrixBuilder.h"

	#include "MathLib/LinAlg/Lis/LisTools.h"
	#include "MathLib/LinAlg/Lis/LisLinearSolver.h"

	typedef AssemblerLib::GlobalSetup<
		AssemblerLib::SerialLisVectorMatrixBuilder,
		AssemblerLib::SerialExecutor> GlobalSetup;

	typedef MathLib::LisLinearSolver LinearSolver;
#else	// LIS
	#include "AssemblerLib/SerialDenseSetup.h"
	#include "MathLib/LinAlg/Dense/DenseTools.h"

	#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
	typedef AssemblerLib::SerialDenseSetup GlobalSetup;

	typedef MathLib::GaussAlgorithm<GlobalSetup::MatrixType, GlobalSetup::VectorType> LinearSolver;
#endif	// LIS

template <typename ShapeType>
struct EigenFixedSizeShapeMatrices
{
	typedef Eigen::Matrix<double, ShapeType::NPOINTS, ShapeType::NPOINTS, Eigen::RowMajor> NodalMatrixType;
	typedef Eigen::Matrix<double, ShapeType::NPOINTS, 1> NodalVectorType;
	typedef Eigen::Matrix<double, ShapeType::DIM, ShapeType::NPOINTS, Eigen::RowMajor> DimNodalMatrixType;
	typedef Eigen::Matrix<double, ShapeType::DIM, ShapeType::DIM, Eigen::RowMajor> DimMatrixType;

	// Dynamic size local matrices are much slower.
	// XXX This is not working because the matrices are not properly
	// initialized.
	//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> NodalMatrixType;
	//typedef Eigen::Matrix<double, Eigen::Dynamic, 1> NodalVectorType;
	//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DimNodalMatrixType;
	//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DimMatrixType;
};

struct LocalAssemblerData
{
	virtual ~LocalAssemblerData() = default;
};

template <typename FemType_, std::size_t _INTEGRATION_ORDER, typename PROPERTY_TYPE>
struct LocalGWAssemblerData : public LocalAssemblerData
{
	typedef typename FemType_::type FemType;
	typedef typename FemType::NodalMatrixType NodalMatrixType;
	typedef typename FemType::NodalVectorType NodalVectorType;

	typedef typename FemType::ShapeMatricesType ShapeMatricesType;

	typedef PROPERTY_TYPE PropertyType;

	static std::size_t constexpr INTEGRATION_ORDER = _INTEGRATION_ORDER;
	static std::size_t constexpr N_INTEGRATION_POINTS =
		BaseLib::pow(INTEGRATION_ORDER, FemType::ShapeFunctionType::DIM);

	typename std::remove_pointer<PROPERTY_TYPE>::type getMaterial()
	{
		return getMaterialProperty(_material);
	}

	void setMaterial(PROPERTY_TYPE const& material)
	{
		_material = material;
	}

	std::array<ShapeMatricesType, N_INTEGRATION_POINTS> _shape_matrices;
private:
	typename std::remove_pointer<PROPERTY_TYPE>::type
	getMaterialProperty(
		typename std::add_pointer<typename std::remove_pointer<PROPERTY_TYPE>::type>::type /*material*/)
	{
		return *_material;
	}

	typename std::remove_pointer<PROPERTY_TYPE>::type
	getMaterialProperty(typename std::remove_pointer<PROPERTY_TYPE>::type /*material*/)
	{
		return _material;
	}
	PROPERTY_TYPE _material;
};

std::unordered_map<std::type_index, std::function<LocalAssemblerData*()>> builder;

template <typename Data>
class ShapeMatricesInitializer
{
public:
	ShapeMatricesInitializer(std::array<typename Data::PropertyType, 4> material_values_) :
		_integration_method(Data::INTEGRATION_ORDER),
		_material_values(material_values_)
	{}

	void operator()(const MeshLib::Element& e,
		Data*& data_ptr)
	{
		std::cerr << std::type_index(typeid(e)).hash_code() << std::endl;
		data_ptr = builder[std::type_index(typeid(e))]();

		Data& data = *data_ptr;
		data.setMaterial(_material_values[e.getValue()]);
		typedef typename Data::FemType::MeshElementType MeshElementType;
		// create FEM Element
		typename Data::FemType fe(*static_cast<const MeshElementType*>(&e));

		for (std::size_t ip(0); ip < Data::N_INTEGRATION_POINTS; ip++) { // ip == number of gauss point
			MathLib::WeightedPoint2D const& wp = _integration_method.getWeightedPoint(ip);

			/*
			static std::size_t const DIM = Data::FemType::ShapeFunctionType::DIM;
			static std::size_t const NODES = Data::FemType::ShapeFunctionType::NPOINTS ;
			data._shape_matrices[ip].init(DIM, NODES);
			*/
			fe.computeShapeFunctions(wp.getCoords(), data._shape_matrices[ip]);
		}
	}

private:
	typename Data::FemType::IntegrationMethod _integration_method;
	std::array<typename Data::PropertyType, 4> _material_values;
};

template <typename Data_>
class LocalGWAssembler
{
public:
	typedef Data_ Data;

	typedef typename Data::NodalVectorType NodalVectorType;
	typedef typename Data::NodalMatrixType NodalMatrixType;

public:
	LocalGWAssembler() :
			_integration_method(Data::INTEGRATION_ORDER)
	{}

	void operator()(NodalMatrixType &localA,
			NodalVectorType & /*rhs*/,
			Data*& data_ptr) const
	{
		localA.setZero();

		Data& data = *data_ptr;
		for (std::size_t ip(0); ip < Data::N_INTEGRATION_POINTS; ip++) { // ip == number of gauss point
			MathLib::WeightedPoint2D const& wp = _integration_method.getWeightedPoint(ip);
			localA += data._shape_matrices[ip].dNdx.transpose() * data.getMaterial() * data._shape_matrices[ip].dNdx * data._shape_matrices[ip].detJ * wp.getWeight();
		}
	}

private:
	typename Data::FemType::IntegrationMethod _integration_method;
};

void prepareBCForSimulation(ProjectData const& project_data,
	std::string const& mesh_name,
	std::vector<std::size_t> & bc_mesh_node_ids,
	std::vector<double> & bc_values)
{
	std::string unique_name;
	// get boundary conditions from ProjectData instance
	std::vector<FEMCondition*> bcs(
		project_data.getConditions(FiniteElement::GROUNDWATER_FLOW, unique_name,
			FEMCondition::BOUNDARY_CONDITION));

	MeshGeoToolsLib::MeshNodeSearcher searcher(*project_data.getMesh(mesh_name));
	for (auto it(bcs.cbegin()); it != bcs.cend(); it++) {
		// fetch geometry obj from condition obj
		GeoLib::GeoObject const* geom_obj((*it)->getGeoObj());
		if (dynamic_cast<GeoLib::Point const*>(geom_obj) != nullptr) {
			GeoLib::Point const& pnt(*dynamic_cast<GeoLib::Point const*>(geom_obj));
			bc_mesh_node_ids.push_back(searcher.getMeshNodeIDForPoint(pnt));
			bc_values.push_back((*it)->getDisValues()[0]);
		} else {
			if (dynamic_cast<GeoLib::Polyline const*>(geom_obj) != nullptr) {
				GeoLib::Polyline const& ply(*dynamic_cast<GeoLib::Polyline const*>(geom_obj));
				std::vector<std::size_t> const& ids(searcher.getMeshNodeIDsAlongPolyline(ply));
				bc_mesh_node_ids.insert(bc_mesh_node_ids.end(), ids.cbegin(), ids.cend());
				for (std::size_t k(0); k<bc_mesh_node_ids.size(); k++)
					bc_values.push_back((*it)->getDisValues()[0]);
			}
		}
	}
}

std::vector<std::vector<std::size_t>>
createDOFMapping(MeshLib::Mesh const& mesh,
	AssemblerLib::MeshComponentMap vec1_composition)
{
	std::vector<std::vector<std::size_t>> map_ele_nodes2vec_entries;
	std::vector<MeshLib::Element*> const& all_eles = mesh.getElements();
	map_ele_nodes2vec_entries.reserve(all_eles.size());
	for (auto e = all_eles.cbegin(); e != all_eles.cend(); ++e) {
		std::size_t const nnodes = (*e)->getNNodes();
		std::size_t const mesh_id = mesh.getID();
		std::vector<MeshLib::Location> vec_items;
		vec_items.reserve(nnodes);
		for (std::size_t j = 0; j < nnodes; j++) {
			vec_items.emplace_back(mesh_id, MeshLib::MeshItemType::Node, (*e)->getNode(j)->getID());
		}

		map_ele_nodes2vec_entries.push_back(
			vec1_composition.getGlobalIndices<AssemblerLib::ComponentOrder::BY_COMPONENT>(
				vec_items));
	}

	return map_ele_nodes2vec_entries;
}

std::tuple<std::string, std::string, std::string, std::size_t>
parseCLI(int argc, char *argv[])
{
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

	TCLAP::ValueArg<std::size_t> time_steps_arg("t", "time_steps",
			"Timesteps to perform [default 1]. Use 1 for steady-state.",
			false, 1, "int");
	cmd.add(time_steps_arg);

	cmd.parse(argc, argv);

	return std::make_tuple(
			mesh_arg.getValue(),
			geometry_arg.getValue(),
			bc_arg.getValue(),
			time_steps_arg.getValue());
}

int main(int argc, char *argv[])
{
	// logog
	LOGOG_INITIALIZE();
	BaseLib::LogogSimpleFormatter *custom_format(new BaseLib::LogogSimpleFormatter);
	logog::Cout *logog_cout(new logog::Cout);
	logog_cout->SetFormatter(*custom_format);

	// Parse CLI arguments.
	std::string mesh_file;
	std::string geometry_file;
	std::string boundary_condition_file;
	std::size_t end_time;
	std::tie(mesh_file, geometry_file, boundary_condition_file, end_time)
		= parseCLI(argc, argv);

	ProjectData project_data;

	// *** read geometry
	{
		FileIO::BoostXmlGmlInterface geo_io(project_data);
		geo_io.readFile(geometry_file);
	}

	// *** read mesh
	project_data.addMesh(FileIO::readMeshFromFile(mesh_file));
	std::string const mesh_name = BaseLib::extractBaseNameWithoutExtension(mesh_file);

	// *** read boundary conditions
	{
		FileIO::BoostXmlCndInterface cnd_io(project_data);
		cnd_io.readFile(boundary_condition_file);
	}

	// *** prepare boundary condition for using in simulation
	std::vector<std::size_t> bc_mesh_node_ids;
	std::vector<double> bc_values;
	prepareBCForSimulation(project_data, mesh_name, bc_mesh_node_ids, bc_values);

	//--------------------------------------------------------------------------
	// Prepare mesh items where data is assigned
	//--------------------------------------------------------------------------
	MeshLib::Mesh const& mesh(*project_data.getMesh(mesh_name));
	const MeshLib::MeshSubset mesh_items_all_nodes(mesh, mesh.getNodes());

	//-------------------------------------------------------------------------
	// Allocate a (global) coefficient matrix, RHS and solution vectors
	//-------------------------------------------------------------------------
	// define a mesh item composition in a vector
	std::vector<MeshLib::MeshSubsets*> vec_comp_dis;
	vec_comp_dis.push_back(new MeshLib::MeshSubsets(&mesh_items_all_nodes));
	AssemblerLib::MeshComponentMap vec1_composition(vec_comp_dis,
			AssemblerLib::ComponentOrder::BY_COMPONENT);

	//--------------------------------------------------------------------------
	// Choose implementation type
	//--------------------------------------------------------------------------
	const GlobalSetup global_setup;

	// allocate a vector and matrix
	typedef GlobalSetup::VectorType GlobalVector;
	typedef GlobalSetup::MatrixType GlobalMatrix;
	std::unique_ptr < GlobalMatrix > A(global_setup.createMatrix(vec1_composition));
	A->setZero();
	std::unique_ptr < GlobalVector > rhs(global_setup.createVector(vec1_composition));
	std::unique_ptr < GlobalVector > x(global_setup.createVector(vec1_composition));

	//--------------------------------------------------------------------------
	// Construct a linear system
	//--------------------------------------------------------------------------
	// create a mapping table from element nodes to entries in the linear system
	std::vector <std::vector<std::size_t>> const map_ele_nodes2vec_entries =
		createDOFMapping(mesh, vec1_composition);

	builder[std::type_index(typeid(MeshLib::Quad))] =
		[](){ return new LocalGWAssemblerData<NumLib::FeQUAD4<EigenFixedSizeShapeMatrices>, 2, double>; };

	// create data structures for properties
	typedef LocalGWAssemblerData<
		NumLib::FeQUAD4<EigenFixedSizeShapeMatrices>,
		2,
		double> LAData;
	std::vector<LAData*> local_assembly_item_vec;
	local_assembly_item_vec.resize(mesh.getNElements());

	std::array<double,4> material_values({{1e-10, 2e-10, 4e-10, 8e-10}});

	//
	// Shape matrices initializer
	//
	typedef LocalGWAssembler<LAData> LA;
	typedef ShapeMatricesInitializer<LocalAssemblerData> SMI;
	SMI shape_matrices_initializer(material_values);

	typedef AssemblerLib::SimpleAssembler<
			MeshLib::Element,
			SMI
			> GlobalInitializer;

	GlobalInitializer global_initializer(shape_matrices_initializer);

	// Call global initializer for each mesh element.
	global_setup.execute(global_initializer, mesh.getElements(), local_assembly_item_vec);

	//
	// Local and global assemblers.
	//
	LA local_gw_assembler;
	typedef typename LA::NodalMatrixType LocalMatrix;
	typedef typename LA::NodalVectorType LocalVector;

	typedef AssemblerLib::VectorMatrixAssembler<
			GlobalMatrix,
			GlobalVector,
			LA,
			LocalMatrix,
			LocalVector > GlobalAssembler;

	GlobalAssembler global_assembler(*A.get(), *rhs.get(), local_gw_assembler,
			AssemblerLib::LocalToGlobalIndexMap(map_ele_nodes2vec_entries));

	// For steady-state problems use 1, which is the default value in CLI.
	for (std::size_t time = 0; time < end_time; ++time)
	{
		A->setZero();

		// Call global assembler for each local assembly item.
		global_setup.execute(global_assembler, local_assembly_item_vec);

		// apply Dirichlet BC
		MathLib::applyKnownSolution(*A, *rhs, bc_mesh_node_ids, bc_values);
		//--------------------------------------------------------------------------
		// solve x=A^-1 rhs
		//--------------------------------------------------------------------------
		LinearSolver ls(*A);
		ls.solve(*rhs, *x);

		{
			std::vector<double> heads;
			heads.reserve(x->size());
			for (std::size_t i = 0; i < x->size(); i++)
				heads.push_back((*x)[i]);

			FileIO::BoostVtuInterface vtu_io;
			vtu_io.setMesh(project_data.getMesh(mesh_name));
			vtu_io.addScalarPointProperty("Head", heads);
			std::string const res_mesh_name(BaseLib::dropFileExtension(mesh_file));
			vtu_io.writeToFile(res_mesh_name+"_with_results.vtu");
		}
	}

	std::remove_if(vec_comp_dis.begin(), vec_comp_dis.end(),
			[](MeshLib::MeshSubsets* v) { delete v; return true; });
	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
