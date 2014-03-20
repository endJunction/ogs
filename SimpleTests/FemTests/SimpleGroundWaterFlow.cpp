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

//#undef LIS
#define LIS

#include <cstdlib>

#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>
#endif

// AssemblerLib
#include "AssemblerLib/GlobalSetup.h"
#include "AssemblerLib/VectorMatrixAssembler.h"
#ifdef LIS
    #include "AssemblerLib/SerialExecutor.h"
    #include "AssemblerLib/SerialLisVectorMatrixBuilder.h"
#else
    #include "AssemblerLib/SerialDenseSetup.h"
#endif

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
#include "GeoObject.h"

// MathLib
#ifdef LIS
    #include "MathLib/LinAlg/Lis/LisTools.h"
    #include "MathLib/LinAlg/Lis/LisLinearSolver.h"
#else
    #include "MathLib/LinAlg/Dense/DenseTools.h"
    #include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
#endif
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

template <typename ElemType>
class LocalFeQuad4AssemblyItem
{
public:
	// definition of vector and matrix types
	typedef Eigen::Matrix<double, ElemType::NPOINTS, ElemType::NPOINTS, Eigen::RowMajor> NodalMatrixType;
	typedef Eigen::Matrix<double, ElemType::NPOINTS, 1> NodalVectorType;
	typedef Eigen::Matrix<double, ElemType::DIM, ElemType::NPOINTS, Eigen::RowMajor> DimNodalMatrixType;
	typedef Eigen::Matrix<double, ElemType::DIM, ElemType::DIM, Eigen::RowMajor> DimMatrixType;

	// type definition of FeQuad4 type
	typedef typename NumLib::FeQUAD4<
		NodalVectorType,
		DimNodalMatrixType,
		DimMatrixType>::type FeQuad4;

	// type definition of ShapeMatricesType
	typedef typename FeQuad4::ShapeMatricesType ShapeMatricesType;

	LocalFeQuad4AssemblyItem() :
		_shape_mat(ShapeMatricesType(3,4)),
		_material(1.0)
	{}

	ShapeMatricesType _shape_mat;
	double _material;
};

template <typename ElemType>
class LocalGWAssembler
{

public:
    typedef Eigen::Matrix<double, ElemType::NPOINTS, ElemType::NPOINTS, Eigen::RowMajor> NodalMatrixType;
    typedef Eigen::Matrix<double, ElemType::NPOINTS, 1> NodalVectorType;
    typedef Eigen::Matrix<double, ElemType::DIM, ElemType::NPOINTS, Eigen::RowMajor> DimNodalMatrixType;
    typedef Eigen::Matrix<double, ElemType::DIM, ElemType::DIM, Eigen::RowMajor> DimMatrixType;

	// Dynamic size local matrices are much slower.
	//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> NodalMatrixType;
	//typedef Eigen::Matrix<double, Eigen::Dynamic, 1> NodalVectorType;
	//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::RowMajor> DimNodalMatrixType;
	//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DimMatrixType;
	typedef typename NumLib::FeQUAD4<
		NodalVectorType,
		DimNodalMatrixType,
		DimMatrixType>::type FeQuad4;

	typedef typename FeQuad4::ShapeMatricesType ShapeMatricesType;

public:
	LocalGWAssembler() :
			_integration_method(2)
	{}

	void operator()(const MeshLib::Element& e, NodalMatrixType &localA,
			NodalVectorType & /*rhs*/,
			ShapeMatricesType& shape)
	{
		localA.setZero();

		// init FeQUAD4
		_fe_quad4.setMeshElement(*static_cast<const MeshLib::Quad*>(&e));

		for (std::size_t ip(0); ip < _integration_method.getNPoints(); ip++) { // ip == number of gauss point

			MathLib::WeightedPoint2D const& wp = _integration_method.getWeightedPoint(ip);
			_fe_quad4.computeShapeFunctions(wp.getCoords(), shape);
			localA += shape.dNdx.transpose() * shape.dNdx * shape.detJ * wp.getWeight();
		}
	}

private:
	typename FeQuad4::IntegrationMethod _integration_method;
	FeQuad4 _fe_quad4;
};


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

	std::vector<FEMCondition*> bcs(
			project_data.getConditions(FiniteElement::GROUNDWATER_FLOW, unique_name,
					FEMCondition::BOUNDARY_CONDITION));

	std::vector < std::size_t > bc_mesh_node_ids;
	std::vector <double> bc_values;
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
#ifdef LIS
	typedef AssemblerLib::GlobalSetup<
		AssemblerLib::SerialLisVectorMatrixBuilder,
		AssemblerLib::SerialExecutor> GlobalSetup;
#else
	typedef AssemblerLib::SerialDenseSetup GlobalSetup;
#endif
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
	std::vector<MeshLib::Element*> const& all_eles = mesh.getElements();

	std::vector < std::vector<std::size_t> > map_ele_nodes2vec_entries;
	map_ele_nodes2vec_entries.reserve(all_eles.size());
	for (auto e = all_eles.cbegin(); e != all_eles.cend(); ++e) {
		std::size_t const nnodes = (*e)->getNNodes();
		std::size_t const mesh_id = mesh.getID();
		std::vector<MeshLib::Location> vec_items;
		vec_items.reserve(nnodes);
		for (std::size_t j = 0; j < nnodes; j++)
			vec_items.emplace_back(mesh_id, MeshLib::MeshItemType::Node, (*e)->getNode(j)->getID());

		map_ele_nodes2vec_entries.push_back(
				vec1_composition.getGlobalIndices<AssemblerLib::ComponentOrder::BY_COMPONENT>(
						vec_items));
	}

	//
	// Local and global assemblers.
	//
	LocalGWAssembler<NumLib::ShapeQuad4> local_gw_assembler;
	typedef typename LocalGWAssembler<NumLib::ShapeQuad4>::NodalMatrixType LocalMatrix;
	typedef typename LocalGWAssembler<NumLib::ShapeQuad4>::NodalVectorType LocalVector;

	typedef AssemblerLib::VectorMatrixAssembler<
			GlobalMatrix,
			GlobalVector,
			MeshLib::Element,
			LocalGWAssembler<NumLib::ShapeQuad4>,
			LocalMatrix,
			LocalVector > GlobalAssembler;

	GlobalAssembler global_assembler(*A.get(), *rhs.get(), local_gw_assembler,
			AssemblerLib::LocalToGlobalIndexMap(map_ele_nodes2vec_entries));

	// create data structures for shape matrices
	typedef typename NumLib::FeQUAD4<
		LocalGWAssembler<NumLib::ShapeQuad4>::NodalVectorType,
		LocalGWAssembler<NumLib::ShapeQuad4>::DimNodalMatrixType,
		LocalGWAssembler<NumLib::ShapeQuad4>::DimMatrixType>::type FeQuad4;

	std::vector<typename FeQuad4::ShapeMatricesType> shape_matrix_vec;
	shape_matrix_vec.reserve(mesh.getNElements());
	std::fill_n(std::back_inserter(shape_matrix_vec), mesh.getNElements(), typename FeQuad4::ShapeMatricesType(3,4));
	
	// Call global assembler for each mesh element.
	global_setup.execute(global_assembler, mesh.getElements(), shape_matrix_vec);

	// apply Dirichlet BC
	MathLib::applyKnownSolution(*A, *rhs, bc_mesh_node_ids, bc_values);
	//--------------------------------------------------------------------------
	// solve x=A^-1 rhs
	//--------------------------------------------------------------------------
#ifdef LIS
	MathLib::LisLinearSolver ls(*A);
#else
	MathLib::GaussAlgorithm<GlobalMatrix, GlobalVector> ls(*A);
#endif
	ls.solve(*rhs, *x);

	if (x->size() > 1000) {
		std::ofstream out("results.txt");
		for (std::size_t i = 0; i < x->size(); i++) {
			out << (*x)[i] << " ";
		}
		out << std::endl;
		out.close();
	} else {
		for (std::size_t i = 0; i < x->size(); i++) {
			std::cout << (*x)[i] << " ";
		}
	}

	std::remove_if(vec_comp_dis.begin(), vec_comp_dis.end(),
			[](MeshLib::MeshSubsets* v) { delete v; return true; });
	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
