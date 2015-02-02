/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
#define PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include "logog/include/logog.hpp"

#include "BaseLib/RunTime.h"

#include "AssemblerLib/LocalAssemblerBuilder.h"
#include "AssemblerLib/VectorMatrixAssembler.h"
#include "AssemblerLib/LocalDataInitializer.h"
#include "AssemblerLib/LocalToGlobalIndexMap.h"

#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/SetMatrixSparsity.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSubset.h"
#include "MeshLib/MeshSubsets.h"
#include "MeshLib/NodeAdjacencyTable.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"

#ifdef USE_PETSC
#include <petscmat.h>
#include "MeshLib/NodePartitionedMesh.h"
#include "MathLib/LinAlg/PETSc/PETScMatrixOption.h"
#endif

#include "BoundaryCondition.h"
#include "GroundwaterFlowFEM.h"
#include "ProcessVariable.h"

namespace ProcessLib
{

template<typename GlobalSetup>
class GroundwaterFlowProcess : public Process
{
    using ConfigTree = boost::property_tree::ptree;

    unsigned const _integration_order = 2;

public:
    GroundwaterFlowProcess(MeshLib::Mesh const& mesh,
            std::vector<ProcessVariable> const& variables,
            ConfigTree const& config)
        : Process(mesh)
    {
        DBUG("Create GroundwaterFlowProcess.");

        // Find the corresponding process variable.
        std::string const name = config.get<std::string>("process_variable");

        auto const& variable = std::find_if(variables.cbegin(), variables.cend(),
                [&name](ProcessVariable const& v) {
                    return v.getName() == name;
                });

        if (variable == variables.end())
            ERR("Expected process variable \'%s\' not found in provided variables list.",
                name.c_str());

        DBUG("Associate hydraulic_head with process variable \'%s\'.",
            name.c_str());
        _hydraulic_head = &*variable;
    }

    void initialize()
    {
        DBUG("Initialize GroundwaterFlowProcess.");

        DBUG("Construct dof mappings.");
        // Create single component dof in every of the mesh's nodes.
        _mesh_subset_all_nodes = new MeshLib::MeshSubset(_mesh, _mesh.getNodes());

        // Collect the mesh subsets in a vector.
        _all_mesh_subsets.push_back(new MeshLib::MeshSubsets(_mesh_subset_all_nodes));

        _local_to_global_index_map.reset(
            new AssemblerLib::LocalToGlobalIndexMap(_all_mesh_subsets));

#ifdef USE_PETSC
        DBUG("Allocate global matrix, vectors, and linear solver.");
        MathLib::PETScMatrixOption mat_opt;
        const MeshLib::NodePartitionedMesh &pmesh
                    = static_cast<const MeshLib::NodePartitionedMesh&>(_mesh);
        mat_opt.d_nz = pmesh.getMaximumNConnectedNodesToNode();
        mat_opt.o_nz = mat_opt.d_nz;
        _A.reset(_global_setup.createMatrix(_local_to_global_index_map->dofSizeGlobal(), mat_opt) );
        _x.reset(_global_setup.createVector(_local_to_global_index_map->dofSizeGlobal()));
        _rhs.reset(_global_setup.createVector(_local_to_global_index_map->dofSizeGlobal()));
        _linearSolver.reset(new typename GlobalSetup::LinearSolver(*_A, "gw_"));
#else
        DBUG("Compute sparsity pattern");
        _node_adjacency_table.createTable(_mesh.getNodes());

        DBUG("Allocate global matrix, vectors, and linear solver.");
        _A.reset(_global_setup.createMatrix(_local_to_global_index_map->dofSize()));
        _x.reset(_global_setup.createVector(_local_to_global_index_map->dofSize()));
        _rhs.reset(_global_setup.createVector(_local_to_global_index_map->dofSize()));
        _linearSolver.reset(new typename GlobalSetup::LinearSolver(*_A));
#endif

        DBUG("Create local assemblers.");
        // Populate the vector of local assemblers.
        _local_assemblers.resize(_mesh.getNElements());
        // Shape matrices initializer
        using LocalDataInitializer = AssemblerLib::LocalDataInitializer<
            GroundwaterFlow::LocalAssemblerDataInterface,
            GroundwaterFlow::LocalAssemblerData,
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType>;

        LocalDataInitializer initializer;

        using LocalAssemblerBuilder =
            AssemblerLib::LocalAssemblerBuilder<
                MeshLib::Element,
                LocalDataInitializer>;

        LocalAssemblerBuilder local_asm_builder(
            initializer, *_local_to_global_index_map);

        DBUG("Calling local assembler builder for all mesh elements.");
        _global_setup.execute(
                local_asm_builder,
                _mesh.getElements(),
                _local_assemblers,
                _hydraulic_conductivity,
                _integration_order);

        DBUG("Create global assembler.");
        _global_assembler.reset(
            new GlobalAssembler(*_A, *_rhs, *_local_to_global_index_map));

        DBUG("Initialize boundary conditions.");
        MeshGeoToolsLib::MeshNodeSearcher& hydraulic_head_mesh_node_searcher =
            MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
                _hydraulic_head->getMesh());

        using BCCI = ProcessVariable::BoundaryConditionCI;
        for (BCCI bc = _hydraulic_head->beginBoundaryConditions();
                bc != _hydraulic_head->endBoundaryConditions(); ++bc)
        {
            (*bc)->initialize(
                hydraulic_head_mesh_node_searcher,
                _dirichlet_bc.global_ids, _dirichlet_bc.values);
        }

    }

    void solve()
    {
        DBUG("Solve GroundwaterFlowProcess.");

        _A->setZero();
        MathLib::setMatrixSparsity(*_A, _node_adjacency_table);
        *_rhs = 0;   // This resets the whole vector.

        // Call global assembler for each local assembly item.
        BaseLib::RunTime assembly_wtimer;
        assembly_wtimer.start();
        //
        _global_setup.execute(*_global_assembler, _local_assemblers);
        //
        _elapsed_ctime += assembly_wtimer.elapsed();

#ifdef USE_PETSC
        std::vector<PetscInt> dbc_pos;
        std::vector<PetscScalar> dbc_value;
        const MeshLib::NodePartitionedMesh &mesh
                    = static_cast<const MeshLib::NodePartitionedMesh&>(_mesh);

        for(std::size_t i=0; i<_dirichlet_bc.global_ids.size(); i++)
        {
             if( mesh.isGhostNode(_dirichlet_bc.global_ids[i]) )
                continue;

             dbc_pos.push_back(static_cast<PetscInt>(mesh.getGlobalNodeID(_dirichlet_bc.global_ids[i])));
             dbc_value.push_back(static_cast<PetscScalar>(_dirichlet_bc.values[i]));
        }
        MathLib::applyKnownSolution(*_A, *_rhs, *_x, dbc_pos, dbc_value);
#else
        // Apply known values from the Dirichlet boundary conditions.
        MathLib::applyKnownSolution(*_A, *_rhs, _dirichlet_bc.global_ids, _dirichlet_bc.values);
#endif
        _linearSolver->solve(*_rhs, *_x);
    }

    void post(std::ostream& os)
    {
        DBUG("Postprocessing GroundwaterFlowProcess.");
        // Postprocessing of the linear system of equations solver results:
        // For example, write _x to _hydraulic_head or convert to velocity etc.

#ifdef USE_PETSC // Note: this is only a test
        int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

        std::vector<PetscScalar>  u(_x->size());
        _x->getGlobalVector(&u[0]);
        if(rank == 0)
        {
            for (std::size_t i = 0; i < u.size(); ++i)
                 os << u[i] << " ";
        }
#else
        for (std::size_t i = 0; i < _x->size(); ++i)
            os << (*_x)[i] << "\n";
#endif
    }

    ~GroundwaterFlowProcess()
    {
        for (auto p : _local_assemblers)
            delete p;

        for (auto p : _all_mesh_subsets)
            delete p;

        delete _mesh_subset_all_nodes;
    }

    /// Explicitly release memory of global vector, matrix and linear solvers
    /// for MPI based parallel computing
    void releaseEquationMemory()
    {
#ifdef USE_PETSC
        INFO("Time elapsed in the assembly using PETSc vector and matrix: %g s.\n",
             _elapsed_ctime);

        delete _A.release();
        delete _rhs.release();
        delete _x.release();
        delete _linearSolver.release();
#endif
    }

private:
    ProcessVariable const* _hydraulic_head = nullptr;

    double const _hydraulic_conductivity = 1e-6;

    MeshLib::MeshSubset const* _mesh_subset_all_nodes = nullptr;
    std::vector<MeshLib::MeshSubsets*> _all_mesh_subsets;

    GlobalSetup _global_setup;
    std::unique_ptr<typename GlobalSetup::LinearSolver> _linearSolver;
    std::unique_ptr<typename GlobalSetup::MatrixType> _A;
    std::unique_ptr<typename GlobalSetup::VectorType> _rhs;
    std::unique_ptr<typename GlobalSetup::VectorType> _x;

    using LocalAssembler = GroundwaterFlow::LocalAssemblerDataInterface<
        typename GlobalSetup::MatrixType, typename GlobalSetup::VectorType>;

    std::vector<LocalAssembler*> _local_assemblers;

    using GlobalAssembler =
        AssemblerLib::VectorMatrixAssembler<
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType>;

    std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap> _local_to_global_index_map;

    std::unique_ptr<GlobalAssembler> _global_assembler;

    /// Global ids in the global matrix/vector where the dirichlet bc is
    /// imposed and their corresponding values.
    struct DirichletBC {
        std::vector<std::size_t> global_ids;
        std::vector<double> values;
    } _dirichlet_bc;

    double _elapsed_ctime = 0.; ///< Clock time.

    MeshLib::NodeAdjacencyTable _node_adjacency_table;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
