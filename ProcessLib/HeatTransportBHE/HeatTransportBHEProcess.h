/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include "ProcessLib/Process.h"
#include "ProcessLib/HeatTransportBHE/LocalAssemblers/HeatTransportBHEProcessAssemblerInterface.h"
#include "HeatTransportBHEProcessData.h"

namespace ProcessLib
{
    namespace HeatTransportBHE
    {
        class HeatTransportBHEProcess final : public Process
        {

        public:
            HeatTransportBHEProcess(
                MeshLib::Mesh& mesh,
                std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
                std::vector<std::unique_ptr<ParameterBase>> const& parameters,
                unsigned const integration_order,
                std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&& process_variables,
                HeatTransportBHEProcessData&& process_data,
                SecondaryVariableCollection&& secondary_variables,
                NumLib::NamedFunctionCaller&& named_function_caller);

            //! \name ODESystem interface
            //! @{

            bool isLinear() const override { return false; }

            void computeSecondaryVariableConcrete(double const t,
                GlobalVector const& x) override;

        private:

            void constructDofTable() override;

            void initializeConcreteProcess(
                NumLib::LocalToGlobalIndexMap const& dof_table,
                MeshLib::Mesh const& mesh,
                unsigned const integration_order) override;

            void assembleConcreteProcess(const double t, GlobalVector const& x,
                GlobalMatrix& M, GlobalMatrix& K,
                GlobalVector& b) override;

            void assembleWithJacobianConcreteProcess(
                const double t, GlobalVector const& x, GlobalVector const& xdot,
                const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
                GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

            HeatTransportBHEProcessData _process_data;

            std::vector<std::unique_ptr<HeatTransportBHELocalAssemblerInterface>>
                _local_assemblers;

            std::vector<std::vector<MeshLib::Element*>> _vec_BHE_elements;
            std::vector<MeshLib::Element*> _vec_soil_elements;
            std::vector<std::vector<MeshLib::Element*>> _vec_BHE_soil_elements;
            std::vector<std::vector<MeshLib::Node*>> _vec_BHE_nodes;
            std::vector<std::unique_ptr<MeshLib::MeshSubset const>> _mesh_subset_BHE_nodes;
            std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_soil_nodes;
            std::vector<int> _vec_BHE_mat_IDs;
        };

    }  // namespace HeatTransportBHE
}  // namespace ProcessLib
