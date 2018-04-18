/**
* \copyright
* Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include <vector>

#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "ProcessLib/HeatTransportBHE/HeatTransportBHEProcessData.h"

#include "IntegrationPointDataBHE.h"
#include "SecondaryData.h"
#include "HeatTransportBHEProcessAssemblerInterface.h"

namespace ProcessLib
{
    namespace HeatTransportBHE
    {
        template <typename ShapeFunction, typename IntegrationMethod,
            int BHE_Dim>
        class HeatTransportBHELocalAssemblerBHE
            : public HeatTransportBHELocalAssemblerInterface
        {
        public:
            using ShapeMatricesType =
                ShapeMatrixPolicyType<ShapeFunction, BHE_Dim>;
            using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
            using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
            using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

            HeatTransportBHELocalAssemblerBHE(
                HeatTransportBHELocalAssemblerBHE const&) = delete;
            HeatTransportBHELocalAssemblerBHE(
                HeatTransportBHELocalAssemblerBHE&&) = delete;

            HeatTransportBHELocalAssemblerBHE(
                MeshLib::Element const& e,
                std::size_t const local_matrix_size,
                std::vector<unsigned> const& dofIndex_to_localIndex,
                bool const is_axially_symmetric,
                unsigned const integration_order,
                HeatTransportBHEProcessData& process_data);

            void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                std::vector<double>& /*local_M_data*/,
                std::vector<double>& /*local_K_data*/,
                std::vector<double>& /*local_b_data*/) override
            {
                OGS_FATAL(
                    "SmallDeformationLocalAssembler: assembly without jacobian is not "
                    "implemented.");
            }

            void assembleWithJacobian(double const t,
                Eigen::VectorXd const& local_u,
                Eigen::VectorXd& local_b,
                Eigen::MatrixXd& local_J) override;

            void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                double const /*t*/,
                double const /*delta_t*/) override
            {
                unsigned const n_integration_points =
                    _integration_method.getNumberOfPoints();
                /*
                for (unsigned ip = 0; ip < n_integration_points; ip++)
                {
                    _ip_data[ip].pushBackState();
                }
                */
            }

            void postTimestepConcrete(std::vector<double> const& /*local_x*/) override;

            Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
                const unsigned integration_point) const override
            {
                auto const& N = _secondary_data.N[integration_point];

                // assumes N is stored contiguously in memory
                return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
            }

        private:

            HeatTransportBHEProcessData& _process_data;
            // FractureProperty const* _fracture_property = nullptr;

            std::vector<IntegrationPointDataBHE,
                Eigen::aligned_allocator<IntegrationPointDataBHE>>
                _ip_data;

            IntegrationMethod _integration_method;
            std::vector<ShapeMatrices, Eigen::aligned_allocator<
                typename ShapeMatricesType::ShapeMatrices>>
                _shape_matrices;
            MeshLib::Element const& _element;
            SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
        };
    }  // namespace HeatTransportBHE
}  // namespace ProcessLib

#include "HeatTransportBHELocalAssemblerBHE_impl.h"
