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

#include "HeatTransportBHEProcessAssemblerInterface.h"
#include "IntegrationPointDataBHE.h"
#include "SecondaryData.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
template <typename ShapeFunction, typename IntegrationMethod, typename BHEType>
class HeatTransportBHELocalAssemblerBHE
    : public HeatTransportBHELocalAssemblerInterface
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, 3 /* GlobalDim */>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    // Using dynamic size, because the number of unknowns in the BHE is runtime
    // value.
    using BheLocalMatrixType =
        typename ShapeMatricesType::template MatrixType<Eigen::Dynamic,
                                                        Eigen::Dynamic>;
    HeatTransportBHELocalAssemblerBHE(
        HeatTransportBHELocalAssemblerBHE const&) = delete;
    HeatTransportBHELocalAssemblerBHE(HeatTransportBHELocalAssemblerBHE&&) =
        delete;

    HeatTransportBHELocalAssemblerBHE(
        MeshLib::Element const& e,
        BHE::BHEAbstract const& bhe,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HeatTransportBHEProcessData& process_data);

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_b_data*/) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

private:
    HeatTransportBHEProcessData& _process_data;

    std::vector<IntegrationPointDataBHE<ShapeMatricesType, BHEType>,
                Eigen::aligned_allocator<
                    IntegrationPointDataBHE<ShapeMatricesType, BHEType>>>
        _ip_data;

    IntegrationMethod _integration_method;

    BHEType const& _bhe;

    std::size_t const element_id;

    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;

    BheLocalMatrixType _R_matrix;

    BheLocalMatrixType _R_s_matrix;

    BheLocalMatrixType _R_pi_s_matrix;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib

#include "HeatTransportBHELocalAssemblerBHE-impl.h"
