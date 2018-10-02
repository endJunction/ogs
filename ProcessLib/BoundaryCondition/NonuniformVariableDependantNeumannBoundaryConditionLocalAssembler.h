/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/PropertyVector.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Parameter/MeshNodeParameter.h"

#include "GenericNonuniformNaturalBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
struct NonuniformVariableDependantNeumannBoundaryConditionData
{
    /* TODO use Parameter */
    MeshLib::PropertyVector<double> const& constant;
    MeshLib::PropertyVector<double> const& prefac1;
    MeshLib::PropertyVector<double> const& prefac2;
    MeshLib::PropertyVector<double> const& prefac3;
    // Used for mapping boundary nodes to bulk nodes.
    NumLib::LocalToGlobalIndexMap const& dof_table_boundary_v2;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class NonuniformVariableDependantNeumannBoundaryConditionLocalAssembler final
    : public GenericNonuniformNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNonuniformNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;
    using NodalVectorType = typename Base::NodalVectorType;
    using NodalMatrixType = typename Base::NodalMatrixType;

public:
    /// The neumann_bc_term factor is directly integrated into the local
    /// element matrix.
    NonuniformVariableDependantNeumannBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        NonuniformVariableDependantNeumannBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          _data(data),
          _local_rhs(local_matrix_size)
    {
    }

    void assemble(std::size_t const mesh_item_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& x, GlobalMatrix& /*K*/,
                  GlobalVector& b, GlobalMatrix* /*Jac*/) override
    {
        _local_rhs.setZero();

        MeshNodeParameter<double> constant_values{"ConstantValues",
                                                  _data.constant};
        MeshNodeParameter<double> prefac1_values{"Prefac1Values",
                                                 _data.prefac1};
        MeshNodeParameter<double> prefac2_values{"Prefac2Values",
                                                 _data.prefac2};
        MeshNodeParameter<double> prefac3_values{"Prefac3Values",
                                                 _data.prefac3};
        // Get element nodes for the interpolation from nodes to
        // integration point.
        NodalVectorType const constant_node_values =
            constant_values.getNodalValuesOnElement(Base::_element, t);
        NodalVectorType const prefac1_node_values =
            prefac1_values.getNodalValuesOnElement(Base::_element, t);
        NodalVectorType const prefac2_node_values =
            prefac2_values.getNodalValuesOnElement(Base::_element, t);
        NodalVectorType const prefac3_node_values =
            prefac3_values.getNodalValuesOnElement(Base::_element, t);
        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();
        auto const indices_v1 =
            NumLib::getIndices(mesh_item_id, dof_table_boundary);
        auto const indices_v2 =
            NumLib::getIndices(mesh_item_id, _data.dof_table_boundary_v2);
        std::vector<double> const local_v1 = x.get(indices_v1);
        std::vector<double> const local_v2 = x.get(indices_v2);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& n_and_weight = Base::_ns_and_weights[ip];
            auto const& N = n_and_weight.N;
            auto const& w = n_and_weight.weight;

            double v1_int_pt = 0.0;
            double v2_int_pt = 0.0;

            NumLib::shapeFunctionInterpolate(local_v1, N, v1_int_pt);
            NumLib::shapeFunctionInterpolate(local_v2, N, v2_int_pt);
            NodalVectorType const neumann_node_values =
                constant_node_values + prefac1_node_values * v1_int_pt +
                prefac2_node_values * v2_int_pt +
                prefac3_node_values * v1_int_pt * v2_int_pt;

            _local_rhs.noalias() += N * neumann_node_values.dot(N) * w;
        }

        b.add(indices_v1, _local_rhs);
    }

private:
    NonuniformVariableDependantNeumannBoundaryConditionData const& _data;
    NodalVectorType _local_rhs;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib
