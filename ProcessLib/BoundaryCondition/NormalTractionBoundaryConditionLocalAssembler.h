/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <iostream>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/FiniteElement/ComputeNormalVector.h"
#include "ProcessLib/Parameter/Parameter.h"

#include "GenericNaturalBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
namespace NormalTractionBoundaryCondition
{

template <typename ShapeMatricesTypeDisplacement, int GlobalDim, int NPoints>
struct IntegrationPointData final
{
    IntegrationPointData(
        typename ShapeMatricesTypeDisplacement::template VectorType<
            NPoints * GlobalDim> const& Nu_times_n_,
        double const integration_weight_)
        : Nu_times_n(Nu_times_n_), integration_weight(integration_weight_)
    {
    }

    // Shape matrix (for displacement) times element's normal vector.
    typename ShapeMatricesTypeDisplacement::template VectorType<
        NPoints * GlobalDim> const Nu_times_n;

    double const integration_weight;
};

class NormalTractionBoundaryConditionLocalAssemblerInterface
{
public:
    virtual void assemble(
        std::size_t const id,
        NumLib::LocalToGlobalIndexMap const& dof_table_boundary, double const t,
        const GlobalVector& /*x*/, GlobalMatrix& /*K*/, GlobalVector& b) = 0;
    virtual ~NormalTractionBoundaryConditionLocalAssemblerInterface() = default;
};

template <typename ShapeFunctionDisplacement, typename IntegrationMethod,
          unsigned GlobalDim>
class NormalTractionBoundaryConditionLocalAssembler final
    : public NormalTractionBoundaryConditionLocalAssemblerInterface
{
public:
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, GlobalDim>;
    using GlobalDimVectorType =
        typename ShapeMatricesTypeDisplacement::GlobalDimVectorType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesTypeDisplacement::GlobalDimNodalMatrixType;

    NormalTractionBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool is_axially_symmetric,
        unsigned const integration_order,
        Parameter<double> const& pressure)
        : _integration_method(integration_order), _pressure(pressure)
    {
        _local_rhs.setZero(local_matrix_size);
        GlobalDimNodalMatrixType const coordinates =
            NumLib::createElementPointsCoordinatesMatrix<
                GlobalDimNodalMatrixType>(e);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        auto const shape_matrices_u =
            initShapeMatrices<ShapeFunctionDisplacement,
                              ShapeMatricesTypeDisplacement, IntegrationMethod,
                              GlobalDim>(e, is_axially_symmetric,
                                         _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            // TODO Extend to rotated 2d meshes and line elements.
            GlobalDimVectorType element_normal =
                NumLib::ComputeNormalVector<GlobalDim>::template calculate<
                    GlobalDimVectorType>(shape_matrices_u[ip].dNdr,
                                         coordinates);
            std::cout << ip << ": normal = \n" << element_normal << "\n";

            typename ShapeMatricesTypeDisplacement::template MatrixType<
                GlobalDim, displacement_size>
                N_u = ShapeMatricesTypeDisplacement::template MatrixType<
                    GlobalDim, displacement_size>::Zero(GlobalDim,
                                                        displacement_size);
            for (int i = 0; i < static_cast<int>(GlobalDim); ++i)
                N_u.template block<1, displacement_size / GlobalDim>(
                       i, i * displacement_size / GlobalDim)
                    .noalias() = shape_matrices_u[ip].N;

            double const integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                shape_matrices_u[ip].integralMeasure *
                shape_matrices_u[ip].detJ;

            _ip_data.emplace_back(N_u.transpose() * element_normal,
                                  integration_weight);
        }
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& /*x*/,
                  GlobalMatrix& /*K*/, GlobalVector& local_rhs)
    {
        _local_rhs.setZero();

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition pos;
        pos.setElementID(id);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);

            auto const& w = _ip_data[ip].integration_weight;
            auto const& Nu_times_n = _ip_data[ip].Nu_times_n;

            _local_rhs.noalias() -=
                Nu_times_n.transpose() * _pressure(t, pos)[0] * w;
        }

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        local_rhs.add(indices, _local_rhs);
    }

private:
    IntegrationMethod const _integration_method;
    Parameter<double> const& _pressure;

    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * GlobalDim;
    std::vector<IntegrationPointData<ShapeMatricesTypeDisplacement, GlobalDim,
                                     ShapeFunctionDisplacement::NPOINTS>,
                Eigen::aligned_allocator<IntegrationPointData<
                    ShapeMatricesTypeDisplacement, GlobalDim,
                    ShapeFunctionDisplacement::NPOINTS>>>
        _ip_data;

    typename ShapeMatricesTypeDisplacement::template VectorType<
        displacement_size>
        _local_rhs;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace NormalTractionBoundaryCondition
}  // namespace ProcessLib
