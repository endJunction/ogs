/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism15.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism6.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra13.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra5.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad8.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad9.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet10.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet4.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri6.h"

#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"

#include <Eigen/Dense>

namespace detail
{
    /// Forwards the Eigen::Matrix type for general N and M.
    /// There is a partial specialization for M = 1 to store the matrix in
    /// column major storage order.
    template <int N, int M>
    struct EigenMatrixType
    {
        using type = Eigen::Matrix<double, N, M, Eigen::RowMajor>;
    };

    /// Specialization for Nx1 matrices which can be stored only in column major
    /// form in Eigen-3.2.5.
    template <int N>
    struct EigenMatrixType<N, 1>
    {
        using type = Eigen::Matrix<double, N, 1, Eigen::ColMajor>;
    };

    /// Specialization for 0xM matrices. Using fixed size Eigen matrices here
    /// would lead to zero sized arrays, which cause compilation errors on
    /// some compilers.
    template <int M>
    struct EigenMatrixType<0, M>
    {
        using type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    };

    template <>
    struct EigenMatrixType<0, 1>
    {
        using type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    };

}   // detail

/// An implementation of ShapeMatrixPolicy using fixed size (compile-time) eigen
/// matrices and vectors.
template <typename ShapeFunction, unsigned GlobalDim>
struct EigenFixedShapeMatrixPolicy
{
    template <int N>
    using VectorType = typename ::detail::EigenMatrixType<N, 1>::type;

    template <int N>
    using RowVectorType = typename ::detail::EigenMatrixType<1, N>::type;

    template <int N, int M>
    using MatrixType = typename ::detail::EigenMatrixType<N, M>::type;

    using NodalMatrixType = MatrixType<ShapeFunction::NPOINTS, ShapeFunction::NPOINTS>;
    using NodalVectorType = VectorType<ShapeFunction::NPOINTS>;
    using NodalRowVectorType = RowVectorType<ShapeFunction::NPOINTS>;
    using DimNodalMatrixType = MatrixType<ShapeFunction::DIM, ShapeFunction::NPOINTS>;
    using DimMatrixType = MatrixType<ShapeFunction::DIM, ShapeFunction::DIM>;
    using GlobalDimNodalMatrixType = MatrixType<GlobalDim, ShapeFunction::NPOINTS>;
    using GlobalDimMatrixType = MatrixType<GlobalDim, GlobalDim>;
    using GlobalDimVectorType = VectorType<GlobalDim>;

    using ShapeMatrices =
        NumLib::ShapeMatrices<
            NodalRowVectorType,
            DimNodalMatrixType,
            DimMatrixType,
            GlobalDimNodalMatrixType>;
};

/// An implementation of ShapeMatrixPolicy using dynamic size eigen matrices and
/// vectors.
template <typename ShapeFunction, unsigned GlobalDim>
struct EigenDynamicShapeMatrixPolicy
{
    // Dynamic size local matrices are much slower in allocation than their
    // fixed counterparts.

    template<int>
    using VectorType =
        Eigen::Matrix<double, Eigen::Dynamic, 1>;

    template<int>
    using RowVectorType =
        Eigen::Matrix<double, 1, Eigen::Dynamic>;

    template<int, int>
    using MatrixType =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    using NodalMatrixType = MatrixType<0,0>;
    using NodalVectorType = VectorType<0>;
    using NodalRowVectorType = RowVectorType<0>;
    using DimNodalMatrixType = MatrixType<0,0>;
    using DimMatrixType = MatrixType<0,0>;
    using GlobalDimNodalMatrixType = MatrixType<0,0>;
    using GlobalDimMatrixType = MatrixType<0,0>;
    using GlobalDimVectorType = VectorType<0>;

    using ShapeMatrices =
        NumLib::ShapeMatrices<
            NodalRowVectorType,
            DimNodalMatrixType,
            DimMatrixType,
            GlobalDimNodalMatrixType>;
};

#ifdef OGS_EIGEN_DYNAMIC_SHAPE_MATRICES
template <typename ShapeFunction, unsigned GlobalDim>
using ShapeMatrixPolicyType = EigenDynamicShapeMatrixPolicy<ShapeFunction, GlobalDim>;

const unsigned OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG = 1;
#else
template <typename ShapeFunction, unsigned GlobalDim>
using ShapeMatrixPolicyType = EigenFixedShapeMatrixPolicy<ShapeFunction, GlobalDim>;

const unsigned OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG = 0;
#endif

extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeHex20  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeHex8   , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapePrism15, 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapePrism6 , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapePyra13 , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapePyra5  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad4  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad8  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad9  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTet10  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTet4   , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTri3   , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTri6   , 3>;

extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 2>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 2>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad4  , 2>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad8  , 2>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad9  , 2>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTri3   , 2>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTri6   , 2>;

extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 1>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 1>;

extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeHex20  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeHex8   , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapePrism15, 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapePrism6 , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapePyra13 , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapePyra5  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad4  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad8  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad9  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTet10  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTet4   , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri3   , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri6   , 3>;

extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 2>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 2>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad4  , 2>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad8  , 2>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad9  , 2>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri3   , 2>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri6   , 2>;

extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 1>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 1>;
//static_assert(std::is_class<ShapeMatrixPolicyType<>::value,
        //"ShapeMatrixPolicyType was not defined.");
