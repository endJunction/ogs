/**
* \copyright
* Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/MechanicsBase.h"

namespace ProcessLib
{
    namespace HeatTransportBHE
    {
        template <typename ShapeMatricesType, typename BMatricesType,
            int GlobalDim>
        struct IntegrationPointDataSoil final
        {
            explicit IntegrationPointDataSoil(
                MaterialLib::Solids::MechanicsBase<GlobalDim>& solid_material)
                : _solid_material(solid_material)
            {
            }

            MaterialLib::Solids::MechanicsBase<GlobalDim>& _solid_material;
            
            double integration_weight;

            typename ShapeMatricesType::NodalRowVectorType N;
            typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx;

            void pushBackState()
            {
            
            }

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        };

    }  // namespace HeatTransportBHE
}  // namespace ProcessLib
