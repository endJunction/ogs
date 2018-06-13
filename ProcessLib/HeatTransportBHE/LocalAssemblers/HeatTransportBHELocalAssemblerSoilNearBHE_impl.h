/**
* \copyright
* Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include "HeatTransportBHELocalAssemblerSoilNearBHE.h"

#include <valarray>
#include <vector>

#include <Eigen/Eigen>

#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MathLib/Point3d.h"

#include "MeshLib/Node.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ProcessLib/LIE/Common/LevelSetFunction.h"
#include "ProcessLib/LIE/Common/Utils.h"

// #include "IntegrationPointDataMatrix.h"
#include "SecondaryData.h"
#include "HeatTransportBHEProcessAssemblerInterface.h"

namespace ProcessLib
{
    namespace HeatTransportBHE
    {
        template <typename ShapeFunction,
            typename IntegrationMethod,
            int GlobalDim>
            HeatTransportBHELocalAssemblerSoilNearBHE<ShapeFunction,
            IntegrationMethod,
            GlobalDim>::
            HeatTransportBHELocalAssemblerSoilNearBHE(
                MeshLib::Element const& e,
                std::size_t const n_variables,
                std::size_t const /*local_matrix_size*/,
                std::vector<unsigned> const& dofIndex_to_localIndex,
                bool const is_axially_symmetric,
                unsigned const integration_order,
                HeatTransportBHEProcessData& process_data)
            : HeatTransportBHELocalAssemblerInterface(
                n_variables * ShapeFunction::NPOINTS * GlobalDim,
                dofIndex_to_localIndex),
            _process_data(process_data),
            _integration_method(integration_order),
            _element(e),
            _is_axially_symmetric(is_axially_symmetric)
        {
            unsigned const n_integration_points =
                _integration_method.getNumberOfPoints();

            _ip_data.reserve(n_integration_points);
            _secondary_data.N.resize(n_integration_points);
        }

        template<typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
            void HeatTransportBHELocalAssemblerSoilNearBHE<ShapeFunction, IntegrationMethod, GlobalDim>::
                assemble(double const t, std::vector<double> const& local_x,
                    std::vector<double>& local_M_data,
                    std::vector<double>& local_K_data,
                    std::vector<double>& local_b_data)
            {
                assert(_element.getDimension() == GlobalDim);

                auto const local_matrix_size = local_x.size();

                assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF_SOIL);

                auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
                    local_M_data, local_matrix_size, local_matrix_size);
                auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
                    local_K_data, local_matrix_size, local_matrix_size);

                unsigned const n_integration_points =
                    _integration_method.getNumberOfPoints();

                SpatialPosition pos;
                pos.setElementID(_element.getID());

                for (unsigned ip = 0; ip < n_integration_points; ip++)
                {
                    pos.setIntegrationPoint(ip);
                    auto const& sm = _shape_matrices[ip];
                    auto const& wp = _integration_method.getWeightedPoint(ip);

                    auto const k_f = _process_data.thermal_conductivity_fluid(t, pos)[0];
                    auto const k_g = _process_data.thermal_conductivity_gas(t, pos)[0];
                    auto const k_s = _process_data.thermal_conductivity_solid(t, pos)[0];

                    auto const heat_capacity_f = _process_data.heat_capacity_fluid(t, pos)[0];
                    auto const heat_capacity_g = _process_data.heat_capacity_gas(t, pos)[0];
                    auto const heat_capacity_s = _process_data.heat_capacity_solid(t, pos)[0];

                    auto const density_f = _process_data.density_fluid(t, pos)[0];
                    auto const density_g = _process_data.density_gas(t, pos)[0];
                    auto const density_s = _process_data.density_solid(t, pos)[0];

                    // for now only using the solid phase parameters

                    // assemble Conductance matrix
                    local_K.noalias() += sm.dNdx.transpose() * k_s * sm.dNdx * sm.detJ *
                        wp.getWeight() * sm.integralMeasure;

                    // assemble Mass matrix
                    local_M.noalias() += sm.N.transpose() * density_s * heat_capacity_s *
                        sm.N * sm.detJ * wp.getWeight() *
                        sm.integralMeasure;
                }

                // debugging
                // std::string sep = "\n----------------------------------------\n";
                // Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
                // std::cout << local_K.format(CleanFmt) << sep;
                // std::cout << local_M.format(CleanFmt) << sep;
            }

        template <typename ShapeFunction,
            typename IntegrationMethod,
            int GlobalDim>
            void HeatTransportBHELocalAssemblerSoilNearBHE<ShapeFunction,
            IntegrationMethod,
            GlobalDim>::
            postTimestepConcrete(std::vector<double> const& /*local_x*/)
        {
            /*
            // Compute average value per element
            const int n = DisplacementDim == 2 ? 4 : 6;
            Eigen::VectorXd ele_stress = Eigen::VectorXd::Zero(n);
            Eigen::VectorXd ele_strain = Eigen::VectorXd::Zero(n);

            unsigned const n_integration_points =
                _integration_method.getNumberOfPoints();
            for (unsigned ip = 0; ip < n_integration_points; ip++)
            {
                auto& ip_data = _ip_data[ip];

                ele_stress += ip_data._sigma;
                ele_strain += ip_data._eps;
            }
            ele_stress /= n_integration_points;
            ele_strain /= n_integration_points;

            (*_process_data._mesh_prop_stress_xx)[_element.getID()] = ele_stress[0];
            (*_process_data._mesh_prop_stress_yy)[_element.getID()] = ele_stress[1];
            (*_process_data._mesh_prop_stress_zz)[_element.getID()] = ele_stress[2];
            (*_process_data._mesh_prop_stress_xy)[_element.getID()] = ele_stress[3];
            if (DisplacementDim == 3)
            {
                (*_process_data._mesh_prop_stress_yz)[_element.getID()] = ele_stress[4];
                (*_process_data._mesh_prop_stress_xz)[_element.getID()] = ele_stress[5];
            }

            (*_process_data._mesh_prop_strain_xx)[_element.getID()] = ele_strain[0];
            (*_process_data._mesh_prop_strain_yy)[_element.getID()] = ele_strain[1];
            (*_process_data._mesh_prop_strain_zz)[_element.getID()] = ele_strain[2];
            (*_process_data._mesh_prop_strain_xy)[_element.getID()] = ele_strain[3];
            if (DisplacementDim == 3)
            {
                (*_process_data._mesh_prop_strain_yz)[_element.getID()] = ele_strain[4];
                (*_process_data._mesh_prop_strain_xz)[_element.getID()] = ele_strain[5];
            }
            */
        }

    }  // namespace HeatTransportBHE
}  // namespace ProcessLib
