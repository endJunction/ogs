/**
* \copyright
* Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include <utility>
#include <vector>

#include "BaseLib/Error.h"

#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "ProcessLib/LocalAssemblerInterface.h"

namespace ProcessLib
{
    namespace HeatTransportBHE
    {
        class HeatTransportBHELocalAssemblerInterface
            : public ProcessLib::LocalAssemblerInterface,
            public NumLib::ExtrapolatableElement
        {
        public:
            HeatTransportBHELocalAssemblerInterface() : _dofIndex_to_localIndex{} {}
            HeatTransportBHELocalAssemblerInterface(
                std::size_t n_local_size, std::vector<unsigned> dofIndex_to_localIndex)
                : _dofIndex_to_localIndex(std::move(dofIndex_to_localIndex))
            {
                _local_T.resize(n_local_size);
                _local_b.resize(_local_T.size());
                _local_A.resize(_local_T.size(), _local_T.size());
            }

            void assembleWithJacobian(double const t,
                std::vector<double> const& local_x_,
                std::vector<double> const& /*local_xdot*/,
                const double /*dxdot_dx*/, const double /*dx_dx*/,
                std::vector<double>& /*local_M_data*/,
                std::vector<double>& /*local_K_data*/,
                std::vector<double>& local_b_data,
                std::vector<double>& local_Jac_data) override
            {
                /*
                auto const local_dof_size = local_x_.size();

                _local_u.setZero();
                for (unsigned i = 0; i < local_dof_size; i++)
                {
                    _local_u[_dofIndex_to_localIndex[i]] = local_x_[i];
                }
                _local_b.setZero();
                _local_J.setZero();

                assembleWithJacobian(t, _local_u, _local_b, _local_J);

                local_b_data.resize(local_dof_size);
                for (unsigned i = 0; i < local_dof_size; i++)
                {
                    local_b_data[i] = _local_b[_dofIndex_to_localIndex[i]];
                }

                local_Jac_data.resize(local_dof_size * local_dof_size);
                for (unsigned i = 0; i < local_dof_size; i++)
                {
                    for (unsigned j = 0; j < local_dof_size; j++)
                    {
                        local_Jac_data[i * local_dof_size + j] = _local_J(
                            _dofIndex_to_localIndex[i], _dofIndex_to_localIndex[j]);
                    }
                }
                */
            }

            virtual void assembleWithJacobian(double const t,
                                              Eigen::VectorXd const& local_T,
                                              Eigen::VectorXd& local_b,
                                              Eigen::MatrixXd& local_A)
            {
                (void)t;
                (void)local_T;
                (void)local_b;
                (void)local_A;
                OGS_FATAL(
                    "HeatTransportBHELocalAssemblerInterface::"
                    "assembleWithJacobian() "
                    "is not implemented");
            }

        private:
            Eigen::VectorXd _local_T;
            Eigen::VectorXd _local_b;
            Eigen::MatrixXd _local_A;

            std::vector<unsigned> const _dofIndex_to_localIndex;
        };

    }  // namespace HeatTransportBHE
}  // namespace ProcessLib
