/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CentralDifferencesJacobianAssembler.h"
#include "BaseLib/Error.h"
#include "LocalAssemblerInterface.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace ProcessLib
{
CentralDifferencesJacobianAssembler::CentralDifferencesJacobianAssembler(
    std::vector<double>&& absolute_epsilons)
    : _absolute_epsilons(std::move(absolute_epsilons))
{
    if (_absolute_epsilons.empty())
        OGS_FATAL("No values for the absolute epsilons have been given.");
}

void CentralDifferencesJacobianAssembler::assembleWithJacobian(
    LocalAssemblerInterface& local_assembler, const double t,
    const std::vector<double>& local_x_data,
    const std::vector<double>& local_xdot_data, const double dxdot_dx,
    const double dx_dx, std::vector<double>& local_M_data,
    std::vector<double>& local_K_data, std::vector<double>& local_b_data,
    std::vector<double>& local_Jac_data)
{
    // TODO do not check in every call.
    if (local_x_data.size() % _absolute_epsilons.size() != 0)
    {
        OGS_FATAL(
            "The number of specified epsilons (%u) and the number of local "
            "d.o.f.s (%u) do not match, i.e., the latter is not divisable by "
            "the former.",
            _absolute_epsilons.size(), local_x_data.size());
    }

    auto const num_r_c =
        static_cast<Eigen::MatrixXd::Index>(local_x_data.size());

    auto const local_x =
        MathLib::toVector<Eigen::VectorXd>(local_x_data, num_r_c);
    auto const local_xdot =
        MathLib::toVector<Eigen::VectorXd>(local_xdot_data, num_r_c);

    auto local_Jac =
        MathLib::createZeroedMatrix(local_Jac_data, num_r_c, num_r_c);

    std::vector<double> local_M_data_;
    std::vector<double> local_K_data_;
    std::vector<double> local_b_data_;
    std::vector<double> local_x_perturbed_data_(local_x_data);


//
//    std::vector<double> M_vect;
//    std::vector<double> M_vect_p;
//    std::vector<double> M_vect_m;
//
//    std::vector<double> K_vect;
//    std::vector<double> K_vect_p;
//    std::vector<double> K_vect_m;
//
//    std::vector<double> b_vect;
//    std::vector<double> b_vect_p;
//    std::vector<double> b_vect_m;

//    local_x_perturbed_data_[9] += 0;
//    local_assembler.assemble(t, local_x_perturbed_data_, M_vect, K_vect, b_vect);
//
//    auto const local_M = MathLib::toMatrix(M_vect, num_r_c, num_r_c);
//    auto const local_K = MathLib::toMatrix(K_vect, num_r_c, num_r_c);
//    auto const local_b = MathLib::toVector<Eigen::VectorXd>(b_vect, num_r_c);
//    auto const local_x_pert = MathLib::toVector<Eigen::VectorXd>(local_x_perturbed_data_, num_r_c);
//    std::cout << "\n=======================\n";
//    std::cout << " local_x:\n" << local_x_pert << "\n";
//    std::cout << "\n=======================\n";
//    std::cout << " local_M_data: \n" << local_M << "\n";
//    std::cout << " local_K_data: \n" << local_K << "\n";
//    std::cout << " local_b_data: \n" << local_b << "\n";
//    std::cout << "\n=======================\n";
//
//    local_x_perturbed_data_[9] += 1e-3;
//    local_assembler.assemble(t, local_x_perturbed_data_, M_vect_p, K_vect_p, b_vect_p);
//
//    auto const local_M_p = MathLib::toMatrix(M_vect_p, num_r_c, num_r_c);
//    auto const local_K_p = MathLib::toMatrix(K_vect_p, num_r_c, num_r_c);
//    auto const local_b_p = MathLib::toVector<Eigen::VectorXd>(b_vect_p, num_r_c);
//    auto const local_x_pert_p = MathLib::toVector<Eigen::VectorXd>(local_x_perturbed_data_, num_r_c);
//
//    std::cout << "\n Plus EPS\n";
//    std::cout << " local_x:\n" << local_x_pert_p << "\n";
//    std::cout << "\n=======================\n";
//    std::cout << " local_M_data: \n" << local_M_p << "\n";
//    std::cout << " local_K_data: \n" << local_K_p << "\n";
//    std::cout << " local_b_data: \n" << local_b_p << "\n";
//    std::cout << "\n=======================\n";
//    local_x_perturbed_data_[9] -= 2*1e-3;
//    local_assembler.assemble(t, local_x_perturbed_data_, M_vect_m, K_vect_m, b_vect_m);
//
//    auto const local_M_m = MathLib::toMatrix(M_vect_m, num_r_c, num_r_c);
//    auto const local_K_m = MathLib::toMatrix(K_vect_m, num_r_c, num_r_c);
//    auto const local_b_m = MathLib::toVector<Eigen::VectorXd>(b_vect_m, num_r_c);
//    auto const local_x_pert_m = MathLib::toVector<Eigen::VectorXd>(local_x_perturbed_data_, num_r_c);
//
//    std::cout << "\n Minus EPS\n";
//    std::cout << " local_x:\n" << local_x_pert_m << "\n";
//    std::cout << "\n=======================\n";
//    std::cout << " local_M_data: \n" << local_M_m << "\n";
//    std::cout << " local_K_data: \n" << local_K_m << "\n";
//    std::cout << " local_b_data: \n" << local_b_m << "\n";
//    std::cout << "\n=======================\n";
//    OGS_FATAL("XXX");

    auto const num_dofs_per_component =
        local_x_data.size() / _absolute_epsilons.size();

    // Residual  res := M xdot + K x - b
    // Computing Jac := dres/dx
    //                = M dxdot/dx + dM/dx xdot + K dx/dx + dK/dx x - db/dx
    //                  (Note: dM/dx and dK/dx actually have the second and
    //                  third index transposed.)
    // The loop computes the dM/dx, dK/dx and db/dx terms, the rest is computed
    // afterwards.

    for (Eigen::MatrixXd::Index i = 0; i < num_r_c; ++i)
    {
        // assume that local_x_data is ordered by component.
        auto const component = i / num_dofs_per_component;
        auto const eps = _absolute_epsilons[component];

        local_x_perturbed_data_[i] += eps;
        local_assembler.assemble(t, local_x_perturbed_data_, local_M_data,
                                 local_K_data, local_b_data);

        local_x_perturbed_data_[i] = local_x_data[i] - eps;
        local_assembler.assemble(t, local_x_perturbed_data_, local_M_data_,
                                 local_K_data_, local_b_data_);

        local_x_perturbed_data_[i] = local_x_data[i];

        if (!local_M_data.empty())
        {
            auto const local_M_p =
                MathLib::toMatrix(local_M_data, num_r_c, num_r_c);
            auto const local_M_m =
                MathLib::toMatrix(local_M_data_, num_r_c, num_r_c);
            local_Jac.col(i).noalias() +=
                // dM/dxi * x_dot
                (local_M_p - local_M_m) * local_xdot / (2.0 * eps);
            local_M_data.clear();
            local_M_data_.clear();
        }
        if (!local_K_data.empty())
        {
            auto const local_K_p =
                MathLib::toMatrix(local_K_data, num_r_c, num_r_c);
            auto const local_K_m =
                MathLib::toMatrix(local_K_data_, num_r_c, num_r_c);
            local_Jac.col(i).noalias() +=
                // dK/dxi * x
                (local_K_p - local_K_m) * local_x / (2.0 * eps);
            local_K_data.clear();
            local_K_data_.clear();
        }
        if (!local_b_data.empty())
        {
            auto const local_b_p =
                MathLib::toVector<Eigen::VectorXd>(local_b_data, num_r_c);
            auto const local_b_m =
                MathLib::toVector<Eigen::VectorXd>(local_b_data_, num_r_c);
            local_Jac.col(i).noalias() -=
                // db/dxi
                (local_b_p - local_b_m) / (2.0 * eps);
            local_b_data.clear();
            local_b_data_.clear();
        }
    }

    // Assemble with unperturbed local x.
    local_assembler.assemble(t, local_x_data, local_M_data, local_K_data,
                             local_b_data);

    // Compute remaining terms of the Jacobian.
    if (dxdot_dx != 0.0 && !local_M_data.empty())
    {
        auto local_M = MathLib::toMatrix(local_M_data, num_r_c, num_r_c);
        local_Jac.noalias() += local_M * dxdot_dx;
    }
    if (dx_dx != 0.0 && !local_K_data.empty())
    {
        auto local_K = MathLib::toMatrix(local_K_data, num_r_c, num_r_c);
        local_Jac.noalias() += local_K * dx_dx;
    }
}

std::unique_ptr<CentralDifferencesJacobianAssembler>
createCentralDifferencesJacobianAssembler(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__jacobian_assembler__type}
    config.checkConfigParameter("type", "CentralDifferences");

    // TODO make non-optional.
    //! \ogs_file_param{prj__processes__process__jacobian_assembler__CentralDifferences__relative_epsilons}
    auto rel_eps = config.getConfigParameterOptional<std::vector<double>>(
        "relative_epsilons");
    //! \ogs_file_param{prj__processes__process__jacobian_assembler__CentralDifferences__component_magnitudes}
    auto comp_mag = config.getConfigParameterOptional<std::vector<double>>(
        "component_magnitudes");

    if (!!rel_eps != !!comp_mag)
    {
        OGS_FATAL(
            "Either both or none of <relative_epsilons> and "
            "<component_magnitudes> have to be specified.");
    }

    std::vector<double> abs_eps;

    if (rel_eps)
    {
        if (rel_eps->size() != comp_mag->size())
        {
            OGS_FATAL(
                "The numbers of components of  <relative_epsilons> and "
                "<component_magnitudes> do not match.");
        }

        abs_eps.resize(rel_eps->size());
        for (std::size_t i = 0; i < rel_eps->size(); ++i)
        {
            abs_eps[i] = (*rel_eps)[i] * (*comp_mag)[i];
        }
    }
    else
    {
        // By default 1e-8 is used as epsilon for all components.
        // TODO: remove this default value.
        abs_eps.emplace_back(1e-8);
    }

    return std::make_unique<CentralDifferencesJacobianAssembler>(
        std::move(abs_eps));
}

}  // ProcessLib
