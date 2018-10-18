/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BHEAbstract.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE  // namespace of borehole heat exchanger
{
/**
 * discharge type of the 2U BHE
 */
enum class BHE_DISCHARGE_TYPE
{
    BHE_DISCHARGE_TYPE_PARALLEL,  // parallel discharge
    BHE_DISCHARGE_TYPE_SERIAL     // serial discharge
};

class BHE_2U final : public BHEAbstract
{
public:
    /**
     * constructor
     */
    BHE_2U(
        BHE::BHE_BOUNDARY_TYPE const bound_type /* type of BHE boundary */,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            bhe_curves /* bhe related curves */,
        BoreholeGeometry const borehole_geometry = {100, 0.013},
        PipeParameters const pipe_geometry =
            {0.016 /* inner radius of the pipline */,
             0.016 /* outer radius of the pipline */,
             0.0029 /* pipe-in wall thickness*/,
             0.0029 /* pipe-out wall thickness*/,
             0.38 /* thermal conductivity of the pipe wall */,
             0.38 /* thermal conductivity of the inner pipe wall */,
             0.38 /* thermal conductivity of the outer pipe wall */},
        RefrigerantParameters const refrigerant_param =
            {
                0.00054741 /* dynamic viscosity of the refrigerant */,
                988.1 /* density of the refrigerant */,
                0.6405 /* thermal conductivity of the refrigerant */,
                4180 /* specific heat capacity of the refrigerant */, 1.0e-4 /* longitudinal dispersivity of the refrigerant in the pipeline */},
        GroutParameters const grout_param =
            {2190 /* density of the grout */, 0.5 /* porosity of the grout */,
             1000 /* specific heat capacity of the grout */,
             2.3 /* thermal conductivity of the grout */},
        ExternallyDefinedRaRb const extern_Ra_Rb =
            {false /* whether Ra and Rb values are used */,
             0.0 /* external defined borehole internal thermal resistance */,
             0.0 /* external defined borehole thermal resistance */},
        ExternallyDefinedThermalResistances const
            extern_def_thermal_resistances =
                {false /* whether user defined R values are used */,
                 0.0 /* external defined borehole thermal resistance */,
                 0.0 /* external defined borehole thermal resistance */,
                 0.0 /* external defined borehole thermal resistance */,
                 0.0 /* external defined borehole thermal resistance */,
                 0.0 /* external defined borehole thermal resistance */},
        double const Q_r = 21.86 /
                           86400 /* total refrigerant flow discharge of BHE */,
        double const omega = 0.04242 /* pipe distance */,
        double const power_in_watt = 0.0 /* injected or extracted power */,
        double const delta_T_val =
            0.0 /* Temperature difference btw inflow and outflow temperature */,
        // double my_ext_Ra = 0.0               /* external defined borehole
        // internal thermal resistance */, double my_ext_Rb = 0.0 /* external
        // defined borehole thermal resistance */, double my_ext_Rfig = 0.0 /*
        // external defined borehole thermal resistance */, double my_ext_Rfog =
        // 0.0           /* external defined borehole thermal resistance */,
        // double my_ext_Rgg1 = 0.0           /* external defined borehole
        // thermal resistance */, double my_ext_Rgg2 = 0.0           /* external
        // defined borehole thermal resistance */, double my_ext_Rgs = 0.0 /*
        // external defined borehole thermal resistance */,
        bool const if_flowrate_curve =
            false /* whether flowrate curve is used*/,
        double const threshold = 0.0 /* Threshold Q value for switching off
                                     the BHE when using Q_Curve_fixed_dT B.C.*/
        ,
        BHE_DISCHARGE_TYPE const type =
            BHE_DISCHARGE_TYPE::BHE_DISCHARGE_TYPE_PARALLEL)
        : BHEAbstract(borehole_geometry, pipe_geometry, refrigerant_param,
                      grout_param, extern_Ra_Rb, extern_def_thermal_resistances,
                      std::move(bhe_curves), bound_type,
                      false /*if_use_ext_Ra_Rb*/, false /*user_defined_R_vals*/,
                      if_flowrate_curve, Q_r, power_in_watt, delta_T_val,
                      threshold),
          _discharge_type(type),
          _omega(omega)
    {
        _u = Eigen::Vector4d::Zero();
        _Nu = Eigen::Vector4d::Zero();

        // get the corresponding curve
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>::
            const_iterator it;
        if (bound_type ==
                BHE_BOUNDARY_TYPE::POWER_IN_WATT_CURVE_FIXED_DT_BOUNDARY ||
            bound_type == BHE_BOUNDARY_TYPE::
                              POWER_IN_WATT_CURVE_FIXED_FLOW_RATE_BOUNDARY ||
            bound_type ==
                BHE_BOUNDARY_TYPE::
                    BUILDING_POWER_IN_WATT_CURVE_FIXED_FLOW_RATE_BOUNDARY)
        {
            it = bhe_curves.find("power_in_watt_curve");
            if (it == bhe_curves.end())
            {
                // curve not found, fatal error
                OGS_FATAL(
                    "Required pow_in_watt_curve cannot be found in the BHE "
                    "parameters!");
            }

            // curve successfully found
            power_in_watt_curve = it->second.get();
        }

        if (if_flowrate_curve)
        {
            use_flowrate_curve = true;

            it = bhe_curves.find("flow_rate_curve");
            if (it == bhe_curves.end())
            {
                OGS_FATAL(
                    "Required flow_rate_curve annot be found in the BHE "
                    "parameters!");
            }

            // curve successfully found
            flowrate_curve = it->second.get();
        }

        constexpr double PI = boost::math::constants::pi<double>();
        double const& D = borehole_geometry.diameter;

        // cross section area calculation
        CSA_i = PI * pipe_geometry.r_inner * pipe_geometry.r_inner;
        CSA_o = PI * pipe_geometry.r_inner * pipe_geometry.r_inner;
        CSA_g1 =
            PI * (0.0625 * D * D -
                  pipe_geometry.r_outer *
                      pipe_geometry.r_outer);  // one fourth of the crosssection
                                               // minus the crossection of pipe
        CSA_g2 =
            PI * (0.0625 * D * D -
                  pipe_geometry.r_outer *
                      pipe_geometry.r_outer);  // one fourth of the crosssection
                                               // minus the crossection of pipe

        // initialization calculation
        initialize();
    };

    static constexpr int number_of_unknowns = 8;

    std::array<double, number_of_unknowns> pipeHeatCapacities() const
    {
        double const& rho_r = refrigerant_param.rho_r;
        double const& heat_cap_r = refrigerant_param.heat_cap_r;
        double const& porosity_g = grout_param.porosity_g;
        double const& rho_g = grout_param.rho_g;
        double const& heat_cap_g = grout_param.heat_cap_g;

        return {{/*i1*/ rho_r * heat_cap_r * CSA_i,
                 /*i2*/ rho_r * heat_cap_r * CSA_i,
                 /*o1*/ rho_r * heat_cap_r * CSA_o,
                 /*o2*/ rho_r * heat_cap_r * CSA_o,
                 /*g1*/ (1.0 - porosity_g) * rho_g * heat_cap_g * CSA_g1,
                 /*g2*/ (1.0 - porosity_g) * rho_g * heat_cap_g * CSA_g1,
                 /*g3*/ (1.0 - porosity_g) * rho_g * heat_cap_g * CSA_g2,
                 /*g4*/ (1.0 - porosity_g) * rho_g * heat_cap_g * CSA_g2}};
    }

    std::array<double, number_of_unknowns> pipeHeatConductions() const
    {
        double const& lambda_r = refrigerant_param.lambda_r;
        double const& rho_r = refrigerant_param.rho_r;
        double const& heat_cap_r = refrigerant_param.heat_cap_r;
        double const& alpha_L = refrigerant_param.alpha_L;
        double const& porosity_g = grout_param.porosity_g;
        double const& lambda_g = grout_param.lambda_g;

        // Here we calculate the laplace coefficients in the governing
        // equations of BHE. These governing equations can be found in
        // 1) Diersch (2013) FEFLOW book on page 952, M.120-122, or
        // 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 10-18.
        return {{// pipe i1, Eq. 18
                 (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()) * CSA_i,
                 // pipe i2, Eq. 18
                 (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()) * CSA_i,
                 // pipe o1,  Eq. 18
                 (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()) * CSA_o,
                 // pipe o2,  Eq. 18
                 (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()) * CSA_o,
                 // pipe g1, Eq. 14
                 (1.0 - porosity_g) * lambda_g * CSA_g1,
                 // pipe g2, Eq. 15
                 (1.0 - porosity_g) * lambda_g * CSA_g1,
                 // pipe g3, Eq. 16
                 (1.0 - porosity_g) * lambda_g * CSA_g2,
                 // pipe g4, Eq. 17
                 (1.0 - porosity_g) * lambda_g * CSA_g2}};
    }

    std::array<Eigen::Vector3d, number_of_unknowns> pipeAdvectionVectors() const
    {
        double const& rho_r = refrigerant_param.rho_r;
        double const& heat_cap_r = refrigerant_param.heat_cap_r;

        return {{// pipe i1, Eq. 10
                 {0, 0, -rho_r * heat_cap_r * _u(0) * CSA_i},
                 // pipe i2, Eq. 11
                 {0, 0, -rho_r * heat_cap_r * _u(0) * CSA_i},
                 // pipe o1, Eq. 12
                 {0, 0, rho_r * heat_cap_r * _u(0) * CSA_o},
                 // pipe o2, Eq. 13
                 {0, 0, rho_r * heat_cap_r * _u(0) * CSA_o},
                 // pipe g1, Eq. 14
                 {0, 0, 0},
                 // pipe g2, Eq. 15
                 {0, 0, 0},
                 // pipe g3, Eq. 16
                 {0, 0, 0},
                 // pipe g4, Eq. 17
                 {0, 0, 0}}};
    }

    template <int NPoints, typename SingleUnknownMatrixType,
              typename RMatrixType, typename RPiSMatrixType,
              typename RSMatrixType>
    void assembleRMatrices(
        int const idx_bhe_unknowns,
        Eigen::MatrixBase<SingleUnknownMatrixType> const& matBHE_loc_R,
        Eigen::MatrixBase<RMatrixType>& R_matrix,
        Eigen::MatrixBase<RPiSMatrixType>& R_pi_s_matrix,
        Eigen::MatrixBase<RSMatrixType>& R_s_matrix) const
    {
        switch (idx_bhe_unknowns)
        {
            case 0:  // R i1 i2
                R_matrix.block(0, 4 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(4 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(NPoints, 5 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(5 * NPoints, NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(0, 0, NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_i1
                R_matrix.block(NPoints, NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_i2
                R_matrix.block(4 * NPoints,
                               4 * NPoints,
                               NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_ig
                R_matrix.block(5 * NPoints,
                               5 * NPoints,
                               NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_ig
                break;
            case 1:  // R o1 o2
                R_matrix.block(2 * NPoints, 6 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(6 * NPoints, 2 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(3 * NPoints, 7 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(7 * NPoints, 3 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(2 * NPoints,
                               2 * NPoints,
                               NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_o1
                R_matrix.block(3 * NPoints,
                               3 * NPoints,
                               NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_o2
                R_matrix.block(6 * NPoints,
                               6 * NPoints,
                               NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_og
                R_matrix.block(7 * NPoints,
                               7 * NPoints,
                               NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_og
                break;
            case 2:  // R g1
                R_matrix.block(4 * NPoints, 6 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(6 * NPoints, 4 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(4 * NPoints, 7 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(7 * NPoints, 4 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(5 * NPoints, 6 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(6 * NPoints, 5 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(5 * NPoints, 7 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(7 * NPoints, 5 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(4 * NPoints,
                               4 * NPoints,
                               NPoints,
                               NPoints) += 2.0 * matBHE_loc_R;  // K_ig
                R_matrix.block(5 * NPoints,
                               5 * NPoints,
                               NPoints,
                               NPoints) += 2.0 * matBHE_loc_R;  // K_ig
                R_matrix.block(6 * NPoints,
                               6 * NPoints,
                               NPoints,
                               NPoints) += 2.0 * matBHE_loc_R;  // K_og
                R_matrix.block(7 * NPoints,
                               7 * NPoints,
                               NPoints,
                               NPoints) += 2.0 * matBHE_loc_R;  // K_og
                break;
            case 3:  // R g2
                R_matrix.block(6 * NPoints, 7 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(7 * NPoints, 6 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(4 * NPoints,
                               4 * NPoints,
                               NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_ig
                R_matrix.block(5 * NPoints,
                               5 * NPoints,
                               NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_ig
                R_matrix.block(6 * NPoints,
                               6 * NPoints,
                               NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_og
                R_matrix.block(7 * NPoints,
                               7 * NPoints,
                               NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_og
                break;
            case 4:  // R s
                R_s_matrix += 1.0 * matBHE_loc_R;

                R_pi_s_matrix.block(4 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_pi_s_matrix.block(5 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_pi_s_matrix.block(6 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_pi_s_matrix.block(7 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(4 * NPoints, 4 * NPoints, NPoints,
                               NPoints) += matBHE_loc_R;  // K_gs
                R_matrix.block(5 * NPoints, 5 * NPoints, NPoints,
                               NPoints) += matBHE_loc_R;  // K_gs
                R_matrix.block(6 * NPoints, 6 * NPoints, NPoints,
                               NPoints) += matBHE_loc_R;  // K_gs
                R_matrix.block(7 * NPoints, 7 * NPoints, NPoints,
                               NPoints) += matBHE_loc_R;  // K_gs
                break;
        }
    }

    /**
     * return the inflow temperature based on outflow temperature and fixed
     * power.
     */
    double getTinByTout(double T_out, double current_time);

    /**
     * required by eigen library,
     * to make sure the dynamically allocated class has
     * aligned "operator new"
     */
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    static constexpr std::pair<int, int> inflow_outflow_bc_component_ids[] = {
        {0, 2}, {3, 1}};

    /// PHI_fig, _PHI_fog, _PHI_gg_1, _PHI_gg_2, _PHI_gs;
    /// Here we store the boundary heat exchange coefficients in the governing
    /// equations of BHE.
    /// These governing equations can be found in
    /// 1) Diersch (2013) FEFLOW book on page 958, M.3, or
    /// 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 90-97.
    std::array<double, number_of_unknowns> boundary_heat_exchange_coefficients;

private:
    void initialize();

    /**
     * calculate thermal resistance
     */
    void calcThermalResistances();

    /**
     * calculate heat transfer coefficient
     */
    void calcHeatTransferCoefficients();

private:
    /**
     * thermal resistances
     */
    double _R_fig, _R_fog;

    /**
     * thermal resistances due to advective flow of refrigerant in the pipes
     */
    double _R_adv_i1, _R_adv_i2, _R_adv_o1, _R_adv_o2;

    /**
     * thermal resistances due to the pipe wall material
     */
    double _R_con_a_i1, _R_con_a_i2, _R_con_a_o1, _R_con_a_o2;

    /**
     * thermal resistances due to the grout transition
     */
    double _R_con_b;

    /**
     * thermal resistances of the grout
     */
    double _R_g;

    /**
     * thermal resistances of the grout soil exchange
     */
    double _R_gs;

    /**
     * thermal resistances due to inter-grout exchange
     */
    double _R_gg_1, _R_gg_2;

    /**
     * discharge type of the 2U BHE
     */
    const BHE_DISCHARGE_TYPE _discharge_type;

    /**
     * cross section area
     */
    double CSA_i, CSA_o, CSA_g1, CSA_g2;
    /**
     * Nusselt number
     */
    Eigen::Vector4d _Nu;
    /**
     * flow velocity inside the pipeline
     */
    Eigen::Vector4d _u;

    /**
     * pipe distance
     */
    double const _omega;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
