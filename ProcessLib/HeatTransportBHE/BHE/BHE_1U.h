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
class BHE_1U final : public BHEAbstract
{
public:
    /**
     * constructor
     */
    BHE_1U(
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
        double const Q_r = std::numeric_limits<double>::quiet_NaN()
        /* total refrigerant flow discharge of BHE; TODO
         * (haibing) Remove the default value */
        ,
        double const omega = 0.06 /* pipe distance */,
        double const power_in_watt = 0.0 /* injected or extracted power */,
        double const delta_T_val =
            0.0 /* Temperature difference btw inflow and outflow temperature */,

        bool if_flowrate_curve = false /* whether flowrate curve is used*/,
        double const threshold = 0.1) /* Threshold Q value for switching off the
                                      BHE when using Q_Curve_fixed_dT B.C.*/
    : BHEAbstract(borehole_geometry, pipe_geometry, refrigerant_param,
                  grout_param, extern_Ra_Rb, extern_def_thermal_resistances,
                  std::move(bhe_curves), bound_type, false /*if_use_ext_Ra_Rb*/,
                  false /*user_defined_R_vals*/, if_flowrate_curve, Q_r,
                  power_in_watt, delta_T_val, threshold),
      _omega(omega)
    {
        _u = Eigen::Vector2d::Zero();
        _Nu = Eigen::Vector2d::Zero();
        // 1U type of BHE has 2 pipes

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
        if (bound_type == BHE_BOUNDARY_TYPE::FIXED_INFLOW_TEMP_CURVE_BOUNDARY)
        {
            it = bhe_curves.find("inflow_temp_curve");
            if (it == bhe_curves.end())
            {
                OGS_FATAL(
                    "Required inflow_temp_curve annot be found in the BHE "
                    "parameters!");
            }

            // curve successfully found
            inflow_temperature_curve = it->second.get();
        }

        updateHeatTransferCoefficients(Q_r);
    };

    static constexpr int number_of_unknowns = 4;

    std::array<double, number_of_unknowns> pipeHeatCapacities() const
    {
        double const& rho_r = refrigerant_param.rho_r;
        double const& heat_cap_r = refrigerant_param.heat_cap_r;
        double const& rho_g = grout_param.rho_g;
        double const& porosity_g = grout_param.porosity_g;
        double const& heat_cap_g = grout_param.heat_cap_g;

        return {{/*i1*/ rho_r * heat_cap_r,
                 /*o1*/ rho_r * heat_cap_r,
                 /*g1*/ (1.0 - porosity_g) * rho_g * heat_cap_g,
                 /*g2*/ (1.0 - porosity_g) * rho_g * heat_cap_g}};
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
        // 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 19-22.
        return {{// pipe i1, Eq. 19
                 (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()),
                 // pipe o1, Eq. 20
                 (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()),
                 // pipe g1, Eq. 21
                 (1.0 - porosity_g) * lambda_g,
                 // pipe g2, Eq. 22
                 (1.0 - porosity_g) * lambda_g}};
    }

    std::array<Eigen::Vector3d, number_of_unknowns> pipeAdvectionVectors() const
    {
        double const& rho_r = refrigerant_param.rho_r;
        double const& heat_cap_r = refrigerant_param.heat_cap_r;

        return {{// pipe i1, Eq. 19
                 {0, 0, -rho_r * heat_cap_r * _u(0)},
                 // pipe o1, Eq. 20
                 {0, 0, rho_r * heat_cap_r * _u(0)},
                 // grout g1, Eq. 21
                 {0, 0, 0},
                 // grout g2, Eq. 22
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
            case 0:  // PHI_fig
                R_matrix.block(0, 2 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(2 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(0, 0, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_i1
                R_matrix.block(2 * NPoints, 2 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_ig
                break;
            case 1:  // PHI_fog
                R_matrix.block(NPoints, 3 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(3 * NPoints, NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(NPoints, NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_o1
                R_matrix.block(3 * NPoints, 3 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_og
                break;
            case 2:  // PHI_gg
                R_matrix.block(2 * NPoints, 3 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(3 * NPoints, 2 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(2 * NPoints, 2 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_ig  // notice we only have
                                         // 1 PHI_gg term here.
                R_matrix.block(3 * NPoints, 3 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_og  // see Diersch 2013 FEFLOW
                                         // book page 954 Table M.2
                break;
            case 3:  // PHI_gs
                R_s_matrix.template block<NPoints, NPoints>(0, 0).noalias() +=
                    1.0 * matBHE_loc_R;

                R_pi_s_matrix.block(2 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_pi_s_matrix.block(3 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(2 * NPoints, 2 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_ig
                R_matrix.block(3 * NPoints, 3 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_og
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
        {0, 1}};

    /// PHI_fig, PHI_fog, PHI_gg, PHI_gs;
    /// Here we store the boundary heat exchange coefficients in the governing
    /// equations of BHE.
    /// These governing equations can be found in
    /// 1) Diersch (2013) FEFLOW book on page 958, M.3, or
    /// 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 90-97.
    std::array<double, number_of_unknowns> boundary_heat_exchange_coefficients;

    static constexpr double pi = boost::math::constants::pi<double>();
    std::array<double, number_of_unknowns> const cross_section_areas = {
        {pi * pipe_param.r_inner * pipe_param.r_inner,
         pi* pipe_param.r_inner* pipe_param.r_inner,
         pi*(0.125 * borehole_geometry.diameter * borehole_geometry.diameter -
             pipe_param.r_outer * pipe_param.r_outer),
         pi*(0.125 * borehole_geometry.diameter * borehole_geometry.diameter -
             pipe_param.r_outer * pipe_param.r_outer)}};

private:
    // TODO (haibing) remove after removing the global Q_r variable. The new
    // corresponding function is updateHeatTransferCoefficients takes the new
    // flow rate as an argument.
    void initialize() {}

    void updateHeatTransferCoefficients(double const flow_rate);

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
     * heating COP curve
     */
    MathLib::PiecewiseLinearInterpolation* const heating_cop_curve = nullptr;

    /**
     * cooling COP curve
     */
    MathLib::PiecewiseLinearInterpolation* const cooling_cop_curve = nullptr;

    /**
     * thermal resistances
     */
    double _R_fig, _R_fog;

    /**
     * thermal resistances due to advective flow of refrigerant in the pipes
     */
    double _R_adv_i1, _R_adv_o1;

    /**
     * thermal resistances due to the pipe wall material
     */
    double _R_con_a_i1, _R_con_a_o1;

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
    double _R_gg;

    /**
     * Nusselt number
     */
    Eigen::Vector2d _Nu;
    /**
     * flow velocity inside the pipeline
     */
    Eigen::Vector2d _u;

    /**
     * pipe distance
     */
    double const _omega;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
