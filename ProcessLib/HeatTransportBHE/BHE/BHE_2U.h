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
class BHE_2U final : public BHEAbstract
{
public:
    /**
     * constructor
     */
    BHE_2U(
        const std::string name /* name of the BHE */,
        BHE::BHE_BOUNDARY_TYPE bound_type /* type of BHE boundary */,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            bhe_curves /* bhe related curves */,
        BoreholeGeometry borehole_geometry = {100, 0.013},
        PipeParameters pipe_geometry =
            {0.016 /* inner radius of the pipline */,
             0.016 /* outer radius of the pipline */,
             0.0029 /* pipe-in wall thickness*/,
             0.0029 /* pipe-out wall thickness*/,
             0.38 /* thermal conductivity of the pipe wall */,
             0.38 /* thermal conductivity of the inner pipe wall */,
             0.38 /* thermal conductivity of the outer pipe wall */},
        RefrigerantParameters refrigerant_param =
            {
                0.00054741 /* dynamic viscosity of the refrigerant */,
                988.1 /* density of the refrigerant */,
                0.6405 /* thermal conductivity of the refrigerant */,
                4180 /* specific heat capacity of the refrigerant */, 1.0e-4 /* longitudinal dispersivity of the refrigerant in the pipeline */},
        GroutParameters grout_param =
            {2190 /* density of the grout */, 0.5 /* porosity of the grout */,
             1000 /* specific heat capacity of the grout */,
             2.3 /* thermal conductivity of the grout */},
        ExternallyDefinedRaRb extern_Ra_Rb =
            {false /* whether Ra and Rb values are used */,
             0.0 /* external defined borehole internal thermal resistance */,
             0.0 /* external defined borehole thermal resistance */},
        ExternallyDefinedThermalResistances extern_def_thermal_resistances =
            {false /* whether user defined R values are used */,
             0.0 /* external defined borehole thermal resistance */,
             0.0 /* external defined borehole thermal resistance */,
             0.0 /* external defined borehole thermal resistance */,
             0.0 /* external defined borehole thermal resistance */,
             0.0 /* external defined borehole thermal resistance */},
        double my_Qr = 21.86 /
                       86400 /* total refrigerant flow discharge of BHE */,
        double my_omega = 0.04242 /* pipe distance */,
        double my_power_in_watt = 0.0 /* injected or extracted power */,
        double my_delta_T_val =
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
        bool if_flowrate_curve = false /* whether flowrate curve is used*/,
        double my_threshold = 0.0 /* Threshold Q value for switching off the BHE
                                     when using Q_Curve_fixed_dT B.C.*/
        ,
        BHE_DISCHARGE_TYPE type =
            BHE_DISCHARGE_TYPE::BHE_DISCHARGE_TYPE_PARALLEL)
        : BHEAbstract(name, BHE_TYPE::TYPE_2U, borehole_geometry, pipe_geometry,
                      refrigerant_param, grout_param, extern_Ra_Rb,
                      extern_def_thermal_resistances, std::move(bhe_curves),
                      bound_type, if_flowrate_curve),
          _discharge_type(type)
    {
        _u = Eigen::Vector4d::Zero();
        _Nu = Eigen::Vector4d::Zero();

        Q_r = my_Qr;

        omega = my_omega;
        power_in_watt_val = my_power_in_watt;
        delta_T_val = my_delta_T_val;
        threshold = my_threshold;
        double const& D = borehole_geometry.diameter;

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

        S_i = PI * 2.0 * pipe_geometry.r_outer;
        S_o = PI * 2.0 * pipe_geometry.r_outer;
        S_g1 = 0.5 * D;
        S_g2 = D;
        S_gs = 0.25 * PI * D;

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

    /**
     * return the number of unknowns needed for 2U BHE
     */
    std::size_t getNumUnknowns() { return 8; }

    void updateFlowRateFromCurve(double current_time)
    {
        if (use_flowrate_curve)
        {
            double Q_r_tmp(0.0);
            Q_r_tmp = flowrate_curve->getValue(current_time);
            updateFlowRate(Q_r_tmp);
        }
    };

    /**
     * calculate thermal resistance
     */
    void calcThermalResistances();

    /**
     * Nusselt number calculation
     */
    void calcNusseltNum();

    /**
     * Renolds number calculation
     */
    void calcRenoldsNum();

    /**
     * Prandtl number calculation
     */
    void calcPrandtlNum();

    /**
     * flow velocity inside the pipeline
     */
    void calcPipeFlowVelocity();

    /**
     * calculate heat transfer coefficient
     */
    void calcHeatTransferCoefficients();

    /**
     * return the coeff of mass matrix,
     * depending on the index of unknown.
     */
    double getMassCoeff(std::size_t idx_unknown);

    /**
     * return the coeff of laplace matrix,
     * depending on the index of unknown.
     */
    void getLaplaceMatrix(std::size_t idx_unknown,
                          Eigen::MatrixXd& mat_laplace);

    /**
     * return the coeff of advection matrix,
     * depending on the index of unknown.
     */
    void getAdvectionVector(std::size_t idx_unknown,
                            Eigen::VectorXd& vec_advection);

    /**
     * return the coeff of boundary heat exchange matrix,
     * depending on the index of unknown.
     */
    double getBoundaryHeatExchangeCoeff(std::size_t idx_unknown);

    /**
     * return the number of grout zones in this BHE.
     */
    std::size_t getNumGroutZones(void) { return 4; };

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
     * heat transfer coefficients
     */
    double _PHI_fig, _PHI_fog, _PHI_gg_1, _PHI_gg_2, _PHI_gs;

    /**
     * Reynolds number
     */
    double _Re;

    /**
     * Prandtl number
     */
    double _Pr;

    /**
     * Nusselt number
     */
    Eigen::Vector4d _Nu;

    /**
     * flow velocity inside the pipeline
     */
    Eigen::Vector4d _u;

    /**
     * discharge type of the 2U BHE
     */
    const BHE_DISCHARGE_TYPE _discharge_type;

    /**
     * pipe distance
     */
    double omega;

    /**
     * specific exchange surfaces S
     */
    double S_i, S_o, S_g1, S_g2, S_gs;
    /**
     * cross section area
     */
    double CSA_i, CSA_o, CSA_g1, CSA_g2;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
