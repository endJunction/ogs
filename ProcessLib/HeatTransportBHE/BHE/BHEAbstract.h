/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

/**
* \file BHEAbstract.h
* 2014/06/04 HS inital implementation
* borehole heat exchanger abstract class
*
* 1) Diersch_2011_CG
* Two very important references to understand this class implementations are: 
* H.-J.G. Diersch, D. Bauer, W. Heidemann, W. R黨aak, P. Sch鋞zl, 
* Finite element modeling of borehole heat exchanger systems: 
* Part 1. Fundamentals, Computers & Geosciences, 
* Volume 37, Issue 8, August 2011, Pages 1122-1135, ISSN 0098-3004, 
* http://dx.doi.org/10.1016/j.cageo.2010.08.003.
*
* 2) FEFLOW_2014_Springer
* FEFLOW: Finite Element Modeling of Flow, Mass and Heat Transport in Porous and Fractured Media
* Diersch, Hans-Joerg, 2014, XXXV, 996 p, Springer. 
* 
*/

#pragma once

#include <iostream>
#include "boost/math/constants/constants.hpp"
#include "BHE_Net_ELE_Abstract.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "GeoLib/Polyline.h"

namespace ProcessLib
{
    namespace HeatTransportBHE
    {
        namespace BHE  // namespace of borehole heat exchanger
        {
            /**
              * list the types of borehole heat exchanger
              */
            enum class BHE_TYPE {
                TYPE_2U,   // two u-tube borehole heat exchanger
                TYPE_1U,   // one u-tube borehole heat exchanger
                TYPE_CXC,  // coaxial pipe with annualar inlet
                TYPE_CXA      // coaxial pipe with centreed inlet
            };

            enum class BHE_BOUNDARY_TYPE {
                FIXED_INFLOW_TEMP_BOUNDARY,
                FIXED_INFLOW_TEMP_CURVE_BOUNDARY,
                POWER_IN_WATT_BOUNDARY,
                POWER_IN_WATT_CURVE_FIXED_DT_BOUNDARY,
                BUILDING_POWER_IN_WATT_CURVE_FIXED_DT_BOUNDARY,
                BUILDING_POWER_IN_WATT_CURVE_FIXED_FLOW_RATE_BOUNDARY,
                POWER_IN_WATT_CURVE_FIXED_FLOW_RATE_BOUNDARY,
                FIXED_TEMP_DIFF_BOUNDARY
            };

            /**
              * list the possible primary variables for the HEAT_TRANSPORT_BHE process
              */
            enum class BHE_PRIMARY_VARS {
                BHE_TEMP_SOIL,
                BHE_TEMP_IN_1,
                BHE_TEMP_IN_2,
                BHE_TEMP_OUT_1,
                BHE_TEMP_OUT_2,
                BHE_TEMP_G_1,
                BHE_TEMP_G_2,
                BHE_TEMP_G_3,
                BHE_TEMP_G_4
            };

            /**
              * discharge type of the 2U BHE
              */
            enum class BHE_DISCHARGE_TYPE {
                BHE_DISCHARGE_TYPE_PARALLEL,   // parallel discharge
                BHE_DISCHARGE_TYPE_SERIAL       // serial discharge
            };



            using namespace boost::math::constants;
            static const double PI = boost::math::constants::pi<double>();

            class BHEAbstract : public BHE_Net_ELE_Abstract
            {
            public:

                struct Borehole_Geometry {
                    /**
                    * length/depth of the BHE
                    * unit is m
                    */
                    double L;

                    /**
                    * diameter of the BHE
                    * unit is m
                    */
                    double D;
                };

                struct Pipe_Parameters {
                    /**
                    * radius of the pipline inner side
                    * unit is m
                    */
                    double r_inner;

                    /**
                    * radius of the pipline outer side
                    * unit is m
                    */
                    double r_outer;

                    /**
                    * pipe-in wall thickness
                    * unit is m
                    */
                    double b_in;

                    /**
                    * pipe-out wall thickness
                    * unit is m
                    */
                    double b_out;

                    /**
                    * thermal conductivity of the pipe wall
                    * unit is kg m sec^-3 K^-1
                    */

                    double lambda_p;
                    /**
                    * thermal conductivity of the inner pipe wall
                    * unit is kg m sec^-3 K^-1
                    */
                    double lambda_p_i;

                    /**
                    * thermal conductivity of the outer pipe wall
                    * unit is kg m sec^-3 K^-1
                    */
                    double lambda_p_o;
                };

                struct Refrigerant_Parameters {
                    /**
                      * dynamics viscosity of the refrigerant
                      * unit is kg m-1 sec-1
                      */
                    double mu_r;

                    /**
                      * density of the refrigerant
                      * unit is kg m-3
                      */
                    double rho_r;

                    /**
                      * thermal conductivity of the refrigerant
                      * unit is kg m sec^-3 K^-1
                      */
                    double lambda_r;

                    /**
                      * specific heat capacity of the refrigerant
                      * unit is m^2 sec^-2 K^-1
                      */
                    double heat_cap_r;

                    /**
                      * longitudinal dispersivity of the
                      * referigerant flow in the pipeline
                      */
                    double alpha_L;
                };

                struct Grout_Parameters {
                    /**
                    * density of the grout
                    * unit is kg m-3
                    */
                    double rho_g;

                    /**
                    * porosity of the grout
                    * unit is [-]
                    */
                    double porosity_g;

                    /**
                    * specific heat capacity of the grout
                    * unit is m^2 sec^-2 K^-1
                    */
                    double heat_cap_g;

                    /**
                    * thermal conductivity of the grout
                    * unit is kg m sec^-3 K^-1
                    */
                    double lambda_g;
                };

                struct Extern_Ra_Rb {
                    /**
                    * whether or not using external given borehole thermal resistance values Ra, Rb
                    */
                    bool use_extern_Ra_Rb;

                    /**
                    * external given borehole internal thermal resistance value
                    */
                    double ext_Ra;

                    /**
                    * external given borehole thermal resistance value
                    */
                    double ext_Rb;
                };

                struct Extern_def_Thermal_Resistances {
                    /**
                    * whether or not using user defined borehole thermal resistance Rfig, Rfog, Rgg, Rgs
                    */
                    bool if_use_defined_therm_resis;

                    /**
                    * external given borehole internal thermal resistance value
                    */
                    double ext_Rfig;

                    /**
                    * external given borehole internal thermal resistance value
                    */
                    double ext_Rfog;

                    /**
                    * external given borehole internal thermal resistance value
                    */
                    double ext_Rgg1;

                    /**
                    * external given borehole internal thermal resistance value
                    */
                    double ext_Rgg2;

                    /**
                    * external given borehole internal thermal resistance value
                    */
                    double ext_Rgs;
                };

                /**
                  * constructor
                  */
                BHEAbstract(BHE_TYPE my_type,
                    const std::string name,
                    Borehole_Geometry borehole_geometry_,
                    Pipe_Parameters pipe_param_,
                    Refrigerant_Parameters refrigerant_param_,
                    Grout_Parameters grout_param_,
                    Extern_Ra_Rb extern_Ra_Rb_,
                    Extern_def_Thermal_Resistances extern_def_thermal_resistances_,
                    std::map<std::string, std::unique_ptr<MathLib::PiecewiseLinearInterpolation >> const& bhe_curves,
                    BHE_BOUNDARY_TYPE my_bound_type = BHE_BOUNDARY_TYPE::FIXED_INFLOW_TEMP_BOUNDARY,
                    bool if_use_ext_Ra_Rb = false,
                    bool user_defined_R_vals = false,
                    bool if_flowrate_curve = false,
                    int n_T_in = 1,
                    int n_T_out = 1)
                    : BHE_Net_ELE_Abstract(name,
                        BHE::BHE_NET_ELE::BHE_NET_BOREHOLE, n_T_in, n_T_out),
                    type(my_type),
                    _name(name),
                    borehole_geometry(borehole_geometry_),
                    pipe_param(pipe_param_),
                    refrigerant_param(refrigerant_param_),
                    grout_param(grout_param_),
                    extern_Ra_Rb(extern_Ra_Rb_),
                    extern_def_thermal_resistances(extern_def_thermal_resistances_),
                    bound_type(my_bound_type),
                    _bhe_curves(bhe_curves),
                    use_flowrate_curve(if_flowrate_curve)
                {};

                /**
                  * destructor
                  */
                virtual ~BHEAbstract() {};

                /**
                  * return the number of unknowns needed for this BHE
                  * abstract function, need to be realized.
                  */
                virtual std::size_t get_n_unknowns() = 0;

                /**
                  * return the number of boundary heat exchange terms for this BHE
                  * abstract function, need to be realized.
                  */
                virtual std::size_t get_n_heat_exchange_terms() = 0;

                /**
                  *
                  */
                virtual void set_T_in_out_global_idx(std::size_t start_idx) = 0;

                /**
                  *
                  */
                virtual void set_T_in_out_bottom_global_idx(std::size_t dof_BHE) = 0;

                /**
                  * return the type of boundary condition on this BHE
                  */
                BHE_BOUNDARY_TYPE get_bound_type() { return bound_type; };

                /**
                  * return the name of the BHE
                  */
                const std::string get_name() { return _name; };

                /**
                  * return the thermal resistance for the inlet pipline
                  * idx is the index, when 2U case,
                  * 0 - the first u-tube
                  * 1 - the second u-tube
                  * needs to be overwritten
                  */
                virtual double get_thermal_resistance_fig(std::size_t idx = 0) = 0;

                /**
                  * return the thermal resistance for the outlet pipline
                  * idx is the index, when 2U case,
                  * 0 - the first u-tube
                  * 1 - the second u-tube
                  * needs to be overwritten
                  */
                virtual double get_thermal_resistance_fog(std::size_t idx = 0) = 0;

                /**
                  * return the thermal resistance
                  */
                virtual double get_thermal_resistance(std::size_t idx = 0) = 0;

                /**
                  * initialization calcultion,
                  * need to be overwritten.
                  */
                virtual void initialize()
                {
                    calc_u();
                    calc_Re();
                    calc_Pr();
                    calc_Nu();
                    calc_thermal_resistances();
                    calc_heat_transfer_coefficients();
                };

                /**
                  * update all parameters based on the new flow rate
                  * not necessarily needs to be overwritten.
                  */
                virtual void update_flow_rate(double new_flow_rate)
                {
                    Q_r = new_flow_rate;
                    calc_u();
                    calc_Re();
                    calc_Pr();
                    calc_Nu();
                    calc_thermal_resistances();
                    calc_heat_transfer_coefficients();
                };

                virtual void update_flowrate_from_curve(double current_time) = 0;

                /**
                  * thermal resistance calculation,
                  * need to be overwritten.
                  */
                virtual void calc_thermal_resistances() = 0;

                /**
                  * Nusselt number calculation,
                  * need to be overwritten.
                  */
                virtual void calc_Nu() = 0;

                /**
                  * Renolds number calculation,
                  * need to be overwritten.
                  */
                virtual void calc_Re() = 0;

                /**
                  * Prandtl number calculation,
                  * need to be overwritten.
                  */
                virtual void calc_Pr() = 0;

                /**
                  * flow velocity inside the pipeline
                  * need to be overwritten.
                  */
                virtual void calc_u() = 0;

                /**
                  * heat transfer coefficient,
                  * need to be overwritten.
                  */
                virtual void calc_heat_transfer_coefficients() = 0;

                /**
                  * return the coeff of mass matrix,
                  * depending on the index of unknown.
                  */
                virtual double get_mass_coeff(std::size_t idx_unknown) = 0;

                /**
                  * return the coeff of laplace matrix,
                  * depending on the index of unknown.
                  */
                virtual void get_laplace_matrix(std::size_t idx_unknown, Eigen::MatrixXd & mat_laplace) = 0;

                /**
                  * return the coeff of advection matrix,
                  * depending on the index of unknown.
                  */
                virtual void get_advection_vector(std::size_t idx_unknown, Eigen::VectorXd & vec_advection) = 0;

                /**
                  * return the coeff of boundary heat exchange matrix,
                  * depending on the index of unknown.
                  */
                virtual double get_boundary_heat_exchange_coeff(std::size_t idx_unknown) = 0;

                /**
                  * return the shift index based on primary variable value
                  */
                virtual int get_loc_shift_by_pv(BHE::BHE_PRIMARY_VARS pv_name) = 0;

                /**
                  * return the number of grout zones in this BHE.
                  */
                virtual std::size_t get_n_grout_zones(void) = 0;

                /**
                  * return the inflow temperature based on outflow temperature and fixed power.
                  */
                virtual double get_Tin_by_Tout(double T_in, double current_time) = 0;

                /**
                  * get the polyline geometry
                  * that is representing this BHE.
                  */
                const GeoLib::Polyline* get_geo_ply() { return _geo_ply; }

                /**
                  * set the polyline geometry
                  * that is representing this BHE.
                  */
                void set_geo_ply(const GeoLib::Polyline* ply) { _geo_ply = ply; }

                /**
                  * set the polyline geometry
                  * that is representing this BHE.
                  */
                void set_ply_eps(double eps) { _ply_eps = eps; }

                /**
                  * set the polyline geometry
                  * that is representing this BHE.
                  */
                double get_ply_eps(void) { return _ply_eps; }

                /**
                  * total refrigerant flow discharge of BHE
                  * unit is m^3/sec
                  */
                double Q_r;

                /**
                  * geometry of the borehole
                  */
                Borehole_Geometry const borehole_geometry;

                /**
                  * geometry of the pipes in the borehole
                  */
                Pipe_Parameters const pipe_param;

                /**
                * parameters of the refrigerant
                */
                Refrigerant_Parameters const refrigerant_param;

                /**
                  * parameters of the grout
                  */
                Grout_Parameters const grout_param;

                /**
                  * Ra Rb values defined by the user externally
                  */
                Extern_Ra_Rb const extern_Ra_Rb;

                /**
                  * thermal resistance values defined by the user externally
                  */
                Extern_def_Thermal_Resistances const extern_def_thermal_resistances;

                /**
                  * power extracted from or injected into the BHE
                  * unit is Watt
                  * if value positive, then injecting power
                  * if value negative, then extracting power
                  */
                double power_in_watt_val;

                /**
                  * temperature difference between inflow and
                  * outflow pipelines
                  */
                double delta_T_val;

                /**
                  * threshold Q value for switching off the BHE
                  * when using the Q_curve_fixed_dT B.C.
                  */
                double threshold;


                /**
                  * map strucutre that contains all the curves related to this BHE
                  */
                std::map<std::string, std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const& _bhe_curves;

                /**
                  * power in watt curve
                  */
                MathLib::PiecewiseLinearInterpolation * _power_in_watt_curve;

                /**
                  * heating COP curve
                  */
                MathLib::PiecewiseLinearInterpolation * _heating_cop_curve;

                /**
                  * cooling COP curve
                  */
                MathLib::PiecewiseLinearInterpolation * _cooling_cop_curve;

                /**
                  * refrigerant flow rate curve
                  */
                MathLib::PiecewiseLinearInterpolation * _flowrate_curve;

                /**
                  * use refrigerant flow rate curve
                  */
                bool use_flowrate_curve;

                /**
                  * for BHEs, the RHS value is zero
                  */
                double get_RHS_value()
                {
                    return 0;
                }

            private:

                /**
                  * the type of the BHE
                  */
                const BHE_TYPE type;

                /**
                  * the type of the boundary condition on this BHE
                  */
                const BHE_BOUNDARY_TYPE bound_type;

                /**
                  * the polyline geometry representing the BHE
                  */
                const GeoLib::Polyline* _geo_ply;

                /**
                  * name of the borehole heat exchanger
                  */
                const std::string _name;

                /**
                  * epsilon value of the BHE polyline
                  */
                double _ply_eps;
            };

        }  // end of namespace BHE
    }  // end of namespace HeatTransportBHE
}  // end of namespace ProcessLib

