/**
 * \author Norbert Grunwald
 * \date   Oct 24, 2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "gtest/gtest.h"

#include <memory>
#include <sstream>
#include "TestMPL.h"
#include <iostream>

#include "Tests/TestTools.h"

#include "MaterialLib/MPL/mpMedium.h"
namespace MPL = MaterialPropertyLib;


TEST(MPL, DensityBilinearTemperaturePressure)
{

    const double beta_T = -1.23e-4;
    const double beta_p = 1.2e-3;
    const double rho_0 = 1234.5;
    const double p_0 = 123456.7;
    const double T_0 = 123.4;

    std::stringstream m;
    m << std::setprecision(16);
    m << "<medium>\n";
    m << "<phases><phase>\n";
    m << "<properties>\n";

    m << "<property>\n";
    m << "<name>density</name>\n";
    m << "<type>BilinearTemperaturePressure</type>\n";
    m << "<reference_density>" << rho_0 << "</reference_density>\n";
    m << "<reference_temperature>" << T_0 << "</reference_temperature>\n";
    m << "<reference_pressure>" << p_0 << "</reference_pressure>\n";
    m << "</property>\n";

    m << "<property>\n";
    m << "<name>thermal_expansivity</name>\n";
    m << "<type>Constant</type>\n";
    m << "<value>" << beta_T << "</value>\n";
    m << "</property>\n";

    m << "<property>\n";
    m << "<name>compressibility</name>\n";
    m << "<type>Constant</type>\n";
    m << "<value>" << beta_p << "</value>\n";
    m << "</property>\n";

    m << "</properties>\n";
    m << "</phase></phases>\n";
    m << "</medium>\n";

    auto const medium = createTestMaterial(m.str());
    auto const& phase = medium.phase(0);

    const double T_min = 100.0;
    const double T_max = 200.0;
    const double dT = 5.0;

    const double p_min = 100000.0;
    const double p_max = 200000.0;
    const double dp = 5000.0;

    for (double T = T_min; T <= T_max; T+=dT)
        for (double p = p_min; p <= p_max; p+=dp)
    {
        MPL::VariableArray variables;
        variables[MPL::Variables::temperature] = T;
        variables[MPL::Variables::phase_pressure] = p;

        auto const rho = MPL::getScalar(phase.property(
                MPL::PropertyEnum::density), variables);
        auto const drhodT = MPL::getScalarDerivative(phase.property(
                MPL::PropertyEnum::density), variables,
                MPL::Variables::temperature);
        auto const drhodp = MPL::getScalarDerivative(phase.property(
                MPL::PropertyEnum::density), variables,
                MPL::Variables::phase_pressure);

        const double rho_expected = rho_0 *
                (1 + beta_p * (p-p_0) + beta_T * (T-T_0));

        const double drhodT_expected = rho_0 * beta_T;
        const double drhodp_expected = rho_0 * beta_p;

        ASSERT_NEAR(rho, rho_expected, 1.e-15);
        ASSERT_NEAR(drhodT, drhodT_expected, 1.e-15);
        ASSERT_NEAR(drhodp, drhodp_expected, 1.e-15);

    }

} // TEST(MPL, SaturationLiakopoulos)



