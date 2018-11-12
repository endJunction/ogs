/**
 * \author Norbert Grunwald
 * \date   Oct 22, 2018
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
#include <iostream>
#include "TestMPL.h"

#include "Tests/TestTools.h"

#include "MaterialLib/MPL/mpMedium.h"
namespace MPL = MaterialPropertyLib;


TEST(MPL, SaturationLiakopoulos)
{

    std::stringstream m;
    m << "<medium>\n";
    m << "<phases></phases>\n";
    m << "<properties><property>\n";
    m << "<name>saturation</name>\n";
    m << "<type>SaturationLiakopoulos</type>\n";
    m << "</property></properties>\n";
    m << "</medium>\n";

    auto medium = createTestMaterial(m.str());

    std::vector<double> p_cap = {
            1.,     2.,     4.,    10.,    20.,    40.,   100.,   200.,   400.,
            1000.,  2000.,  4000., 10000., 11000., 12000., 13000., 14000., 15000.,
            16000., 17000., 18000., 19000., 20000., 21000., 22000., 23000., 24000.,
            25000.};

    std::vector<double> sL_expected = {
            9.99999999980278E-01, 9.99999999893874E-01, 9.99999999428925E-01,
            9.99999994717360E-01, 9.99999971573611E-01, 9.99999847034881E-01,
            9.99998585017872E-01, 9.99992385845569E-01, 9.99959027505329E-01,
            9.99620989750676E-01, 9.97960509527072E-01, 9.89025306316687E-01,
            8.98480153023156E-01, 8.72047653500129E-01, 8.41949781074817E-01,
            8.08047390932890E-01, 7.70207837433721E-01, 7.28304201637326E-01,
            6.82214658686562E-01, 6.31821952151409E-01, 5.77012951742652E-01,
            5.17678277081630E-01, 4.53711974565622E-01, 3.85011237456148E-01,
            3.11476161552399E-01, 2.33009530459143E-01, 2.00000000000000E-01,
            2.00000000000000E-01};

    std::vector<double> dsLdpc_expected = {
            -4.78830438000000E-11,-1.28831623708299E-10,-3.46627656684535E-10,
            -1.28257198533247E-09,-3.45082138226893E-09,-9.28460027858602E-09,
            -3.43543510815759E-08,-9.24320277084059E-08,-2.48692799523419E-07,
            -9.20198984332438E-07,-2.47583945961004E-06,-6.66136469842884E-06,
            -2.46480036475078E-05,-2.82414092788214E-05,-3.19775105440375E-05,
            -3.58493645810797E-05,-3.98508851067620E-05,-4.39766819229822E-05,
            -4.82219393859309E-05,-5.25823224924466E-05,-5.70539030257785E-05,
            -6.16331005828163E-05,-6.63166348476062E-05,-7.11014865038198E-05,
            -7.59848648803150E-05,-8.09641809129672E-05, 0.00000000000000E+00,
             0.00000000000000E+00};

    std::vector<double> d2sLdpc2_expected = {
            -6.83721982420200E-11,-9.19793377465406E-11,-1.23737407744961E-10,
            -1.83138453785623E-10,-2.46371392587090E-10,-3.31437018444824E-10,
            -4.90545779093823E-10,-6.59918461824164E-10,-8.87771121098725E-10,
            -1.31395212972828E-09,-1.76762558218859E-09,-2.37794066322163E-09,
            -3.51948844082765E-09,-3.66599166447537E-09,-3.80505727548593E-09,
            -3.93763905271720E-09,-4.06450563171039E-09,-4.18628694118842E-09,
            -4.30350670307317E-09,-4.41660578158615E-09,-4.52595934058384E-09,
            -4.63188970116860E-09,-4.73467614494485E-09,-4.83456250375259E-09,
            -4.93176311648190E-09,-5.02646756198373E-09, 0.00000000000000E+00,
             0.00000000000000E+00};


    for (size_t i= 0; i < p_cap.size(); i++)
    {
        MPL::VariableArray variables;
        variables[MPL::Variables::capillary_pressure] = p_cap[i];

        auto const sL = MPL::getScalar(medium.property(
                MPL::PropertyEnum::saturation), variables);
        auto const dsLdpc = MPL::getScalarDerivative(medium.property(
                MPL::PropertyEnum::saturation), variables,
                MPL::Variables::capillary_pressure);
        auto const d2sLdpc2 = MPL::getScalarDerivative(medium.property(
                MPL::PropertyEnum::saturation), variables,
                MPL::Variables::capillary_pressure,
                MPL::Variables::capillary_pressure);

        ASSERT_NEAR(sL, sL_expected[i], 1.e-12);
        ASSERT_NEAR(dsLdpc, dsLdpc_expected[i], 1.e-12);
        ASSERT_NEAR(d2sLdpc2, d2sLdpc2_expected[i], 1.e-12);

    }

} // TEST(MPL, SaturationLiakopoulos)

TEST(MPL, SaturationFredlund)
{

    std::stringstream m;
    m << "<medium>\n";
    m << "<phases></phases>\n";
    m << "<properties><property>\n";
    m << "<name>saturation</name>\n";
    m << "<type>SaturationFredlund</type>\n";
    m << "<s_max>1.0</s_max>\n";
    m << "<parameter_a>5706.0</parameter_a>\n";
    m << "<parameter_n>0.606</parameter_n>\n";
    m << "<parameter_m>2.617</parameter_m>\n";
    m << "<parameter_psi_r>10000.0</parameter_psi_r>\n";
    m << "</property></properties>\n";
    m << "</medium>\n";

    auto medium = createTestMaterial(m.str());

    std::vector<double> p_cap = {
            1.0E+00,2.0E+00,4.0E+00,1.0E+01,2.0E+01,4.0E+01,1.0E+02,2.0E+02,
            4.0E+02,1.0E+03,2.0E+03,3.0E+03,4.0E+03,1.0E+04,2.0E+04,4.0E+04,
            1.0E+05,2.0E+05,4.0E+05,1.0E+06,2.0E+06,4.0E+06,1.0E+07,2.0E+07,
            4.0E+07,1.0E+08,2.0E+08,4.0E+08,1.0E+09};

    std::vector<double> sL_expected = {
            9.94918837186489E-01,9.92280224352111E-01,9.88283636872069E-01,
            9.79714860181653E-01,9.69370863744307E-01,9.53951226425036E-01,
            9.21877577488932E-01,8.84882004188767E-01,8.33011555813153E-01,
            7.36362190071533E-01,6.40953910644821E-01,5.77662176082160E-01,
            5.30464750781053E-01,3.77735196082819E-01,2.73089529423746E-01,
            1.88206961420621E-01,1.09856647093762E-01,7.20872687740453E-02,
            4.73762521146771E-02,2.75673128759096E-02,1.85303824241269E-02,
            1.25654997371142E-02,7.56270206457136E-03,5.12710336399544E-03,
            3.42198310961213E-03,1.89745305811559E-03,1.10989392526772E-03,
            5.34660710851223E-04,0.00000000000000E+00};

    std::vector<double> dsLdpc_expected = {
            -3.06879846171774E-03,-2.32657249392803E-03,-1.76008152232360E-03,
            -1.21065057316003E-03,-9.06311319200425E-04,-6.72545905853371E-04,
            -4.43847748285339E-04,-3.16176436487250E-04,-2.18102551901694E-04,
            -1.24167786653234E-04,-7.50814580120344E-05,-5.37768774652272E-05,
            -4.15830475825951E-05,-1.61544190560502E-05,-6.90725091094361E-06,
            -2.65511354671603E-06,-6.63271884244165E-07,-2.19489331271572E-07,
            -7.11322809603052E-08,-1.60072502407637E-08,-5.24376682296463E-09,
            -1.74638757432949E-09,-4.19871104245975E-10,-1.45832018197477E-10,
            -5.14801089480583E-11,-1.32788832948586E-11,-4.82816311421309E-12,
            -1.77112987745904E-12,-4.74891885765972E-13};

    std::vector<double> d2sLdpc2_expected = {
            1.22236872836409E-03,4.66237252439269E-04,1.78035813782725E-04,
            4.99931754583531E-05,1.91762852038913E-05,7.37487149661781E-06,
            2.09189497658792E-06,8.05849829297494E-07,3.08014166573882E-07,
            8.34368222842290E-08,2.94451723556196E-08,1.54771110937140E-08,
            9.61665294500674E-09,1.84533887160621E-09,4.51316784990743E-10,
            9.59893445719224E-11,1.04076122334185E-11,1.77318596539052E-12,
            2.90048047400246E-13,2.59325395531261E-14,4.19209569374042E-15,
            6.86842066843188E-16,6.45835478778439E-17,1.10357291856394E-17,
            1.91949614169229E-18,1.94818144671141E-19,3.50673661936609E-20,
            6.38306591325196E-21,6.80328096397525E-22};


    for (size_t i= 0; i < p_cap.size(); i++)
    {
        MPL::VariableArray variables;
        variables[MPL::Variables::capillary_pressure] = p_cap[i];

        auto const sL = MPL::getScalar(medium.property(
                MPL::PropertyEnum::saturation), variables);
        auto const dsLdpc = MPL::getScalarDerivative(medium.property(
                MPL::PropertyEnum::saturation), variables,
                MPL::Variables::capillary_pressure);
        auto const d2sLdpc2 = MPL::getScalarDerivative(medium.property(
                MPL::PropertyEnum::saturation), variables,
                MPL::Variables::capillary_pressure,
                MPL::Variables::capillary_pressure);

        ASSERT_NEAR(sL, sL_expected[i], 1.e-12);
        ASSERT_NEAR(dsLdpc, dsLdpc_expected[i], 1.e-12);
        ASSERT_NEAR(d2sLdpc2, d2sLdpc2_expected[i], 1.e-12);

    }

} // TEST(MPL, SaturationFredlund)

TEST(MPL, SaturationBrooksCorey)
{

    std::stringstream m;
    m << "<medium>\n";
    m << "<phases></phases>\n";
    m << "<properties><property>\n";
    m << "<name>saturation</name>\n";
    m << "<type>SaturationBrooksCorey</type>\n";
    m << "<residual_liquid_saturation>0.1</residual_liquid_saturation>\n";
    m << "<residual_gas_saturation>0.2</residual_gas_saturation>\n";
    m << "<lambda>2.0</lambda>\n";
    m << "<entry_pressure>5000.0</entry_pressure>\n";
    m << "</property></properties>\n";
    m << "</medium>\n";

    auto medium = createTestMaterial(m.str());


    std::vector<double> p_cap = {
            5.0E+03,6.0E+03,7.0E+03,8.0E+03,9.0E+03,1.0E+04,1.2E+04,1.4E+04,
            1.6E+04,1.8E+04,2.0E+04,2.5E+04,3.0E+04,3.5E+04,4.0E+04,4.5E+04,
            5.0E+04,1.0E+05,1.5E+05,2.0E+05,2.5E+05,3.0E+05,3.5E+05,4.0E+05,
            4.5E+05,5.0E+05};

    std::vector<double> sL_expected = {
            8.00000000000000E-01,5.86111111111111E-01,4.57142857142857E-01,
            3.73437500000000E-01,3.16049382716049E-01,2.75000000000000E-01,
            2.21527777777777E-01,1.89285714285714E-01,1.68359375000000E-01,
            1.54012345679012E-01,1.43750000000000E-01,1.28000000000000E-01,
            1.19444444444444E-01,1.14285714285714E-01,1.10937500000000E-01,
            1.08641975308642E-01,1.07000000000000E-01,1.01750000000000E-01,
            1.00777777777777E-01,1.00437500000000E-01,1.00280000000000E-01,
            1.00194444444444E-01,1.00142857142857E-01,1.00109375000000E-01,
            1.00086419753086E-01,1.00070000000000E-01};

    std::vector<double> dsLdpc_expected = {
            -2.80000000000000E-04,-1.62037037037037E-04,-1.02040816326530E-04,
            -6.83593750000000E-05,-4.80109739368998E-05,-3.50000000000000E-05,
            -2.02546296296296E-05,-1.27551020408163E-05,-8.54492187500000E-06,
            -6.00137174211248E-06,-4.37500000000000E-06,-2.24000000000000E-06,
            -1.29629629629629E-06,-8.16326530612244E-07,-5.46875000000000E-07,
            -3.84087791495198E-07,-2.80000000000000E-07,-3.50000000000000E-08,
            -1.03703703703703E-08,-4.37500000000000E-09,-2.24000000000000E-09,
            -1.29629629629629E-09,-8.16326530612244E-10,-5.46875000000000E-10,
            -3.84087791495199E-10,-2.80000000000000E-10};

    std::vector<double> d2sLdpc2_expected = {
            1.68000000000000E-07,8.10185185185185E-08,4.37317784256559E-08,
            2.56347656250000E-08,1.60036579789666E-08,1.05000000000000E-08,
            5.06365740740740E-09,2.73323615160349E-09,1.60217285156250E-09,
            1.00022862368541E-09,6.56250000000000E-10,2.68800000000000E-10,
            1.29629629629629E-10,6.99708454810495E-11,4.10156250000000E-11,
            2.56058527663465E-11,1.68000000000000E-11,1.05000000000000E-12,
            2.07407407407407E-13,6.56250000000000E-14,2.68800000000000E-14,
            1.29629629629629E-14,6.99708454810495E-15,4.10156250000000E-15,
            2.56058527663466E-15,1.68000000000000E-15};


    for (size_t i= 0; i < p_cap.size(); i++)
    {
        MPL::VariableArray variables;
        variables[MPL::Variables::capillary_pressure] = p_cap[i];

        auto const sL = MPL::getScalar(medium.property(
                MPL::PropertyEnum::saturation), variables);
        auto const dsLdpc = MPL::getScalarDerivative(medium.property(
                MPL::PropertyEnum::saturation), variables,
                MPL::Variables::capillary_pressure);
        auto const d2sLdpc2 = MPL::getScalarDerivative(medium.property(
                MPL::PropertyEnum::saturation), variables,
                MPL::Variables::capillary_pressure,
                MPL::Variables::capillary_pressure);

        ASSERT_NEAR(sL, sL_expected[i], 1.e-12);
        ASSERT_NEAR(dsLdpc, dsLdpc_expected[i], 1.e-12);
        ASSERT_NEAR(d2sLdpc2, d2sLdpc2_expected[i], 1.e-12);

    }

} // TEST(MPL, SaturationBrooksCorey)
