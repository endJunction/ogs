/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>
#include <autocheck/autocheck.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <iostream>
#include <limits>

#include "Tests/TestTools.h"
#include "Tests/AutoCheckTools.h"

#include "NumLib/Fem/Integration/IntegrationGaussLobattoRegular.h"

template <int Order>
struct Polynomial;

namespace autocheck
{
template <int Order, typename Gen = generator<double>>
struct RandomPolynomialGenerator
{

    randomTupleGenerator<double, Order + 1, Gen>
        _a_generator;
    randomTupleGenerator<double, Order, Gen> _b_generator;

    using result_type = Polynomial<Order>;

    result_type operator()(std::size_t size = 0)
    {
        result_type p{_a_generator(size), _b_generator(size)};
        return std::move(p);
    }
};

}  // namespace autocheck

using namespace NumLib;
using namespace boost::math;
namespace ac = autocheck;

TEST(NumLibFemIntegrationGaussLobatto, NumberOfPointsRegular)
{
    // check position indices
    // dim = 1
    {
        using IntegrationMethod = IntegrationGaussLobattoRegular<1>;
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(1, 0)[0]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(1, 1)[0]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(1, 2)[0]);

        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(2, 0)[0]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(2, 1)[0]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(2, 2)[0]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(2, 3)[0]);
    }
    // dim = 2
    {
        using IntegrationMethod = IntegrationGaussLobattoRegular<2>;
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(1, 0)[0]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(1, 0)[1]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(1, 1)[0]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(1, 1)[1]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(1, 2)[0]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(1, 2)[1]);

        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(1, 3)[0]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(1, 3)[1]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(1, 4)[0]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(1, 4)[1]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(1, 5)[0]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(1, 5)[1]);

        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(1, 6)[0]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(1, 6)[1]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(1, 7)[0]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(1, 7)[1]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(1, 8)[0]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(1, 8)[1]);

        // 4 integration point per dimension
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(2, 0)[0]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(2, 0)[1]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(2, 1)[0]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(2, 1)[1]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(2, 2)[0]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(2, 2)[1]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(2, 3)[0]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(2, 3)[1]);

        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(2, 4)[0]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(2, 4)[1]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(2, 5)[0]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(2, 5)[1]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(2, 6)[0]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(2, 6)[1]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(2, 7)[0]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(2, 7)[1]);

        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(2, 8)[0]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(2, 8)[1]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(2, 9)[0]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(2, 9)[1]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(2, 10)[0]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(2, 10)[1]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(2, 11)[0]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(2, 11)[1]);

        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(2, 12)[0]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(2, 12)[1]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(2, 13)[0]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(2, 13)[1]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(2, 14)[0]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(2, 14)[1]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(2, 15)[0]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(2, 15)[1]);

        // 5 integration point per dimension
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(3, 0)[0]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(3, 0)[1]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(3, 1)[0]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(3, 1)[1]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(3, 2)[0]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(3, 2)[1]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(3, 3)[0]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(3, 3)[1]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(3, 4)[0]);
        ASSERT_EQ(4, IntegrationMethod::getPositionIndices(3, 4)[1]);

        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(3, 5)[0]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(3, 5)[1]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(3, 6)[0]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(3, 6)[1]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(3, 7)[0]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(3, 7)[1]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(3, 8)[0]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(3, 8)[1]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(3, 9)[0]);
        ASSERT_EQ(4, IntegrationMethod::getPositionIndices(3, 9)[1]);

        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(3, 10)[0]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(3, 10)[1]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(3, 11)[0]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(3, 11)[1]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(3, 12)[0]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(3, 12)[1]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(3, 13)[0]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(3, 13)[1]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(3, 14)[0]);
        ASSERT_EQ(4, IntegrationMethod::getPositionIndices(3, 14)[1]);

        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(3, 15)[0]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(3, 15)[1]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(3, 16)[0]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(3, 16)[1]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(3, 17)[0]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(3, 17)[1]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(3, 18)[0]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(3, 18)[1]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(3, 19)[0]);
        ASSERT_EQ(4, IntegrationMethod::getPositionIndices(3, 19)[1]);

        ASSERT_EQ(4, IntegrationMethod::getPositionIndices(3, 20)[0]);
        ASSERT_EQ(0, IntegrationMethod::getPositionIndices(3, 20)[1]);
        ASSERT_EQ(4, IntegrationMethod::getPositionIndices(3, 21)[0]);
        ASSERT_EQ(1, IntegrationMethod::getPositionIndices(3, 21)[1]);
        ASSERT_EQ(4, IntegrationMethod::getPositionIndices(3, 22)[0]);
        ASSERT_EQ(2, IntegrationMethod::getPositionIndices(3, 22)[1]);
        ASSERT_EQ(4, IntegrationMethod::getPositionIndices(3, 23)[0]);
        ASSERT_EQ(3, IntegrationMethod::getPositionIndices(3, 23)[1]);
        ASSERT_EQ(4, IntegrationMethod::getPositionIndices(3, 24)[0]);
        ASSERT_EQ(4, IntegrationMethod::getPositionIndices(3, 24)[1]);
    }
}

template <int Order>
struct Polynomial
{
    double operator()(double const x) const
    {
        double result = a_i[0];
        for (int i = 1; i < Order + 1; ++i)
        {
            result += a_i[i] * std::pow(x - b_i[i-1], i);
        }
        return result;
    }

    double integrate(double const x_low, double const x_high) const
    {
        double result = a_i[0] * (x_high - x_low);
        for (int i = 1; i < Order + 1; ++i)
        {
            result += a_i[i] / (i + 1) *
                      (std::pow(x_high - b_i[i - 1], i + 1) -
                       std::pow(x_low - b_i[i - 1], i + 1));
        }
        return result;
    }

    std::array<double, Order+1> a_i;
    std::array<double, Order> b_i;
};

template <int Order>
std::ostream& operator<<(std::ostream& os, Polynomial<Order> const& p)
{
    return os << "a_i : " << p.a_i << ", b_i = " << p.b_i << "\n";
}

TEST(GaussLobatto, PolyIntegration)
{
    {
        Polynomial<0> p;
        p.a_i = {{1}};
        p.b_i = {{}};
        ASSERT_EQ(1, p(0.));
        ASSERT_EQ(1, p(1.));

        ASSERT_EQ(1, p.integrate(0., 1.));
        ASSERT_EQ(2, p.integrate(-1., 1.));
        ASSERT_EQ(3, p.integrate(-1., 2.));
    }

    {
        Polynomial<2> p;
        p.a_i = {{1, 2, 3}};
        p.b_i = {{4, 5}};
        EXPECT_EQ(68., p(0.));
        EXPECT_EQ(43, p(1.));

        EXPECT_EQ(55, p.integrate(0., 1.));
        EXPECT_EQ(138, p.integrate(-1., 1.));
        EXPECT_EQ(171, p.integrate(-1., 2.));
    }
}

template <typename IntegrationMethod, int Order>
double integratePoly(IntegrationMethod const& integration_method,
                     Polynomial<Order> const& p)
{
    double sum = 0;

    auto const n_integration_points = integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const w = integration_method.getWeightedPoint(ip).getWeight();
        auto const x = integration_method.getWeightedPoint(ip).getCoords()[0];
        sum += p(x) * w;
    }

    return sum;
}

template <typename IntegrationMethod>
struct IntegrationTester
{
    template <int Order>
    bool operator()(Polynomial<Order> const& p) const
    {
        double const computed_value = integratePoly(_integration_method, p);
        double const true_value = p.integrate(-1, 1);
        if (true_value == 0)
        {
            return std::abs(computed_value - true_value) < abs_error;
        }
        return std::abs((computed_value - true_value) / true_value) < rel_error;
    }

    IntegrationMethod const& _integration_method;
    double const abs_error;
    double const rel_error;

};

struct NumLibFemIntegrationGaussLobatto1D : public ::testing::Test
{
    using IntegrationMethod = IntegrationGaussLobattoRegular<1>;
    ac::gtest_reporter gtest_reporter;
};

TEST_F(NumLibFemIntegrationGaussLobatto1D, P0_1)
{
    ac::check<Polynomial<0>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(1), 1e-15,
                                             1e-11},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<0>()),
        gtest_reporter);
}

TEST_F(NumLibFemIntegrationGaussLobatto1D, P1_1)
{
    ac::check<Polynomial<1>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(1), 1e-15,
                                             1e-12},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<1>()),
        gtest_reporter);
}

TEST_F(NumLibFemIntegrationGaussLobatto1D, P2_1)
{
    ac::check<Polynomial<2>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(1), 1e-15,
                                             1e-11},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<2>()),
        gtest_reporter);
}

TEST_F(NumLibFemIntegrationGaussLobatto1D, P3_1)
{
    ac::check<Polynomial<3>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(1), 1e-15,
                                             1e-12},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<3>()),
        gtest_reporter);
}

TEST_F(NumLibFemIntegrationGaussLobatto1D, P4_2)
{
    ac::check<Polynomial<4>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(2), 1e-15,
                                             1e-11},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<4>()),
        gtest_reporter);
}

TEST_F(NumLibFemIntegrationGaussLobatto1D, P5_2)
{
    ac::check<Polynomial<5>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(2), 1e-15,
                                             1e-12},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<5>()),
        gtest_reporter);
}

TEST_F(NumLibFemIntegrationGaussLobatto1D, P5_3)
{
    ac::check<Polynomial<5>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(3), 1e-15,
                                             1e-12},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<5>()),
        gtest_reporter);
}

TEST_F(NumLibFemIntegrationGaussLobatto1D, P6_3)
{
    ac::check<Polynomial<6>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(3), 1e-15,
                                             1e-11},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<6>()),
        gtest_reporter);
}

TEST_F(NumLibFemIntegrationGaussLobatto1D, P7_3)
{
    ac::check<Polynomial<7>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(3), 1e-15,
                                             1e-12},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<7>()),
        gtest_reporter);
}

TEST_F(NumLibFemIntegrationGaussLobatto1D, P8_4)
{
    ac::check<Polynomial<8>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(4), 1e-15,
                                             1e-12},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<8>()),
        gtest_reporter);
}

TEST_F(NumLibFemIntegrationGaussLobatto1D, P9_4)
{
    ac::check<Polynomial<9>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(4), 1e-15,
                                             1e-12},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<9>()),
        gtest_reporter);
}

TEST_F(NumLibFemIntegrationGaussLobatto1D, P10_5)
{
    ac::check<Polynomial<10>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(5), 1e-15,
                                             1e-12},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<10>()),
        gtest_reporter);
}

TEST_F(NumLibFemIntegrationGaussLobatto1D, P11_5)
{
    ac::check<Polynomial<11>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(5), 1e-15,
                                             1e-12},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<11>()),
        gtest_reporter);
}

TEST_F(NumLibFemIntegrationGaussLobatto1D, P12_6)
{
    ac::check<Polynomial<12>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(6), 1e-15,
                                             1e-12},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<12>()),
        gtest_reporter);
}

TEST_F(NumLibFemIntegrationGaussLobatto1D, P13_6)
{
    ac::check<Polynomial<13>>(
        IntegrationTester<IntegrationMethod>{IntegrationMethod(6), 1e-15,
                                             1e-12},
        1000, ac::make_arbitrary(ac::RandomPolynomialGenerator<13>()),
        gtest_reporter);
}
