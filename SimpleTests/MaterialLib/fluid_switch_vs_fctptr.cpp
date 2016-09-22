#include <celero/Celero.h>

#include "Applications/ApplicationsLib/LogogSetup.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/ConstantFluidProperty.h"

using namespace MaterialLib::Fluid;

class LiquidDensitySwitchCase final : public FluidProperty
{
public:
    explicit LiquidDensitySwitchCase(std::array<double, 5>& parameters)
        : _beta(parameters[0]),
          _rho0(parameters[1]),
          _temperature0(parameters[2]),
          _p0(parameters[3]),
          _bulk_modulus(parameters[4])
    {
    }

    std::string getName() const override
    {
        return "Temperature or pressure dependent liquid density";
    }

    double getValue(const ArrayType& var_vals) const override
    {
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        const double p = var_vals[static_cast<int>(PropertyVariableType::pl)];
        return _rho0 / (1. + _beta * (T - _temperature0)) /
               (1. - (p - _p0) / _bulk_modulus);
    }

    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var) const override
    {
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        const double p = var_vals[static_cast<int>(PropertyVariableType::pl)];
        switch (var)
        {
            case PropertyVariableType::T:
                return dLiquidDensity_dT(T, p);
            case PropertyVariableType::pl:
                return dLiquidDensity_dp(T, p);
            default:
                return 0.;
        }
    }

private:
    const double _beta, _rho0, _temperature0, _p0, _bulk_modulus;

    double dLiquidDensity_dT(const double T, const double p) const
    {
        const double fac_T = 1. + _beta * (T - _temperature0);
        return -_beta * _rho0 / (fac_T * fac_T) /
               (1. - (p - _p0) / _bulk_modulus);
    }
    double dLiquidDensity_dp(const double T, const double p) const
    {
        const double fac_p = 1. - (p - _p0) / _bulk_modulus;
        return _rho0 / (1. + _beta * (T - _temperature0)) /
               (fac_p * fac_p * _bulk_modulus);
    }
};

class LiquidDensityFunctionPointer final : public FluidProperty
{
public:
    explicit LiquidDensityFunctionPointer(std::array<double, 5>& parameters)
        : _beta(parameters[0]),
          _rho0(parameters[1]),
          _temperature0(parameters[2]),
          _p0(parameters[3]),
          _bulk_modulus(parameters[4])
    {
    }

    std::string getName() const override
    {
        return "Temperature or pressure dependent liquid density";
    }

    double getValue(const ArrayType& var_vals) const override
    {
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        const double p = var_vals[static_cast<int>(PropertyVariableType::pl)];
        return _rho0 / (1. + _beta * (T - _temperature0)) /
               (1. - (p - _p0) / _bulk_modulus);
    }

    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var) const override
    {
        const int func_id = static_cast<int>(var);
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        const double p = var_vals[static_cast<int>(PropertyVariableType::pl)];
        return (this->*_derivative_functions[func_id])(T, p);
    }

private:
    const double _beta, _rho0, _temperature0, _p0, _bulk_modulus;

    double dLiquidDensity_dT(const double T, const double p) const
    {
        const double fac_T = 1. + _beta * (T - _temperature0);
        return -_beta * _rho0 / (fac_T * fac_T) /
               (1. - (p - _p0) / _bulk_modulus);
    }

    double dLiquidDensity_dp(const double T, const double p) const
    {
        const double fac_p = 1. - (p - _p0) / _bulk_modulus;
        return _rho0 / (1. + _beta * (T - _temperature0)) /
               (fac_p * fac_p * _bulk_modulus);
    }

    typedef double (LiquidDensityFunctionPointer::*DerivativeFunctionPointer)(
        const double, const double) const;

    static DerivativeFunctionPointer
        _derivative_functions[PropertyVariableNumber];
};

LiquidDensityFunctionPointer::DerivativeFunctionPointer
    LiquidDensityFunctionPointer::_derivative_functions
        [PropertyVariableNumber] = {
            &LiquidDensityFunctionPointer::dLiquidDensity_dT,
            &LiquidDensityFunctionPointer::dLiquidDensity_dp, nullptr};

/// BENCHMARKS ////////////////////////////////////////////////////////////////

struct BaseFixture : public celero::TestFixture
{
    const double T0 = 273.15;
    const double p0 = 1.e+5;
    const double rho0 = 999.8;
    const double K = 2.15e+9;
    const double beta = 2.e-4;
    std::array<double, 5> parameters{{T0, p0, rho0, K, beta}};

    const std::array<double, 3> vars = {{273.15 + 60.0, 1.e+6, 0.}};
    const double T = vars[0];
    const double p = vars[1];

    double const eps = 1e-10;

    FluidProperty const* rhoConstant;
    FluidProperty const* rhoSwitch;
    FluidProperty const* rhoFctPtr;

    void setUp(int64_t experimentValue) override
    {
        rhoConstant = new ConstantFluidProperty{rho0};
        rhoSwitch = new LiquidDensitySwitchCase{parameters};
        rhoFctPtr = new LiquidDensityFunctionPointer{parameters};
    }

    void tearDown() override {
        delete rhoConstant;
        delete rhoSwitch;
        delete rhoFctPtr;
    }
    /*
std::vector<std::pair<int64_t, uint64_t>> getExperimentValues()
    const override
{
    std::vector<std::pair<int64_t, uint64_t>> vectorSizes;
    for (int i = 0; i < 5; ++i)
        vectorSizes.push_back(
            std::pair<int64_t, uint64_t>(std::pow(16, i), 0));
    return vectorSizes;
}
    */
};

BASELINE_F(FluidDensity2, Constant, BaseFixture, 100, 1000000)
{
    celero::DoNotOptimizeAway(
        rhoConstant->getdValue(vars, PropertyVariableType::T));
    celero::DoNotOptimizeAway(
        rhoConstant->getdValue(vars, PropertyVariableType::pl));
}
BENCHMARK_F(FluidDensity2, FctPtr, BaseFixture, 100, 1000000)
{
    celero::DoNotOptimizeAway(
        rhoFctPtr->getdValue(vars, PropertyVariableType::T));
    celero::DoNotOptimizeAway(
        rhoFctPtr->getdValue(vars, PropertyVariableType::pl));
}

BENCHMARK_F(FluidDensity2, switchCase, BaseFixture, 100, 1000000)
{
    celero::DoNotOptimizeAway(
        rhoSwitch->getdValue(vars, PropertyVariableType::T));
    celero::DoNotOptimizeAway(
        rhoSwitch->getdValue(vars, PropertyVariableType::pl));
}

BASELINE_F(FluidDensity1, Constant, BaseFixture, 100, 1000000)
{
    celero::DoNotOptimizeAway(
        rhoConstant->getdValue(vars, PropertyVariableType::T));
}

BENCHMARK_F(FluidDensity1, FctPtrT, BaseFixture, 100, 1000000)
{
    celero::DoNotOptimizeAway(
        rhoFctPtr->getdValue(vars, PropertyVariableType::T));
}
BENCHMARK_F(FluidDensity1, FctPtrP, BaseFixture, 100, 1000000)
{
    celero::DoNotOptimizeAway(
        rhoFctPtr->getdValue(vars, PropertyVariableType::pl));
}

BENCHMARK_F(FluidDensity1, switchCaseT, BaseFixture, 100, 1000000)
{
    celero::DoNotOptimizeAway(
        rhoSwitch->getdValue(vars, PropertyVariableType::T));
}
BENCHMARK_F(FluidDensity1, switchCaseP, BaseFixture, 100, 1000000)
{
    celero::DoNotOptimizeAway(
        rhoSwitch->getdValue(vars, PropertyVariableType::pl));
}
/*
BENCHMARK_F(t, t20, BaseFixture, 30, 10000000)
{
    auto intersections = polygon->getAllIntersectionPoints(l);
    celero::DoNotOptimizeAway(intersections.size());
}
*/

int main(int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    celero::Run(argc, argv);
}
