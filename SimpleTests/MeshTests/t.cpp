#include <celero/Celero.h>

#include "Applications/ApplicationsLib/LogogSetup.h"
#include "GeoLib/Polygon.h"

struct BaseFixture : public celero::TestFixture
{
    std::vector<GeoLib::Point*> points;

    GeoLib::Point a{0.5, -1, 0};
    GeoLib::Point b{0.5, 1, 0};
    GeoLib::LineSegment l{&a, &b};
    GeoLib::Polygon* polygon;

    void setUp(int64_t experimentValue) override
    {
        for (double i = 0.; i < static_cast<double>(experimentValue); ++i)
        {
            points.emplace_back(new GeoLib::Point{0, 0, i});
            points.emplace_back(new GeoLib::Point{1, 0, i});
            points.emplace_back(new GeoLib::Point{1, 1, i});
        }
        points.emplace_back(new GeoLib::Point{0, 0, 0});
        polygon = new GeoLib::Polygon{points, false};
    }

    void tearDown() override
    {
        for (auto* p_ptr : points)
            delete p_ptr;
        delete polygon;
    }

    std::vector<std::pair<int64_t, uint64_t>> getExperimentValues()
        const override
    {
        std::vector<std::pair<int64_t, uint64_t>> vectorSizes;
        for (int i = 0; i < 5; ++i)
            vectorSizes.push_back(
                std::pair<int64_t, uint64_t>(std::pow(16, i), 0));
        return vectorSizes;
    }
};

BASELINE_F(t, t, BaseFixture, 30, 10000000)
{
    celero::DoNotOptimizeAway(polygon->getAllIntersectionPoints(l));
}

BENCHMARK_F(t, t10, BaseFixture, 30, 10000000)
{
    auto intersections = polygon->getAllIntersectionPoints(l);
    celero::DoNotOptimizeAway(intersections);
}

BENCHMARK_F(t, t20, BaseFixture, 30, 10000000)
{
    auto intersections = polygon->getAllIntersectionPoints(l);
    celero::DoNotOptimizeAway(intersections.size());
}
int main(int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    celero::Run(argc, argv);
}

