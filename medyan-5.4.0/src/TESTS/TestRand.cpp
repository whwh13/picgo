#include <cstdint>

#include <catch2/catch.hpp>

#include "Rand.h"

namespace medyan {

TEST_CASE("Safe exponential distribution test", "[Distribution]") {
    using std::isinf;
    using std::uint64_t;
    using namespace rand;

    // Our hand-made "random" engine.
    struct Engine {
        uint64_t value = 0;
        uint64_t called = 0;
        uint64_t operator()() {
            ++called;
            const auto ret = value;
            if(value >= UINT32_MAX) value = 0;
            else ++value;
            return ret;
        }
        static constexpr uint64_t min() { return 0; }
        static constexpr uint64_t max() { return UINT32_MAX; }
    };

    Engine g;

    SECTION("lambda == 0") {
        const float lambdaPrev = 1.0f;
        std::exponential_distribution< float > d(lambdaPrev);

        g.value = 0;
        const auto res0 = safeExpDist(d, 0, g);
        CHECK(isinf(res0));
        CHECK(d.lambda() == lambdaPrev);
        CHECK(g.called == 0);

        g.value = Engine::max() / 2;
        const auto res1 = safeExpDist(d, 0, g);
        CHECK(isinf(res1));
        CHECK(d.lambda() == lambdaPrev);
        CHECK(g.called == 0);

        g.value = Engine::max();
        const auto res2 = safeExpDist(d, 0, g);
        CHECK(isinf(res2));
        CHECK(d.lambda() == lambdaPrev);
        CHECK(g.called == 0);
    }

    SECTION("lambda == inf") {
        const float lambdaPrev = 1.0f;
        std::exponential_distribution< float > d(lambdaPrev);

        g.value = 0;
        const auto res0 = safeExpDist(d, std::numeric_limits< float >::infinity(), g);
        CHECK(res0 == 0);
        CHECK(d.lambda() == lambdaPrev);
        CHECK(g.called == 0);
    }

    SECTION("0 < lambda < inf") {
        std::exponential_distribution< float > d(2.0f);
        const float lambda = 1.0f;

        uint64_t cnt = 0, cntInf = 0;
        for(uint64_t i = 0xffff'fff0; i < 0xffff'ffff; ++i) {
            // These are the values that can create some infinities in certain implementations

            g.value = i;
            const auto res = safeExpDist(d, lambda, g);
            ++cnt;
            if(isinf(res)) ++cntInf;
        }

        REQUIRE(cnt > 0);
        CHECK(cntInf == 0);
        CHECK(d.lambda() == lambda);
        CHECK(g.called >= cnt);
    }

}

TEST_CASE("Auxiliary random functions", "[Distribution]") {
    using namespace medyan;

    // Check that the random unit vector function indeed produces unit vectors.
    for(int rep = 0; rep < 5; ++rep) {
        const auto v1 = Rand::randUnitVector3<double>();
        CHECK(magnitude2(v1) == Approx(1.0));
        const auto v2 = Rand::randUnitVector3<float>();
        CHECK(magnitude2(v2) == Approx(1.0f));
    }
}

} // namespace medyan
