#include <catch2/catch.hpp>

#include "Mechanics/ForceField/Membrane/MembraneMeshless.hpp"
#include "TESTS/Mechanics/ForceField/TestFFCommon.hpp"

namespace medyan {

static auto defaultMembraneMeshlessFF10() {
    // The 10 vertices are neighbors of each other.
    // Used when test coordinates contain only meshless vertices.
    constexpr int nv = 10;

    MembraneMeshlessFF ff;
    ff.verticesInfo.resize(nv);
    for(int i = 0; i < nv; ++i) {
        for(int j = i+1; j < nv; ++j) {
            ff.vertexPairsInfo.push_back({
                i*5, j*5,
                i,   j,
            });
        }
    }
    return ff;
}

TEST_CASE("Force field: Membrane meshless", "[ForceField]") {
    // The test case checks whether energy and force are consistent

    using namespace std;
    using namespace test_ff_common;
    using VF = std::vector< floatingpoint >;

    // Prepare data
    //---------------------------------
    struct TestInput {
        std::string        name;
        MembraneMeshlessFF ff;
        VF                 coord;
    };


    // Multiple test sets
    std::vector< TestInput > testInputs;

    // Input: standard symmetric case
    testInputs.push_back({
        "generic",
        // FF parameters.
        [] {
            MembraneMeshlessFF ff = defaultMembraneMeshlessFF10();
            return ff;
        }(),
        // Coordinates.
        {
            0,     0,     0,      1.0, 0.0,

            1,     0,     -0.2,   1.1, 0.1,
            0.62,  0.78,  -0.2,   -0.2, 0.1,
            -0.22, 0.97,  -0.3,   0.1, 0.3,
            -0.90, 0.43,  -0.2,   1.9, 0.0,
            -0.90, -0.43, -0.1,   4.2, 0.1,
            -0.22, -0.97, -0.24,  3.0, 0.1,
            0.62,  -0.78, -0.2,   -0.2, 0.1,

            2,     0,     0,      1.0, 0.7,
            3,     0.1,   0.2,    0.0, 1.5,
        },
    });


    // Prepare test parameters
    //---------------------------------
    const floatingpoint moveMag      = std::is_same< floatingpoint, float >::value ? 5e-3 : 1e-5;

    const floatingpoint diffDeRelEps = std::is_same< floatingpoint, float >::value ? 6e-2 : 5e-4;

    const size_t repsTot     = 10;
    const size_t repsPassReq = 9;


    // Run tests for each set of input
    //---------------------------------
    for(auto& ti : testInputs) {

        // Prepare functions
        //---------------------------------
        const auto calcEnergy = [&](const VF& c) {
            auto tempC = c;
            auto ret = ti.ff.computeEnergy(tempC.data());
            return ret;
        };
        const auto calcForce = [&](const VF& c, VF& f) {
            auto tempC = c;
            ti.ff.computeForces(tempC.data(), f.data());
        };


        // Run the test
        //---------------------------------
        size_t repsPass = 0;
        for(size_t rep = 1; rep <= repsTot; ++rep) {
            const auto res = testEnergyForceConsistency(ti.coord, calcEnergy, calcForce, moveMag, diffDeRelEps);
            if(res.passed)
                ++repsPass;
            else
                WARN(ti.name << " (Rep " << rep << '/' << repsTot << " fail) E: " << calcEnergy(ti.coord) << " Actual de: " << res.deActual << " Expected de: " << res.deExpected);
        }
        CHECK(repsPass >= repsPassReq);
    }
}

} // namespace medyan
