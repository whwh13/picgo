#include <catch2/catch.hpp>

#include "Mechanics/ForceField/Membrane/FixedVertexAttachmentStretching.hpp"
#include "TESTS/Mechanics/ForceField/TestFFCommon.hpp"

namespace medyan {

TEST_CASE("Force field: Fixed vertex attachment stretching", "[ForceField]") {
    // The test case checks whether energy and force are consistent

    using namespace std;
    using namespace test_ff_common;
    using FF = FixedVertexAttachmentStretching;
    using VF = std::vector< floatingpoint >;

    // Prepare data
    //---------------------------------
    struct TestInput {
        std::string                  name;
        std::vector<FF::Interaction> interactions;
        VF                           coord;
    };


    // Multiple test sets
    std::vector< TestInput > testInputs;

    // Input: standard symmetric case
    testInputs.push_back({
        "generic",
        // Interactions.
        {
            { {  2.0, 4.0,  8.0 }, 0, 1.0 },
            { { -2.1, 3.0, -9.0 }, 3, 50.0 },
        },
        // Coordinates.
        {
            0.0, 1.0, 2.0,
            3.0, 4.0, 5.0,
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
            FF ff;
            ff.interactions = ti.interactions;
            auto ret = ff.computeEnergy(tempC.data());
            return ret;
        };
        const auto calcForce = [&](const VF& c, VF& f) {
            auto tempC = c;
            FF ff;
            ff.interactions = ti.interactions;
            ff.computeForces(tempC.data(), f.data());
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
