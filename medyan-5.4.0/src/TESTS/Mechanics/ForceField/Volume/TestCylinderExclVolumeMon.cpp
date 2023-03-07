#include <numeric> // iota
#include <type_traits> // is_same
#include <vector>

#include <catch2/catch.hpp>

#include "Mechanics/ForceField/Volume/CylinderVolumeMon.hpp"
#include "TESTS/Mechanics/ForceField/TestFFCommon.hpp"


TEST_CASE("Force field: Cylinder excluded volume by monomer", "[ForceField]") {
    // The test case checks whether energy and force are consistent

    using namespace std;
    using namespace medyan;
    using namespace test_ff_common;
    using VF = std::vector< floatingpoint >;

    // Prepare data
    //---------------------------------
    struct TestInput {
        std::string      name;
        std::vector< CylinderVolumeEachInteractionInfo > interactions;
        VF               coord;
    };


    // Multiple test sets
    std::vector< TestInput > testInputs;

    // Input: standard symmetric case
    testInputs.push_back({
        "standard-symmetric",
        {
            CylinderVolumeEachInteractionInfo {
                // Cylinder 1.
                {
                    // coord index.
                    0, 3,
                    // monomer index range.
                    0, 20,
                    // monomer interval and equilibrium length.
                    5, 1.0,
                },
                // Cylinder 2.
                {
                    // coord index.
                    6, 9,
                    // monomer index range.
                    -9, 11,
                    // monomer interval and equilibrium length.
                    4, 1.0,
                },
                // kvol.
                0.25,
            },
        },
        // coordinates.
        {
            10, 0, 0,
            -10, 0, 0,

            0, 10, 5,
            0, -10, 5,
        },
    });

    // Input: general cases
    testInputs.push_back({
        "general",
        {
            CylinderVolumeEachInteractionInfo {
                // Cylinder 1.
                {
                    // coord index.
                    0, 3,
                    // monomer index range.
                    0, 20,
                    // monomer interval and equilibrium length.
                    3, 1.0,
                },
                // Cylinder 2.
                {
                    // coord index.
                    6, 9,
                    // monomer index range.
                    -9, 11,
                    // monomer interval and equilibrium length.
                    4, 1.0,
                },
                // kvol.
                100.0,
            },
            CylinderVolumeEachInteractionInfo {
                // Cylinder 1.
                {
                    // coord index.
                    0, 3,
                    // monomer index range.
                    0, 20,
                    // monomer interval and equilibrium length.
                    3, 1.0,
                },
                // Cylinder 2.
                {
                    // coord index.
                    12, 15,
                    // monomer index range.
                    -25, 25,
                    // monomer interval and equilibrium length.
                    10, 2.0,
                },
                // kvol.
                2e3,
            },
        },
        // coordinates
        {
            10, 0, 0,
            -10, 0, 0,

            0, 10, 5,
            3, -10, 5,

            3, 70, 30,
            2, -20, 0,
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
            floatingpoint ret = 0;
            for(auto& i : ti.interactions) {
                ret += energy(i, c.data());
            }
            return ret;
        };
        const auto calcForce = [&](const VF& c, VF& f) {
            for(auto& i : ti.interactions) {
                force(i, c.data(), f.data());
            }
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
