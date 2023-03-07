#include <numeric> // iota
#include <type_traits> // is_same
#include <vector>

#include "catch2/catch.hpp"

#include "Mechanics/ForceField/Boundary/BoundaryCylinderRepulsionExp.h"
#include "Structure/BoundaryElementImpl.h"
#include "TESTS/Mechanics/ForceField/TestFFCommon.hpp"

using namespace medyan::test_ff_common;

TEST_CASE("Force field: Boundary bead repulsion", "[ForceField]") {
    // The test case checks whether energy and force are consistent

    using namespace std;
    using namespace medyan;
    using VF = std::vector< floatingpoint >;

    {
        // Prepare data
        //---------------------------------
        PlaneBoundaryElement element1({1.0, 2.0, 3.0}, {0.0, 0.0, 1.0}, 10, 2.7);
        REQUIRE(BoundaryElement::getBoundaryElements().size() == 1);

        // Test inputs
        std::vector< int > numNeighbors { 3 };
        std::vector< floatingpoint > coords {
            // On plane
            5.0, -100.0, 3.0,

            // Above plane
            0.0, 0.0, 13.0,

            // Below plane
            -1.0, -1e8, -7.0,
        };
        std::vector< floatingpoint > kreps { 1e5, 1e4, 1e3 };
        std::vector< floatingpoint > slens { 1.0, 1.5, 2.0 };
        std::vector< int > beadSet { 0, 3, 6 };


        // Prepare test parameters
        //---------------------------------
        const floatingpoint moveMag      = std::is_same< floatingpoint, float >::value ? 5e-3 : 1e-5;

        const floatingpoint diffDeRelEps = std::is_same< floatingpoint, float >::value ? 6e-2 : 5e-4;

        const size_t repsTot     = 10;
        const size_t repsPassReq = 9;


        // Prepare functions
        //---------------------------------
        const auto calcEnergy = [&](const VF& c) {
            auto tempC = c;
            auto ret = BoundaryCylinderRepulsionExp{}.energy(
                tempC.data(),
                beadSet.data(),
                kreps.data(),
                slens.data(),
                numNeighbors.data()
            );
            return ret;
        };
        const auto calcForce = [&](const VF& c, VF& f) {
            auto tempC = c;
            BoundaryCylinderRepulsionExp{}.forces(
                tempC.data(),
                f.data(),
                beadSet.data(),
                kreps.data(),
                slens.data(),
                numNeighbors.data()
            );
        };


        // Run the test
        //---------------------------------
        size_t repsPass = 0;
        for(size_t rep = 1; rep <= repsTot; ++rep) {
            const auto res = testEnergyForceConsistency(coords, calcEnergy, calcForce, moveMag, diffDeRelEps);
            if(res.passed)
                ++repsPass;
            else
                WARN("(Rep " << rep << '/' << repsTot << " fail) E: " << calcEnergy(coords) << " Actual de: " << res.deActual << " Expected de: " << res.deExpected);
        }
        REQUIRE(repsPass >= repsPassReq);
    }

    // Post condition
    //---------------------------------
    REQUIRE(BoundaryElement::getBoundaryElements().size() == 0);

}
