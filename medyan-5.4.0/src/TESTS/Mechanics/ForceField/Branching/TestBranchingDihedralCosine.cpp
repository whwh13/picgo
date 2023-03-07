#include <numeric> // iota
#include <type_traits> // is_same
#include <vector>

#include "catch2/catch.hpp"

#include "Mechanics/ForceField/Branching/BranchingDihedralCosine.h"
#include "TESTS/Mechanics/ForceField/TestFFCommon.hpp"

using namespace medyan::test_ff_common;

TEST_CASE("Force field: Branching Dihedral Cosine", "[ForceField]") {
    // The test case checks whether energy and force are consistent

    using namespace std;
    using namespace medyan;
    using VF = vector< floatingpoint >;

    // Prepare data
    //---------------------------------
    VF kdih;
    VF pos;
    VF coords;

    // small dihedral angle
    coords.insert(coords.end(), {
        0.0, -1.0, 0.0,
        0.0, 1.0, 0.0,
        1.0, 0.0, 0.0,
        1.5, 1.0, 0.1
    });
    pos.push_back(0.6);
    kdih.push_back(1.0);

    // 90 deg dihedral angle (<90)
    coords.insert(coords.end(), {
        0.0, -1.0, 0.0,
        0.0, 1.0, 0.0,
        1.0, 0.0, 0.0,
        1.5, 0.1, 1.0
    });
    pos.push_back(0.5);
    kdih.push_back(1.0);

    // 90 deg dihedral angle (>90)
    coords.insert(coords.end(), {
        0.0, -1.0, 0.0,
        0.0, 1.0, 0.0,
        1.0, 0.0, 0.0,
        1.5, -0.1, 1.0
    });
    pos.push_back(0.5);
    kdih.push_back(1.0);

    // large dihedral angle
    coords.insert(coords.end(), {
        0.0, -1.0, 0.0,
        0.0, 1.0, 0.0,
        1.0, 0.0, 0.0,
        1.5, -1.0, 0.1
    });
    pos.push_back(0.4);
    kdih.push_back(1.0);

    // special case dihedral angle
    coords.insert(coords.end(), {
        0.0, -1.0, 0.0,
        0.0, 1.0, 0.0,
        1.0, 0.0, 0.0,
        1.5, -1.0, 0.1
    });
    pos.push_back(1.0);
    kdih.push_back(1.0);

    const auto nint = kdih.size();
    vector< unsigned > beadSet(coords.size() / 3);
    std::iota(beadSet.begin(), beadSet.end(), 0u);
    for(auto& x : beadSet) x *= 3;

    // Prepare functions
    //---------------------------------
    const auto calcEnergy = [&](const VF& c) {
        return BranchingDihedralCosine {}.energy(
            c.data(), nint, beadSet.data(), kdih.data(), pos.data()
        );
    };
    const auto calcForce = [&](const VF& c, VF& f) {
        BranchingDihedralCosine {}.forces(
            c.data(), f.data(), nint, beadSet.data(), kdih.data(), pos.data(), nullptr
        );
    };

    // Prepare test parameters
    //---------------------------------
    const floatingpoint moveMag      = std::is_same< floatingpoint, float >::value ? 5e-3 : 1e-5;

    const floatingpoint diffDeRelEps = std::is_same< floatingpoint, float >::value ? 6e-2 : 5e-4;

    const size_t repsTot     = 10;
    const size_t repsPassReq = 9;

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
