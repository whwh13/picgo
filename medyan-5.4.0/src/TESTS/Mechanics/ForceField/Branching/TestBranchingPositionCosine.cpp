#include <numeric> // iota
#include <type_traits> // is_same
#include <vector>

#include "catch2/catch.hpp"

#include "Mechanics/ForceField/Branching/BranchingPositionCosine.h"
#include "TESTS/Mechanics/ForceField/TestFFCommon.hpp"

using namespace medyan::test_ff_common;

TEST_CASE("Force field: Branching Position Cosine", "[ForceField]") {
    // The test case checks whether energy and force are consistent

    using namespace std;
    using namespace medyan;
    using V3 = Vec< 3, floatingpoint >;
    using VA3 = VecArray< 3, floatingpoint >;

    // Prepare data
    //---------------------------------
    vector< floatingpoint > kpos;
    vector< floatingpoint > pos;
    VA3 coords;

    //
    coords.push_back(V3 { 0.0, -1.0, 0.0 });
    coords.push_back(V3 { 0.0, 1.0, 0.0 });
    coords.push_back(V3 { 1.0, 0.0, 0.0 });
    pos.push_back(0.6);
    kpos.push_back(1.0);

    // 90 deg dihedral angle (<90)
    coords.push_back(V3 { 0.0, -1.0, 0.0 });
    coords.push_back(V3 { 0.0, 1.0, 0.0 });
    coords.push_back(V3 { 1.0, 0.0, 0.0 });
    pos.push_back(0.5);
    kpos.push_back(1.0);

    // 90 deg dihedral angle (>90)
    coords.push_back(V3 { 0.0, -1.0, 0.0 });
    coords.push_back(V3 { 0.0, 1.0, 0.0 });
    coords.push_back(V3 { 1.0, 0.0, 0.0 });
    pos.push_back(0.5);
    kpos.push_back(1.0);

    // large dihedral angle
    coords.push_back(V3 { 0.0, -1.0, 0.0 });
    coords.push_back(V3 { 0.0, 1.0, 0.0 });
    coords.push_back(V3 { 1.0, 0.0, 0.0 });
    pos.push_back(0.4);
    kpos.push_back(1.0);

    // special case dihedral angle
    coords.push_back(V3 { 0.0, -1.0, 0.0 });
    coords.push_back(V3 { 0.0, 1.0, 0.0 });
    coords.push_back(V3 { 1.0, 0.0, 0.0 });
    pos.push_back(1.0);
    kpos.push_back(1.0);

    const auto nint = kpos.size();
    vector< unsigned > beadSet(coords.size());
    std::iota(beadSet.begin(), beadSet.end(), 0u);

    // Prepare functions
    //---------------------------------
    const auto calcEnergy = [&](const VA3& c) {
        return BranchingPositionCosine {}.energy(
            c.data(), beadSet.data(), kpos.data(), pos.data()
        );
    };
    const auto calcForce = [&](const VA3& c, VA3& f) {
        BranchingPositionCosine {}.forces(
            c.data(), f.data(), beadSet.data(), kpos.data(), pos.data(), nullptr
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
