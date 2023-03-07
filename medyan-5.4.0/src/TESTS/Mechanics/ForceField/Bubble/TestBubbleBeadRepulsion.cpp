#include <catch2/catch.hpp>

#include "Mechanics/ForceField/Bubble/BubbleCylinderRepulsion.h"
#include "TESTS/Mechanics/ForceField/TestFFCommon.hpp"

namespace medyan {

TEST_CASE("Bubble bead repulsion", "[ForceField]") {
    using namespace std;
    using VF = vector<floatingpoint>;

    struct TestInput {
        VF coord;
        Size ndof;
        vector<BubbleCylinderRepulsion::PairInteraction> pairs;
    };

    // Add 2 bubbles, 1 movable and 1 fixed.
    // They are close to 3 beads, leading to 6 interactions in total.

    // Bubble information.
    const floatingpoint krep = 1;
    const floatingpoint slen = 1;
    const floatingpoint radii[] = { 0.2, 1.0 };

    TestInput input {
        {
            // Beads.
            -1, 0.1, 0,
            0, -0.1, 0,
            1, 0.2, 0,
            // Bubble.
            0, -1, 1,
            0, 1, 1,     // fixed.
        },
        12, // Number of dof.
        [&] {
            vector<BubbleCylinderRepulsion::PairInteraction> pairs;
            for(int i = 0; i < 2; ++i) {
                for(int j = 0; j < 3; ++j) {
                    pairs.push_back({
                        i*3 + 9, j*3,
                        krep, slen, radii[i]
                    });
                }
            }
            return pairs;
        }(),
    };


    // Prepare functions
    //---------------------------------
    const auto fillDep = [&](VF& c) {
        for(int i = input.ndof; i < c.size(); ++i) {
            c[i] = input.coord[i];
        }
    };
    const auto backProp = [&](const VF& c, VF& f) {
        for(int i = input.ndof; i < f.size(); ++i) {
            f[i] = 0;
        }
    };
    const auto calcEnergy = [&](VF& c) {
        fillDep(c);

        floatingpoint energy = 0;
        for(auto& pair : input.pairs) {
            energy += BubbleBeadRepulsionExp{}.energy(
                c.data(),
                pair.bubbleCoordIndex, pair.beadCoordIndex,
                pair.krep, pair.slen, pair.radius
            );
        }
        return energy;
    };
    const auto calcForce = [&](VF& c, VF& f) {
        fillDep(c);

        for(auto& pair : input.pairs) {
            BubbleBeadRepulsionExp{}.forces(
                c.data(), f.data(),
                pair.bubbleCoordIndex, pair.beadCoordIndex,
                pair.krep, pair.slen, pair.radius
            );
        }
        backProp(c, f);
    };


    {
        using namespace test_ff_common;

        // Prepare test parameters
        //---------------------------------
        const floatingpoint moveMag      = std::is_same< floatingpoint, float >::value ? 5e-3 : 1e-5;
        const floatingpoint diffDeRelEps = std::is_same< floatingpoint, float >::value ? 6e-2 : 5e-4;

        const Size repsTot     = 10;
        const Size repsPassReq = 9;

        // Run the test
        //---------------------------------
        Size repsPass = 0;
        for(Index rep = 1; rep <= repsTot; ++rep) {
            const auto res = testEnergyForceConsistency(input.coord, (int)input.ndof, calcEnergy, calcForce, moveMag, diffDeRelEps);
            if(res.passed) {
                ++repsPass;
            }
            else {
                WARN("(Rep " << rep << '/' << repsTot << " fail) E: " << calcEnergy(input.coord) << " Actual de: " << res.deActual << " Expected de: " << res.deExpected);
            }
        }
        CHECK(repsPass >= repsPassReq);
    }
}

} // namespace medyan
