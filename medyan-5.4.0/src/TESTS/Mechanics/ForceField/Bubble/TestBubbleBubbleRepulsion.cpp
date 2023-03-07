#include <catch2/catch.hpp>

#include "Mechanics/ForceField/Bubble/BubbleBubbleRepulsion.h"
#include "TESTS/Mechanics/ForceField/TestFFCommon.hpp"

namespace medyan {

TEST_CASE("Bubble bubble repulsion", "[ForceField]") {
    using namespace std;
    using VF = vector<floatingpoint>;

    struct TestInput {
        VF coord;
        Size ndof;
        vector<BubbleBubbleRepulsion::PairInteraction> pairs;
    };

    // Add 4 bubbles, in which 2 are fixed.
    // 10 interactions in total.

    // Bubble information.
    const floatingpoint krep = 1;
    const floatingpoint slen = 1;
    const floatingpoint radii[] = { 0.2, 3.0, 0.2, 0.2 };

    TestInput input {
        {
            1, 0.5, 0,
            0.5, 0.5, 10,
            0, 0, 0,     // fixed.
            0, 1, 0,     // fixed.
        },
        6, // Number of dof.
        [&] {
            vector<BubbleBubbleRepulsion::PairInteraction> pairs;
            for(int i = 0; i < 4; ++i) {
                for(int j = i+1; j < 4; ++j) {
                    pairs.push_back({
                        i*3, j*3,
                        krep, slen, radii[i], radii[j]
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
            energy += BubbleBubbleRepulsionExp{}.energy(
                c.data(),
                pair.coordIndex1, pair.coordIndex2,
                pair.radius1, pair.radius2,
                pair.krep, pair.slen
            );
        }
        return energy;
    };
    const auto calcForce = [&](VF& c, VF& f) {
        fillDep(c);

        for(auto& pair : input.pairs) {
            BubbleBubbleRepulsionExp{}.forces(
                c.data(), f.data(),
                pair.coordIndex1, pair.coordIndex2,
                pair.radius1, pair.radius2,
                pair.krep, pair.slen
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
