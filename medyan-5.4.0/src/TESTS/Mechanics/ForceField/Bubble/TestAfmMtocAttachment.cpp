#include <string_view>

#include <catch2/catch.hpp>

#include "Mechanics/ForceField/Bubble/AFMAttachment.h"
#include "Mechanics/ForceField/Bubble/MTOCAttachment.h"
#include "TESTS/Mechanics/ForceField/TestFFCommon.hpp"

namespace medyan {

TEST_CASE("AFM and MTOC attachment", "[ForceField]") {
    using namespace std;
    using VF = vector<floatingpoint>;

    struct TestInput {
        VF coord;
        Size ndof;
        // MTOC attachment interactions take the same form.
        vector<AFMAttachment::PairInteraction> pairs;
    };

    // Add 2 bubbles, 1 movable and 1 fixed.
    // They are each attached to 2 unique beads, leading to 4 interactions in total.

    // Bubble information.
    const floatingpoint kstr = 1;
    const floatingpoint radii[] = { 0.2, 0.6 };

    TestInput input {
        {
            // Beads.
            -1, 0.1, 0,
            -0.5, -0.1, 0,
            0.5, 0.2, 0,
            1, 0.1, 0,
            // Bubble.
            0, -1, 1,
            0, 1, 1,     // fixed.
        },
        13, // Number of dof.
        {
            { 12, 0, kstr, radii[0] },
            { 12, 3, kstr, radii[0] },
            { 15, 6, kstr, radii[1] },
            { 15, 9, kstr, radii[1] },
        },
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
    const auto makeCalcEnergy = [&](auto&& ffImpl) {
        return [&](VF& c) {
            fillDep(c);

            floatingpoint energy = 0;
            for(auto& pair : input.pairs) {
                energy += ffImpl.energy(
                    c.data(),
                    pair.bubbleCoordIndex, pair.beadCoordIndex,
                    pair.kstr, pair.radius
                );
            }
            return energy;
        };
    };
    const auto makeCalcForce = [&](auto&& ffImpl) {
        return [&](VF& c, VF& f) {
            fillDep(c);

            for(auto& pair : input.pairs) {
                ffImpl.forces(
                    c.data(), f.data(),
                    pair.bubbleCoordIndex, pair.beadCoordIndex,
                    pair.kstr, pair.radius
                );
            }
            backProp(c, f);
        };
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
        const auto runTestFor = [&](auto&& ffImpl, std::string_view name) {
            INFO("Testing " << name);
            const auto calcEnergy = makeCalcEnergy(ffImpl);
            const auto calcForce  = makeCalcForce(ffImpl);

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
        };

        runTestFor(AFMAttachmentHarmonic{}, "AFM harmonic");
        runTestFor(MTOCAttachmentHarmonic{}, "MTOC harmonic");
    }
}

} // namespace medyan
