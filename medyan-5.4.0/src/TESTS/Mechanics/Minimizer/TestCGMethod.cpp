#include <algorithm>
#include <array>
#include <string>

#include <catch2/catch.hpp>

#include "Mechanics/Minimizer/CGMethod.hpp"

TEST_CASE("Conjugate gradient minimizer", "[Minimizer]") {
    using namespace std;
    using namespace medyan;

    // Prepare tests on minimization problems.
    //---------------------------------

    // Harmonic in every dimension
    const auto problemHarmonic10 = [](auto&& cg, ConjugateGradientDescentSearch descentSearchMethod, string cgName) {
        INFO("CG name is " << cgName);

        constexpr int n = 10;
        const vector<floatingpoint> optCoord { -500.0, -400.0, -300.0, -200.0, -100.0, 0.0, 100.0, 200.0, 300.0, 400.0 };
        const vector<floatingpoint> weight { 3, 1, 0.1, 0.5, 1, 10, 20, 50, 100, 0.1 };
        const auto energyFunc = [&](const floatingpoint* coord) {
            floatingpoint en = 0;
            for(int i = 0; i < n; ++i) {
                en += 0.5 * weight[i] * (coord[i] - optCoord[i]) * (coord[i] - optCoord[i]);
            }
            return en;
        };
        const auto energyFuncIndividual = [&](const floatingpoint* coord) {
            EnergyReport res;
            res.total = energyFunc(coord);
            return res;
        };
        const auto forceFunc = [&](const floatingpoint* coord, floatingpoint* force, int numVal) {
            for(int i = 0; i < n; ++i) {
                force[i] = -weight[i] * (coord[i] - optCoord[i]);
            }
        };
        // Recovery function is no-op.
        const auto funcRecoveryOnError = [](std::vector<floatingpoint>&, const std::vector<floatingpoint>&) {};

        // Various initial conditions
        vector<vector<floatingpoint>> inits {
            vector<floatingpoint>(800, n),
            vector<floatingpoint>(optCoord.begin(), optCoord.end()),  // already at minimum.
            vector<floatingpoint>{ 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100 },
        };

        // Minimization parameters
        ConjugateGradientParams cgParams;
        cgParams.descentSearchMethod = descentSearchMethod;
        cgParams.gradTol = 0.1;
        cgParams.maxDist = 10;
        cgParams.lambdaMax = 1;
        cgParams.tryToRecoverInLineSearchError = false;

        const auto eachTest = [&](vector<floatingpoint>& init, string name) {
            vector<floatingpoint> force(0.0, n);
            auto minRes = cg.minimize(
                cgParams,
                init, n,
                energyFunc,
                forceFunc,
                energyFuncIndividual,
                funcRecoveryOnError,
                &force
            );
            REQUIRE(minRes.success());

            // Check force below tolerance.
            const int numForceAboveTol = std::count_if(force.begin(), force.end(), [&](floatingpoint d) { return abs(d) > cgParams.gradTol; });
            {
                INFO("In " << name << ", " << numForceAboveTol << " forces are above tolerance " << cgParams.gradTol);
                CHECK(numForceAboveTol == 0);
            }
        };

        // Run tests
        eachTest(inits[0], "input_0");
        eachTest(inits[1], "input_1_already_minimized");
        eachTest(inits[2], "input_2");
    };

    // Run tests on various cg flavors.
    //---------------------------------
    CGMethod cgMethod;

    problemHarmonic10(cgMethod, ConjugateGradientDescentSearch::steepest,       "steepest");
    problemHarmonic10(cgMethod, ConjugateGradientDescentSearch::fletcherRieves, "fletcher_rieves");
    problemHarmonic10(cgMethod, ConjugateGradientDescentSearch::polakRibiere,   "polak_ribiere");

}
