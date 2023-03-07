#ifndef MEDYAN_Mechanics_ForceField_ForceFieldDiagnose_hpp
#define MEDYAN_Mechanics_ForceField_ForceFieldDiagnose_hpp

// Defines functions that diagnose force fields in its current state.
// Note that this is a diagnosis tool for force fields not minimization algorithm, but some information in the current minimization state would be required.

#include <algorithm> // fill
#include <tuple>

#include "Mechanics/ForceField/ForceFieldManager.h"
#include "Util/Math/PolynomialRegression.hpp"

namespace medyan {

// Sometimes, the line search algorithm fails to find a lower energy state given a search direction, along which the energy is promised to decrease by the minimization algorithm. This diagnosis helps to find whether any force field is culprit to this problem.
inline void diagnoseForLineSearch(
    ForceFieldManager& ffm,
    const std::vector<floatingpoint>& coord,
    const std::vector<floatingpoint>& searchDir,
    int numDof,
    floatingpoint currentLineSearchLambdaTol
) {
    // Local storage for coordinates.
    auto coordSearch = coord;
    std::vector<floatingpoint> forceTemp(coord.size());

    // Function that computes all energy values for a range of lambda values around 0.
    // energyFunc: coord -> float
    const auto gatherEnergy = [&](auto&& energyFunc, floatingpoint maxLambdaOneSide, int numSamplesOneSide) {
        std::vector<floatingpoint> lambdas;
        std::vector<floatingpoint> energyRes;
        lambdas.reserve(2 * numSamplesOneSide + 1);
        energyRes.reserve(2 * numSamplesOneSide + 1);

        const auto dl = maxLambdaOneSide / numSamplesOneSide;
        for(int i = -numSamplesOneSide; i <= numSamplesOneSide; ++i) {
            // Compute new coordinates
            const auto lambda = i * dl;
            for(int j = 0; j < numDof; ++j) {
                coordSearch[j] = coord[j] + lambda * searchDir[j];
            }

            // Compute actual energy and gather results.
            lambdas.push_back(lambda);
            energyRes.push_back(energyFunc(coordSearch));
        }

        return std::tuple { lambdas, energyRes };
    };

    // Given force function, compute the expected energy derivative along the search direction at the current coord.
    // forceFunc: (coord, force&) -> void
    const auto energyDerivativeAlongSearchDir = [&](auto&& forceFunc) {
        fill(forceTemp.begin(), forceTemp.end(), 0);
        forceFunc(coord, forceTemp);

        floatingpoint dotRes = 0;
        for(int i = 0; i < numDof; ++i) {
            dotRes -= forceTemp[i] * searchDir[i];
        }

        return dotRes;
    };

    // Diagnosis functions.
    //---------------------------------

    log::info("quad reports 2nd and 1st coefficients in a quad fit. lin reports 1st coefficient in a linear fit.");
    // All force fields together.
    {
        log::info("All force fields included.");
        floatingpoint lambdaMax = currentLineSearchLambdaTol;
        for(int i = 0; i < 4; ++i) {
            // Get energy derivative from sampling.
            int numSamplesOneSide = 50;
            const auto [lambdas, energies] = gatherEnergy(
                [&ffm](std::vector<floatingpoint>& curCoord) {
                    return ffm.computeEnergy(curCoord.data());
                },
                lambdaMax,
                numSamplesOneSide
            );
            const auto reg2 = polynomialRegression<2>(lambdas, energies);
            const auto reg1 = polynomialRegression<1>(lambdas, energies);

            // Print some results.
            log::info("    lambdaMax={:.3e} quad=({:.3e}, {:.3e}) lin={:.3e}", lambdaMax, reg2[2], reg2[1], reg1[1]);

            // Reduce lambdaMax range and try again
            lambdaMax *= 0.2;
        }

        // Get energy derivative from force.
        const auto deExpected = energyDerivativeAlongSearchDir(
            [&](const std::vector<floatingpoint>& curCoord, std::vector<floatingpoint>& force) {
                coordSearch = curCoord;
                ffm.computeForces(coordSearch.data(), force.data(), coordSearch.size());
            }
        );
        log::info("    deExpected={:.3e}", deExpected);
    }

    // Each individual force field.
    for(auto& pff : ffm.getForceFields()) {
        log::info("Force field: {}", pff->getName());
        floatingpoint lambdaMax = currentLineSearchLambdaTol;
        for(int i = 0; i < 4; ++i) {
            // Get energy derivative from sampling.
            int numSamplesOneSide = 50;
            const auto [lambdas, energies] = gatherEnergy(
                [&](std::vector<floatingpoint>& curCoord) {
                    // Precompute dependent variables.
                    ffm.beforeComputeEnergy(curCoord.data());

                    // Then compute energies.
                    return pff->computeEnergy(curCoord.data());
                },
                lambdaMax,
                numSamplesOneSide
            );
            const auto reg2 = polynomialRegression<2>(lambdas, energies);
            const auto reg1 = polynomialRegression<1>(lambdas, energies);

            // Print some results.
            log::info("    lambdaMax={:.3e} quad=({:.3e}, {:.3e}) lin={:.3e}", lambdaMax, reg2[2], reg2[1], reg1[1]);

            // Reduce lambdaMax range and try again
            lambdaMax *= 0.2;
        }

        // Get energy derivative from force.
        const auto deExpected = energyDerivativeAlongSearchDir(
            [&](const std::vector<floatingpoint>& curCoord, std::vector<floatingpoint>& force) {
                coordSearch = curCoord;

                // Precompute dependent coordinates.
                ffm.beforeComputeForce(coordSearch.data());

                // Then compute forces on all possible variables.
                pff->computeForces(coordSearch.data(), force.data());

                // Then propagate forces onto only independent variables.
                for(auto& pff2 : ffm.getForceFields()) {
                    pff2->propagateDependentForces(coordSearch.data(), force.data());
                }
            }
        );
        log::info("    deExpect={:.3e}", deExpected);

    }


}

} // namespace medyan

#endif
