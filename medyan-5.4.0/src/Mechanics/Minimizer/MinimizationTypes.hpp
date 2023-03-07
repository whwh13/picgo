#ifndef MEDYAN_Mechanics_Minimizer_MinimizationTypes_Hpp
#define MEDYAN_Mechanics_Minimizer_MinimizationTypes_Hpp

#include <cstdint>
#include <numeric>

#include "Mechanics/ForceField/Types.hpp"

namespace medyan {
struct MinimizationResult {
    static constexpr std::uint_fast32_t errorLineSearch = 1 << 0;

    EnergyReport energiesBefore;
    EnergyReport energiesAfter;

    std::uint_fast32_t errorFlags = 0;

    // Stats.
    int           numEnergyCall = 0;
    int           energyCallLimit = 0;
    int           numForceCall = 0;
    int           forceCallLimit = 0;

    int           numCGMethodCall = 0;

    floatingpoint gradInfNorm = std::numeric_limits<floatingpoint>::infinity();
    floatingpoint gradTol = 0;
    floatingpoint energyRelChange = std::numeric_limits<floatingpoint>::infinity();
    floatingpoint energyRelTol = 0;


    // Modifiers.
    void updateWithMinimizationResult(const MinimizationResult& res) {
        if(numCGMethodCall == 0) {
            energiesBefore = res.energiesBefore;
        }
        energiesAfter   = res.energiesAfter;

        errorFlags      |= res.errorFlags;

        numEnergyCall   += res.numEnergyCall;
        energyCallLimit += res.energyCallLimit;
        numForceCall    += res.numForceCall;
        forceCallLimit  += res.forceCallLimit;

        ++numCGMethodCall;

        gradInfNorm     = res.gradInfNorm;
        gradTol         = res.gradTol;
        energyRelChange = res.energyRelChange;
        energyRelTol    = res.energyRelTol;
    }

    // Accessors.
    bool success() const { return errorFlags == 0; }
    bool callLimitExceeded() const {
        return
            (energyCallLimit != 0 && numEnergyCall > energyCallLimit) ||
            (forceCallLimit  != 0 && numForceCall  > forceCallLimit);
    }
    bool gradInfNormConverged() const { return gradInfNorm < gradTol; }
    bool energyRelChangeConverged() const { return energyRelChange < energyRelTol; }
    bool converged() const { return gradInfNormConverged() || energyRelChangeConverged(); }
};

} // namespace medyan

#endif
