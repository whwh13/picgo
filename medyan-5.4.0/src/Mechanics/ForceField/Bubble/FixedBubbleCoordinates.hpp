#ifndef MEDYAN_Mechanics_ForceField_Bubble_FixedBubbleCoordinates_hpp
#define MEDYAN_Mechanics_ForceField_Bubble_FixedBubbleCoordinates_hpp

#include <algorithm>

#include "Mechanics/ForceField/ForceField.h"
#include "Structure/SubSystem.h"

namespace medyan {

// This class serves solely as the purpose of appending fixed bubble coordinates.
// Does not actually perform energy or force calculations.
class FixedBubbleCoordinates : public ForceField {
public:
    // The fixed bubbles are indexed the same as the loop index of those bubbles.
    // The coordinates are ordered in (x1, y1, z1, x2, y2, z2, ...) contiguously.
    std::vector<floatingpoint> fixedBubbleCoords;
    Index fixedBubbleStartingCoordIndex = 0;

    virtual std::string getName() override { return "FixedBubbleCoordinates"; }

    virtual void vectorize(const FFCoordinateStartingIndex& si) override {
        fixedBubbleStartingCoordIndex = si.fixedBubble;

        fixedBubbleCoords.clear();
        for(auto& bb : si.ps->bubbles) {
            if(bb.fixed) {
                fixedBubbleCoords.insert(fixedBubbleCoords.end(), bb.coord.begin(), bb.coord.end());
            }
        }
    }

    virtual void computeDependentCoordinates(floatingpoint* coord) const override {
        std::copy(fixedBubbleCoords.begin(), fixedBubbleCoords.end(), coord + fixedBubbleStartingCoordIndex);
    }
    // propagateDependentForces function is a no-op.
    // pushForwardIndependentTangentVector function is a no-op.

    virtual floatingpoint computeEnergy(floatingpoint* coord) override { return 0; }
    virtual void computeForces(floatingpoint *coord, floatingpoint *force) override {}
};

} // namespace medyan

#endif
