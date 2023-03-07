#ifndef MEDYAN_Mechanics_ForceField_Special_BeadMoveConstraintFF_hpp
#define MEDYAN_Mechanics_ForceField_Special_BeadMoveConstraintFF_hpp

#include <vector>

#include "MathFunctions.h"
#include "Mechanics/ForceField/ForceField.h"
#include "Structure/Bead.h"
#include "Structure/DofSerializer.hpp"

namespace medyan {

struct BeadMoveConstraintFF : ForceField {
    // Force constant.
    floatingpoint k = 1;
    // Maximum move range perminimization.
    floatingpoint rmax = 160;

    std::vector<int> beadSet;
    std::vector<Vec<3, floatingpoint>> initPos;

    // Constructor sets the force constants.
    BeadMoveConstraintFF(floatingpoint k, floatingpoint rmax) : k(k), rmax(rmax) {}

    virtual std::string getName() override { return "BeadMoveConstraint"; }
    virtual void vectorize(const FFCoordinateStartingIndex& si) override {
        beadSet.clear();
        beadSet.reserve(Bead::getBeads().size());
        initPos.clear();
        initPos.reserve(Bead::getBeads().size());
        for(auto pb : Bead::getBeads()) {
            beadSet.push_back(findBeadCoordIndex(*pb, si));
            initPos.push_back(pb->coord);
        }
    }

    virtual floatingpoint computeEnergy(floatingpoint* coord) override {
        using namespace mathfunc;

        floatingpoint energy = 0;
        for(int i = 0; i < beadSet.size(); ++i) {
            const auto dist2 = distance2(
                makeRefVec<3, floatingpoint>(coord + beadSet[i]),
                initPos[i]
            );
            energy += fene(FeneDistSq{ dist2 }, k, rmax);
        }

        return energy;
    }

    virtual void computeForces(floatingpoint* coord, floatingpoint* force) override {
        using namespace mathfunc;

        for(int i = 0; i < beadSet.size(); ++i) {
            const auto dist = makeRefVec<3, floatingpoint>(coord + beadSet[i]) - initPos[i];
            const auto dist2 = magnitude2(dist);
            makeRefVec<3, floatingpoint>(force + beadSet[i])
                -= dFeneCoeff(FeneDistSq{ dist2 }, k, rmax) * dist;
        }
    }
};

} // namespace medyan

#endif
