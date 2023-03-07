#ifndef MEDYAN_Mechanics_ForceField_Membrane_FixedVertexCoordinates_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_FixedVertexCoordinates_hpp

#include <algorithm>

#include "Mechanics/ForceField/ForceField.h"
#include "Structure/SubSystem.h"

namespace medyan {

// This class serves solely as the purpose of appending fixed vertex coordinates.
// Does not actually perform energy or force calculations.
class FixedVertexCoordinates : public ForceField {
public:
    // The pinned vertex coord indices are cached in the mesh vertex attributes.
    // The coordinates are ordered in (x1, y1, z1, x2, y2, z2, ...) contiguously.
    std::vector<FP> fixedVertexCoords;
    Index fixedVertexStartingCoordIndex = 0;

    virtual std::string getName() override { return "FixedVertexCoordinates"; }

    virtual void vectorize(const FFCoordinateStartingIndex& si) override {
        fixedVertexStartingCoordIndex = si.fixedVertex;

        fixedVertexCoords.clear();
        for(auto& m : si.ps->membranes) {
            auto& mesh = m.getMesh();
            for(auto& v : mesh.getVertices()) {
                const auto& vObj = v.attr.vertex(*si.ps);
                if(vObj.pinned) {
                    fixedVertexCoords.insert(fixedVertexCoords.end(), vObj.coord.begin(), vObj.coord.end());
                }
            }
        }
    }

    virtual void computeDependentCoordinates(FP* coord) const override {
        std::copy(fixedVertexCoords.begin(), fixedVertexCoords.end(), coord + fixedVertexStartingCoordIndex);
    }
    // propagateDependentForces function is a no-op.
    // pushForwardIndependentTangentVector function is a no-op.

    virtual FP computeEnergy(FP* coord) override { return 0; }
    virtual void computeForces(FP *coord, FP *force) override {}
};

} // namespace medyan

#endif
