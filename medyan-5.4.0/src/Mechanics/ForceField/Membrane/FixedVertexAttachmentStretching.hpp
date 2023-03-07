#ifndef MEDYAN_Mechanics_ForceField_Membrane_FixedVertexAttachmentStretching_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_FixedVertexAttachmentStretching_hpp

#include "Mechanics/ForceField/ForceField.h"
#include "Structure/SubSystem.h"

namespace medyan {

struct FixedVertexAttachmentStretching : public ForceField {
    struct Interaction {
        Vec<3, FP> ac;     // Attachment coordinate.
        Index vci = 0;     // Starting index of the vertex coordinates.
        FP k = 0;          // Stretching coefficient.
    };

    // Cached data.
    std::vector<Interaction> interactions;

    // ForceField interface.
    std::string getName() override { return "FixedVertexAttachmentStretching"; }
    virtual void vectorize(const FFCoordinateStartingIndex& si, const SimulConfig& conf) override {
        interactions.clear();
        interactions.reserve(si.ps->fixedVertexAttachments.size());

        for(const auto& obja : si.ps->fixedVertexAttachments) {
            const auto& objv = si.ps->vertices[obja.vertexSysIndex];
            const auto& mesh = si.ps->membranes[objv.getParentSysIndex()].getMesh();
            const Index vci = mesh.attribute(mesh.vertexIndex(objv.getTopoIndex())).cachedCoordIndex;
            interactions.push_back({
                obja.coord,
                vci,
                obja.kStretch,
            });
        }
    }

    virtual FP computeEnergy(FP* coord) override {
        FP en = 0;
        for(auto& in : interactions) {
            en += ienergy(coord, in);
        }
        return en;
    }
    virtual void computeForces(FP* coord, FP* force) override {
        for(auto& in : interactions) {
            iforce(coord, force, in);
        }
    }


    // Interaction kernel.
    static FP ienergy(const FP* coord, const Interaction& in) {
        auto vc = makeRefVec<3>(coord + in.vci);
        return in.k / 2 * distance2(in.ac, vc);
    }
    static void iforce(const FP* coord, FP* force, const Interaction& in) {
        auto vc = makeRefVec<3>(coord + in.vci);
        auto vf = makeRefVec<3>(force + in.vci);
        auto f = in.k * (vc - in.ac);
        vf -= f;
    }
};

} // namespace medyan

#endif
