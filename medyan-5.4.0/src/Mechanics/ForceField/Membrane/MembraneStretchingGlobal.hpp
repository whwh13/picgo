#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneStretchingGlobal_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneStretchingGlobal_hpp

#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Membrane/MembraneStretchingImpl.hpp"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/Membrane.hpp"

namespace medyan {

struct MembraneStretchingGlobal : public ForceField {

    SubSystem* ps = nullptr;

    virtual void vectorize(const FFCoordinateStartingIndex& si) override {
        ps = si.ps;
    }

    virtual FP computeEnergy(FP* coord) override {
        double en = 0;

        for(auto& m: ps->membranes) {
            const auto& mesh = m.getMesh();

            if(mesh.metaAttribute().vertexSystem == MembraneMeshVertexSystem::general
                && !mesh.metaAttribute().hasLipidReservoir
            ) {
                double totarea = 0;
                for(auto& v : mesh.getVertices()) {
                    totarea += v.attr.gVertex.astar;
                }
                totarea /= 3;

                double enMem = MembraneStretchingHarmonic{}.energy(
                    totarea,
                    m.mMembrane.kArea,
                    m.mMembrane.eqArea
                );

                if(!std::isfinite(enMem)) {
                    return inffp;
                }

                en += enMem;
            }
        }

        return en;
    }

    virtual void computeForces(FP* coord, FP* force) override {

        for (auto& m: ps->membranes) {

            const auto& mesh = m.getMesh();
            assertValidIndexCacheForFF(mesh);

            if(mesh.metaAttribute().vertexSystem == MembraneMeshVertexSystem::general
                && !mesh.metaAttribute().hasLipidReservoir
            ) {

                const auto eqarea = m.mMembrane.eqArea;
                const auto karea = m.mMembrane.kArea;

                double totarea = 0;
                for(auto& v : mesh.getVertices()) {
                    totarea += v.attr.gVertex.astar;
                }
                totarea /= 3;

                for(auto& v : mesh.getVertices()) {
                    const auto& va = v.attr;

                    MembraneStretchingHarmonic{}.forces(
                        force + va.cachedCoordIndex,
                        totarea,
                        va.gVertex.dAstar, // Do not divide by 3 here!
                        karea,
                        eqarea
                    );
                }
            }
        } // end for (m : membranes)
    }

    virtual std::string getName() override { return "MembraneStretchingGlobal"; }

};

} // namespace medyan

#endif
