#ifndef MEDYAN_Mechanics_ForceField_VolumeConservation_VolConsrvMembrane_hpp
#define MEDYAN_Mechanics_ForceField_VolumeConservation_VolConsrvMembrane_hpp

#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/VolumeConservation/VolConsrvMembraneHarmonic.hpp"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/FuncMembraneGeo.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"

namespace medyan {

// The force field of elasticity of volume enclosed by the membrane
struct VolumeConservationMembrane : public ForceField {

    VolumeConservationMembraneHarmonic impl;
    SubSystem* ps = nullptr;

    virtual void vectorize(const FFCoordinateStartingIndex& si) override {
        ps = si.ps;
    }

    virtual std::string getName() override { return "VolumeConservation"; }

    virtual FP computeEnergy(FP* coord) override {
        using namespace std;

        double en = 0;

        for(auto& m: ps->membranes) {

            const auto& mesh = m.getMesh();

            const double kBulk = m.mMembrane.kVolume;
            const double eqVolume = m.mMembrane.eqVolume;

            double volume = m.mMembrane.volumeOffset;
            for(const auto& t : mesh.getTriangles())
                volume += t.attr.gTriangle.coneVolume;

            const double enMem = impl.energy(volume, kBulk, eqVolume);

            if(!isfinite(enMem)) {
                log::error("In {} energy calculation, a membrane has energy {}", getName(), enMem);

                return inffp;
            }
            else {
                en += enMem;
            }
            
        }

        return en;
    }

    virtual void computeForces(FP* coord, FP* force) override {
        using namespace std;
        using MT = Membrane::MeshType;

        for (auto& m: ps->membranes) {

            const auto& mesh = m.getMesh();
            medyan::assertValidIndexCacheForFF(mesh);

            const double kBulk = m.mMembrane.kVolume;
            const double eqVolume = m.mMembrane.eqVolume;

            double volume = m.mMembrane.volumeOffset;
            for(const auto& t : mesh.getTriangles()) volume += t.attr.gTriangle.coneVolume;

            const Size numVertices = mesh.getVertices().size();
            for(MT::VertexIndex vi {0}; vi < numVertices; ++vi) {
                const auto& dVolume = mesh.attribute(vi).gVertex.dVolume;

                impl.forces(force + mesh.attribute(vi).cachedCoordIndex, volume, dVolume, kBulk, eqVolume);
            }
        }

    }

};

} // namespace medyan

#endif
