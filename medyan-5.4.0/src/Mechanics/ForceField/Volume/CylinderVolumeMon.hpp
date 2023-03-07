#ifndef MEDYAN_Mechanics_ForceField_Volume_CylinderVolumeMon_hpp
#define MEDYAN_Mechanics_ForceField_Volume_CylinderVolumeMon_hpp

#include <memory>
#include <vector>

#include "MathFunctions.h"
#include "Mechanics/ForceField/ForceField.h"
#include "Structure/DofSerializer.hpp"
#include "Structure/NeighborListImpl.h"
#include "SubSystem.h"

namespace medyan {

struct CylinderVolumeEachInteractionInfo {
    struct CylinderInfo {
        Index c1Index = 0;
        Index c2Index = 0;
        Index monIndexAtMinus = 0;
        Index monIndexPastPlus = 0;
        Size  monInterval = 1;
        floatingpoint monEqLen = 0;
    };
    CylinderInfo cyl1;
    CylinderInfo cyl2;
    // The force constant of volume exclusion.
    // When it is multiplied by the interaction kernel (1/d^n), create a quantity with energy per square unit equilibrium length.
    floatingpoint kvol = 0;
};

inline floatingpoint energy(const CylinderVolumeEachInteractionInfo& info, const floatingpoint* coords) {
    using namespace mathfunc;
    const auto c11 = makeRefVec<3>(coords + info.cyl1.c1Index);
    const auto c12 = makeRefVec<3>(coords + info.cyl1.c2Index);
    const auto c21 = makeRefVec<3>(coords + info.cyl2.c1Index);
    const auto c22 = makeRefVec<3>(coords + info.cyl2.c2Index);
    const auto numMon1 = info.cyl1.monIndexPastPlus - info.cyl1.monIndexAtMinus;
    const auto numMon2 = info.cyl2.monIndexPastPlus - info.cyl2.monIndexAtMinus;
    const auto m1Start = ceildiv(info.cyl1.monIndexAtMinus, info.cyl1.monInterval) * info.cyl1.monInterval;
    const auto m2Start = ceildiv(info.cyl2.monIndexAtMinus, info.cyl2.monInterval) * info.cyl2.monInterval;
    const auto segEqLen1 = info.cyl1.monInterval * info.cyl1.monEqLen;
    const auto segEqLen2 = info.cyl2.monInterval * info.cyl2.monEqLen;
    const auto coeff = info.kvol * segEqLen1 * segEqLen2;

    floatingpoint en = 0;

    for(Index m1 = m1Start; m1 < info.cyl1.monIndexPastPlus; m1 += info.cyl1.monInterval) {
        const auto loc1 = (m1 - info.cyl1.monIndexAtMinus + 0.5) / numMon1;
        const auto coord1 = c11 + loc1 * (c12 - c11);
        for(Index m2 = m2Start; m2 < info.cyl2.monIndexPastPlus; m2 += info.cyl2.monInterval) {
            const auto loc2 = (m2 - info.cyl2.monIndexAtMinus + 0.5) / numMon2;
            const auto coord2 = c21 + loc2 * (c22 - c21);
            const auto dist2 = distance2(coord1, coord2);
            en += coeff / (dist2 * dist2);
        }
    }

    return en;
}

inline void force(const CylinderVolumeEachInteractionInfo& info, const floatingpoint* coords, floatingpoint* forces) {
    using namespace mathfunc;
    const auto c11 = makeRefVec<3>(coords + info.cyl1.c1Index);
    const auto c12 = makeRefVec<3>(coords + info.cyl1.c2Index);
    const auto c21 = makeRefVec<3>(coords + info.cyl2.c1Index);
    const auto c22 = makeRefVec<3>(coords + info.cyl2.c2Index);
    auto f11 = makeRefVec<3>(forces + info.cyl1.c1Index);
    auto f12 = makeRefVec<3>(forces + info.cyl1.c2Index);
    auto f21 = makeRefVec<3>(forces + info.cyl2.c1Index);
    auto f22 = makeRefVec<3>(forces + info.cyl2.c2Index);
    const auto numMon1 = info.cyl1.monIndexPastPlus - info.cyl1.monIndexAtMinus;
    const auto numMon2 = info.cyl2.monIndexPastPlus - info.cyl2.monIndexAtMinus;
    const auto m1Start = ceildiv(info.cyl1.monIndexAtMinus, info.cyl1.monInterval) * info.cyl1.monInterval;
    const auto m2Start = ceildiv(info.cyl2.monIndexAtMinus, info.cyl2.monInterval) * info.cyl2.monInterval;
    const auto segEqLen1 = info.cyl1.monInterval * info.cyl1.monEqLen;
    const auto segEqLen2 = info.cyl2.monInterval * info.cyl2.monEqLen;
    const auto coeff = info.kvol * segEqLen1 * segEqLen2;
    const auto coeff4 = coeff * 4;

    for(Index m1 = m1Start; m1 < info.cyl1.monIndexPastPlus; m1 += info.cyl1.monInterval) {
        const auto loc1 = (m1 - info.cyl1.monIndexAtMinus + 0.5) / numMon1;
        const auto coord1 = c11 + loc1 * (c12 - c11);
        for(Index m2 = m2Start; m2 < info.cyl2.monIndexPastPlus; m2 += info.cyl2.monInterval) {
            const auto loc2 = (m2 - info.cyl2.monIndexAtMinus + 0.5) / numMon2;
            const auto coord2 = c21 + loc2 * (c22 - c21);
            const auto r12 = coord2 - coord1;
            const auto dist2 = magnitude2(r12);
            const auto grad1 = (coeff4 / (dist2 * dist2 * dist2)) * r12;
            // grad2 = -grad1
            f11 -= (1 - loc1) * grad1;
            f12 -= loc1 * grad1;
            f21 -= -(1 - loc2) * grad1;
            f22 -= -loc2 * grad1;
        }
    }
}

struct CylinderVolumeMon : ForceField {

    // Context pointer and neighbor lists, initialized at construction.
    std::unique_ptr< CylinderCylinderNL > neighborList;
    std::vector<short> hnlIdVec;

    // Volume exclusion vectorization information.
    SubSystem* ps = nullptr;
    std::vector< CylinderVolumeEachInteractionInfo > interactions;
    
    // Initialize the volume force field.
    // If phnl is nullptr, the original hybrid neighbor list is initialized internally and used. Otherwise, the hybrid neighbor list pointed to by phnl is used.
    CylinderVolumeMon(HybridCylinderCylinderNL* phnl, const SimulConfig& conf) {
        if(phnl) {
            // Set hybrid neighbor list.
            auto numFilTypes = conf.geoParams.cylinderNumMon.size();
            for(int i = 0; i < numFilTypes; ++i) {
                for(int j = 0; j < numFilTypes; ++j) {
                    hnlIdVec.push_back(phnl->setneighborsearchparameters(i, j, true, false, conf.mechParams.VolumeCutoff, 0.0));
                }
            }
        }
        else {
            // Set original neighbor list.
            neighborList = std::make_unique<CylinderCylinderNL>(conf.mechParams.VolumeCutoff);
        }
    }

    virtual void vectorize(const FFCoordinateStartingIndex& si) override {
        // TODO replace the following with function parameters.
        const auto& geoParams = SysParams::Geometry();
        const auto& mechParams = SysParams::Mechanics();

        ps = si.ps;
        interactions.clear();

        const auto judgeNeighbor = [](const Cylinder& c1, const Cylinder& c2) {
            return !(c1.getBranchingCylinder() == &c2 || c2.getBranchingCylinder() == &c1);
        };
        const auto judgeAndAddInteraction = [&, this](const Cylinder& c1, const Cylinder& c2) {
            if(!judgeNeighbor(c1, c2)) return;

            auto& info = interactions.emplace_back();

            // Set indices.
            info.cyl1.c1Index = findBeadCoordIndex(*c1.getFirstBead(), si);
            info.cyl1.c2Index = findBeadCoordIndex(*c1.getSecondBead(), si);
            info.cyl2.c1Index = findBeadCoordIndex(*c2.getFirstBead(), si);
            info.cyl2.c2Index = findBeadCoordIndex(*c2.getSecondBead(), si);

            // Set monomers.
            info.cyl1.monIndexAtMinus = c1.getFirstBead()->monomerSerial;
            info.cyl1.monIndexPastPlus = c1.getSecondBead()->monomerSerial;
            info.cyl1.monEqLen = geoParams.monomerSize[c1.getType()];
            info.cyl2.monIndexAtMinus = c2.getFirstBead()->monomerSerial;
            info.cyl2.monIndexPastPlus = c2.getSecondBead()->monomerSerial;
            info.cyl2.monEqLen = geoParams.monomerSize[c2.getType()];

            info.cyl1.monInterval = mechParams.volumeExclusionMonomerInterval[c1.getType()];
            info.cyl2.monInterval = mechParams.volumeExclusionMonomerInterval[c2.getType()];

            // Set force constant.
            info.kvol = std::max(
                c1.getMCylinder()->getExVolConst(),
                c2.getMCylinder()->getExVolConst()
            );
        };

        for(auto pc : Cylinder::getCylinders()) {
            if(neighborList) {
                // Use original neighbor list.
                auto neighbors = neighborList->getNeighbors(pc);
                for(auto pnc : neighbors) {
                    judgeAndAddInteraction(*pc, *pnc);
                }
            }
            else {
                // Use hybrid neighbor list.
                auto phnl = ps->getHNeighborList();
                for(auto id : hnlIdVec) {
                    auto neighbors = phnl->getNeighborsstencil(id, pc);
                    for(auto pnc : neighbors) {
                        judgeAndAddInteraction(*pc, *pnc);
                    }
                }
            }
        }
    }

    virtual std::string getName() override { return "CylinderExcludedVolumeMonomer"; }

    virtual floatingpoint computeEnergy(floatingpoint *coords) override {
        floatingpoint en = 0;
        for(auto& i : interactions) {
            en += energy(i, coords);
        }
        return en;
    }
    virtual void computeForces(floatingpoint *coords, floatingpoint *forces) override {
        for(auto& i : interactions) {
            force(i, coords, forces);
        }
    }

    virtual vector<NeighborList*> getNeighborLists() override {
        return { neighborList.get() };
    }

};

} // namespace medyan

#endif
