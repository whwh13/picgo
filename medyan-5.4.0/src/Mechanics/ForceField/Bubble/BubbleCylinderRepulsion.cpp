
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BubbleCylinderRepulsion.h"

#include <algorithm> // max

#include "Bubble.h"
#include "Bead.h"
#include "Structure/Filament.h"
#include "Structure/DofSerializer.hpp"

namespace medyan {

void BubbleCylinderRepulsion::vectorize(const FFCoordinateStartingIndex& si) {
    ps_ = si.ps;

    //count interactions
    pairInteractions_.clear();
    for(auto& bb : ps_->bubbles) {
        for(auto pb : ps_->opBubbleBeadNL.value().getNeighbors(bb.sysIndex)) {
            pairInteractions_.push_back({
                findBubbleCoordIndex(bb, si),
                findBeadCoordIndex(*pb, si),
                bb.getRepulsionConst(),
                bb.getScreeningLength(),
                bb.getRadius(),
            });
        }
    }
}


floatingpoint BubbleCylinderRepulsion::computeEnergy(floatingpoint* coord) {
    floatingpoint energy = 0;
    for(auto& pair : pairInteractions_) {
        const auto u = _FFType.energy(
            coord,
            pair.bubbleCoordIndex, pair.beadCoordIndex,
            pair.krep, pair.slen, pair.radius
        );
        energy += u;
    }

    return energy;
}

void BubbleCylinderRepulsion::computeForces(floatingpoint *coord, floatingpoint *force) {
    for(auto& pair : pairInteractions_) {
        _FFType.forces(
            coord, force,
            pair.bubbleCoordIndex, pair.beadCoordIndex,
            pair.krep, pair.slen, pair.radius
        );
    }
}

namespace {

void bubbleCylinderRepulsionLoadForce(
    const BubbleBeadRepulsionExp& interaction,
    floatingpoint          radius,
    floatingpoint          kRep,
    floatingpoint          screenLen,
    const Bead&            bo,
    Bead&                  bd,
    BubbleCylinderRepulsion::LoadForceEnd end,
    const Vec< 3, floatingpoint >& bubbleCoord
) {
    using LoadForceEnd = BubbleCylinderRepulsion::LoadForceEnd;

    auto& loadForces = (end == LoadForceEnd::Plus ? bd.loadForcesP : bd.loadForcesM);
    auto& lfi        = (end == LoadForceEnd::Plus ? bd.lfip        : bd.lfim       );

    // Direction of polymerization
    const auto dir = normalizedVector(bd.coordinate() - bo.coordinate());

    // Array of coordinate values to update
    const auto monSize = SysParams::Geometry().monomerSize   [bd.getType()];
    const auto cylSize = SysParams::Geometry().cylinderNumMon[bd.getType()];

    for (int i = 0; i < cylSize; i++) {

        const auto newCoord = bd.coordinate() + (i * monSize) * dir;

        // Projection magnitude ratio on the direction of the cylinder
        // (Effective monomer size) = (monomer size) * proj
        const auto proj = std::max< floatingpoint >(dot(normalizedVector(bubbleCoord - newCoord), dir), 0.0);
        const auto loadForce = interaction.loadForces(bubbleCoord, bd.coord, radius, kRep, screenLen);

        // The load force stored in bead also considers effective monomer size.
        loadForces[i] += proj * loadForce;
    }

    //reset lfi
    lfi = 0;

} // void bubbleCylinderRepulsionLoadForce(...)

} // namespace (anonymous)

void BubbleCylinderRepulsion::computeLoadForces() {
    
    for (auto& bb : ps_->bubbles) {
        
        //total number of neighbor cylinders
        int cmax = ps_->opBubbleBeadNL.value().getNeighbors(bb.sysIndex).size();
        for(int ni = 0; ni < cmax; ni++){            
            floatingpoint kRep = bb.getRepulsionConst();
            floatingpoint screenLength = bb.getScreeningLength();
            
            floatingpoint radius = bb.getRadius();
            
            Bead* b = ps_->opBubbleBeadNL.value().getNeighbors(bb.sysIndex)[ni];
            auto& fil = *static_cast<Filament*>(b->getParent());

            {
                auto pc = fil.getPlusEndCylinder();
                if(pc->getSecondBead() == b) {
                    bubbleCylinderRepulsionLoadForce(
                        _FFType, radius, kRep, screenLength,
                        *pc->getFirstBead(), *b, LoadForceEnd::Plus,
                        bb.coord
                    );
                }
            }
            {
                auto pc = fil.getMinusEndCylinder();
                if(pc->getFirstBead() == b) {
                    bubbleCylinderRepulsionLoadForce(
                        _FFType, radius, kRep, screenLength,
                        *pc->getSecondBead(), *b, LoadForceEnd::Minus,
                        bb.coord
                    );
                }
            }
        }
    }
}
void BubbleCylinderRepulsion::computeLoadForce(SubSystem& sys, Cylinder* c, LoadForceEnd end) const {
    for (auto& bb : sys.bubbles) {
        
        //total number of neighbor cylinders
        int cmax = sys.opBubbleBeadNL.value().getNeighbors(bb.sysIndex).size();
        for(int ni = 0; ni < cmax; ni++){
            
            floatingpoint kRep = bb.getRepulsionConst();
            floatingpoint screenLength = bb.getScreeningLength();
            
            floatingpoint radius = bb.getRadius();
            
            Bead* pb = sys.opBubbleBeadNL.value().getNeighbors(bb.sysIndex)[ni];
            if(end == LoadForceEnd::Plus ? (c->getSecondBead() == pb) : (c->getFirstBead() == pb)) {
                bubbleCylinderRepulsionLoadForce(
                    _FFType, radius, kRep, screenLength,
                    (end == LoadForceEnd::Plus ? *c->getFirstBead() : *c->getSecondBead()),
                    *pb,
                    end,
                    bb.coord
                );
                break;
            }
        } // End loop neighbor of bb
    } // End loop bubbles
}

} // namespace medyan
