
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

#ifndef MEDYAN_MTOCAttachment_h
#define MEDYAN_MTOCAttachment_h

#include <vector>

#include "common.h"
#include "SysParams.h"
#include "Mechanics/ForceField/Bubble/MTOCAttachmentHarmonic.h"
#include "Mechanics/ForceField/ForceField.h"
#include "Structure/SubSystem.h"

namespace medyan {

/// Represents an attachment potential of a MTOC.
class MTOCAttachment : public ForceField {
public:
    struct PairInteraction {
        Index bubbleCoordIndex = 0;
        Index beadCoordIndex = 0;
        floatingpoint kstr = 0.0;
        floatingpoint radius = 0.0;
    };

    MTOCAttachmentHarmonic impl;
private:

    SubSystem* ps_ = nullptr;
    std::vector<PairInteraction> pairInteractions_;

public:
    virtual std::string getName() override {return "MTOCAttachment";}

    virtual void vectorize(const FFCoordinateStartingIndex&) override;

    virtual floatingpoint computeEnergy(floatingpoint *coord) override;
    virtual void computeForces(floatingpoint *coord, floatingpoint *f) override;
    //virtual void computeForcesAux(double *coord, double *f);
    
    virtual void computeAuxParams(SubSystem& sys) override {
        for(auto& mtoc : sys.mtocs) {
            auto& bb = mtoc.getBubble(sys);
            // Index 0:2 for bubble, 3:5 for bead.
            floatingpoint coordset[6] {};
            floatingpoint forceset[6] {};

            coordset[0] = bb.coord[0];
            coordset[1] = bb.coord[1];
            coordset[2] = bb.coord[2];

            for(auto pf : mtoc.getFilaments()) {
                auto& bead = *pf->getMinusEndCylinder()->getFirstBead();
                coordset[3] = bead.coord[0];
                coordset[4] = bead.coord[1];
                coordset[5] = bead.coord[2];

                impl.forces(
                    coordset, forceset,
                    0, 3,
                    mtoc.attachmentStretchingK, bb.getRadius()
                );
            }

            mtoc.attachmentForce[0] = forceset[0];
            mtoc.attachmentForce[1] = forceset[1];
            mtoc.attachmentForce[2] = forceset[2];
        }
    }
};

} // namespace medyan

#endif
