
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_AFMAttachment_h
#define MEDYAN_AFMAttachment_h

#include <vector>

#include "common.h"
#include "Mechanics/ForceField/Bubble/AFMAttachmentHarmonic.h"
#include "Mechanics/ForceField/ForceField.h"
#include "Structure/SubSystem.h"
#include "SysParams.h"

namespace medyan {

/// Represents an attachment potential of a AFM.
class AFMAttachment : public ForceField {
public:
    struct PairInteraction {
        Index bubbleCoordIndex = 0;
        Index beadCoordIndex = 0;
        floatingpoint kstr = 0.0;
        floatingpoint radius = 0.0;
    };

    AFMAttachmentHarmonic impl;
private:

    SubSystem* ps_ = nullptr;
    std::vector<PairInteraction> pairInteractions_;
    
public:
    virtual std::string getName() override { return "AFMAttachment"; }
    
    virtual void vectorize(const FFCoordinateStartingIndex&) override;
    
    virtual floatingpoint computeEnergy(floatingpoint *coord) override;
    virtual void computeForces(floatingpoint *coord, floatingpoint *f) override;
    
    virtual void computeAuxParams(SubSystem& sys) override {
        for(auto& afm : sys.afms) {
            auto& bb = afm.getBubble(sys);
            // Index 0:2 for bubble, 3:5 for bead.
            floatingpoint coordset[6] {};
            floatingpoint forceset[6] {};

            coordset[0] = bb.coord[0];
            coordset[1] = bb.coord[1];
            coordset[2] = bb.coord[2];

            for(auto pf : afm.getFilaments()) {
                auto& bead = *pf->getMinusEndCylinder()->getFirstBead();
                coordset[3] = bead.coord[0];
                coordset[4] = bead.coord[1];
                coordset[5] = bead.coord[2];

                impl.forces(
                    coordset, forceset,
                    0, 3,
                    afm.attachmentStretchingK, bb.getRadius()
                );
            }

            afm.attachmentForce[0] = forceset[0];
            afm.attachmentForce[1] = forceset[1];
            afm.attachmentForce[2] = forceset[2];
        }
    }
};

} // namespace medyan

#endif


