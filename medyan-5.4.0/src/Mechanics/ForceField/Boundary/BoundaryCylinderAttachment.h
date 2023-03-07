
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

#ifndef MEDYAN_BoundaryCylinderAttachment_h
#define MEDYAN_BoundaryCylinderAttachment_h

#include <vector>

#include "common.h"

#include "Mechanics/ForceField/ForceField.h"
#include "SysParams.h"
#include "Util/Math/Vec.hpp"

namespace medyan {
//FORWARD DECLARATIONS
class BoundaryElement;
class Bead;

/// Represents an attractive interaction between a cylinder and its pin point near a boundary.
/// @note - This BoundaryInteractions implementation really does not involve Boundaries at all.
///         The pinned position of the beads is preset by the special protocol in which Bead
///         elements are chosen with the criteria of being a certain distance away from the Boundary.
///         In reality, any Bead could be pinned at a certain position, but for now associating
///         this potential with a Boundary makes more sense.

template <class BAttachmentInteractionType>
class BoundaryCylinderAttachment : public ForceField {
    
private:
    BAttachmentInteractionType _FFType;
    
    std::vector<int> beadSet;
    
    ///Array describing the constants in calculation
    std::vector<FP> kattr;
    std::vector< Vec< 3, FP > > pins; ///< coordinates of pins for each bead

public:
    
    ///Array describing indexed set of interactions
    ///For filaments, this is a 1-point potential (only the active bead is updated)
    const static int n = 1;
    
    virtual void vectorize(const FFCoordinateStartingIndex&, const SimulConfig&) override;
    
    virtual FP computeEnergy(FP *coord) override;
    //@{
    /// Tepulsive force calculation
    virtual void computeForces(FP *coord, FP *f) override;
    
    virtual std::string getName() override { return "BoundaryCylinderAttachment"; }
};

} // namespace medyan

#endif
