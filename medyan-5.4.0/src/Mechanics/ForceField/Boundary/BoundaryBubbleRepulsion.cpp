
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

#include "BoundaryBubbleRepulsion.h"

#include "BoundaryBubbleRepulsionExp.h"
#include "BoundaryElement.h"

#include "Bubble.h"
#include "Bead.h"
#include "Structure/DofSerializer.hpp"

namespace medyan {

template <class BRepulsionInteractionType>
void BoundaryBubbleRepulsion<BRepulsionInteractionType>::vectorize(const FFCoordinateStartingIndex& si, const SimulConfig&) {
    ps_ = si.ps;

    //count interactions
    int nint = 0;
    for (auto be: BoundaryElement::getBoundaryElements())
    {
        for(auto bbIndex : ps_->opBoundaryBubbleNL.value().getNeighbors(be)) {
            ++nint;
        }
    }

    beadSet.assign(n * nint, 0);
    krep.assign(nint, 0);
    slen.assign(nint, 0);
    auto beList = BoundaryElement::getBoundaryElements();

    int nbe = BoundaryElement::getBoundaryElements().size();
    int i = 0;
    int ni = 0;
    int bindex = 0;

    nneighbors.assign(nbe, 0);//stores number of interactions per boundary element.

    int cumnn=0;
    for (i = 0; i < nbe; i++) {

        auto be = BoundaryElement::getBoundaryElements()[i];//beList[i];

        nneighbors[i] = 0;
        auto idx = 0;

        for (auto bbIndex : ps_->opBoundaryBubbleNL.value().getNeighbors(be)) {
            bindex = findBubbleCoordIndex(ps_->bubbles[bbIndex], si);
            beadSet[cumnn+idx] = bindex;
            krep[cumnn+idx] = be->getRepulsionConst();
            slen[cumnn+idx] = be->getScreeningLength();
            idx++;

        }
        nneighbors[i]=idx;
        cumnn+=idx;
    }


}


template <class BRepulsionInteractionType>
floatingpoint BoundaryBubbleRepulsion<BRepulsionInteractionType>::computeEnergy(floatingpoint *coord) {
    
    return _FFType.energy(coord, beadSet.data(), krep.data(), slen.data(), nneighbors.data());
    
}

template <class BRepulsionInteractionType>
void BoundaryBubbleRepulsion<BRepulsionInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
    
    _FFType.forces(coord, f, beadSet.data(), krep.data(), slen.data(), nneighbors.data());
}


///Template specializations
template floatingpoint BoundaryBubbleRepulsion<BoundaryBubbleRepulsionExp>::computeEnergy(floatingpoint *coord);
template void BoundaryBubbleRepulsion<BoundaryBubbleRepulsionExp>::computeForces(floatingpoint *coord, floatingpoint *f);
template void BoundaryBubbleRepulsion<BoundaryBubbleRepulsionExp>::vectorize(const FFCoordinateStartingIndex&, const SimulConfig&);

} // namespace medyan
