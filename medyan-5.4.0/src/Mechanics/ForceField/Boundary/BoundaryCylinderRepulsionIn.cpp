
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

#include "BoundaryCylinderRepulsionIn.h"

#include <algorithm> // max

#include "BoundaryCylinderRepulsionExpIn.h"
#include "BoundaryElement.h"
#include "BoundaryElementImpl.h"

#include "Bead.h"
#include "Cylinder.h"

#include "MathFunctions.h"
#include "cross_check.h"
#include "CUDAcommon.h"

namespace medyan {
using namespace mathfunc;

template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsionIn<BRepulsionInteractionType>::vectorize(const FFCoordinateStartingIndex& si, const SimulConfig& conf) {
    
    //count interactions
    int nint = 0;
    for (auto be: BoundaryElement::getBoundaryElements())
    {
        
        for(auto &c : _neighborList->getNeighbors(be))
        {
            if(c->isMinusEnd()) nint++;
            nint++;
        }
    }
    CUDAcommon::tmin.numinteractions[9] += nint;
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
        auto nn = _neighborList->getNeighbors(be).size();
        
        nneighbors[i] = 0;
        auto idx = 0;
        
        for (ni = 0; ni < nn; ni++) {
            if(_neighborList->getNeighbors(be)[ni]->isMinusEnd())
            {
                bindex = _neighborList->getNeighbors(be)[ni]->getFirstBead()->getIndex() * 3 + si.bead;
                beadSet[cumnn+idx] = bindex;
                krep[cumnn+idx] = be->getRepulsionConst();
                slen[cumnn+idx] = be->getScreeningLength();
                idx++;
            }
            bindex = _neighborList->getNeighbors(be)[ni]->getSecondBead()->getIndex() * 3 + si.bead;
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
floatingpoint BoundaryCylinderRepulsionIn<BRepulsionInteractionType>::computeEnergy(floatingpoint *coord) {
    floatingpoint U_ii=0.0;
    
    U_ii = _FFType.energy(coord, beadSet.data(), krep.data(), slen.data(), nneighbors.data());

    return U_ii;
}

template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsionIn<BRepulsionInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
    _FFType.forces(coord, f, beadSet.data(), krep.data(), slen.data(), nneighbors.data());
#ifdef DETAILEDOUTPUT
    floatingpoint maxF = 0.0;
    floatingpoint mag = 0.0;
    for(int i = 0; i < CGMethod::N/3; i++) {
        mag = 0.0;
        for(int j = 0; j < 3; j++)
            mag += f[3 * i + j]*f[3 * i + j];
        mag = sqrt(mag);
        //        std::cout<<"SL "<<i<<" "<<mag*mag<<" "<<forceAux[3 * i]<<" "<<forceAux[3 * i + 1]<<" "<<forceAux[3 * i +
        //                2]<<endl;
        if(mag > maxF) maxF = mag;
    }
    std::cout<<"max "<<getName()<<" "<<maxF<<endl;
#endif
}

namespace {
    
    template< typename InteractionType >
    void boundaryCylinderRepulsionLoadForce(
                                            const InteractionType& interaction, floatingpoint kRep, floatingpoint screenLen,
                                            const Bead& bo, Bead& bd, typename BoundaryCylinderRepulsionIn< InteractionType >::LoadForceEnd end,
                                            BoundaryElement* be
                                            ) {
        using LoadForceEnd = typename BoundaryCylinderRepulsionIn< InteractionType >::LoadForceEnd;
        
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
            const auto proj = std::max< floatingpoint >(-dot(vector2Vec< 3 >(be->normal(vec2Vector(newCoord))), dir), 0.0);
            const auto loadForce = interaction.loadForces(be->distance(vec2Vector(newCoord)), kRep, screenLen);
            
            // The load force stored in bead also considers effective monomer size.
            loadForces[i] += proj * loadForce;
        }
        
        //reset lfi
        lfi = 0;
        
    } // void boundaryRepulsionLoadForce(...)
    
} // namespace (anonymous)

template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsionIn<BRepulsionInteractionType>::computeLoadForces() {
    //    std::cout<<"BOUNDARY REPULSION LOAD FORCES DOES NOT USE VECTORIZED FORCES/COORDINATES"<<endl;
    for (auto be: BoundaryElement::getBoundaryElements()) {
        
        for(auto &c : _neighborList->getNeighbors(be)) {
            
            floatingpoint kRep = be->getRepulsionConst();
            floatingpoint screenLength = be->getScreeningLength();
            
            
            //potential acts on second cylinder bead unless this is a minus end
            if(c->isPlusEnd()) {
                boundaryCylinderRepulsionLoadForce(
                                                   _FFType, kRep, screenLength,
                                                   *c->getFirstBead(), *c->getSecondBead(), LoadForceEnd::Plus,
                                                   be
                                                   );
            }
            
            if(c->isMinusEnd()) {
                boundaryCylinderRepulsionLoadForce(
                                                   _FFType, kRep, screenLength,
                                                   *c->getSecondBead(), *c->getFirstBead(), LoadForceEnd::Minus,
                                                   be
                                                   );
            }
            
        }
        
    }
}
template< typename InteractionType >
void BoundaryCylinderRepulsionIn< InteractionType >::computeLoadForce(SubSystem& sys, Cylinder* c, LoadForceEnd end) const {
    for (auto be : BoundaryElement::getBoundaryElements()) {
        
        for(auto cyl : _neighborList->getNeighbors(be)) if(c == cyl) {
            
            floatingpoint kRep = be->getRepulsionConst();
            floatingpoint screenLength = be->getScreeningLength();
            
            boundaryCylinderRepulsionLoadForce(
                                               _FFType, kRep, screenLength,
                                               (end == LoadForceEnd::Plus ? *c->getFirstBead() : *c->getSecondBead()),
                                               (end == LoadForceEnd::Plus ? *c->getSecondBead() : *c->getFirstBead()),
                                               end,
                                               be
                                               );
            
            break;
            
        }
        
    }
}

// Explicit template instantiations
template class BoundaryCylinderRepulsionIn< BoundaryCylinderRepulsionExpIn >;

} // namespace medyan
