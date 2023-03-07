
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

#ifndef MEDYAN_AFM_h
#define MEDYAN_AFM_h

#include "common.h"
#include "Chemistry/EmissionAbsorption.hpp"
#include "Structure/BoundaryElementImpl.h"
#include "Structure/Bubble.h"
#include "Util/StableVector.hpp"

namespace medyan {
//FORWARD DECLARATIONS
class Filament;

///A class to represent the structure of a AFM attachemnt system
/*!
 *  The AFM class provides a base structure for a mechanical AFM, which includes
 *  a Bubble (or two bubbles) representing the AFM and [Filaments](@ref Filament) which are
 *  attached to it at the minus end by some sort of potential.
 *
 *  This class only has functions to set and get the constituent Bubble and Filament
 *  and is only mechanically relevant for now, but may be extended to have chemical
 *  properties in the future.
 */
class AFM {
public:
    using CoordinateType = Bubble::CoordinateType;
    using DyRateType = EmissionAbsorptionContainer::DyRateType;

private:
    StableVectorIndex<Bubble> bubbleSysIndex_ { -1 }; // A bubble that physically represents the AFM.
    vector<Filament*> _filaments; ///< An ordered vector of filaments in the AFM bubble
    PlaneBoundaryElement* _pbe;
    
public:

    StableVectorIndex<AFM> sysIndex {};

    // Current compartment index. Updated by position update.
    Index compartmentIndex = 0;

    // Mechanical constants.
    floatingpoint attachmentStretchingK = 0;

    // The pulling force from all filament attachments.
    CoordinateType attachmentForce {};

    //----------------------------------
    // Chemsitry data.
    //----------------------------------
    std::vector< EmissionAbsorptionContainer > vecEmiAbs;

    // Dynamic rate parameters. Only works when dyRateType is set to "force".
    floatingpoint emiForce1 = 0; // Minimum force required to activate emission.
    floatingpoint emiForce2 = 0; // Minimum force that maximizes emission rate.

    ///Constructor
    AFM() = default;
    
    //@{
    ///Setters
    template< typename Context >
    void setBubbleSysIndex(Context& sys, StableVectorIndex<Bubble> index) {
        bubbleSysIndex_ = index;
        sys.bubbles[index].setAFMIndex(sysIndex);
    }
    void setPlaneBoundaryElement(PlaneBoundaryElement* pbe) {
        _pbe = pbe;
    }
    
    void addFilament(Filament* f) {_filaments.push_back(f);}
    //@}
    
    //@{
    ///Getters
    template< typename Context >
    Bubble& getBubble(Context& sys) const { return sys.bubbles[bubbleSysIndex_]; }
    auto getPlaneBoundaryElement() const { return _pbe; }
    const vector<Filament*>& getFilaments() {return _filaments;}
    //@}

    // Setup chemistry.
    template< typename Context >
    void setChemistry(Context& sys, const std::vector<ReactionEmissionAbsorptionSetupInit>& vecEmiAbsInit, const ChemistryData& chemData) {
        auto& grid = *sys.getCompartmentGrid();
        compartmentIndex = grid.getCompartmentIndex(sys.bubbles[bubbleSysIndex_].coord);
        vecEmiAbs = setEmiAbs(sys, compartmentIndex, vecEmiAbsInit, chemData);
    }
    template< typename Context >
    void clearChemistry(Context& sys) {
        vecEmiAbs.clear();
    }

    // Update position.
    template< typename Context >
    void updatePosition(Context& sys) {
        using namespace std;

        if(bubbleSysIndex_.value < 0) {
            return;
        }

        auto& grid = *sys.getCompartmentGrid();
        auto newci = grid.getCompartmentIndex(sys.bubbles[bubbleSysIndex_].coord);
        if(newci != compartmentIndex) {
            // Update the diffusing species in emi-abs reactions.
            for(auto& ea : vecEmiAbs) {
                updateCompartmentEmiAbs(sys, ea, newci);
            }

            // Record new compartment index.
            compartmentIndex = newci;
        }
    }

    // Update reaction rates.
    void updateReactionRates() {
        for(auto& ea : vecEmiAbs) {
            if(ea.dyRateType == DyRateType::force) {
                ea.prEmi->setRateMulFactor(getEmiRateFactorByForce(), ReactionBase::mechanochemical);
                ea.prEmi->updatePropensity();
            }
        }
    }

    // Auxiliary function that computes the force scaling rate based on the force.
    floatingpoint getEmiRateFactorByForce() const {
        auto force = magnitude(attachmentForce);
        if(force < emiForce1) {
            return 0;
        }
        if(force < emiForce2) {
            return (force - emiForce1) / (emiForce2 - emiForce1);
        }
        return 1;
    }

    void printSelf() const;
    
    //GetType implementation just returns zero (no AFM types yet)
    int getType() {return 0;}
    
};

} // namespace medyan

#endif

