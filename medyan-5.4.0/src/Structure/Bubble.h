
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

#ifndef MEDYAN_Structure_Bubble_h
#define MEDYAN_Structure_Bubble_h

#include <optional>

#include "common.h"
#include "MathFunctions.h"
#include "SysParams.h"
#include "Util/Math/Vec.hpp"
#include "Util/StableVector.hpp"

namespace medyan {

// Forward declarations.
struct MTOC;
struct AFM;

/// Represents a dummy point potential that is involved in mechanical equilibration. This
/// object has no chemical reactions or properties associated with it.

/*!
 *   Bubbles are artificial objects in minimization to physically represent hard spheres.
 *   They contain coordinates as well as force constants for any number of potentials. A 
 *   bubble contains a bead which is the center of the bubble, as well as a radius representing
 *   the physical size of the bubble.
 */

class Bubble {

public:
    using CoordinateType = Vec<3, floatingpoint>;

private:
    double birthTime_ = 0.0;

    int type_;     ///< The type of bubble
    

    floatingpoint _radius;       ///< The radius of this bubble
    floatingpoint _kRepuls;      ///< Repulsion constant for bubble-bubble and bubble-cylinder interactions
    floatingpoint _screenLength; ///< Screening length for a repulsive potential
    floatingpoint _MTOCBendingK; ///< use for MTOC-MT bending force field
    floatingpoint _AFMBendingK; ///< use for AFM-filament bending force field
    
    
    

    // If the bubble is used in MTOC or AFM, store the corresponding indices.
    std::optional<StableVectorIndex<MTOC>> mtocSysIndex_;
    std::optional<StableVectorIndex<AFM>> afmSysIndex_;

public:
    CoordinateType coord {};
    CoordinateType force {};

    // If the bubble is fixed, the coordinates will not be free coordinates in energy minimization.
    bool           fixed = false;

    // Stable index. This will not change during its lifetime.
    // Can be used anywhere.
    // Should never be updated.
    StableVector<Bubble>::Index sysIndex {};
    // Looping index. This value may change for each bubble, but will be contiguous for all movable bubbles.
    // Used as the sequence in mechanical vectorization.
    // Updated during DOF serialization.
    Index          loopIndex = 0;

    /// Default constructor, which only sets birth time.
    Bubble() : birthTime_(tau()) {}

    Bubble(const Bubble&) = default;

    // Mechanical property setters.
    void setMechanicalProperties(const MechParams& mechParams) {
        _kRepuls = mechParams.BubbleK[type_];
        _radius = mechParams.BubbleRadius[type_];
        _screenLength = mechParams.BubbleScreenLength[type_];
        _MTOCBendingK = mechParams.MTOCBendingK.empty() ? 0 : mechParams.MTOCBendingK[type_];
        _AFMBendingK = mechParams.AFMBendingK.empty() ? 0 : mechParams.AFMBendingK[type_];
    }

    //@{
    /// Getters

    floatingpoint getRadius()          const {return _radius;}
    floatingpoint getRepulsionConst()  const {return _kRepuls;}
    floatingpoint getScreeningLength() const {return _screenLength;}
	floatingpoint getMTOCBendingK()    const {return _MTOCBendingK;}
    floatingpoint getAFMBendingK()     const {return _AFMBendingK;}

    auto getId() const { return sysIndex.value; }
    auto getBirthTime() const { return birthTime_; }

    void setType(int type) { type_ = type; }
    int getType() const { return type_; }
    //@}

    void setMTOCIndex(StableVectorIndex<MTOC> mtocSysIndex) { mtocSysIndex_ = mtocSysIndex; }
    bool isMTOC() const { return mtocSysIndex_.has_value(); }
    auto getMTOCIndex() const { return mtocSysIndex_.value(); }

    void setAFMIndex(StableVectorIndex<AFM> afmSysIndex) { afmSysIndex_ = afmSysIndex; }
    bool isAFM() const { return afmSysIndex_.has_value(); }
    auto getAFMIndex() const { return afmSysIndex_.value(); }

    /// Print bubble information
    void printSelf()const;

    // This function can only be called if isAFM() is true.
    template< typename Context >
    void updatePositionManually(Context& sys) {
        //if reaching the desire position
        if(iter > SysParams::Chemistry().StepTotal) {
            iter = 1;
            currentStep++;
        }
        //All position updates will be finished in 1 second
        //Step displacement is 1 /StepTotal
        if(tau() > (currentStep * SysParams::Chemistry().StepTime + iter * 1 / SysParams::Chemistry().StepTotal)){
            floatingpoint step;
            
            if(currentStep > SysParams::Chemistry().IterChange){
                step = SysParams::Chemistry().AFMStep2;
            }
            else{
                step = SysParams::Chemistry().AFMStep1;
            }

            coord[2] += step;

            // Update boundary element coordinate.
            sys.afms[getAFMIndex()].getPlaneBoundaryElement()->updateCoords(mathfunc::vec2Vector(coord));

            iter++;
        }
    }
    double iter = 1;
    int currentStep = 1;

};

} // namespace medyan

#endif
