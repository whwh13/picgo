
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

#ifndef MEDYAN_ForceField_h
#define MEDYAN_ForceField_h

#include <vector>

#include "common.h"
#include "Mechanics/ForceField/Types.hpp"
#include "Structure/SubSystem.h"

namespace medyan {
//FORWARD DECLARATIONS
class Cylinder;
class NeighborList;
class HybridNeighborList;

/// An abstract class to represent various force field calculations
/*!
 *  ForceField is used for force calculations between elements in the SubSystem.
 *  Specific implementations of the ForceField class will have different potentials.
 */
class ForceField {
    
public:
    using LoadForceEnd = ForceFieldTypes::LoadForceEnd;

    virtual ~ForceField() = default;

    /// Get the name of this forcefield
    virtual std::string getName() = 0;
    
    /// Produce a vectorized version of interaction data
    /// Could include constants, positions, neighbors list data, etc
    // This version is deprecated in favor of the version taking simulation configuration.
    virtual void vectorize(const FFCoordinateStartingIndex&) {
        throw std::runtime_error("ForceField::vectorize(const FFCoordinateStartingIndex&) is deprecated. Use ForceField::vectorize(const FFCoordinateStartingIndex&, const SimulConfig&) instead.");
    }
    // Produce a vectorized version of interaction data.
    virtual void vectorize(const FFCoordinateStartingIndex& si, const SimulConfig&) { vectorize(si); }

    /// Compute total energy of this forcefield in the system
    /// @return  the energy value if valid. If an inf or NaN value has been
    /// calculated, return -1.
    virtual floatingpoint computeEnergy(floatingpoint *coord) = 0;

    /// Compute forces of this forcefield in the system. Update Bead
    /// forces accordingly.
    virtual void computeForces(floatingpoint *coord, floatingpoint *f) = 0;

    // Some force fields have the knowledge of computing specific dependent variables necessary for this or other force fields to use.
    virtual void computeDependentCoordinates(floatingpoint* coord) const {}
    // Propagate the forces accumulated on dependent coordinates onto independent coordinates using the chain rule.
    virtual void propagateDependentForces(const floatingpoint* coord, floatingpoint* force) const {}
    // Push forward a tangent vector of independent variable space to all variable space.
    // Note: will record on the same vector, so the original vector must have full size.
    virtual void pushForwardIndependentTangentVector(const floatingpoint* coord, floatingpoint* vec) const {}
    
    ///Compute all load forces on beads in this system.
    ///Updates all Bead's load force components for Reaction updating.
    virtual void computeLoadForces() {}
    // Compute all auxiliary force-related parameters for the system.
    // For example, some forces on certain elements can be used for changing chemical rates.
    // Note:
    // - This function may not access the vectorization information, because it can be called outside energy minimization.
    virtual void computeAuxParams(SubSystem& sys) {}
    // Compute load force for a specific cylinder only.
    virtual void computeLoadForce(SubSystem& sys, Cylinder* pc, LoadForceEnd end) const { }
    
    /// In the case of a calculation error, print the culprit of the FF error.
    /// Typically, will just print the Trackable element where the error came from.
    virtual void whoIsCulprit() {}

    // Force buffer accessor
    const auto& getForceBuffer() const { return forceBuffer_; }

    /// Get all neighbor lists associated with a ForceField
    virtual std::vector<NeighborList*> getNeighborLists() { return {}; }

    // Assign stretchforces for Linker and Motor. Can be extended to other FFs as well.
    // Notes:
    // - Requires valid vectorization.
    // - Shuold only be called once per minimization.
    virtual void assignforcemags() {}

protected:
    std::vector< floatingpoint > forceBuffer_;
};

} // namespace medyan

#endif
