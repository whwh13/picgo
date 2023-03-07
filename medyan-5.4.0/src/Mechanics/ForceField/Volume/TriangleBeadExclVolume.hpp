#ifndef MEDYAN_Mechanics_ForceField_Volume_TriangleBeadExclVolume_hpp
#define MEDYAN_Mechanics_ForceField_Volume_TriangleBeadExclVolume_hpp

#include <memory> // unique_ptr
#include <string>
#include <vector>

#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Volume/TriangleBeadExclVolRepulsion.hpp"
#include "Structure/NeighborListImpl.h"
#include "Structure/Bead.h"
#include "Structure/SurfaceMesh/Triangle.hpp"
#include "SysParams.h"
#include "Util/Io/Log.hpp"

namespace medyan {
//FORWARD DECLARATIONS
class Cylinder;


/// Represents an excuded volume interaction between a triangle and a cylinder (bead).
class TriangleBeadExclVolume : public ForceField {
    
private:
    TriangleBeadExclVolRepulsion interaction;
    std::unique_ptr<TriangleFilBeadNL> neighborList_;  // Neighbor list of triangle-bead


public:

    // (temp) stores the bead starting index in the coordinate array
    Index       beadStartingIndex = 0;
    SubSystem*  ps = nullptr;

    ///Constructor
    TriangleBeadExclVolume(SubSystem& sys) :
        neighborList_(std::make_unique< TriangleFilBeadNL >(
            SysParams::Mechanics().triangleBeadVolume.cutoff,
            SysParams::Mechanics().triangleBeadVolume.cutoffMech,
            &sys
        ))
    { }

    virtual void vectorize(const FFCoordinateStartingIndex& si) override {
        beadStartingIndex = si.bead;
        ps = si.ps;
    }

    virtual FP computeEnergy(FP* coord) override;
    virtual void computeForces(FP* coord, FP* force) override;

    virtual void computeLoadForces() override;
    virtual void computeLoadForce(SubSystem& sys, Cylinder* c, LoadForceEnd end) const override;

    /// Get the neighbor list for this interaction
    virtual std::vector<NeighborList*> getNeighborLists() override {
        return { neighborList_.get() };
    }

    virtual std::string getName() override {return "TriangleBeadExcludedVolume";}

};

} // namespace medyan

#endif
