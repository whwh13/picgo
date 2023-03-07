
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

#include "BoundarySurfaceImpl.h"

#include "BoundaryElementImpl.h"

#include "SubSystem.h"

#include "SysParams.h"
#include "MathFunctions.h"

namespace medyan {
using namespace mathfunc;

Plane::Plane(SubSystem* s, vector<floatingpoint> coords, vector<floatingpoint> normal ) :
    BoundarySurface(s, 3), _coords(coords), _normal(normal) {
    
    //Create a plane boundary element
    _boundaryElements.emplace_back(s->addTrackable<PlaneBoundaryElement>
                                   (coords, normal,
                                   SysParams::Boundaries().BoundaryK,
                                   SysParams::Boundaries().BScreenLength));
}

Sphere::Sphere(SubSystem* s, vector<floatingpoint> coords, floatingpoint radius)
    : BoundarySurface(s, 3), _coords(coords) {
    
    //Create a sphere boundary element
    _boundaryElements.emplace_back(s->addTrackable<SphereBoundaryElement>
                                   (coords, radius,
                                   SysParams::Boundaries().BoundaryK,
                                   SysParams::Boundaries().BScreenLength));
    
}

CylinderZ::CylinderZ(SubSystem* s, vector<floatingpoint> coords, floatingpoint radius, floatingpoint height)
    : BoundarySurface(s, 3), _coords(coords) {
    
    //Create a cylindricalZ boundary element
    _boundaryElements.emplace_back(s->addTrackable<CylindricalZBoundaryElement>
                                   (coords, radius, height,
                                   SysParams::Boundaries().BoundaryK,
                                   SysParams::Boundaries().BScreenLength));
}

CylinderXYZ::CylinderXYZ(SubSystem* s, vector<floatingpoint> coords, floatingpoint radius, floatingpoint height)
: BoundarySurface(s, 3), _coords(coords) {
    
    //Create a cylindricalZ boundary element
    _boundaryElements.emplace_back(s->addTrackable<CylindricalXYZBoundaryElement>
                                   (coords, radius, height,
                                    SysParams::Boundaries().BoundaryK,
                                    SysParams::Boundaries().BScreenLength));
}

HalfSphereZ::HalfSphereZ(SubSystem* s, vector<floatingpoint> coords, floatingpoint radius, bool up)
    : BoundarySurface(s, 3), _coords(coords) {
    
    //Create a half sphere Z boundary element
    _boundaryElements.emplace_back(s->addTrackable<HalfSphereZBoundaryElement>
                                   (coords, radius, up,
                                   SysParams::Boundaries().BoundaryK,
                                   SysParams::Boundaries().BScreenLength));
    
}

} // namespace medyan
