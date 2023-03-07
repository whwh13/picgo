
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

#ifndef MEDYAN_BoundarySurfaceImpl_h
#define MEDYAN_BoundarySurfaceImpl_h

#include <cmath>

#include "common.h"

#include "BoundarySurface.h"

namespace medyan {
/// A simple implementation of the BoundarySurface class.
class Plane: public BoundarySurface {
    
private:
    vector<floatingpoint> _coords; ///< Coordinates of center
    vector<floatingpoint> _normal; ///< Normal vector
    
public:
    
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of plane
    ///@param normal - normal vector to plane
    Plane(SubSystem* s, vector<floatingpoint> coords, vector<floatingpoint> normal);
};


/// A simple implementation of the BoundarySurface class.
/// @note this represents a full sphere, not just a half
class Sphere: public BoundarySurface {
    
private:
    vector<floatingpoint> _coords; ///< Center of sphere
    
public:
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of plane
    ///@param normal - normal vector to sphere
    Sphere(SubSystem* s, vector<floatingpoint> coords, floatingpoint radius);
};


/// A simple implementation of the BoundarySurface class.
class CylinderZ: public BoundarySurface {
    
private:
    vector<floatingpoint> _coords; ///< Center of cylinder
    
public:
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of cylinder
    ///@param normal - normal vector to sphere
    CylinderZ(SubSystem* s, vector<floatingpoint> coords, floatingpoint radius, floatingpoint height);
    
};

/// A simple implementation of the BoundarySurface class.
class HalfSphereZ: public BoundarySurface {
    
private:
    vector<floatingpoint> _coords; ///< Center of half-sphere
    
public:
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of half sphere
    ///@param normal - normal vector to sphere
    HalfSphereZ(SubSystem* s, vector<floatingpoint> coords, floatingpoint radius, bool up);
    
};

/// A simple implementation of the BoundarySurface class.
class CylinderXYZ: public BoundarySurface {
    
private:
    vector<floatingpoint> _coords; ///< Center of cylinder
    
public:
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of cylinder
    ///@param normal - normal vector to sphere
    CylinderXYZ(SubSystem* s, vector<floatingpoint> coords, floatingpoint radius, floatingpoint height);
    
};

} // namespace medyan

#endif
