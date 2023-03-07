
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

#include "BoundaryImpl.h"

#include "BoundarySurfaceImpl.h"
#include "BoundaryElement.h"

#include "Compartment.h"

#include "Controller/GController.h"
#include "SysParams.h"

namespace medyan {
BoundaryCubic::BoundaryCubic(SubSystem* s, vector<BoundaryMove> move)

    : Boundary(s, 3, BoundaryShape::Cube, move){
    
    //Get full system size (want planes to be slightly inside compartment grid)
    floatingpoint zeroX = 25;
    floatingpoint zeroY = 25;
    floatingpoint zeroZ = 25;
    
    floatingpoint sysX = GController::getSize()[0] - zeroX;
    floatingpoint sysY = GController::getSize()[1] - zeroY;
    floatingpoint sysZ = GController::getSize()[2] - zeroZ;
    
    //Create boundary surfaces, add to vector
    //X normal planes
    _boundarySurfaces.emplace_back(new Plane(s, {zeroX, sysY / 2, sysZ / 2}, {1, 0, 0}));
    _boundarySurfaces.emplace_back(new Plane(s, {sysX, sysY / 2, sysZ / 2}, {-1, 0, 0}));
    
    //Y normal planes
    _boundarySurfaces.emplace_back(new Plane(s, {sysX / 2, zeroY, sysZ / 2}, {0, 1, 0}));
    _boundarySurfaces.emplace_back(new Plane(s, {sysX / 2, sysY, sysZ / 2}, {0, -1, 0}));
    
    //Z normal planes
    _boundarySurfaces.emplace_back(new Plane(s, {sysX / 2, sysY / 2, zeroZ}, {0, 0, 1}));
    _boundarySurfaces.emplace_back(new Plane(s, {sysX / 2, sysY / 2, sysZ}, {0, 0, -1}));

    volume();
    
}

bool BoundaryCubic::within(Compartment* C) {
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        //go half a compartment dist in the direction of normal
        auto coordinate = mathfunc::vec2Vector(C->coordinates());
        
        //initial check of coord
        if(be->distance(coordinate) > 0) continue;
        
        //if not, see if any part is in bounds
        auto normal = be->normal(coordinate);
        
        auto effCoordinate =
        {coordinate[0] + SysParams::Geometry().compartmentSizeX * normal[0] / 2,
         coordinate[1] + SysParams::Geometry().compartmentSizeY * normal[1] / 2,
         coordinate[2] + SysParams::Geometry().compartmentSizeZ * normal[2] / 2};
        
        //check if this is within boundary
        if(be->distance(effCoordinate) <= 0)
            return false;
    
    }
    return true;
}


bool BoundaryCubic::within(const vector<floatingpoint>& coordinates) {
    
    // check if all planes return positive distance
    // (means in front of plane, relative to normal)
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        if(be->distance(coordinates) <= 0)
            return false;
    }
    return true;
}

floatingpoint BoundaryCubic::distance(const vector<floatingpoint>& coordinates) {
    
    // loop through, get smallest distance
    floatingpoint smallestDist = numeric_limits<floatingpoint>::infinity();
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        floatingpoint dist = be->distance(coordinates);
        
        if(dist < smallestDist) smallestDist = dist;
        
    }
    return smallestDist;
}

floatingpoint BoundaryCubic::lowerdistance(const vector<floatingpoint>& coordinates) {
    
    // loop through, get smallest distance
    floatingpoint smallestDist = numeric_limits<floatingpoint>::infinity();
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        floatingpoint dist = be->lowerdistance(coordinates);
        
        if(dist < smallestDist) smallestDist = dist;
        
    }
    return smallestDist;
}

//Qin, the same as lowerdistance for now
floatingpoint BoundaryCubic::sidedistance(const vector<floatingpoint>& coordinates) {
    
    // loop through, get smallest distance
    floatingpoint smallestDist = numeric_limits<floatingpoint>::infinity();
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        floatingpoint dist = be->lowerdistance(coordinates);
        
        if(dist < smallestDist) smallestDist = dist;
        
    }
    return smallestDist;
}

vector<floatingpoint> BoundaryCubic::normal(vector<floatingpoint>& coordinates) {
    
    // loop through, get smallest distance
    BoundaryElement* closestPlane = nullptr;
    floatingpoint smallestDist = numeric_limits<floatingpoint>::infinity();
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        floatingpoint dist = be->distance(coordinates);
        
        if(dist < smallestDist) {
            smallestDist = dist;
            closestPlane = be;
        }
        
    }
    //return normal of plane
    return closestPlane->normal(coordinates);
}

floatingpoint BoundaryCubic::getboundaryelementcoord(int bidx) {

    for(auto &bs : _boundarySurfaces) {

        auto be = bs->boundaryElements()[0].get();

        if(be->normal({0,0,0})[0] > 0 && bidx == 0)
            return be->_coords[0];

        else if (be->normal({0,0,0})[0] < 0 && bidx == 1)
            return be->_coords[0];

        else if (be->normal({0,0,0})[1] > 0 && bidx == 2)
            return be->_coords[1];

        else if (be->normal({0,0,0})[1] < 0 && bidx == 3)
            return be->_coords[1];

        else if (be->normal({0,0,0})[2] > 0&& bidx == 4)
            return be->_coords[2];

        else if (be->normal({0,0,0})[2] < 0&&  bidx == 5 )
            return be->_coords[2];

    }
    return std::numeric_limits<floatingpoint>::quiet_NaN();
};

void BoundaryCubic::move(vector<floatingpoint> dist) {

    //do nothing
    if(_move.size() ==0 ) return;
    else if(_move.size() == 1 && _move[0] == BoundaryMove::All) {

        for(auto &bs : _boundarySurfaces) {

            auto be = bs->boundaryElements()[0].get();

            //right
            if(be->normal({0,0,0})[0] < 0)
                be->updateCoords({be->_coords[0] + dist[1], be->_coords[1],
                                  be->_coords[2]});
                //left
            else if (be->normal({0,0,0})[0] > 0)
                be->updateCoords({be->_coords[0] + dist[0], be->_coords[1],
                                  be->_coords[2]});
                //top
            else if (be->normal({0,0,0})[1] < 0)
                be->updateCoords({be->_coords[0], be->_coords[1] + dist[3],
                                  be->_coords[2]});
                //bottom
            else if (be->normal({0,0,0})[1] > 0)
                be->updateCoords({be->_coords[0], be->_coords[1] + dist[2],
                                  be->_coords[2]});
                //front
            else if (be->normal({0,0,0})[2] < 0)
                be->updateCoords({be->_coords[0], be->_coords[1], be->_coords[2] +
                                                                  dist[5]});
                //back
            else if (be->normal({0,0,0})[2] > 0)
                be->updateCoords({be->_coords[0], be->_coords[1], be->_coords[2] +
                                                                  dist[4]});

        }
    }
    else {
        for (auto bm:_move) {
            if (bm == BoundaryMove::None) return;
                //move the left plane
            else if (bm == BoundaryMove::Left) {

                for (auto &bs : _boundarySurfaces) {

                    auto be = bs->boundaryElements()[0].get();

                    if (be->normal({0, 0, 0})[0] > 0) {

                        be->updateCoords(
                                {be->_coords[0] + dist[0], be->_coords[1], be->_coords[2]});
                        return;
                    }
                }
            }
            else if (bm == BoundaryMove::Right) {

                for (auto &bs : _boundarySurfaces) {

                    auto be = bs->boundaryElements()[0].get();

                    if (be->normal({0, 0, 0})[0] < 0) {

                        be->updateCoords(
                                {be->_coords[0] + dist[1], be->_coords[1], be->_coords[2] });
                        return;
                    }
                }
            }
            else if (bm == BoundaryMove::Front) {

                for (auto &bs : _boundarySurfaces) {

                    auto be = bs->boundaryElements()[0].get();

                    if (be->normal({0, 0, 0})[1] > 0) {

                        be->updateCoords(
                                {be->_coords[0], be->_coords[1] + dist[2], be->_coords[2]});
                        return;
                    }
                }
            }
            else if (bm == BoundaryMove::Back) {

                for (auto &bs : _boundarySurfaces) {

                    auto be = bs->boundaryElements()[0].get();

                    if (be->normal({0, 0, 0})[1] < 0) {

                        be->updateCoords(
                                {be->_coords[0], be->_coords[1] + dist[3], be->_coords[2]});
                        return;
                    }
                }
            }
            else if (bm == BoundaryMove::Bottom) {

                for (auto &bs : _boundarySurfaces) {

                    auto be = bs->boundaryElements()[0].get();

                    if (be->normal({0, 0, 0})[2] > 0) {

                        be->updateCoords(
                                {be->_coords[0], be->_coords[1], be->_coords[2] + dist[4]});
                        return;
                    }
                }
            }
            else if (bm == BoundaryMove::Top) {

                for (auto &bs : _boundarySurfaces) {

                    auto be = bs->boundaryElements()[0].get();

                    if (be->normal({0, 0, 0})[2] < 0) {

                        be->updateCoords(
                                {be->_coords[0], be->_coords[1], be->_coords[2] + dist[5]});
                        return;
                    }
                }
            }
        }
    }
}

void BoundaryCubic::volume(){
    // array
    floatingpoint spanarray[3]={0.0,0.0,0.0}; //stores span along x y and z axes.
    for(auto &bs : _boundarySurfaces) {

        auto be = bs->boundaryElements()[0].get();

        auto coord = be->_coords;

        //right
        if(be->normal({0,0,0})[0] > 0)
            spanarray[0] -= coord[0];
        //left
        else if (be->normal({0,0,0})[0] < 0)
            spanarray[0] += coord[0];
        //top
        else if (be->normal({0,0,0})[1] > 0)
            spanarray[1] -= coord[1];
        //bottom
        else if (be->normal({0,0,0})[1] < 0)
            spanarray[1] += coord[1];
        //front
        else if (be->normal({0,0,0})[2] > 0)
            spanarray[2] -= coord[2];
        //back
        else if (be->normal({0,0,0})[2] < 0)
            spanarray[2] += coord[2];

    }

    Boundary::systemvolume = (50 + spanarray[0]) * (50 + spanarray[1]) * (50 +
            spanarray[2]);
}


BoundarySpherical::BoundarySpherical(SubSystem* s, floatingpoint diameter, vector<BoundaryMove> move)

    : Boundary(s, 3, BoundaryShape::Sphere, move) {
    
    floatingpoint sysX = GController::getSize()[0];
    floatingpoint sysY = GController::getSize()[1];
    floatingpoint sysZ = GController::getSize()[2];
        
    _boundarySurfaces.emplace_back(
    new Sphere(s, {sysX / 2, sysY / 2, sysZ / 2}, diameter / 2));
    volume();
}

bool BoundarySpherical::within(Compartment* C) {
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        //go half a compartment dist in the direction of normal
        auto coordinate = mathfunc::vec2Vector(C->coordinates());
        
        //initial check of coord
        if(be->distance(coordinate) > 0) continue;
        
        //if not, see if any part is in bounds
        auto normal = be->normal(coordinate);
        
        auto effCoordinate =
        {coordinate[0] + SysParams::Geometry().compartmentSizeX * normal[0] / 2,
         coordinate[1] + SysParams::Geometry().compartmentSizeY * normal[1] / 2,
         coordinate[2] + SysParams::Geometry().compartmentSizeZ * normal[2] / 2};
        
        //check if this is within boundary
        if(be->distance(effCoordinate) <= 0)
            return false;
        
    }
    return true;
}

bool BoundarySpherical::within(const vector<floatingpoint>& coordinates) {
    
    //check if the boundary element returns a positive distance
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    
    return be->distance(coordinates) > 0;
    
}

floatingpoint BoundarySpherical::distance(const vector<floatingpoint>& coordinates) {
    
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    
    floatingpoint dist = be->distance(coordinates);
    
    if(dist > 0) return dist;
    else return numeric_limits<floatingpoint>::infinity();
}

//the same as distance
floatingpoint BoundarySpherical::lowerdistance(const vector<floatingpoint>& coordinates) {
    
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    
    floatingpoint dist = be->distance(coordinates);
    
    if(dist > 0) return dist;
    else return numeric_limits<floatingpoint>::infinity();
}

floatingpoint BoundarySpherical::sidedistance(const vector<floatingpoint>& coordinates) {
    
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    
    floatingpoint dist = be->distance(coordinates);
    
    if(dist > 0) return dist;
    else return numeric_limits<floatingpoint>::infinity();
}


vector<floatingpoint> BoundarySpherical::normal(vector<floatingpoint>& coordinates) {
    
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    return be->normal(coordinates);
}

void BoundarySpherical::volume(){
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    floatingpoint r[1];
    be->elementeqn(r);
    Boundary::systemvolume = 4/3 * 22/7 * r[0] * r[0] * r[0];
}


BoundaryCapsule::BoundaryCapsule(SubSystem* s, floatingpoint diameter, vector<BoundaryMove> move)

    : Boundary(s, 3, BoundaryShape::Capsule, move) {
    
    floatingpoint sysX = GController::getSize()[0];
    floatingpoint sysY = GController::getSize()[1];
    floatingpoint sysZ = GController::getSize()[2];
    floatingpoint height = sysZ - diameter;
    
    _boundarySurfaces.emplace_back(
    new CylinderZ(s, {sysX / 2, sysY / 2, sysZ / 2}, diameter / 2, height));
    _boundarySurfaces.emplace_back(
    new HalfSphereZ(s, {sysX / 2, sysY / 2, sysZ / 2 + height / 2}, diameter / 2, false));
    _boundarySurfaces.emplace_back(
    new HalfSphereZ(s, {sysX / 2, sysY / 2, sysZ / 2 - height / 2}, diameter / 2, true));
    volume();
}

bool BoundaryCapsule::within(Compartment* C) {
    
    //just calls regular within for now
    return within(mathfunc::vec2Vector(C->coordinates()));
}


bool BoundaryCapsule::within(const vector<floatingpoint>& coordinates) {
    
    //check if the boundary elements return a positive distance
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        floatingpoint dist = be->distance(coordinates);
        if(dist <= 0) return false;
    }
    return true;
}

floatingpoint BoundaryCapsule::distance(const vector<floatingpoint>& coordinates) {
    
    // loop through, get smallest distance
    floatingpoint smallestDist = numeric_limits<floatingpoint>::infinity();
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        floatingpoint dist = be->distance(coordinates);
        
        if(dist < smallestDist) smallestDist = dist;
        
    }
    return smallestDist;
}

//Do not use it for now
floatingpoint BoundaryCapsule::lowerdistance(const vector<floatingpoint>& coordinates) {
    
    // loop through, get smallest distance
    floatingpoint smallestDist = numeric_limits<floatingpoint>::infinity();
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        floatingpoint dist = be->distance(coordinates);
        
        if(dist < smallestDist) smallestDist = dist;
        
    }
    return smallestDist;
}

floatingpoint BoundaryCapsule::sidedistance(const vector<floatingpoint>& coordinates) {
    
    // loop through, get smallest distance
    floatingpoint smallestDist = numeric_limits<floatingpoint>::infinity();
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        floatingpoint dist = be->distance(coordinates);
        
        if(dist < smallestDist) smallestDist = dist;
        
    }
    return smallestDist;
}

void BoundaryCapsule::volume(){
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    floatingpoint dim[2];
    be->elementeqn(dim);
    floatingpoint r = dim[0];
    floatingpoint h = dim[1];
    floatingpoint pi = 22/7;
    Boundary::systemvolume = (h + 4/3 * r) * pi * r * r;
}

BoundaryCylinder::BoundaryCylinder(SubSystem* s, floatingpoint diameter, vector<BoundaryMove> move)

    : Boundary(s, 3, BoundaryShape::Cylinder, move) {
    
    floatingpoint sysX = GController::getSize()[0];
    floatingpoint sysY = GController::getSize()[1];
    floatingpoint sysZ = GController::getSize()[2];
    
    floatingpoint height = sysZ;
    
    _boundarySurfaces.emplace_back(
            new CylinderXYZ(s, {sysX / 2, sysY / 2, sysZ / 2}, diameter / 2, height));
    volume();

}

bool BoundaryCylinder::within(Compartment* C) {

    // project compartment to a 2D cylinderical coordinate
    floatingpoint comX = GController::getCompartmentSize()[0];
    floatingpoint comY = GController::getCompartmentSize()[1];
    auto r = SysParams::Boundaries().diameter / 2;
    auto x = C->coordinates()[0] - r;
    auto y = C->coordinates()[1] - r;

    auto r1 = sqrt((x - comX / 2) * (x - comX / 2) + (y - comY / 2) * (y - comY / 2));
    auto r2 = sqrt((x + comX / 2) * (x + comX / 2) + (y - comY / 2) * (y - comY / 2));
    auto r3 = sqrt((x - comX / 2) * (x - comX / 2) + (y + comY / 2) * (y + comY / 2));
    auto r4 = sqrt((x + comX / 2) * (x + comX / 2) + (y + comY / 2) * (y + comY / 2));
    
    //cout << "x= " << C->coordinates()[0] << " y = " << C->coordinates()[1] << endl;
    
    if (r1 < r || r2 < r || r3 < r || r4 < r) return true;
    else return false;
    
    //return within(C->coordinates());
}


bool BoundaryCylinder::within(const vector<floatingpoint>& coordinates) {
    
    //check if the boundary elements return a positive distance
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();

        floatingpoint dist = be->distance(coordinates);
        if(dist <= 0) return false;
    }
    return true;
}

floatingpoint BoundaryCylinder::distance(const vector<floatingpoint>& coordinates) {
    
    // loop through, get smallest distance
    floatingpoint smallestDist = numeric_limits<floatingpoint>::infinity();
    
    for(auto &bs : _boundarySurfaces) {
        
        auto be = bs->boundaryElements()[0].get();
        
        floatingpoint dist = be->distance(coordinates);
        
        if(dist < smallestDist) smallestDist = dist;
        
    }
    return smallestDist;
}

//Qin
//lower distance should only return the distance between beads and the lower boundary
floatingpoint BoundaryCylinder::lowerdistance(const vector<floatingpoint>& coordinates) {
    
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    
    floatingpoint dist = be->lowerdistance(coordinates);
    
    if(dist > 0) return dist;
    else return numeric_limits<floatingpoint>::infinity();
}

floatingpoint BoundaryCylinder::sidedistance(const vector<floatingpoint>& coordinates) {
    
    auto be = _boundarySurfaces[0]->boundaryElements()[0].get();
    
    floatingpoint dist = be->sidedistance(coordinates);
    
    if(dist > 0) return dist;
    else return numeric_limits<floatingpoint>::infinity();
}

void BoundaryCylinder::volume(){
    double radius = GController::getSize()[0];
    double sysZ = GController::getSize()[2];

    Boundary::systemvolume = radius * radius * 3.14159 * sysZ;
}

} // namespace medyan
