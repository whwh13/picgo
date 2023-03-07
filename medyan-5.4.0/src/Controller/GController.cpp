
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

#include "Controller/GController.h"

#include "SubSystem.h"
#include "CompartmentGrid.h"
#include "BoundaryImpl.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include "Rand.h"

namespace medyan {

using namespace mathfunc;
unsigned int GController::getCompartmentID(const vector<floatingpoint> &coords)
{
    //Check if out of bounds
    unsigned int index = 0;
    unsigned int i = 0;
    for(auto x: coords)
    {
        //Flatten the coordinates to 1D, get integer index
        if(i == 0) {
            if(x < 0 || x >= (_compartmentSize[0] * _grid[0]))
                throw OutOfBoundsException();

            index += int(x / _compartmentSize[0]);
        }
        else if(i == 1) {
            if(x < 0 || x >= (_compartmentSize[1] * _grid[1]))
                throw OutOfBoundsException();

            index += int(x / _compartmentSize[1]) * _grid[0];
        }
        else {
            if(x < 0 || x >= (_compartmentSize[2] * _grid[2]))
                throw OutOfBoundsException();

            index += int(x / _compartmentSize[2]) * _grid[0] * _grid[1];
        }
        i++;
    }

    return index;
}

Compartment& GController::getCompartment(const int index){
    try {
        return _compartmentGrid->getCompartment(index);
    }
    catch (const std::exception& e) {
        LOG(ERROR) << "Bad compartment access at index = " << index;
        throw;
    }
}

Compartment& GController::getCompartment(const vector<size_t> &indices)
{
    size_t index = 0;
    size_t i = 0;
    for(auto x: indices)
    {
        //Flatten the indices to 1D
        if(i == 0) {
            if(x >= _grid[0])
                throw OutOfBoundsException();

            index += x;
        }
        else if(i == 1) {
            if(x >= _grid[1])
                throw OutOfBoundsException();

            index += x * _grid[0];
        }
        else {
            if(x >= _grid[2])
                throw OutOfBoundsException();

            index += x * _grid[0] * _grid[1];
        }

        i++;
    }
    try {
        return _compartmentGrid->getCompartment(index);
    }
    catch (const std::exception& e) {
        LOG(ERROR) << "Bad compartment access at index = " << index;
        cout << "Indices = " << indices[0] << " " << indices[1] << " " << indices[2] << endl;
        throw;
    }
}


void GController::generateConnections()
{
    for(size_t i=0U; i<_grid[0]; ++i) {

        for(size_t j=0U; j<_grid[1]; ++j) {

            for(size_t k=0U; k<_grid[2]; ++k)
            {
                vector<size_t> indices{i,j,k};
                Compartment *target = &getCompartment(indices);
                
                medyan::Vec<3, floatingpoint> coordinates =
                   {indices[0] * _compartmentSize[0] + _compartmentSize[0] / 2,
                    indices[1] * _compartmentSize[1] + _compartmentSize[1] / 2,
                    indices[2] * _compartmentSize[2] + _compartmentSize[2] / 2};
                target->setCoordinates(coordinates);
                //Go through all 27 neighbors
                int stencilcount = 0;
                for(int ii: {-1,0,1}){
                    for(int jj: {-1,0,1}){
                        for(int kk: {-1,0,1}){
                            //Consider the target bin itself as a neighbor.
                            stencilcount++;
                            int iprime = i+ii;
                            int jprime = j+jj;
                            int kprime = k+kk;

                            if(iprime<0 || iprime==int(_grid[0]) || jprime<0 ||
                               jprime==int(_grid[1]) || kprime<0 ||
                               kprime==int(_grid[2]))
                                continue;
                            vector<size_t> currentIndices{size_t(iprime), size_t
                                    (jprime), size_t(kprime)};
                            Compartment *neighbor = &getCompartment(currentIndices);
                            //27 enclosing neighbors
                            target->addenclosingNeighbour(neighbor, stencilcount -1);
                            if(ii ==0 && jj == 0 && kk ==0) continue;
                            if(jj==0 && kk>=0 && ii <=0)
                                target->adduniquepermuteNeighbour(neighbor, stencilcount -1);
                            else if(jj==0 && kk==1 && ii == 1)
                                target->adduniquepermuteNeighbour(neighbor, stencilcount -1);
                            else if(jj == 1)
                                target->adduniquepermuteNeighbour(neighbor, stencilcount -1);

                            if(ii != 0 && jprime == j && kprime == k)
                                target->setNeighbor(neighbor->getId(), (ii < 0? 0: 1));
                            else if(jj != 0 && iprime == i && kprime == k)
                                target->setNeighbor(neighbor->getId(), (jj < 0? 2: 3));
                            else if(kk != 0 && iprime == i && jprime == j)
                                target->setNeighbor(neighbor->getId(), (kk < 0? 4: 5));
                        }
                    }
                }
            }
        }
    }

}

CompartmentGrid* GController::initializeGrid(const GeoParams& geoParams) {

    //Initial parameters of system

    _compartmentSize = {geoParams.compartmentSizeX,
                        geoParams.compartmentSizeY,
                        geoParams.compartmentSizeZ};

    _grid = {geoParams.NX,
             geoParams.NY,
             geoParams.NZ};

    _size = {_compartmentSize[0] * _grid[0],
             _compartmentSize[1] * _grid[1],
             _compartmentSize[2] * _grid[2]};

    _centerGrid = {_compartmentSize[0] * _grid[0] / 2,
                   _compartmentSize[1] * _grid[1] / 2,
                   _compartmentSize[2] * _grid[2] / 2};


    //Check that grid and compartmentSize match nDim
    if(
        _grid[0] != 0 && _grid[1] != 0 && _grid[2]!=0 &&
        _compartmentSize[0] != 0 &&
        _compartmentSize[1] != 0 &&
        _compartmentSize[2] != 0){
    }
    else {
        log::error("Grid parameters are invalid. Exiting.");
        throw std::runtime_error("Invalid grid parameters");
    }

    //Set the instance of this grid with given parameters
    _subSystem->compartmentGrid = std::make_unique<CompartmentGrid>(geoParams);
    _compartmentGrid = _subSystem->getCompartmentGrid();

    //Create connections based on dimensionality
    generateConnections();

    return _compartmentGrid;
}

Boundary* GController::initializeBoundary(BoundParams::BoundaryType& BTypes) {

    BoundParams::BoundaryType type;
    vector<BoundaryMove> move;
    for(auto bm:BTypes.boundaryMove){
        if(bm == "NONE") move.push_back(BoundaryMove::None);
        else if(bm == "LEFT") {

            move.push_back(BoundaryMove::Left);
        }
        else if(bm == "RIGHT") {

            move.push_back(BoundaryMove::Right);
        }
        else if(bm == "FRONT") {

            move.push_back(BoundaryMove::Front);
        }
        else if(bm == "BACK") {

            move.push_back(BoundaryMove::Back);
        }
        else if(bm == "BOTTOM") {

            move.push_back(BoundaryMove::Bottom);
        }
        else if(bm == "TOP") {

            move.push_back(BoundaryMove::Top);
        }
        else if(bm == "ALL") {

            move.push_back(BoundaryMove::All);
        }
            //if nothing is specified, don't move boundaries
        else if(bm == "") {
            move.push_back(BoundaryMove::None);
        }
        else {
            log::error("Given boundary movement {} not yet implemented. Exiting.", bm);
            throw std::runtime_error("Invalid boundary movement");
        }
    }

    if(BTypes.boundaryShape == "CUBIC")
        _boundary = new BoundaryCubic(_subSystem, move);

    else if(BTypes.boundaryShape == "SPHERICAL") {

        if(move.size() > 0) {
            if(move[0] != BoundaryMove::None){

                log::error("Moving boundaries for a spherical shape not yet implemented. Exiting.");
                throw std::runtime_error("Invalid boundary movement");
            }
        }

        _boundary = new BoundarySpherical(_subSystem,
                    SysParams::Boundaries().diameter, move);
    }

    else if(BTypes.boundaryShape == "CAPSULE") {

        if(move.size() > 0) {
            if(move[0] != BoundaryMove::None){

                log::error("Moving boundaries for a capsule shape not yet implemented. Exiting.");
                throw std::runtime_error("Invalid boundary movement");
            }
        }
        _boundary = new BoundaryCapsule(_subSystem,
                    SysParams::Boundaries().diameter, move);
    }

    else if(BTypes.boundaryShape == "CYLINDER") {

        if(move.size() > 0) {
            if(move[0] != BoundaryMove::None){

                log::error("Moving boundaries for a cylinder shape not yet implemented. Exiting.");
                throw std::runtime_error("Invalid boundary movement");
            }
        }
        _boundary = new BoundaryCylinder(_subSystem,
                                        SysParams::Boundaries().diameter, move);
    }

    else{
        log::error("Given boundary shape not yet implemented. Exiting.");
        throw std::runtime_error("Invalid boundary shape");
    }

    _subSystem->addBoundary(_boundary);
    return _boundary;
}

void GController::setActiveCompartments() {

    //initialize all compartments equivalent to cproto
    for(auto& C : _compartmentGrid->getCompartments())
        if(_boundary->within(C.get())) C->setAsActive();
}

void GController::findCompartments(
    const medyan::Vec<3, floatingpoint>& coords,
    Compartment* ccheck, floatingpoint dist,
    vector<Compartment*>& compartments
) {

    //base case : if c and ccheck are not within range, return
    if(distance2(coords, ccheck->coordinates()) > (dist * dist) ) return;

    //recursive case, c and ccheck are in range. call for all neighbors
    else {
        //if not already in list, add it
        auto it = find(compartments.begin(), compartments.end(), ccheck);

        if( it == compartments.end()) {
            //add the compartment
            compartments.push_back(ccheck);

            //recursively call for all neighbors
            for(auto cnindex : ccheck->getNeighborIndices()) if(cnindex != -1) {
                findCompartments(coords, &_compartmentGrid->getCompartment(cnindex), dist, compartments);
            }
        }
    }
}

Compartment& GController::getRandomCompartment() {

    //return a compartment that is activated
    while(true) {

        //create a random coordinate
        vector<floatingpoint> coord =
        {_grid[0] * _compartmentSize[0] * Rand::randfloatingpoint(0,0.999),
         _grid[1] * _compartmentSize[1] * Rand::randfloatingpoint(0,0.999),
        _grid[2] * _compartmentSize[2] * Rand::randfloatingpoint(0,0.999)};
        
        Compartment& c = getCompartment(coord);
        if(c.isActivated()) return c;
    }
}



vector<int>    GController::_grid = {};
vector<floatingpoint> GController::_size = {};
vector<floatingpoint> GController::_compartmentSize = {};
vector<floatingpoint> GController::_centerGrid = {};

} // namespace medyan
