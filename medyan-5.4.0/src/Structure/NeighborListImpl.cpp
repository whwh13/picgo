
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

#include "NeighborListImpl.h"

#include "Bead.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bubble.h"
#include "BoundaryElement.h"
#include "Structure/SurfaceMesh/FuncMembraneGeo.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Triangle.hpp"

#include "Controller/GController.h"
#include "MathFunctions.h"
#include "CUDAcommon.h"
#include "NeighborListImplCUDA.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
#include "Util/Math/TriangleArithmetics.hpp"

namespace medyan {
using namespace mathfunc;
#ifdef NLSTENCILLIST
void CylinderCylinderNL::updateallcylinderstobin() {
    for(auto cyl:Cylinder::getCylinders())
        updatebin(cyl);
}

void CylinderCylinderNL::assignallcylinderstobin() {
    for(auto cyl:Cylinder::getCylinders())
        assignbin(cyl);
}

void CylinderCylinderNL::assignbin(Cylinder* cyl){
    Bin* _bin;
    try {_bin = getBin(cyl->coordinate);}
    catch (exception& e) {
        cout << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    _bin->addCylinder(cyl);
    cyl->_binvec.push_back(_bin);
}

void CylinderCylinderNL::unassignbin(Cylinder* cyl, Bin* bin){
    bin->removeCylinder(cyl);
}

void CylinderCylinderNL::updatebin(Cylinder *cyl){
    Bin* _bin;
//    std::cout<<coordinate[0]<<" "<<coordinate[1]<<" "<<coordinate[2]<<endl;
    try {_bin = getBin(cyl->coordinate);}
    catch (exception& e) {
        cout << e.what();
        cyl->printSelf();
        exit(EXIT_FAILURE);
    }

    if(_bin != cyl->_binvec.at(_ID)) {
        auto oldBin = cyl->_binvec.at(_ID);
        auto newBin = _bin;

        //remove from old compartment, add to new
        oldBin->removeCylinder(cyl);
        cyl->_binvec.at(_ID) = newBin;
        _bin->addCylinder(cyl);
    }
}

void CylinderCylinderNL::generateConnections() {
    for(size_t i=0U; i<_grid[0]; ++i) {

        for(size_t j=0U; j<_grid[1]; ++j) {

            for(size_t k=0U; k<_grid[2]; ++k) {
                vector<size_t> indices{i,j,k};
                Bin *target = getBin(indices);//defined in this file.

                medyan::Vec<3, floatingpoint> coordinates =
                        {indices[0] * _binSize[0] + _binSize[0] / 2,
                         indices[1] * _binSize[1] + _binSize[1] / 2,
                         indices[2] * _binSize[2] + _binSize[2] / 2};
                target->setCoordinates(coordinates);
                int stencilcount = 0;

                //Go through all neighbors to get the neighbors list
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
                                Bin *neighbor = getBin(currentIndices);
                                target->addNeighbour(neighbor);
                                target->stencilID.push_back(stencilcount-1);
                        }
                    }
                }
            }
        }
    }

}

void CylinderCylinderNL::initializeBinGrid() {

//    //Initial parameters of system
    floatingpoint searchdist = 1.125 * (_rMax);
    _binSize = {searchdist, searchdist, searchdist};
    {
        _size.push_back(int(SysParams::Geometry().NX * SysParams::Geometry()
                .compartmentSizeX));
        if( (_size[0]) % int(_binSize[0]) ==0)
            _grid.push_back(_size[0]/_binSize[0]);
        else
            _grid.push_back(_size[0]/_binSize[0] + 1);

        _size.push_back(int(SysParams::Geometry().NY * SysParams::Geometry()
                .compartmentSizeY));
        if( (_size[1]) % int(_binSize[1]) ==0)
            _grid.push_back(_size[1]/_binSize[1]);
        else
            _grid.push_back(_size[1]/_binSize[1] + 1);

        _size.push_back(int(SysParams::Geometry().NZ * SysParams::Geometry()
                .compartmentSizeZ));
        if( (_size[2]) % int(_binSize[2]) ==0)
            _grid.push_back(_size[2]/_binSize[2]);
        else
            _grid.push_back(_size[2]/_binSize[2] + 1);
    }

    //Check that grid and compartmentSize match nDim
    if(
        _grid[0] != 0 && _grid[1] != 0 && _grid[2]!=0 &&
        _binSize[0] != 0 &&
        _binSize[1] != 0 &&
        _binSize[2] != 0){
    }
    else {
        cout << "Bin parameters for CylinderCylinderNeighborLists are invalid. Exiting." <<
             endl;
        exit(EXIT_FAILURE);
    }
    int size = 1;
    for(auto x: _grid) {
        if(x != 0) size*=x;
    }
    //Set the instance of this grid with given parameters
    _binGrid = new BinGrid(size, _ID, _binSize);
    //Create connections based on dimensionality
    generateConnections();
#ifdef CUDAACCL_NLS
    binGridv = new bin[size];
    auto binvec = _binGrid->getBins();
    for(int i = 0; i < size; i ++){
        auto bcoord = binvec[i]->coordinates();
        binGridv[i].binID = binvec[i]->_ID;
        for(int dim =0;dim<3;dim++)
            binGridv[i].bincoord[dim] = bcoord.at(dim);
        auto binneighbors = binvec[i]->getNeighbours();
        int count = 0;
        for(auto n:binneighbors) {
            binGridv[i].neighbors[count] = n->_ID;
            binGridv[i].binstencilID[count] = binvec[i]->stencilID[count];
            count++;
        }
    }
    if(stream_NL == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream_NL));
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_binGrid, size * sizeof(bin)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_binGrid, binGridv, size * sizeof(bin),
                                       cudaMemcpyHostToDevice));
#endif
}

//You need a vector of all grids so you can loop through and update respective coordinates.
Bin* CylinderCylinderNL::getBin(const vector<floatingpoint> &coords) {
    //Check if out of bounds
    size_t index = 0;
    size_t i = 0;
    for(auto x: coords)
    {
        //Flatten the coordinates to 1D, get integer index
        if(i == 0) {
            if(x < 0 || x >= (_binSize[0] * _grid[0]))
                throw OutOfBoundsException();

            index += int(x / _binSize[0]);
        }
        else if(i == 1) {
            if(x < 0 || x >= (_binSize[1] * _grid[1]))
                throw OutOfBoundsException();

            index += int(x / _binSize[1]) * _grid[0];
        }
        else {
            if(x < 0 || x >= (_binSize[2] * _grid[2]))
                throw OutOfBoundsException();

            index += int(x / _binSize[2]) * _grid[0] * _grid[1];
        }
        i++;
    }

    try {
        return _binGrid->getBin(index);
    }
    catch (exception& e){
        cout << "Bad bin access at..." << endl;
        cout << "Bin index = " << index << endl;
        cout << "Coords = " << coords[0] << " " << coords[1] << " " << coords[2] << endl;
        throw NaNCoordinateException();
    }
}

Bin* CylinderCylinderNL::getBin(const vector<size_t> &indices) {
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
        return _binGrid->getBin(index);
    }
    catch (exception& e){
        cout << "Bad Bin access at..." << endl;
        cout << "Bin index = " << index << endl;
        cout << "Indices = " << indices[0] << " " << indices[1] << " " << indices[2] << endl;
        throw NaNCoordinateException();
    }
}

void CylinderCylinderNL::updateNeighborsbin(Cylinder* cylinder, bool runtime){
    //clear existing
    _list4mbin[cylinder].clear();
    auto binvec = cylinder->_binvec;//The different bins that this cylinder belongs to.
    if(binvec.size()<=_ID)
        assignbin(cylinder);
    binvec = cylinder->_binvec;
    auto parentbin =  binvec.at(_ID);
    vector<Bin*> _neighboringBins = binvec.at(_ID)//Get the bin that belongs to the
                    // current binGrid of interest for this NL.
            ->getNeighbours();
    //
//    int ncyls2 = 0;
//    int tcyl2 = 0;
    int nbincount = 0;
    auto nbinstencil = parentbin->stencilID;// A standard templated numbering of
    // neighboring bins is implemented i.e. based on position w.r.t. bin of interest,
    // neighboring bins are given a particular ID.nbinstencil stores the set of such
    // neighbors that is close to bin of interest. Bins close to the boundary will have
    // < 27 elements in the stencilID vector.
    for(auto &bin : _neighboringBins) {
//        bin->vectorize();
        bool isbinneeded = _binGrid->iswithincutoff(cylinder->coordinate,
                                                     parentbin->coordinates(),
                                                     nbinstencil.at(nbincount), _rMax);
        nbincount++;
        if(isbinneeded) {
            for (auto &ncylinder : bin->getCylinders()) {
                bool checkstatus = false;
                if ((cylinder->getType() == NLcyltypes[0] &&
                     ncylinder->getType() == NLcyltypes[1])||
                     (cylinder->getType() == NLcyltypes[1] &&
                     ncylinder->getType() == NLcyltypes[0])) {
                    checkstatus = true;
                }
                if (checkstatus) {
                    //Don't add the same cylinder!
                    if (cylinder == ncylinder) continue;

                    //Dont add if ID is more than cylinder for half-list
                    if (!_full && cylinder->getId() <= ncylinder->getId())
                        continue;

                    //Don't add if belonging to same parent
                    if (cylinder->getParent() == ncylinder->getParent()) {

                        //if not cross filament, check if not neighboring
                        auto dist = fabs(cylinder->getPosition() -
                                         ncylinder->getPosition());
                        if (dist <= ChemParams::minCylinderDistanceSameFilament) continue;
                    }
                    //Dont add if not within range
                    floatingpoint dist = twoPointDistance(cylinder->coordinate,
                                                   ncylinder->coordinate);
                    if (dist > _rMax || dist < _rMin) continue;
                    //If we got through all of this, add it!
                    _list4mbin[cylinder].push_back(ncylinder);
//                    ncyls2++;

                    //if runtime, add to other list as well if full
                    if (runtime && _full) {
                        _list4mbin[ncylinder].push_back(cylinder);
//                        ncyls2++;
                    }
                }
            }
        }
    }
}

vector<Cylinder*> CylinderCylinderNL::getNeighborsstencil(Cylinder* cylinder) {


    return _list4mbin[cylinder];
}
#endif

void CylinderCylinderNL::updateNeighbors(Cylinder* cylinder, bool runtime) {

    //clear existing
    _list[cylinder].clear();

    //Find surrounding compartments (For now its conservative)
    vector<Compartment*> compartments;
    auto searchDist = SysParams::Geometry().largestCompartmentSide;

    GController::findCompartments(mathfunc::vector2Vec<3>(cylinder->coordinate),
                                  cylinder->getCompartment(),
                                  searchDist + _rMax, compartments);
//    std::cout<<" neighboring cmps "<<compartments.size()<<endl;
//    int tcyl = 0;
    for(auto &comp : compartments) {
//        tcyl += comp->getCylinders().size();
        for(auto &ncylinder : comp->getCylinders()) {
            //Don't add the same cylinder!
            if(cylinder == ncylinder) continue;

            //Dont add if ID is more than cylinder for half-list
            if(!_full && cylinder->getId() <= ncylinder->getId())
                continue;

            //Don't add if belonging to same parent
            if(cylinder->getParent() == ncylinder->getParent()) {

                //if not cross filament, check if not neighboring
                auto dist = fabs(cylinder->getPosition() -
                                 ncylinder->getPosition());
                if(dist <= ChemParams::minCylinderDistanceSameFilament) continue;
            }
            //Dont add if not within range
            floatingpoint distsq = twoPointDistancesquared(cylinder->coordinate,
                                           ncylinder->coordinate);
            if(distsq > (_rMax * _rMax) || distsq < (_rMin * _rMin)) continue;

            //If we got through all of this, add it!
            _list[cylinder].push_back(ncylinder);

            //if runtime, add to other list as well if full
            if(runtime && _full)
                _list[ncylinder].push_back(cylinder);
        }
    }
//    std::cout<<"Total cyls "<<tcyl<<endl;
}

void CylinderCylinderNL::addNeighbor(Neighbor* n) {

    //return if not a cylinder!
    Cylinder* cylinder;
    if(!(cylinder = dynamic_cast<Cylinder*>(n))) return;

    //update neighbors
#ifdef NLORIGINAL
    updateNeighbors(cylinder, true);
#endif
#ifdef NLSTENCILLIST
    updateNeighborsbin(cylinder, true);
#endif
}

void CylinderCylinderNL::removeNeighbor(Neighbor* n) {

    Cylinder* cylinder;
    if(!(cylinder = dynamic_cast<Cylinder*>(n))) return;
#ifdef NLORIGINAL
    _list.erase(cylinder);

    //remove from other lists
    for(auto it = _list.begin(); it != _list.end(); it++) {

        auto cit = find(it->second.begin(), it->second.end(), cylinder);
        if(cit != it->second.end()) it->second.erase(cit);
    }
#endif
#ifdef NLSTENCILLIST
//    std::cout<<"Removing neighbors of "<<cylinder<<" from NL "<<_ID<<endl;
    //Remove from NeighborList
    _list4mbin.erase(cylinder);
    //Remove from bin
    Bin* bin = cylinder->_binvec.at(_ID);
    unassignbin(cylinder, bin);
    //remove from other lists
//    std::cout<<"Removed from cylinders ";
    for(auto it = _list4mbin.begin(); it != _list4mbin.end(); it++) {
        auto cit = find(it->second.begin(), it->second.end(), cylinder);
        {
            if (cit != it->second.end()) {
                it->second.erase(cit);
//                std::cout<<it->first<<" ";
            }
        }
    }
//    std::cout<<endl;
#endif
}

void CylinderCylinderNL::reset() {
//    std::cout<<"Total number of bins "<< _binGrid->getBins().size()<<endl;
//    for(auto bin:_binGrid->getBins()){
//        std::cout<<bin->getCylinders().size()<<" ";
//    }
//    std::cout<<endl;
#ifdef CUDAACCL_NL
    _list.clear();

//    pair_cIndex_cmp.clear();
    nint = 0;
    pair_cIndex_cnIndex.clear();
    //1. Get total interactions
    for(auto cylinder: Cylinder::getCylinders()) {
        auto searchDist = SysParams::Geometry().largestCompartmentSide;
        vector<Compartment *> compartments;
        GController::findCompartments(mathfunc::vector2Vec<3>(cylinder->coordinate),
                                      cylinder->getCompartment(),
                                      searchDist + _rMax, compartments);
        for (auto c:compartments) {
            nint += c->getCylinders().size();
        }
    }
    //2. Assign optimal blocks and threads
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the maximum occupancy for a full device launch
    if(nint>0) {
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
                                                       CylinderCylinderNLCUDA, 0, 0);
        blocksnthreads.clear();
        blocksnthreads.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreads.push_back(blockSize);
        std::cout<<"NL blocks and threads "<<blocksnthreads.at(0)<<" "<<blocksnthreads.at(1)<<endl;
    }

#endif
    //loop through all neighbor keys
#ifdef CUDA_TIMETRACK
    chrono::high_resolution_clock::time_point mins, mine;
#endif
#ifdef NLORIGINAL//serial
#ifdef CUDA_TIMETRACK
    mins = chrono::high_resolution_clock::now();
#endif
    _list.clear();

    for(auto cylinder: Cylinder::getCylinders()) {
        updateNeighbors(cylinder);
    }

#endif


#ifdef NLSTENCILLIST
    #ifdef CUDA_TIMETRACK
    mins = chrono::high_resolution_clock::now();
    #endif
/*    chrono::high_resolution_clock::time_point mins, mine;
    mins = chrono::high_resolution_clock::now();*/
//    tot = 0;
    _list4mbin.clear();
    //check and reassign cylinders to different bins if needed.
    updateallcylinderstobin();

    for(auto cylinder: Cylinder::getCylinders()) {
        updateNeighborsbin(cylinder);
//        tot += _list4mbin[cylinder].size();
    }
//    std::cout<<"reset NLSTENCILLIST size "<<" "<<tot<<endl;
    /*mine= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_sten(mine - mins);
    std::cout<<"NLSTEN reset time "<<elapsed_sten.count()<<endl;*/
    #ifdef CUDA_TIMETRACK
    std::cout<<"reset NLSTENCILLIST size "<<" "<<tot<<endl;
    mine= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_sten(mine - mins);
    std::cout<<"NLSTEN reset time "<<elapsed_sten.count()<<endl;
    #endif
#endif

//        std::cout<<cylinder->_dcIndex<<" "<<_list.size()<<" "<<vec_numpairs<<" "<<_full<<endl;
#ifdef CUDAACCL_NL
        for(auto cylinder: Cylinder::getCylinders()) {
        //Find surrounding compartments (For now its conservative)
        vector<Compartment*> compartments;
        auto searchDist = SysParams::Geometry().largestCompartmentSide;

        GController::findCompartments(mathfunc::vector2Vec<3>(cylinder->coordinate),
                                      cylinder->getCompartment(),
                                      searchDist + _rMax, compartments);
        for (auto c:compartments) {
            for(auto ncyl:c->getCylinders()){
                pair_cIndex_cnIndex.push_back(cylinder->getStableIndex());
                pair_cIndex_cnIndex.push_back(ncyl->getStableIndex());
            }
//            pair_cIndex_cmp.push_back(cylinder->_dcIndex);
//            pair_cIndex_cmp.push_back(GController::getCompartmentID(c->coordinates()));
        }
    }
    //get total number of pairs.
    int vec_numpairs =0;//TODO remove later
//    for(auto cylinder: Cylinder::getCylinders()) {
//        vec_numpairs += _list[cylinder].size();
//    }
    std::cout<<pair_cIndex_cnIndex.size()<<" "<<nint<<endl;

/*    int *cpu_pair_cIndex_cnIndex;
    cpu_pair_cIndex_cnIndex = new int[pair_cIndex_cnIndex.size()];
    for (auto i = 0; i < pair_cIndex_cnIndex.size(); i++)
        cpu_pair_cIndex_cnIndex[i] = pair_cIndex_cnIndex.at(i);


    int cpu_pair_cIndex_cmp[pair_cIndex_cmp.size()];
    for (auto i = 0; i < pair_cIndex_cmp.size(); i++)
        cpu_pair_cIndex_cmp[i] = pair_cIndex_cmp.at(i);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pair_cIndex_cmp,  pair_cIndex_cmp.size() * sizeof(int)),
                                "cuda data transfer", " NeighborListImpl.h");
    CUDAcommon::handleerror(cudaMemcpy(gpu_pair_cIndex_cmp, cpu_pair_cIndex_cmp, pair_cIndex_cmp.size()
                                           *sizeof(int), cudaMemcpyHostToDevice));*/
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pair_cIndex_cnIndex,  pair_cIndex_cnIndex.size() * sizeof(int)),
                            "cuda data transfer", " NeighborListImpl.h");

    CUDAcommon::handleerror(cudaMemcpy(gpu_pair_cIndex_cnIndex, pair_cIndex_cnIndex.data(), pair_cIndex_cnIndex.size()
                                                                                 *sizeof(int), cudaMemcpyHostToDevice));
//    delete cpu_pair_cIndex_cnIndex;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_NL, 2 * nint * sizeof(int)));

    int cpu_params2[2];
    cpu_params2[0] = nint;
    cpu_params2[1] = int(_full);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params2, 2 * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_params2, cpu_params2, 2 * sizeof(int), cudaMemcpyHostToDevice));
//    if(gpu_numpairs == NULL) {
//    numpairs[0] = 0;
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_numpairs, sizeof(int)));
//    CUDAcommon::handleerror(cudaMemcpy(gpu_numpairs, numpairs, 1 * sizeof(int), cudaMemcpyHostToDevice));
//    }

//        int *a;
//        a = new int[1];
//        a[0] = 0;
//        int *gpu_a;
//        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_a,  sizeof(int)),
//                                "cuda data transfer", " NeighborListImpl.h");
//        CUDAcommon::handleerror(cudaMemcpy(gpu_a, a, sizeof(int), cudaMemcpyHostToDevice));
//        testfunction<<<1,1>>>(gpu_a);

        //Gather all necessary variables.
    auto cylcylnlvars = CUDAcommon::getCylCylNLvars();
    floatingpoint* coord_com = cylcylnlvars.gpu_coord_com;
    int *beadSet = cylcylnlvars.gpu_beadSet;
    int *cylID = cylcylnlvars.gpu_cylID;
    int *filID = cylcylnlvars.gpu_filID;
    int *cmpIDlist = cylcylnlvars.gpu_cmpID;
    int *fvecposition = cylcylnlvars.gpu_fvecpos;
//    int *cylvecpospercmp = cylcylnlvars.gpu_cylvecpospercmp;
//    if(gpu_params == NULL) {
        floatingpoint cpu_params[2];
        cpu_params[0] = floatingpoint(_rMin);
        cpu_params[1] = floatingpoint(_rMax);
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params,  2 * sizeof(floatingpoint)),
                                "cuda data transfer", " NeighborListImpl.h");
        CUDAcommon::handleerror(cudaMemcpy(gpu_params, cpu_params, 2 * sizeof(floatingpoint), cudaMemcpyHostToDevice));
//    }

//    std::cout<<_rMin<<" "<<_rMax<<endl;

    resetintvariableCUDA<<<1,1>>>(gpu_numpairs);
    CylinderCylinderNLCUDA<<<blocksnthreads[0],blocksnthreads[1]>>> (coord_com, beadSet, cylID, filID, cmpIDlist,
            fvecposition, gpu_pair_cIndex_cnIndex, gpu_params, gpu_numpairs, gpu_NL, gpu_params2);
//    CUDAcommon::handleerror(cudaDeviceSynchronize());
    //TODO make this Async and synchronize stream before binding manager call.
    cudaMemcpy(numpairs, gpu_numpairs,  sizeof(int), cudaMemcpyDeviceToHost);

    std::cout<<"Number of neighbors "<<numpairs[0]<<" "<<vec_numpairs<<" Full "<<_full<<endl;

    if(true){
        //copy forces back to CUDA
//        cudaMemcpy(numpairs, gpu_numpairs,  sizeof(int), cudaMemcpyDeviceToHost);
//        std::cout<<"Number of neighbors "<<numpairs[0]<<" "<<vec_numpairs<<" Full "<<_full<<endl;
        int *NL;
        NL = new int[2 * numpairs[0]];
    cudaMemcpy(NL, gpu_NL, 2 * numpairs[0] * sizeof(int), cudaMemcpyDeviceToHost);

//    for(auto ii = 0; ii < numpairs[0]; ii++)
//        std::cout<<NL[2*ii]<<" "<<NL[2*ii+1]<<endl;
//    std::cout<<"-----------------------------------------------------"<<endl;
//    for(auto c:Cylinder::getCylinders()){
//        for(auto cn:_list[c])
//            std::cout<<c->_dcIndex<<" "<<cn->_dcIndex<<endl;
//    }
//    std::cout<<endl;

        //set neighborlist.
        for(auto c:Cylinder::getCylinders()){
            _list[c].clear();
        }
        for(auto id = 0; id < numpairs[0]; id++){
            Cylinder *cylinder = NULL;
            Cylinder *ncylinder = NULL;
            cylinder = Cylinder::getCylinders()[NL[2 * id]];
            ncylinder = Cylinder::getCylinders()[NL[2 * id + 1]];
//            std::cout<<cylinder->_dcIndex<<" "<<NL[2 * id]<<" "<<ncylinder->_dcIndex<<" "<<NL[2 * id +1]<<endl;
//            for(auto c:Cylinder::getCylinders()){
//                if(c->_dcIndex == NL[2 * id])
//                    cylinder = c;
//                else if(c->_dcIndex == NL[2 * id +1])
//                    ncylinder = c;
//                if(cylinder != NULL && ncylinder != NULL)
//                    _list[cylinder].push_back(ncylinder);
//            }

            if(cylinder == NULL || ncylinder == NULL || cylinder->getStableIndex() != NL[2 * id] || ncylinder->getStableIndex() !=
               NL[2 * id + 1]) {
                cout << "Error. Could not find neighborlist from CUDA in Cylinder Database. Check Cylinder IDs Exiting."
                        "." << endl;
                exit(EXIT_FAILURE);
            }
            else{
                _list[cylinder].push_back(ncylinder);
            }
        }
        if(cudacpyforces) {
            CUDAcommon::handleerror(cudaFree(gpu_NL), "cudaFree", "NeighborListImpl.cu");
            CUDAcommon::handleerror(cudaFree(gpu_numpairs), "cudaFree", "NeighborListImpl.cu");
        }
        delete NL;
    }
    CUDAcommon::handleerror(cudaFree(gpu_pair_cIndex_cnIndex),"cudaFree","NeighborListImpl.cu");
    CUDAcommon::handleerror(cudaFree(gpu_params2),"cudaFree","NeighborListImpl.cu");

    CUDAcommon::handleerror(cudaFree(gpu_params),"cudaFree","NeighborListImpl.cu");
#endif
#ifdef CUDAACCL_NLS
    int *gpu_array;
    int *gpu_stage;
    int *gpu_nvec;
    int Nva[1];
    int Nv = 64;
    int stage[2];
    stage[0] =-1;stage[1] = 0;
    Nva[0]= Nv;
    vector<int> testarray(Nv);
    for(int i = 0; i < Nv; i++) {
        testarray.at(i) = (rand() % Nv);
        std::cout<<testarray.at(i)<<" ";
    }
//    vector<floatingpoint> testarray = { 19.0, 128.0, 9.0, 101.0, 14.0, 98.0, 45.0, 32.0, 4.0, 2.0, 50.0, 30.0, 95.0, 38.0, 17.0, 40.0};
    std::cout<<endl;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_array, Nv*sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_array, testarray.data(), Nv* sizeof(int),
                                       cudaMemcpyHostToDevice));
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_stage, 2*sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_stage,stage , 2 * sizeof(int),
                                       cudaMemcpyHostToDevice));
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_nvec, sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_nvec, Nva, sizeof(int),
                                       cudaMemcpyHostToDevice));

//    resetbitonic<<<1,1>>>(gpu_stage);

    for(int i = 0;i<21; i ++){
        incrementbitonic<<<1,1>>>(gpu_stage);
        bitonicsort<<<4,8>>>(gpu_array, gpu_stage, gpu_nvec);
        int cpu_nvec[Nv];
        CUDAcommon::handleerror(cudaMemcpy(cpu_nvec, gpu_array, Nv * sizeof(int),
                                           cudaMemcpyDeviceToHost),"NL", "NLImpl.cu");
        for(int i = 0; i < Nv; i++)
            std::cout<<cpu_nvec[i]<<" ";
        std::cout<<endl;
    }
    int cpu_nvec[Nv];
    CUDAcommon::handleerror(cudaMemcpy(cpu_nvec, gpu_array, Nv * sizeof(int),
                                            cudaMemcpyDeviceToHost),"NL", "NLImpl.cu");
    for(int i = 0; i < Nv; i++)
        std::cout<<cpu_nvec[i]<<" ";
    std::cout<<endl;
#endif
}

vector<Cylinder*> CylinderCylinderNL::getNeighbors(Cylinder* cylinder) {
    return _list[cylinder];
}

//BOUNDARYELEMENT - CYLINDER

void BoundaryCylinderNL::updateNeighbors(BoundaryElement* be) {

    //clear existing
    _list[be].clear();

    //loop through beads, add as neighbor
    for (auto &c : Cylinder::getCylinders()) {


        floatingpoint dist = be->distance(c->coordinate);

        //If within range, add it
        if(dist < _rMax) _list[be].push_back(c);
    }
}

void BoundaryCylinderNL::addNeighbor(Neighbor* n) {

    //return if not a boundary element!
    BoundaryElement* be;
    if(!(be = dynamic_cast<BoundaryElement*>(n))) return;

    //update neighbors
    updateNeighbors(be);
}

void BoundaryCylinderNL::removeNeighbor(Neighbor* n) {

    BoundaryElement* be;
    if(!(be = dynamic_cast<BoundaryElement*>(n))) return;

    _list.erase(be);
}

void BoundaryCylinderNL::addDynamicNeighbor(DynamicNeighbor* n) {

    //return if not a cylinder!
    Cylinder* c;

    if(!(c = dynamic_cast<Cylinder*>(n))) return;

    for(auto it = _list.begin(); it != _list.end(); it++) {

        //if within range, add it
        if(it->first->distance(c->coordinate) < _rMax)
            it->second.push_back(c);
    }
}

void BoundaryCylinderNL::removeDynamicNeighbor(DynamicNeighbor* n) {

    //return if not a cylinder!
    Cylinder* c;

    if(!(c = dynamic_cast<Cylinder*>(n))) return;

    for(auto it = _list.begin(); it != _list.end(); it++) {

        auto cit = find(it->second.begin(), it->second.end(), c);
        if(cit != it->second.end()) it->second.erase(cit);
    }
}

void BoundaryCylinderNL::reset() {

    _list.clear();

    //loop through all neighbor keys
    for(auto boundary: BoundaryElement::getBoundaryElements()) {
        updateNeighbors(boundary);
    }
}

vector<Cylinder*> BoundaryCylinderNL::getNeighbors(BoundaryElement* be) {
    return _list[be];
}




/// Triangle - Beads (filament)

void TriangleFilBeadNL::addNeighbor(Neighbor* n) {
    using MT = Membrane::MeshType;

    if(Triangle* t = dynamic_cast<Triangle*>(n)) {
        const auto& mesh = t->getParent(*ps).getMesh();
        const MT::TriangleIndex ti { t->getTopoIndex() };
        const auto hei0 = mesh.halfEdge(ti);
        const auto hei1 = mesh.next(hei0);
        const auto hei2 = mesh.next(hei1);
        const Vec< 3, floatingpoint > v0 (mesh.attribute(mesh.target(hei0)).getCoordinate(*ps));
        const Vec< 3, floatingpoint > v1 (mesh.attribute(mesh.target(hei1)).getCoordinate(*ps));
        const Vec< 3, floatingpoint > v2 (mesh.attribute(mesh.target(hei2)).getCoordinate(*ps));

        for(auto b : Bead::getBeads()) {
            const auto dist = trianglePointDistance(
                v0, v1, v2,
                b->coordinate()
            );

            if(dist < _rMax) {
                listBT_[b].push_back(t->sysIndex);
                listTB_[t->sysIndex].push_back(b);
            }
            if(dist < rMaxMech_) {
                listBTMech_[b].push_back(t->sysIndex);
                listTBMech_[t->sysIndex].push_back(b);
            }
        }
    }
    else if(Bead* b = dynamic_cast<Bead*>(n)) {

        for(auto t : ps->triangles) {
            const auto& mesh = t.getParent(*ps).getMesh();
            const MT::TriangleIndex ti { t.getTopoIndex() };
            const auto hei0 = mesh.halfEdge(ti);
            const auto hei1 = mesh.next(hei0);
            const auto hei2 = mesh.next(hei1);
            const Vec< 3, floatingpoint > v0 (mesh.attribute(mesh.target(hei0)).getCoordinate(*ps));
            const Vec< 3, floatingpoint > v1 (mesh.attribute(mesh.target(hei1)).getCoordinate(*ps));
            const Vec< 3, floatingpoint > v2 (mesh.attribute(mesh.target(hei2)).getCoordinate(*ps));

            const auto dist = trianglePointDistance(
                v0, v1, v2,
                b->coordinate()
            );

            if(dist < _rMax) {
                listBT_[b].push_back(t.sysIndex);
                listTB_[t.sysIndex].push_back(b);
            }
            if(dist < rMaxMech_) {
                listBTMech_[b].push_back(t.sysIndex);
                listTBMech_[t.sysIndex].push_back(b);
            }
        } // End loop triangles
    }
}

void TriangleFilBeadNL::removeNeighbor(Neighbor* n) {
    
    if(Triangle* t = dynamic_cast<Triangle*>(n)) {
        removeNeighbor_(t->sysIndex, listTB_, listBT_);
        removeNeighbor_(t->sysIndex, listTBMech_, listBTMech_);
    }
    else if(Bead* b = dynamic_cast<Bead*>(n)) {
        removeNeighbor_(b, listBT_, listTB_);
        removeNeighbor_(b, listBTMech_, listTBMech_);
    }
}

void TriangleFilBeadNL::reset() {
    
    listBT_.clear();
    listTB_.clear();
    listBTMech_.clear();
    listTBMech_.clear();

    for(auto t : ps->triangles) {

        auto& mesh = t.getParent(*ps).getMesh();
        const auto vis = medyan::vertexIndices(
            mesh,
            Membrane::MeshType::TriangleIndex { t.getTopoIndex() }
        );

        for(auto b : Bead::getBeads()) {
            const auto dist = trianglePointDistance(
                mesh.attribute(vis[0]).getCoordinate(*ps),
                mesh.attribute(vis[1]).getCoordinate(*ps),
                mesh.attribute(vis[2]).getCoordinate(*ps),
                b->coordinate()
            );

            if(dist < _rMax) {
                listBT_[b].push_back(t.sysIndex);
                listTB_[t.sysIndex].push_back(b);
            }
            if(dist < rMaxMech_) {
                listBTMech_[b].push_back(t.sysIndex);
                listTBMech_[t.sysIndex].push_back(b);
            }
        }
    }
}

} // namespace medyan
