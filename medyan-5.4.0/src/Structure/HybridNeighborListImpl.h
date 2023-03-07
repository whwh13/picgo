
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

#ifndef MEDYAN_HybridNeighborListImpl_h
#define MEDYAN_HybridNeighborListImpl_h
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
#include <unordered_map>

#include <vector>

#include "common.h"

#include "HybridNeighborList.h"
#include "DynamicNeighbor.h"
#include "BinGrid.h"
#include "SysParams.h"

namespace medyan {
//FORWARD DECLARATIONS
class Cylinder;

//Hybrid Manager to control calculate all neighborlists in a single search.
//For every pair-wise distance between two filament types, there are fullstatus and uniquestatus requirements.
//Unique status determines whether
class HybridCylinderCylinderNL : public HybridNeighborList {

private:
    static const bool CROSSCHECK_NL_SWITCH = false;

    vector<unordered_map<Cylinder*, vector<Cylinder*>>> _list4mbinvec;
//    unordered_map<Cylinder*, vector<Cylinder*>> _list4mbin;
    ///Helper function to update neighbors
    ///@param runtime - specifying whether the cylinder is being
    ///created/destroyed at runtime vs at a full neighbor list update.
    void updateNeighbors(Cylinder* cylinder, bool runtime = false);
    static short totalhybridNL;
public:
    //In the implementation of Hybrid NeighborList, each compartment has a single HybridNeighborList pointer.
    // So, _ID is always 0.

    ///< The neighbors list, as a hash map
    void initializeBinGrid();//Segments reaction volume into a collection of bins whose dimension is determined by the _largestrMax.
    void generateConnections();//Assigns 3-d index and coordinate to each bin. Also determines neighbors of each bin.

    vector<int> _grid; ///< Number of bins in each dimension
    vector<floatingpoint> _binSize; ///< Bin size in each dimension
    vector<int> _size;       ///< Size of entire grid spanned in each dimension
//    short NLcyltypes[2] = {0,0};// The two types of cylinders that engage in this neighbors
    // List
    BinGrid* _binGrid;
    Bin* getBin(const vector<floatingpoint> &coords);//Get bin based on Coordinates of bin cneter-of-mass
    Bin* getBin(const vector<size_t> &indices);
    void assignallcylinderstobin();//Assigns all cylinders to bins in the grid
    void assignbin(Cylinder* cyl);//Associates cylinder with a bin based on coordinates.
    void unassignbin(Cylinder* cyl, Bin* bin);//Disassociates cylinder from the bin.
    void updateallcylinderstobin();//updates bin associations of all cylinders
    void updatebin(Cylinder* cyl);//updates bin association of a cylinder.
    void updateNeighborsbin(Cylinder* cylinder, bool runtime = false);
    vector<Cylinder*> getNeighborsstencil(short HNLID, Cylinder* cylinder);//Each unique neighborList in the HybridNeighborList has an associated ID HNLID.
    //This ID is assigned when the parameters are set.

    //For any pair wise cylinder map requested, additional parameters such as fullstatus and uniquestatus are required.

    short setneighborsearchparameters(short ftype1, short ftype2,bool uniquestatus, bool fullstatus,
                                      float rMax, float rMin){
        short returnHNLID = -1;//Each unique neighborList in the HybridNeighborList has an associated ID HNLID.
        //This ID is assigned when the parameters are set.
        vector<short> ftypepairs;
        float localrMinsq = rMin * rMin;
        float localrMaxsq = rMax * rMax;
        if(ftype1 < ftype2)
            ftypepairs = {ftype1,ftype2};
        else
            ftypepairs = {ftype2, ftype1};
/*        if(uniquestatus){
            unordered_map<Cylinder*, vector<Cylinder*>> _nlstbin;
            _list4mbinvec.push_back(_nlstbin);
            _rMaxsqvec.push_back(localrMaxsq);
            _rMinsqvec.push_back(localrMinsq);
            _filamentIDvec.push_back(ftypepairs);
            HNLIDvec.push_back(totalhybridNL); // assign ID
            returnHNLID = totalhybridNL;
            totalhybridNL++;//increment total NLs by 1
            uniquestatusvec.push_back(uniquestatus);
            _fullstatusvec.push_back(fullstatus);
            _smallestrMinsq = min(localrMinsq,_smallestrMinsq);
            _largestrMaxsq = max(localrMaxsq, _largestrMaxsq);
            float localmaxcylsize = max(SysParams::Geometry().cylinderSize[ftype1],
                                        SysParams::Geometry().cylinderSize[ftype2]);
            _maxcylsize= max(localmaxcylsize,_maxcylsize);
        }*/
//        else{
            //search through the vectors to find if the non-unique NL
            //If it is not found, add.
            bool isfound = false;
            unordered_map<Cylinder*, vector<Cylinder*>> _nlstbin;
            for(short idx = 0; idx <totaluniquefIDpairs; idx++) {
                vector<short> fIDpair = _filamentIDvec[idx];
                if (isfound|| fIDpair[0] != ftypepairs[0] || fIDpair[1] != ftypepairs[1])
                    continue;
                for (int idx2 = 0; idx2 < uniquestatusvec[idx].size(); idx2++) {
                    if (uniquestatus) {
                        isfound = true;
                        _rMaxsqvec[idx].push_back(localrMaxsq);
                        _rMinsqvec[idx].push_back(localrMinsq);
                        _fullstatusvec[idx].push_back(fullstatus);
                        _smallestrMinsq = min(localrMinsq,_smallestrMinsq);
                        _largestrMaxsq = max(localrMaxsq, _largestrMaxsq);
                        float localmaxcylsize = max(SysParams::Geometry().cylinderSize[ftype1],
                                                    SysParams::Geometry().cylinderSize[ftype2]);
                        _maxcylsize= max(localmaxcylsize,_maxcylsize);
                        HNLIDvec[idx].push_back(totalhybridNL); // assign ID
                        returnHNLID = totalhybridNL;
                        _list4mbinvec.push_back(_nlstbin);
                        totalhybridNL++;
                        break;
                    } else {
                        if (uniquestatusvec[idx][idx2] != uniquestatus) continue;
                        isfound = true;
                        returnHNLID = HNLIDvec[idx][idx2];
                            //insert and break
                        _rMaxsqvec[idx][idx2] = max(localrMaxsq, _rMaxsqvec[idx][idx2]);
                        _rMinsqvec[idx][idx2] = min(localrMinsq, _rMinsqvec[idx][idx2]);
                        _smallestrMinsq = min(localrMinsq, _smallestrMinsq);
                        _largestrMaxsq = max(localrMaxsq, _largestrMaxsq);
                        break;
                    }
                }
                int countbounds = uniquestatusvec[idx].size();
                if(!isfound && !uniquestatus &&  countbounds > 0){
                    HNLIDvec[idx].push_back(totalhybridNL);
                    returnHNLID = totalhybridNL;
                    totalhybridNL++;
                    //set and break
                    _rMaxsqvec[idx].push_back(localrMaxsq);
                    _rMinsqvec[idx].push_back(localrMinsq);
                    uniquestatusvec[idx].push_back(uniquestatus);
                    _fullstatusvec[idx].push_back(fullstatus);
                    _smallestrMinsq = min(localrMinsq, _smallestrMinsq);
                    _largestrMaxsq = max(localrMaxsq, _largestrMaxsq);
                    _list4mbinvec.push_back(_nlstbin);
                     isfound = true;
                    break;
                }
            }
            if(isfound == false){
                vector<float> localrMinsqvec = {localrMinsq};
                vector<float> localrMaxsqvec = {localrMaxsq};
                vector<short> totalhybridNLvec = {totalhybridNL};
                vector<bool> localuniquestatusvec ={uniquestatus};
                vector<bool> localfullstatusvec ={fullstatus};
                unordered_map<Cylinder*, vector<Cylinder*>> _nlstbin;
//                if(uniquestatus){
                _list4mbinvec.push_back(_nlstbin);
                _rMaxsqvec.push_back(localrMaxsqvec);
                _rMinsqvec.push_back(localrMinsqvec);
                _filamentIDvec.push_back(ftypepairs);
                HNLIDvec.push_back(totalhybridNLvec); // assign ID
                returnHNLID = totalhybridNL;
                totalhybridNL++;//increment total NLs by 1
                uniquestatusvec.push_back(localuniquestatusvec);
                _fullstatusvec.push_back(localfullstatusvec);
                _smallestrMinsq = min(localrMinsq,_smallestrMinsq);
                _largestrMaxsq = max(localrMaxsq, _largestrMaxsq);
                float localmaxcylsize = max(SysParams::Geometry().cylinderSize[ftype1],
                                            SysParams::Geometry().cylinderSize[ftype2]);
                _maxcylsize= max(localmaxcylsize,_maxcylsize);
                totaluniquefIDpairs++;
/*                }
                else {
                    _list4mbinvec.push_back(_nlstbin);
                    _rMaxsqvec.push_back(localrMaxsq);
                    _rMaxsqvec.push_back(localrMinsq);
                    _filamentIDvec.push_back(ftypepairs);
                    HNLIDvec.push_back(totalhybridNL); // assign ID
                    returnHNLID = totalhybridNL;
                    totalhybridNL++;//increment total NLs by 1
                    uniquestatusvec.push_back(uniquestatus);
                    _fullstatusvec.push_back(fullstatus);
                    _smallestrMinsq = min(rMin * rMin, _smallestrMinsq);
                    _largestrMaxsq = max(rMax * rMax, _largestrMaxsq);
                    float localmaxcylsize = max(SysParams::Geometry().cylinderSize[ftype1],
                                                SysParams::Geometry().cylinderSize[ftype2]);
                    _maxcylsize = max(localmaxcylsize, _maxcylsize);
                }*/
            }
//        }
        return returnHNLID;
    }


    //constructor
    HybridCylinderCylinderNL() {}

    //Creates bingrid and assigns cylinders to each bin.
    virtual void initializeHybridNeighborList(){
        initializeBinGrid();
        assignallcylinderstobin();
    }
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);

    //@{
    /// The implementation of these functions calls the static version,
    /// all cylinders are dynamic
    virtual void addDynamicNeighbor(DynamicNeighbor* n) {addNeighbor(n);}
    virtual void removeDynamicNeighbor(DynamicNeighbor* n) {removeNeighbor(n);}
    //@}
    virtual void reset();

    /// Get all cylinder neighbors
    vector<Cylinder*> getNeighbors(Cylinder* cylinder);

};
#endif

} // namespace medyan

#endif
