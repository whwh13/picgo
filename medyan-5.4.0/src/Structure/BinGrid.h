
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

#ifndef MEDYAN_BinGrid_h
#define MEDYAN_BinGrid_h

#include <type_traits> // enable_if

#include "common.h"
#include "Composite.h"
#include "Bin.h"
#include "MathFunctions.h"

namespace medyan {
/*!
* BinGrid is a vector of numBins Bins where numBins is provided during instantiation.
* Each Bin in BinGrid is indexed in NeighborListImpl.
* BinGrid pointer belongs to the respective BinGrid.
 * PrototypeBin exists so we can copy properties to all other bins in bingrid.
*/
class BinGrid : public Composite {
private:
//    Bin _prototype_bin; ///< Prototype bin, to be configured
    ///< before initialization
    short _bgType = -1;

public:
    vector<floatingpoint> _binSize; //size in nm of bins that are part of bingrid
    vector<vector<floatingpoint>> binedgevectors; //vectors that represent edges of each bin.
    vector<vector<floatingpoint>> binvertexoffset; //offset vectors to get coordinates of each bin
    // from the bin center coordinates.
    vector<vector<floatingpoint>> binplanenormal={{1.0,0.0,0.0},
                                           {0.0,1.0,0.0},
                                           {0.0,0.0,1.0},
                                           {0.0,0.0,-1.0},
                                           {0.0,-1.0,0.0},
                                           {-1.0,0.0,0.0}};//plane normal vectors
    //vector tracking vertexIDs that form each of the 12 edges of bins.
    vector<vector<int>> edgevertexID ={{0,1},{0,2},{1,3},{2,3},{0,4},{1,5},{2,
                                          6},{3,7},{4,5},{4,6},{5,7},{6,7}};
    //To determine if a bin is essential to obtain all neighbors of a cylinder, shortest
    // distance, shortest distance of the point to each compartment is needed. Shortest
    // distance can be determined based on distance to the edge or vertex depending on
    // the neighbor.
    vector<int> shrtdistVorEorPID = {0,0,1,1,0,2,2,3,3,4,1,5,2,-1,3,6,4,7,4,8,5,9,5,
                                    10,6,11,7};
    vector<int> shrtdistType = {2,1,2,1,3,1,2,1,2,1,3,1,3,0,3,1,3,1,2,1,2,1,3,1,2,1,2};
    //0 self, 1 edge sharing, 2 vertex sharing, 3 face sharing.

    /// Constructor, creates a number of Compartment instances
    BinGrid(int numBins, short bgType, vector<floatingpoint> binSize):
            _bgType(bgType), _binSize(binSize) {

        //vectors corresponding to each of the 12 edges.
        binedgevectors = {{0,0,binSize[2]},
                          {0, binSize[1], 0},
                          {0, binSize[1], 0},
                          {0,0,binSize[2]},
                          {binSize[0],0,0},
                          {binSize[0],0,0},
                          {binSize[0],0,0},
                          {binSize[0],0,0},
                          {0,0,binSize[2]},
                          {0, binSize[1], 0},
                          {0, binSize[1], 0},
                          {0,0,binSize[2]}};
        //offset to get coordinates of 8 vertices of any bin.
        binvertexoffset = {{-binSize[0]/2, -binSize[1]/2, -binSize[2]/2},
                           {-binSize[0]/2, -binSize[1]/2, binSize[2]/2},
                           {-binSize[0]/2, binSize[1]/2, -binSize[2]/2},
                           {-binSize[0]/2, binSize[1]/2, binSize[2]/2},
                           {binSize[0]/2, -binSize[1]/2, -binSize[2]/2},
                           {binSize[0]/2, -binSize[1]/2, binSize[2]/2},
                           {binSize[0]/2, binSize[1]/2, -binSize[2]/2},
                           {binSize[0]/2, binSize[1]/2, binSize[2]/2}};


        //add children
        for(size_t i=0; i<numBins; ++i)
            addChild(unique_ptr<Component>(new Bin(bgType, i)));
    }

    /// Get bins that this grid holds
    vector<Bin*> getBins() {

        vector<Bin*> bins;

        for (auto &c : children())
            bins.push_back((Bin*)c.get());

        return bins;
    }

    void updatecindices(){
        for (auto &c : children())
            ((Bin*)c.get())->updatecindices();
    }

//    /// Get a bin at certain index
    Bin* getBin(int index) {
        return (Bin*)(children()[index].get());
    }

    /// Get name of this bin grid
    virtual string getFullName() const {return string("BinGrid");};

    /// Get the protobin from this grid, in order to configure and then initialize
//    Bin& getProtoBin() {return _prototype_bin;}
//    const Bin& getProtoBin() const {return _prototype_bin;}

    /// Print properties of this grid
    virtual void printSelf() const override {
        cout << getFullName() << endl;
        cout << "Number of Bin objects: " << numberOfChildren() << endl;
        cout << "Type: " << _bgType <<endl;
        for(auto &c : children())
            c->printSelf();
    }
    ///GetType implementation just returns zero (no CompartmentGrid types yet)
    virtual int getType() {return _bgType;}

    vector<floatingpoint> getbinsize(){ return _binSize;}

    bool iswithincutoff(vector<floatingpoint> cylcoord, vector<floatingpoint> bincoord, int
                        NbinstencilID, floatingpoint cutoff){
        int type = shrtdistType[NbinstencilID];
        if(type ==0) return true;
        else if(type ==1){
            int edgeID = shrtdistVorEorPID[NbinstencilID];
            vector<floatingpoint> edgevec = {binedgevectors.at(edgeID)[0], binedgevectors.at
                    (edgeID)[1], binedgevectors.at(edgeID)[2]};

            int vertexID = edgevertexID[edgeID][0];
            vector<floatingpoint> vertexoffset = { binvertexoffset.at(vertexID)[0],
                                            binvertexoffset.at(vertexID)[1],
                                            binvertexoffset.at(vertexID)[2]};
            vector<floatingpoint> vcoord = {bincoord[0] + vertexoffset[0],
                                     bincoord[1] + vertexoffset[1],
                                     bincoord[2] + vertexoffset[2]};
            vector<floatingpoint> pointvec = {cylcoord[0] - vcoord[0],cylcoord[1] - vcoord[1],
                                       cylcoord[2] - vcoord[2]};
            vector<floatingpoint> normal = mathfunc::normalizeVector(edgevec);

            floatingpoint proj = mathfunc::scalarprojection(normal,pointvec);
            floatingpoint d = mathfunc::magnitude(pointvec);
            floatingpoint dist = sqrt(d*d - proj*proj);
/*            std::cout<<"edgeID "<<edgeID<<"  vertexID "<<vertexID<<" vcoord "
                    ""<<vcoord[0]<<" "<<vcoord[1]<<" "<<vcoord[2]<<endl;*/
            if(dist <= cutoff) return true;
            else return false;
        }
        else if(type == 2){
            int vertexID = shrtdistVorEorPID[NbinstencilID];
            vector<floatingpoint> vertexoffset = { binvertexoffset.at(vertexID)[0],
                                            binvertexoffset.at(vertexID)[1],
                                            binvertexoffset.at(vertexID)[2]};
            vector<floatingpoint> vcoord = {bincoord[0] + vertexoffset[0],
                                     bincoord[1] + vertexoffset[1],
                                     bincoord[2] + vertexoffset[2]};
            floatingpoint dist = mathfunc::twoPointDistance(vcoord, cylcoord);
            if(dist <= cutoff) return true;
            else return false;
        }
        else if(type == 3){
            int planeID = shrtdistVorEorPID[NbinstencilID];
            int vertexID = 7;
            if(planeID<=2) vertexID = 0;
            vector<floatingpoint> vertexoffset = { binvertexoffset.at(vertexID)[0],
                                            binvertexoffset.at(vertexID)[1],
                                            binvertexoffset.at(vertexID)[2]};
            vector<floatingpoint> vcoord = {bincoord[0] + vertexoffset[0],
                                     bincoord[1] + vertexoffset[1],
                                     bincoord[2] + vertexoffset[2]};
            vector<floatingpoint> planeeqn = {binplanenormal.at(planeID)[0],
                                       binplanenormal.at(planeID)[1],
                                       binplanenormal.at(planeID)[2]};
            floatingpoint d = -mathfunc::scalarprojection(planeeqn,vcoord);
            planeeqn.push_back(d);
            floatingpoint dist = abs(mathfunc::getdistancefromplane(cylcoord.data(),planeeqn.data
                    ()));
//            std::cout<<planeeqn[0]<<" "<<planeeqn[1]<<" "<<planeeqn[2]<<" "
//                    ""<<planeeqn[3]<<" "<<dist<<endl;
            if(dist <= cutoff)
                return true;
            else return false;
        }
        return false;
    }

    template< typename VecType, std::enable_if_t<VecType::vec_size == 3>* = nullptr >
    bool iswithincutoff(const VecType& cylcoord, vector<floatingpoint> bincoord, int
    NbinstencilID, floatingpoint cutoff){
        int type = shrtdistType[NbinstencilID];
        if(type ==0) return true;
        else if(type ==1){
            int edgeID = shrtdistVorEorPID[NbinstencilID];
            vector<floatingpoint> edgevec = {binedgevectors.at(edgeID)[0], binedgevectors.at
                    (edgeID)[1], binedgevectors.at(edgeID)[2]};

            int vertexID = edgevertexID[edgeID][0];
            vector<floatingpoint> vertexoffset = { binvertexoffset.at(vertexID)[0],
                                            binvertexoffset.at(vertexID)[1],
                                            binvertexoffset.at(vertexID)[2]};
            vector<floatingpoint> vcoord = {bincoord[0] + vertexoffset[0],
                                     bincoord[1] + vertexoffset[1],
                                     bincoord[2] + vertexoffset[2]};
            vector<floatingpoint> pointvec = {cylcoord[0] - vcoord[0],cylcoord[1] - vcoord[1],
                                       cylcoord[2] - vcoord[2]};
            vector<floatingpoint> normal = mathfunc::normalizeVector(edgevec);

            floatingpoint proj = mathfunc::scalarprojection(normal,pointvec);
            floatingpoint d = mathfunc::magnitude(pointvec);
            floatingpoint dist = sqrt(d*d - proj*proj);
/*            std::cout<<"edgeID "<<edgeID<<"  vertexID "<<vertexID<<" vcoord "
                    ""<<vcoord[0]<<" "<<vcoord[1]<<" "<<vcoord[2]<<endl;*/
            if(dist <= cutoff) return true;
            else return false;
        }
        else if(type == 2){
            int vertexID = shrtdistVorEorPID[NbinstencilID];
            vector<floatingpoint> vertexoffset = { binvertexoffset.at(vertexID)[0],
                                            binvertexoffset.at(vertexID)[1],
                                            binvertexoffset.at(vertexID)[2]};
            medyan::Vec< 3, floatingpoint > vcoord {bincoord[0] + vertexoffset[0],
                                                      bincoord[1] + vertexoffset[1],
                                                      bincoord[2] + vertexoffset[2]};
            floatingpoint dist = medyan::distance(vcoord, cylcoord);
            if(dist <= cutoff) return true;
            else return false;
        }
        else if(type == 3){
            int planeID = shrtdistVorEorPID[NbinstencilID];
            int vertexID = 7;
            if(planeID<=2) vertexID = 0;
            vector<floatingpoint> vertexoffset = { binvertexoffset.at(vertexID)[0],
                                            binvertexoffset.at(vertexID)[1],
                                            binvertexoffset.at(vertexID)[2]};
            vector<floatingpoint> vcoord = {bincoord[0] + vertexoffset[0],
                                     bincoord[1] + vertexoffset[1],
                                     bincoord[2] + vertexoffset[2]};
            vector<floatingpoint> planeeqn = {binplanenormal.at(planeID)[0],
                                       binplanenormal.at(planeID)[1],
                                       binplanenormal.at(planeID)[2]};
            floatingpoint d = -mathfunc::scalarprojection(planeeqn,vcoord);
            planeeqn.push_back(d);
            floatingpoint dist = abs(mathfunc::getdistancefromplane(cylcoord.value.data(), planeeqn.data
                    ()));
//            std::cout<<planeeqn[0]<<" "<<planeeqn[1]<<" "<<planeeqn[2]<<" "
//                    ""<<planeeqn[3]<<" "<<dist<<endl;
            if(dist <= cutoff)
                return true;
            else return false;
        }
        return true;
    }

};

} // namespace medyan

#endif //MEDYAN_BINGRID_H
