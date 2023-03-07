
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

#ifndef MEDYAN_Bin_h
#define MEDYAN_Bin_h
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>

#include "common.h"
#include "Composite.h"

namespace medyan {
//#include "Cylinder.h"
/*!
 * Bin is a spatial voxel subset of reaction volume. Bin is exclusively used to calculate
 * stenciled cell list corresponding to each NeighborsList.
 * Bin stores the size of spatial cube that it represents along with the center of mass
 * of it.
 * Bin keeps trackable Cylinders in it to easily access all elements at once.
 * Bins that correspond to a particular NeighborList are stored in a single BinGrid.
*/
class Cylinder;
class Bin : public Composite {

private:
	#ifdef DEBUGCONSTANTSEED
	using bincylinderdatatype = unordered_set<Cylinder*, HashbyId<Cylinder*>, customEqualId<Cylinder*>>;

	#else
	using bincylinderdatatype = unordered_set<Cylinder*>; ///< Set of cylinders
    // that are in this bin
	#endif

protected:

    bincylinderdatatype _cylinders; ///< Set of cylinders
    // that are in this bin

    vector<Bin*> _neighbours; ///< Neighbors of the bin
    vector<Bin*> _uniquepermutationneighbours;
    ///OTHER BIN PROPERTIES
    vector<floatingpoint> _coords;  ///< Coordinates of this bin
    short _binGridType = -1;
    vector<int> cindicesvector;

public:
    short _ID = 0;
    vector<int> stencilID; /// template stencil has neighboring bins numbered from 0-26.
/// stencilID vector stores which of those 27 bins are neighbors to the current bin.
/// 13 refers to self

    // Default constructor, only takes in number of dimensions
    Bin(short binGridType, short ID) : _neighbours(), _binGridType(binGridType){
        _ID = ID;
    }

    // Constructor which clones another bin
    Bin(const Bin &C) : _neighbours()
    {
                // Should eventually clone beads, cylinders, boundary elements.... not clear yet
    }

    // Assignment operator
    Bin& operator=(const Bin &other);

    /// Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~Bin() noexcept
    {
        removeFromNeighboursList();

        // Should eventually delete beads, cylinders, boundary elements....not yet clear
    }

    ///Setter and getter for coordinates
    virtual void setCoordinates(vector<floatingpoint> coords) {_coords = coords;}
    virtual const vector<floatingpoint>& coordinates() {return _coords;}

    /// Remove this bin from the neighbor list
    virtual void removeFromNeighboursList() {

        for(auto &n : _neighbours)
        {
            n->removeNeighbour(this);
        }
    }

    virtual string getFullName() const override {return "Bin";};

    ///Add a cylinder to this bin
    void addCylinder(Cylinder* c);

    ///Remove a cylinder from this bin
    ///@note does nothing if cylinder is not in bin already
    void removeCylinder(Cylinder* c);

    ///get the cylinders in this bin
    bincylinderdatatype& getCylinders() {return _cylinders;};

    /// Add a neighboring bin to this bins list of neighbors
    void addNeighbour(Bin *comp) {
        auto nit = find(_neighbours.begin(),_neighbours.end(), comp);
        if(nit==_neighbours.end())
            _neighbours.push_back(comp);
        else
            throw runtime_error(
                    "Bin::addNeighbour(): Bin is already a neighbour");
    }

    /// Remove a neighboring bin
    void removeNeighbour(Bin *comp) {
        auto nit = find(_neighbours.begin(),_neighbours.end(), comp);
        if(nit!=_neighbours.end())
            _neighbours.erase(nit);
    }

    //Add unique permutation neighbor
    void adduniquepermutationNeighbour(Bin *comp) {
        auto nit = find(_uniquepermutationneighbours.begin(),
                _uniquepermutationneighbours.end(), comp);
        if(nit==_uniquepermutationneighbours.end())
            _uniquepermutationneighbours.push_back(comp);
        else
            throw runtime_error(
                    "Bin::addNeighbour(): Bin is already a unique permutation neighbour");
    }

    void removeuniquepermutationNeighbour(Bin *comp) {
        auto nit = find(_uniquepermutationneighbours.begin(),
                _uniquepermutationneighbours.end(), comp);
        if(nit!=_uniquepermutationneighbours.end())
            _uniquepermutationneighbours.erase(nit);
    }

    /// Clone a bin
    /// @note - this does not clone the neighbors, just reactions and species
//    virtual Bin* clone() {
//        Bin *B = new Bin(*this);
//        return B;
//    }

    /// Gives the number of neighbors to this bin
    size_t numberOfNeighbours() const {return _neighbours.size();}

    /// Get the vector list of neighbors to this bin
    vector<Bin*> getNeighbours() {return _neighbours;}
//    const vector<Bin*>& getNeighbours() const {return _neighbours;}

    vector<Bin*> getuniquepermutationNeighbours() {return _uniquepermutationneighbours;}

    /// Print properties of this bin
    virtual void printSelf() const override {
        cout << this->getFullName() << "\n"
             << "Number of neighbors: " << numberOfNeighbours() << endl;
    }

    //GetType implementation just returns zero (no Bin types yet)
    virtual int getType() override {return 0;}

    void updatecindices();

    vector<int> getcindices(){
        updatecindices();
        return cindicesvector;}

};

} // namespace medyan

#endif
