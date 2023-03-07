
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

#ifndef MEDYAN_NeighborListImpl_h
#define MEDYAN_NeighborListImpl_h

#include <algorithm> // remove
#include <array>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include "common.h"
#include "NeighborList.h"
#include "DynamicNeighbor.h"
#include "BinGrid.h"
#include "Structure/Bead.h"
#include "Structure/BoundaryElement.h"
#include "Structure/Bubble.h"
#include "Structure/CellList.hpp"
#include "Structure/Cylinder.h"
#include "Structure/Filament.h"
#include "Structure/SurfaceMesh/Triangle.hpp"
#include "SysParams.h"
#include "Util/StableVector.hpp"

namespace medyan {
//FORWARD DECLARATIONS
class SubSystem;


struct NeighborListCellList3D {

    floatingpoint               boxsize = 0;
    std::array<int, 3>          ncells {};
    CellListManager< int, int > list;

    void init(floatingpoint boxsize, std::array<int, 3> ncells) {
        this->boxsize = boxsize;
        this->ncells = ncells;
        const int ncell = ncells[0] * ncells[1] * ncells[2];

        list.clearAll();
        for(int i = 0; i < ncell; ++i) {
            auto headindex = list.addHead(i);
            // Important: the head index is the same as the cell index.
            assert(headindex == i);
        }
    }
    void init(floatingpoint boxsize, const Eigen::Vector3d& domain) {
        ncells[0] = static_cast<int>(std::ceil(domain[0] / boxsize));
        ncells[1] = static_cast<int>(std::ceil(domain[1] / boxsize));
        ncells[2] = static_cast<int>(std::ceil(domain[2] / boxsize));
        init(boxsize, ncells);
    }

    // Cell traversal utilities.
    //----------------------------------

    bool isValidCellIndex3(std::array<int, 3> index3) const {
        return index3[0] >= 0 && index3[0] < ncells[0] &&
               index3[1] >= 0 && index3[1] < ncells[1] &&
               index3[2] >= 0 && index3[2] < ncells[2];
    }

    // Assuming that x is the fastest varying index.
    std::array<int, 3> getCellIndex3(int index) const {
        int ix = index % ncells[0];
        index /= ncells[0];
        int iy = index % ncells[1];
        index /= ncells[1];
        int iz = index;
        return {ix, iy, iz};
    }
    int getCellIndex(std::array<int, 3> index3) const {
        return index3[0] + ncells[0] * (index3[1] + ncells[1] * index3[2]);
    }

    std::array<int, 3> resolveCellIndex3(const Eigen::Vector3d& coord) const {
        std::array<int, 3> index3;
        index3[0] = static_cast<int>(std::floor(coord[0] / boxsize));
        index3[1] = static_cast<int>(std::floor(coord[1] / boxsize));
        index3[2] = static_cast<int>(std::floor(coord[2] / boxsize));
        return index3;
    }
    int resolveCellIndex(const Eigen::Vector3d& coord) const {
        return getCellIndex(resolveCellIndex3(coord));
    }

    // Accessors.
    auto numCells() const { return ncells[0] * ncells[1] * ncells[2]; }

    // Loop over all neighbors of the specified cell.
    //
    // Note:
    // - Func: cellIndex -> void.
    // - The current cell is also included.
    // - If halfDir is true, will search for only half of the 26 neighbors. This leads to 14 cells to be visited including the current cell.
    // - Out-of-range cells will not be used.
    template< bool halfDir, typename Func >
    void forEachNeighbor(int cellIndex, Func&& func) const {
        auto index3 = getCellIndex3(cellIndex);
        for(int dz = (halfDir ? 0 : -1); dz <= 1; ++dz) {
            for(int dy = (halfDir && dz == 0 ? 0 : -1); dy <= 1; ++dy) {
                for(int dx = (halfDir && dz == 0 && dy == 0 ? 0 : -1); dx <= 1; ++dx) {
                    auto nindex3 = index3;
                    nindex3[0] += dx;
                    nindex3[1] += dy;
                    nindex3[2] += dz;
                    if(isValidCellIndex3(nindex3)) {
                        func(getCellIndex(nindex3));
                    }
                }
            }
        }
    }

    // Loop over all neighboring pairs of elements.
    //
    // Note:
    // - Func: int (1st element), int (2nd element) -> void.
    template< typename Func >
    void forAllPairs(Func&& func, bool includeSelf = false) const {
        const int nc = numCells();
        for(int i = 0; i < nc; ++i) {
            const auto doForNeighbor = [&, this](int j) {
                for(auto ei : list.getElements(i)) {
                    for(auto ej : list.getElements(j)) {
                        if(includeSelf || ei != ej) {
                            func(ei, ej);
                        }
                    }
                }
            };
            forEachNeighbor<true>(i, doForNeighbor);
        }
    }

    // Element manipulation.
    //----------------------------------

    // Add an element to the cell list.
    // Parameters:
    // - elementHandle:  The unique representation of the element in int.
    // - coord:          The coordinate of the element.
    // Returns the index of the element in the cell list.
    auto add(int elementHandle, const Eigen::Vector3d& coord) {
        auto cellIndex = resolveCellIndex(coord);
        return list.addElement(elementHandle, cellIndex);
    }

    // Note:
    // - Does not reset the cell list.
    void clearElements() {
        list.clearElements();
    }
};


/// An implementation of NeighborList for Cylinder-Cylinder interactions
/// This can be a half or full list depending on the usage.
class CylinderCylinderNL : public NeighborList {

private:
    unordered_map<Cylinder*, vector<Cylinder*>> _list;
    ///< The neighbors list, as a hash map

    bool _full; ///<Specifying whether this is a full or half list
    unordered_map<Cylinder*, vector<Cylinder*>> _list4mbin;
#ifdef CUDAACCL_NL
    vector<int> blocksnthreads;
    int nint;
//    vector<int> pair_cIndex_cmp;
    vector<int> pair_cIndex_cnIndex;
    int *gpu_pair_cIndex_cnIndex;
//    int *gpu_pair_cIndex_cmp;
    floatingpoint *gpu_params = NULL;
    int *gpu_NL = NULL;
    int *gpu_numpairs = NULL;
    int *gpu_params2;
    int numpairs[1];
#endif

    ///Helper function to update neighbors
    ///@param runtime - specifying whether the cylinder is being
    ///created/destroyed at runtime vs at a full neighbor list update.
    void updateNeighbors(Cylinder* cylinder, bool runtime = false);

public:
#ifdef CUDAACCL_NL
    bool cudacpyforces = false;

    int getNLsize() {
        return numpairs[0];
    }
    int* getNLsizeCUDA(){
        return gpu_numpairs;
    }
    int* getNLCUDA(){
        return gpu_NL;
    }
#endif
    short _ID; //ID helps link binGridType to NeighborList.
#ifdef NLSTENCILLIST
    /* New implementation of NeighborList which  uses stencils to ensure that ALL
    relevant neghbors of a bindingSite/cylinder can be obtained from just parsing through
    the nearest 27 bins. In stencil neighborList implementation, space is divided into
    sub-volumes referred to as bins. The bin sizes are determined based on the bindingdistance
    and cylinder length.
    */
    void initializeBinGrid();//Initializes bins based on Grid dimensions
    void generateConnections() ;//Assigns coordinate for each bin. Generate a set of neighbors for each bin.

    vector<int> _grid; ///< Number of bins in each dimension
    vector<floatingpoint> _binSize; ///< Bin size in each dimension
    vector<int> _size;       ///< Size of entire grid spanned in each dimension
    floatingpoint maxcylsize; // maximum of the two cylinder sizes that are part of this
    // CylinderCylinderNL.
    short NLcyltypes[2] = {0,0};// The two types of cylinders that engage in this neighbors
    // List
    BinGrid* _binGrid;
    Bin* getBin(const vector<floatingpoint> &coords);// returns Bin pointer corresponding to a bin center coordinate
    Bin* getBin(const vector<size_t> &indices);// returns Bin pointer that coresponds to an integer coordinate or bins.
    //Integer coordiante/index is given by coordinate/bin_dimension
    void assignallcylinderstobin();//Assigns all cylinders to respective bins
    void assignbin(Cylinder* cyl);//Associates a Bin pointer to the Cylinder based on coordinate.
    void unassignbin(Cylinder* cyl, Bin* bin);//Removes cylinder from the Bin
    void updateallcylinderstobin();//Checks cylinder coordinates and reassigns bins
    void updatebin(Cylinder* cyl);
    void updateNeighborsbin(Cylinder* cylinder, bool runtime = false);
    vector<Cylinder*> getNeighborsstencil(Cylinder* cylinder);
//    void setbinvars(){
//        initializeBinGrid();
//        assignallcylinderstobin();
//        NLcyltypes[0] = 0;
//        NLcyltypes[1] = 0;
//        _ID = SysParams::numcylcylNL;
//        SysParams::numcylcylNL++;
//        //Determine binSize based on the longer of the two cylinders involved in the NL.
//        std::cout<<"Cylinder size "<<SysParams::Geometry()
//                .cylinderSize[NLcyltypes[0]]<<endl;
//        maxcylsize = max(SysParams::Geometry().cylinderSize[NLcyltypes[0]],
//                         SysParams::Geometry().cylinderSize[NLcyltypes[1]]);
//    }
#endif
#ifdef CUDAACCL_NLS
    struct bin *binGridv;
    struct bin *gpu_binGrid;
    cudaStream_t  stream_NL;
#endif
    //While Excluded volume neighborlist is not a full list, linker and motor
    // neighborlists are.
    CylinderCylinderNL(float rMax, float rMin = 0.0, bool full = false, short ID = 0)
            : NeighborList(rMax, rMin), _full(full) {
#ifdef NLSTENCILLIST
        //Right now only two cylinders of same type can be considered for NL.
        NLcyltypes[0] = 0;
        NLcyltypes[1] = 0;
        maxcylsize = max(SysParams::Geometry().cylinderSize[NLcyltypes[0]],
                         SysParams::Geometry().cylinderSize[NLcyltypes[1]]);
//        std::cout<<"Cylinder size "<<SysParams::Geometry()
//                .cylinderSize[NLcyltypes[0]]<<endl;
        initializeBinGrid();
        assignallcylinderstobin();

        _ID = SysParams::numcylcylNL;
//        std::cout<<"NL ID "<<SysParams::numcylcylNL<<endl;
        SysParams::numcylcylNL++;
        //Determine binSize based on the longer of the two cylinders involved in the NL.

#endif
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


/// An implementation of NeighborList for BoundaryElement-Cylinder interactions
class BoundaryCylinderNL : public NeighborList {

private:
    unordered_map<BoundaryElement*, vector<Cylinder*>> _list;
    ///< The neighbors list, as a hash map

    ///Helper function to update neighbors
    void updateNeighbors(BoundaryElement* be);

public:
    BoundaryCylinderNL(float rMax): NeighborList(rMax) {}

    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);

    virtual void addDynamicNeighbor(DynamicNeighbor* n);
    virtual void removeDynamicNeighbor(DynamicNeighbor* n);

    virtual void reset();

    /// Get all Cylinder neighbors of a boundary element
    vector<Cylinder*> getNeighbors(BoundaryElement* be);
};


/// An implementation of NeighborList for BoundaryElement-Bubble interactions
class BoundaryBubbleNL {
public:
    using BoundaryElementElement = BoundaryElement*;
    using BubbleElement = StableVectorIndex<Bubble>;
private:
    floatingpoint rMax_ = 0;
    std::unordered_map<BoundaryElementElement, std::vector<BubbleElement>> list_;
    ///< The neighbors list, as a hash map

    ///Helper function to update neighbors
    template< typename Context >
    void updateNeighbors(Context& sys, BoundaryElementElement be) {
        auto& listbe = list_[be];
        // Clear existing.
        listbe.clear();
        // Loop through all bubbles and add as neighbor.
        for(auto& b : sys.bubbles) {
            floatingpoint dist = be->distance(mathfunc::vec2Vector(b.coord));
            // If within cutoff, add as neighbor.
            if(dist < rMax_) {
                listbe.push_back(b.sysIndex);
            }
        }
    }

public:
    BoundaryBubbleNL(floatingpoint rMax): rMax_(rMax) {}

    template< typename Context >
    void addNeighbor(Context& sys, BoundaryElementElement be) {
        updateNeighbors(sys, be);
    }
    template< typename Context >
    void removeNeighbor(Context& sys, BoundaryElementElement be) {
        list_.erase(be);
    }
    template< typename Context >
    void addDynamicNeighbor(Context& sys, BubbleElement b) {
        for(auto& [be, neighbors] : list_) {
            if(be->distance(mathfunc::vec2Vector(sys.bubbles[b].coord)) < rMax_) {
                neighbors.push_back(b);
            }
        }
    }
    template< typename Context >
    void removeDynamicNeighbor(Context& sys, BubbleElement b) {
        for(auto& [be, neighbors] : list_) {
            auto it = std::find(neighbors.begin(), neighbors.end(), b);
            if(it != neighbors.end()) {
                neighbors.erase(it);
            }
        }
    }
    template< typename Context >
    void reset(Context& sys) {
        list_.clear();
        for(auto boundary : BoundaryElement::getBoundaryElements()) {
            updateNeighbors(sys, boundary);
        }
    }

    /// Get all Bubble neighbors of a boundary element
    auto& getNeighbors(BoundaryElementElement be) const {
        return list_.at(be);
    }
};

/// An implementation of NeighborList for Bubble-Bubble interactions
/// @note - This is currently implemented as a half list only
class BubbleBubbleNL {
public:
    using BubbleElement = StableVectorIndex<Bubble>;

private:
    floatingpoint rMax_ = 0;
    // This is a half list.
    std::unordered_map<BubbleElement, std::vector<BubbleElement>, StableVectorIndexHash<Bubble>> list_;
    ///< The neighbors list, as a hash map

    ///Helper function to update neighbors
    template< typename Context >
    void updateNeighbors(Context& sys, BubbleElement bb) {
        auto& listbb = list_[bb];
        // Clear existing.
        listbb.clear();
        // Loop through all bubbles and add as neighbor.
        for(auto& b : sys.bubbles) {
            auto dist2 = distance2(b.coord, sys.bubbles[bb].coord);
            // If within cutoff, add as neighbor.
            if(dist2 < rMax_ * rMax_) {
                listbb.push_back(b.sysIndex);
            }
        }
    }

public:
    BubbleBubbleNL(floatingpoint rMax): rMax_(rMax) {}

    template< typename Context >
    void addDynamicNeighbor(Context& sys, BubbleElement bb) {
        updateNeighbors(sys, bb);
    }
    template< typename Context >
    void removeDynamicNeighbor(Context& sys, BubbleElement bb) {
        list_.erase(bb);
        // Remove from other lists.
        for(auto& [b, neighbors] : list_) {
            auto it = std::find(neighbors.begin(), neighbors.end(), bb);
            if(it != neighbors.end()) {
                neighbors.erase(it);
            }
        }
    }

    template< typename Context >
    void reset(Context& sys) {
        list_.clear();
        for(auto& b : sys.bubbles) {
            updateNeighbors(sys, b.sysIndex);
        }
    }

    /// Get all Bubble neighbors of a bubble
    auto& getNeighbors(BubbleElement bb) const {
        return list_.at(bb);
    }
};

/// An implementation of NeighborList for Bubble-Cylinder interactions
class BubbleBeadNL {
public:
    using BubbleElement = StableVectorIndex<Bubble>;
    using BeadElement = Bead*;

private:
    floatingpoint rMax_ = 0;
    std::unordered_map<BubbleElement, std::vector<BeadElement>, StableVectorIndexHash<Bubble>> list_;
    ///< The neighbors list, as a hash map

    ///Helper function to update neighbors
    template< typename Context >
    void updateNeighbors(Context& sys, BubbleElement bbIndex) {
        auto& listbb = list_[bbIndex];
        // Clear existing.
        listbb.clear();
        // Loop through all cylinders and add as neighbor.
        for(auto& b: Bead::getBeads()) {
            if(shouldBeNeighbors(sys, bbIndex, b)) {
                listbb.push_back(b);
            }
        }
    }

public:
    BubbleBeadNL(float rMax): rMax_(rMax) {}

    template< typename Context >
    void addDynamicNeighbor(Context& sys, BubbleElement bb) {
        updateNeighbors(sys, bb);
    }
    template< typename Context >
    void addDynamicNeighbor(Context& sys, BeadElement b) {
        for(auto& [bbIndex, neighbors] : list_) {
            if(shouldBeNeighbors(sys, bbIndex, b)) {
                neighbors.push_back(b);
            }
        }
    }
    template< typename Context >
    void removeDynamicNeighbor(Context& sys, BubbleElement bbIndex) {
        list_.erase(bbIndex);
    }
    template< typename Context >
    void removeDynamicNeighbor(Context& sys, BeadElement b) {
        for(auto& [bbIndex, neighbors] : list_) {
            auto it = std::find(neighbors.begin(), neighbors.end(), b);
            if(it != neighbors.end()) {
                neighbors.erase(it);
            }
        }
    }

    template< typename Context >
    void reset(Context& sys) {
        list_.clear();
        for(auto& bb : sys.bubbles) {
            updateNeighbors(sys, bb.sysIndex);
        }
    }

    /// Get all Cylinder neighbors of a bubble
    auto& getNeighbors(BubbleElement bbIndex) const {
        return list_.at(bbIndex);
    }

    // Filters pairs and find whether they should be neighbors.
    template< typename Context >
    bool shouldBeNeighbors(Context& sys, BubbleElement bbIndex, BeadElement pb) const {
        auto& bb = sys.bubbles[bbIndex];
        auto& fil = *static_cast<Filament*>(pb->getParent());
        // If the bubble is used in MTOC or AFM, and the filament containing the bead is attached to the MTOC or AFM, then do not add as neighbor.
        if(bb.isMTOC()) {
            auto& mtocFils = sys.mtocs[bb.getMTOCIndex()].getFilaments();
            if(std::find(mtocFils.begin(), mtocFils.end(), &fil) != mtocFils.end()) {
                return false;
            }
        }
        if(bb.isAFM()) {
            auto& afmFils = sys.afms[bb.getAFMIndex()].getFilaments();
            if(std::find(afmFils.begin(), afmFils.end(), &fil) != afmFils.end()) {
                return false;
            }
        }
        const auto dist2 = distance2(pb->coord, bb.coord);
        // If within cutoff, add as neighbor.
        return dist2 < rMax_ * rMax_;
    }
};


class TriangleFilBeadNL: public NeighborList {
public:
    using TriangleElement = medyan::StableVector<Triangle>::Index;
    using BeadElement     = Bead*;

    using TriangleElementHash = StableVectorIndexHash<Triangle>;

    SubSystem* ps = nullptr;
private:
    double     rMaxMech_ = 0.0;

    std::unordered_map< BeadElement, std::vector<TriangleElement> > listBT_, listBTMech_;
    std::unordered_map< TriangleElement, std::vector<BeadElement>, TriangleElementHash > listTB_, listTBMech_;

    template< typename A, typename B, typename HashA, typename HashB >
    void removeNeighbor_(
        A a,
        std::unordered_map< A, std::vector<B>, HashA >& listAB,
        std::unordered_map< B, std::vector<A>, HashB >& listBA
    ) {
        if(listAB.find(a) != listAB.end()) {
            for(auto b : listAB[a]) {
                auto& as = listBA[b];
                as.erase(std::remove(as.begin(), as.end(), a), as.end());
            }
            listAB.erase(a);
        }
    }

public:
    TriangleFilBeadNL(double rMax, double rMaxMech, SubSystem* ps): NeighborList(rMax), rMaxMech_(rMaxMech), ps(ps) {}

    virtual void addNeighbor(Neighbor* n) override;
    virtual void removeNeighbor(Neighbor* n) override;

    //@{
    /// The implementation of these functions calls the static version,
    /// all Triangles and Cylinders are dynamic
    virtual void addDynamicNeighbor(DynamicNeighbor* n) override { addNeighbor(n); }
    virtual void removeDynamicNeighbor(DynamicNeighbor* n) override { removeNeighbor(n); }
    //@}

    virtual void reset() override;

    /// Get all Cylinder neighbors of a triangle
    bool hasNeighbor(BeadElement     b) const { return listBT_.find(b) != listBT_.end(); }
    bool hasNeighbor(TriangleElement t) const { return listTB_.find(t) != listTB_.end(); }
    const auto& getNeighbors(BeadElement     b) const { return listBT_.at(b); }
    const auto& getNeighbors(TriangleElement t) const { return listTB_.at(t); }

    bool hasNeighborMech(BeadElement     b) const { return listBTMech_.find(b) != listBTMech_.end(); }
    bool hasNeighborMech(TriangleElement t) const { return listTBMech_.find(t) != listTBMech_.end(); }
    const auto& getNeighborsMech(BeadElement     b) const { return listBTMech_.at(b); }
    const auto& getNeighborsMech(TriangleElement t) const { return listTBMech_.at(t); }
};

} // namespace medyan

#endif
