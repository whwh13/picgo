
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

#ifndef MEDYAN_SubSystem_h
#define MEDYAN_SubSystem_h

#include <functional>
#include <optional>
#include <vector>
#include <unordered_set>

#include "common.h"

#include "Trackable.h"
#include "Database.h"
#include "Movable.h"
#include "Reactable.h"

#include "CUDAcommon.h"
#include "HybridNeighborList.h"
#include "HybridNeighborListImpl.h"
#include "SysParams.h"
#include "Chemistry/ChemSim.h"
#include "Chemistry/DissipationTracker.h"
#include "Controller/GController.h"
#include "Mechanics/ForceField/Types.hpp"
#include "Mechanics/ForceField/Membrane/MembraneBendingTypes.hpp"
#include "Mechanics/Minimizer/MinimizationTypes.hpp"
#include "Structure/Bead.h"
#include "Structure/Boundary.h"
#include "Structure/BranchingPoint.h"
#include "Structure/Bubble.h"
#include "Structure/CompartmentGrid.h"
#include "Structure/Cylinder.h"
#include "Structure/DynamicNeighbor.h"
#include "Structure/Filament.h"
#include "Structure/Linker.h"
#include "Structure/MotorGhost.h"
#include "Structure/NeighborListImpl.h"
#include "Structure/Special/AFM.h"
#include "Structure/Special/MTOC.h"
#include "Structure/SurfaceMesh/Edge.hpp"
#include "Structure/SurfaceMesh/FixedVertexAttachment.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Triangle.hpp"
#include "Structure/SurfaceMesh/Vertex.hpp"
#include "Util/StableVector.hpp"

#ifdef SIMDBINDINGSEARCH
    #include "Util/DistModule/dist_common.h"
#endif

namespace medyan {
    
// Forward declarations.
template< typename MemType > class MembraneRegion;


/// Manages all [Movables](@ref Movable) and [Reactables](@ref Reactable). Also holds all
/// [NeighborLists](@ref NeighborList) associated with chemical or mechanical interactions,
/// as well as the CompartmentGrid which contains all chemical structural information, and the
/// system Boundary.

/*! This is a class which handles all changes and information regarding the simulated system.
 *  This class operates as a top manager and provides connections between smaller parts
 *  of the system. All creation and changes go through this class and will be redirected
 *  to lower levels. See the databases for more documentation on the explicit creation of
 *  subsystem objects at initialization and during runtime.
 *
 *  This class has functions to add or remove Trackable elements from the system, as well
 *  as update Movable and Reactable instances in the system. It also can update the
 *  NeighborList container that it holds for the system.
 */
class SubSystem {

public:
    ///Default constructor does nothing
    ///constructor creates and stores an instance of HybridCylinderCylinderNL

    SubSystem() {
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        _HneighborList = new HybridCylinderCylinderNL();
#endif
}

    ~SubSystem() {}

    /// Add a Trackable to the SubSystem
    template<class T, typename ...Args>
    T *addTrackable(Args &&...args) {
	    minsT = chrono::high_resolution_clock::now();

        //create instance
        T *t = new T(forward<Args>(args)...);
        t->addToSubSystem();

        //if movable or reactable, add
//        if (t->_movable) addMovable((Movable *) t);

//        if (t->_reactable) addReactable((Reactable *) t);

        //if neighbor, add
        if (t->_dneighbor) {
	        minsN = chrono::high_resolution_clock::now();
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
            _HneighborList->addDynamicNeighbor((DynamicNeighbor *) t);
#endif
            for (auto nlist : _neighborLists)
                nlist->addDynamicNeighbor((DynamicNeighbor *) t);

            // Temporary workaround for several neighbor lists. TODO: Remove this.
            if(auto b = dynamic_cast<Bead*>(t); b) {
                if(opBubbleBeadNL.has_value()) {
                    opBubbleBeadNL->addDynamicNeighbor(*this, b);
                }
            }
        } else if (t->_neighbor) {
        	minsN = chrono::high_resolution_clock::now();
        	for (auto nlist : _neighborLists)
        		nlist->addNeighbor((Neighbor *) t);
            // Temporary workaround for several neighbor lists. TODO: Remove this.
            if(auto be = dynamic_cast<BoundaryElement*>(t); be) {
                if(opBoundaryBubbleNL.has_value()) {
                    opBoundaryBubbleNL->addNeighbor(*this, be);
                }
            }
        	mineN = chrono::high_resolution_clock::now();
        	chrono::duration<floatingpoint> elapsed_time(mineN - minsN);
        	timeneighbor += elapsed_time.count();
        }
	    mineT = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_timeT(mineT - minsT);
	    timetrackable += elapsed_timeT.count();
        return t;
    }

    /// Remove a trackable from the SubSystem
    /// Deleting the actual object should be executed by the actual
    /// callback and/or controlling function.
    template<class T>
    void removeTrackable(T *t) {
        //remove from subsystem
        t->removeFromSubSystem();

        //if movable or reactable, remove
//        if (t->_movable) removeMovable((Movable *) t);

//        if (t->_reactable) removeReactable((Reactable *) t);

        //if neighbor, remove
        if (t->_dneighbor) {
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
            _HneighborList->removeDynamicNeighbor((DynamicNeighbor *) t);
#endif
            for (auto nlist : _neighborLists)
                nlist->removeDynamicNeighbor((DynamicNeighbor *) t);

            // Temporary workaround for several neighbor lists. TODO: Remove this.
            if(auto b = dynamic_cast<Bead*>(t); b) {
                if(opBubbleBeadNL.has_value()) {
                    opBubbleBeadNL->removeDynamicNeighbor(*this, b);
                }
            }
        } else if (t->_neighbor) {
            for (auto nlist : _neighborLists)
                nlist->removeNeighbor((Neighbor *) t);
            // Temporary workaround for several neighbor lists. TODO: Remove this.
            if(auto be = dynamic_cast<BoundaryElement*>(t); be) {
                if(opBoundaryBubbleNL.has_value()) {
                    opBoundaryBubbleNL->removeNeighbor(*this, be);
                }
            }
        }
    }


    // Get current simulation time.
    auto tau() const { return ::tau(); }

    /// Get the subsystem boundary
    Boundary *getBoundary() { return _boundary; }

    /// Add a boundary to this subsystem
    void addBoundary(Boundary *boundary) { _boundary = boundary; }

    // Region in membrane
    medyan::MembraneRegion< medyan::Membrane >* getRegionInMembrane() const { return _regionInMembrane; }
    void setRegionInMembrane(medyan::MembraneRegion< medyan::Membrane >* r) { _regionInMembrane = r; }

    /// Add a neighbor list to the subsystem
    void addNeighborList(NeighborList *nl) { _neighborLists.push_back(nl); }

    /// Reset all neighbor lists in subsystem
    void resetNeighborLists();
    //create vectors of cylinder information.
    void vectorizeCylinder(medyan::SimulConfig&);


    //@{
    /// CompartmentGrid management
    CompartmentGrid* getCompartmentGrid() const { return compartmentGrid.get(); }
    //@]

    /// Update the binding managers of the system
    void updateBindingManagers(medyan::SimulConfig&);

#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    //getter for HneighborList
    HybridCylinderCylinderNL* getHNeighborList(){return _HneighborList;}

    void initializeHNeighborList(){_HneighborList->initializeHybridNeighborList();}
#endif

    const auto& getCylinderLoadForceFunc() const { return _cylinderLoadForceFunc; }
    template< typename Func >
    void setCylinderLoadForceFunc(Func&& f) { _cylinderLoadForceFunc = std::forward<Func>(f); }

    const auto& getNeighborLists() const { return _neighborLists; }

    static floatingpoint SIMDtime;
    static floatingpoint SIMDtimeV2;
    static floatingpoint HYBDtime;
	static floatingpoint timeneighbor;
	static floatingpoint timedneighbor;
	static floatingpoint timetrackable;


    //----------------------------------
    // The chemistry simulation engine.
    //----------------------------------
    std::unique_ptr<ChemSim> pChemSim;


    //----------------------------------
    // The global compartment grid.
    //----------------------------------
    std::unique_ptr<CompartmentGrid> compartmentGrid;


    //----------------------------------
    // All trackable elements in the system.
    //----------------------------------
    StableVector< Membrane > membranes;
    StableVector< Vertex >   vertices;
    StableVector< Edge >     edges;
    StableVector< Triangle > triangles;
    StableVector< MeshlessSpinVertex > meshlessSpinVertices;
    StableVector< FixedVertexAttachment > fixedVertexAttachments;

    StableVector< Bubble > bubbles;
    StableVector< MTOC >   mtocs;
    StableVector< AFM >    afms;


    //----------------------------------
    // Neighbor lists.
    //----------------------------------
    std::vector<NeighborList*> _neighborLists; ///< All neighborlists in the system
    // Used only in Hybrid binding Manager cases
    #if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        HybridCylinderCylinderNL* _HneighborList;
    #endif

    // Neighbor list between meshless vertices.
    NeighborListCellList3D meshlessSpinVertexCellList;

    // Neighbor list involving bubbles.
    std::optional<BoundaryBubbleNL> opBoundaryBubbleNL;
    std::optional<BubbleBubbleNL> opBubbleBubbleNL;
    std::optional<BubbleBeadNL> opBubbleBeadNL;

    //----------------------------------
    // Protein curvature mismatch parameters and cached states.
    //----------------------------------

    // Protein curvature mismatch parameters. Indexed curvature-mismatch setups.
    // The data is initiailized by MController initiailization.
    std::vector< ProteinCurvatureMismatchParams > proteinCurvatureMismatchParams;

    // Actual curvature mismatch data stored in the system.
    // The data is initialized by vectorization of membrane bending force field, and the actual data is populated by membrane bending energy computation.
    AllVertexCurvatureMismatchParams allVertexCurvatureMismatchParams;

    // Maps the index of a membrane diffusing species to the index of the curvature mismatch setup.
    // If a protein does not have a curvature mismatch setup, the index is set to -1.
    // The data is initialized by MController initialization.
    std::vector<int> indexMapMembraneDiffusingSpeciesToCurvatureMismatchSetup;

    // Maps the index of a desorption reaction to the index of the curvature mismatch setup.
    // If the reactant of a desorption reaction does not have a curvature mismatch setup, the index is set to -1.
    // The data is initialized by MController initialization.
    std::vector<int> indexMapDesorptionReactionToCurvatureMismatchSetup;

    //----------------------------------
    // Simulation state reports.
    //----------------------------------
    std::unique_ptr<DissipationTracker> pdt;
    MinimizationResult prevMinResult;


private:
    static const bool CROSSCHECK_SWITCH = false;
	chrono::high_resolution_clock::time_point minsN, mineN, minsT,mineT;
    Boundary* _boundary; ///< Boundary pointer
    medyan::MembraneRegion< medyan::Membrane >* _regionInMembrane; // The region inside membrane. The value is set by the controller


    //Cylinder vector
    floatingpoint* cylsqmagnitudevector = nullptr;
    static bool initialize;

    // The observer pointer of force field manager used in MController
    std::function< void(Cylinder*, ForceFieldTypes::LoadForceEnd) > _cylinderLoadForceFunc;

    chrono::high_resolution_clock::time_point minsSIMD, mineSIMD, minsHYBD, mineHYBD;

public:
    // Some initialization procedures.
    void initializeProteinCurvatureMismatch(const SimulConfig& conf) {
        const auto np = conf.chemistryData.speciesMembraneDiffusing.size();
        const auto nd = conf.chemistryData.reactionsAdsorptionDesorption.size();
        const auto ncm = conf.mechParams.proteinCurvatureMismatchSetups.size();
        // Reset indices.
        proteinCurvatureMismatchParams.clear();
        indexMapMembraneDiffusingSpeciesToCurvatureMismatchSetup.assign(np, -1);
        indexMapDesorptionReactionToCurvatureMismatchSetup.assign(nd, -1);
        // Set protein curvature mismatch parameters.
        for(int si = 0; si < ncm; ++si) {
            const auto& setup = conf.mechParams.proteinCurvatureMismatchSetups[si];
            int cmSpeciesIndex = std::find_if(
                conf.chemistryData.speciesMembraneDiffusing.begin(), conf.chemistryData.speciesMembraneDiffusing.end(),
                [&](const auto& s) { return s.name == setup.speciesName; }
            ) - conf.chemistryData.speciesMembraneDiffusing.begin();

            proteinCurvatureMismatchParams.push_back({
                cmSpeciesIndex,
                setup.kBending,
                setup.eqCurv,
                conf.chemistryData.speciesMembraneDiffusing[cmSpeciesIndex].area,
            });

            indexMapMembraneDiffusingSpeciesToCurvatureMismatchSetup[cmSpeciesIndex] = si;
        }
        // Map desorption reaction index to CM index.
        for(int ri = 0; ri < nd; ++ri) {
            // Find species index. Must exist.
            int speciesIndex = std::find_if(
                conf.chemistryData.speciesMembraneDiffusing.begin(), conf.chemistryData.speciesMembraneDiffusing.end(),
                [&](const auto& s) { return s.name == conf.chemistryData.reactionsAdsorptionDesorption[ri].speciesName2D; }
            ) - conf.chemistryData.speciesMembraneDiffusing.begin();

            indexMapDesorptionReactionToCurvatureMismatchSetup[ri] = indexMapMembraneDiffusingSpeciesToCurvatureMismatchSetup[speciesIndex];
        }
    }
};

} // namespace medyan

#endif
