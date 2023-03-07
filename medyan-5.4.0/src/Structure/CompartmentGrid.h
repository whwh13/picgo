
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

#ifndef MEDYAN_CompartmentGrid_h
#define MEDYAN_CompartmentGrid_h

#include <cmath>

#include "common.h"
#include "Rand.h"
#include "SysParams.h"
#include "Structure/CellList.hpp"
#include "Structure/Compartment.h"
#include "Util/StableVector.hpp"

namespace medyan {

// Forward declarations.
class Cylinder;
class ChemSim;
class Edge;
class Triangle;
class Vertex;


/// A simple n-dimensional grid of Compartment objects.

/*!
 *  The CompartmentGrid class is a singleton grid of Compartment objects, each of them 
 *  seperately holding internal and diffusion reactions, species information, as well as
 *  spatial information. This class is n-dimensional, and the dimension is specified at 
 *  runtime.
 *  All compartments within CompartmentGrid are indexed by the GController, and
 *  this class is responsible for assignment of compartment neighbors upon 
 *  initialization. All initialization of the CompartmentGrid should be done through 
 *  GController.initializeGrid().
 *
 *  The _prototype_compartment is used such that all reactions and species can be added 
 *  to the prototype, and this configuration can be copied to all compartments in the grid.
 */
class CompartmentGrid {
private:
    Compartment _prototype_compartment; ///< Prototype compartment, to be configured
                                        ///< before initialization
    
    SpeciesPtrContainerVector _bulkSpecies;    ///< Bulk species in this grid
    ReactionPtrContainerVector _bulkReactions; ///< Bulk reactions in this grid
    
public:
    // Geometric parameters.
    //-------------------------------------------------------------------------

    // Number of compartments in each dimension.
    std::array< Size, 3 >          shape {};
    // Side lengths of a compartment. Indexed by axis (x, y, z).
    std::array< floatingpoint, 3 > compartmentLengths {};
    // Areas of a compartment. Indexed by planes (yz, zx, xy).
    std::array< floatingpoint, 3 > compartmentAreas {};
    floatingpoint                  compartmentVolume = 0;
    // Side lengths of the entire grid. Indexed by axis (x, y, z).
    std::array< floatingpoint, 3 > gridLengths {};
    floatingpoint                  gridVolume = 0;
    // Fraction of span when choosing random coordinates in the entire grid.
    std::array<std::array<floatingpoint, 3>, 2> fracGridSpan {{{0, 0, 0}, {1, 1, 1}}};

    // All compartments (except for the proto compartment) managed by this class.
    std::vector<std::unique_ptr<Compartment>> compartmentList;

    medyan::CellListManager< Cylinder*, Compartment* > cylinderCellList;
    medyan::CellListManager< medyan::StableVector<medyan::Triangle>::Index, Compartment* > triangleCellList;
    medyan::CellListManager< medyan::StableVector<medyan::Edge>::Index,     Compartment* > edgeCellList;
    medyan::CellListManager< medyan::StableVector<medyan::Vertex>::Index,   Compartment* > vertexCellList;

    /// Constructor, creates a number of Compartment instances
    CompartmentGrid(const GeoParams& geoParams) {
        // Set geometry.
        shape = {
            geoParams.NX,
            geoParams.NY,
            geoParams.NZ,
        };
        compartmentLengths = {
            geoParams.compartmentSizeX,
            geoParams.compartmentSizeY,
            geoParams.compartmentSizeZ,
        };
        gridLengths = {
            geoParams.NX * geoParams.compartmentSizeX,
            geoParams.NY * geoParams.compartmentSizeY,
            geoParams.NZ * geoParams.compartmentSizeZ,
        };
        compartmentAreas = {
            compartmentLengths[1] * compartmentLengths[2],
            compartmentLengths[2] * compartmentLengths[0],
            compartmentLengths[0] * compartmentLengths[1],
        };
        compartmentVolume = compartmentLengths[0] * compartmentLengths[1] * compartmentLengths[2];
        gridVolume = gridLengths[0] * gridLengths[1] * gridLengths[2];
        fracGridSpan = geoParams.fracGridSpan;

        const auto numCompartments = shape[0] * shape[1] * shape[2];
        // Add children.
        for(Index i = 0; i < numCompartments; ++i) {
            compartmentList.push_back(std::make_unique<Compartment>());
            auto c = compartmentList.back().get();
            c->_ID = i;

            cylinderCellList.addHead(c, c->cylinderCell);
            c->cylinderCell.manager = &cylinderCellList;
            triangleCellList.addHead(c, c->triangleCell);
            c->triangleCell.manager = &triangleCellList;
            edgeCellList.addHead(c, c->edgeCell);
            c->edgeCell.manager = &edgeCellList;
            vertexCellList.addHead(c, c->vertexCell);
            c->vertexCell.manager = &vertexCellList;
        }
    }

    ~CompartmentGrid() {
        for(auto& C : compartmentList) {
            C->clearReactions();
        }
    }
    
    /// Get compartments that this grid holds
    const auto& getCompartments() const {
        return compartmentList;
    }
    
    // Get a compartment at certain index, index3 or coordinate.
    Compartment& getCompartment(Index index) const {
        
        return *compartmentList[index];
    }
    Compartment& getCompartment(const std::array<Index, 3>& index3) const {
        return getCompartment(getCompartmentIndex(index3));
    }
    Compartment& getCompartment(Index ix, Index iy, Index iz) const {
        return getCompartment(getCompartmentIndex(ix, iy, iz));
    }
    // Find the compartment with a certain coordinate.
    Compartment& getCompartment(const Vec<3, floatingpoint>& coord) const {
        return getCompartment(getCompartmentIndex(coord));
    }

    // Find the compartment 1D index with 3D index.
    // Note:
    // - The smaller index is the fastest changing one (fortran ordering).
    // - The index is not checked for out of bounds.
    Index getCompartmentIndex(Index ix, Index iy, Index iz) const {
        return ix + shape[0] * (iy + shape[1] * iz);
    }
    Index getCompartmentIndex(const std::array<Index, 3>& index3) const {
        return getCompartmentIndex(index3[0], index3[1], index3[2]);
    }

    // Find the compartment 3D index at a certain coordinate.
    std::array<Index, 3> getCompartmentIndex3(const Vec<3, floatingpoint>& coord) const {
        const auto i = static_cast<Index>(std::floor(coord[0] / compartmentLengths[0]));
        const auto j = static_cast<Index>(std::floor(coord[1] / compartmentLengths[1]));
        const auto k = static_cast<Index>(std::floor(coord[2] / compartmentLengths[2]));
        if(
            i < 0 || i >= shape[0] ||
            j < 0 || j >= shape[1] ||
            k < 0 || k >= shape[2]
        ) {
            throw std::runtime_error("Coordinate out of range of the grid.");
        }
        return {i, j, k};
    }

    // Find the compartment flat index at a certain coordinate.
    Index getCompartmentIndex(const Vec<3, floatingpoint>& coord) const {
        return getCompartmentIndex(getCompartmentIndex3(coord));
    }

    /// Set all compartments as active. Used at initialization
    void setAllAsActive() const {
        for (auto& c : getCompartments()) c->setAsActive();
    }
    
    /// Get name of this compartment grid
    std::string getFullName() const {return string("CompartmentGrid");};
    
    /// Get the protocompartment from this grid, in order to configure and then initialize
    Compartment& getProtoCompartment() {return _prototype_compartment;}
    const Compartment& getProtoCompartment() const {return _prototype_compartment;}
    
    /// Add reactions to all compartments in the grid
    /// @param - chem, a ChemSim object that controls the reaction algorithm
    void addChemSimReactions(medyan::ChemSim* chem);

    /// Print properties of this grid
    void printSelf() const {
        cout << getFullName() << endl;
        cout << "Number of Compartment objects: " << compartmentList.size() << endl;
        for(auto &c : compartmentList)
            c->printSelf();
    }

    const auto& getSpeciesBulk() const {
        return _bulkSpecies;
    }
    /// Add a bulk species to this grid
    template<typename ...Args>
    Species* addSpeciesBulk (Args&& ...args) {
        _bulkSpecies.addSpecies<Species>(forward<Args>(args)...);
        return _bulkSpecies.findSpeciesByIndex(_bulkSpecies.size() - 1);
    }
    
    /// Remove bulk species
    void removeSpeciesBulk(const string& name) {_bulkSpecies.removeSpecies(name);}
    
    /// Bulk species finder functions
    Species* findSpeciesBulkByName(const string& name) const {
        return _bulkSpecies.findSpeciesByName(name);
    }
    Species* findSpeciesBulkByMolecule(int molecule) const {
        return _bulkSpecies.findSpeciesByMolecule(molecule);
    }
    
    /// Add a bulk reaction to this compartment grid
    template<unsigned short M, unsigned short N, typename ...Args>
    ReactionBase* addBulkReaction (Args&& ...args)
    {
        ReactionBase *r = _bulkReactions.addReaction<M,N>(forward<Args>(args)...);
        return r;
    }
    
    /// Add an bulk reaction to this compartment grid
    /// @param species, rate - specifying the species and rate that should be assigned
    template<template <unsigned short M, unsigned short N> class RXN, unsigned short M, unsigned short N>
    ReactionBase* addBulk(initializer_list<Species*> species, float rate) {
        ReactionBase *r = _bulkReactions.add<RXN,M,N>(species,rate);
        return r;
    }
    
    /// Add a unique bulk reaction pointer to this compartment
    ReactionBase* addBulkReactionUnique (unique_ptr<ReactionBase> &&reaction) {
        ReactionBase *r = _bulkReactions.addReactionUnique(move(reaction));
        return r;
    }
    
    /// Remove all bulk reactions that have a given species
    /// @param s - species whose reactions should be removed
    void removeBulkReactions (Species* s) {_bulkReactions.removeReactions(s);}

    /// Remove a bulk reaction
    void removeBulkReaction(ReactionBase *r) {_bulkReactions.removeReaction(r);}
    
    
    /// Count the number of diffusing species with a given name
    species_copy_t countDiffusingSpecies(const string& name);
    /// Count the number of bulk species with a given name
    species_copy_t  countBulkSpecies(const string& name);

    //----------------------------------
    // Auxiliary generators.
    //----------------------------------

    // Get random coordinates in a given compartment.
    auto getRandomCoordinatesIn(const Compartment& c) const {
        auto& cc = c.coordinates();
        return Vec<3, floatingpoint> {
            cc[0] + compartmentLengths[0] * Rand::randfloatingpoint(-1,1) / 2,
            cc[1] + compartmentLengths[1] * Rand::randfloatingpoint(-1,1) / 2,
            cc[2] + compartmentLengths[2] * Rand::randfloatingpoint(-1,1) / 2,
        };
    }
    auto getRandomCenterCoordinatesIn(const Compartment& c) const {
        return getRandomCoordinatesIn(c);
    }

    // Get random coordinates in this grid.
    auto getRandomCoordinates() const {
        return Vec<3, floatingpoint> {
            Rand::randfloatingpoint(fracGridSpan[0][0], fracGridSpan[1][0]) * gridLengths[0],
            Rand::randfloatingpoint(fracGridSpan[0][1], fracGridSpan[1][1]) * gridLengths[1],
            Rand::randfloatingpoint(fracGridSpan[0][2], fracGridSpan[1][2]) * gridLengths[2],
        };
    }
    auto getRandomCenterCoordinates() const {
        return getRandomCoordinates();
    }
};

} // namespace medyan

#endif
