
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

#ifndef MEDYAN_GController_h
#define MEDYAN_GController_h

#include <memory>
#include <vector>

#include "common.h"
#include "MathFunctions.h"
#include "Structure/CompartmentGrid.h"
#include "Structure/SubSystem.h"
#include "SysParams.h"

namespace medyan {

/// An exception to be thrown when an index/coordinate is out of bounds of the grid
class OutOfBoundsException : public exception {
    
    virtual const char* what() const throw() {
        return "An element is out of the bounds of the grid. Try adjusting minimization parameters.";
    }
};

/// An exception to be thrown when an index/coordinate is NaN
struct NaNCoordinateException : std::exception {
    
    virtual const char* what() const noexcept override {
        return "A element coordinate is NaN. Try adjusting minimization parameters.";
    }
    
};

//FORWARD DECLARATIONS
class Boundary;
class Compartment;

/// Used to control the geometry of the CompartmentGrid, as well as the geometry of
/// the entire system
/*!
 *  The GeometryController class is used by the SubSystem to control the geometry of 
 *  the simulation, which includes the CompartmentGrid geometry as well as any 
 *  [Boundaries] (@ref Boundary) that are in place. It has functionality to initialize
 *  a CompartmentGrid based on the given dimensionality as well as find the correct
 *  Compartment based on a set of coordinates.
 *  @note - The geometry controller currently only supports three-dimensional grids.
 */
class GController {

friend class Output;
    
private:
    static vector<int> _grid; ///< Size of each dimension, in compartment lengths
    static vector<floatingpoint> _compartmentSize; ///< Compartment size in each dimension
    static vector<floatingpoint> _centerGrid; ///< The center of the grid
    static vector<floatingpoint> _size;       ///< The size of the full grid in each dimension

    inline static CompartmentGrid* _compartmentGrid = nullptr; ///< The compartment grid
    
    Boundary* _boundary;   ///< The boundary that this controls
    
    SubSystem* _subSystem; ///< SubSystem ptr
    
    ///Generate all neighbors lists for each compartment in the CompartmentGrid
    void generateConnections();

public:
    ///Constructor sets SubSystem
    GController(SubSystem* ps) : _subSystem(ps) {}
    
    /// Initialize and return the grid based on input parameters
    /// Used at system initialization.
    CompartmentGrid* initializeGrid(const GeoParams& geoParams);
    
    /// Initialize and return a boundary based on input parameters
    /// Used at system initialization.
    Boundary* initializeBoundary(BoundParams::BoundaryType& BType);
    
    /// Set compartments in compartment grid as active based on boundary.
    /// Used at system initialization.
    void setActiveCompartments();

    /// Get the SubSystem ptr
    const SubSystem* getSubSystem() const {return _subSystem;}
    
    //@{
    /// Get a compartment based on coordinates or indices
    static Compartment& getCompartment(const vector<size_t> &indices);
    static Compartment& getCompartment(const vector<floatingpoint> &coords) {
        return getCompartment(mathfunc::vector2Vec<3>(coords));
    }
    static Compartment& getCompartment(const medyan::Vec<3, floatingpoint>& coords) {
        // Check if out of bounds.
        int index = 0;
        for(int i = coords.size() - 1; i >= 0; --i) {
            if(coords[i] < 0 || coords[i] >= (_compartmentSize[i] * _grid[i])) {
                throw OutOfBoundsException();
            }
            index *= _grid[i];
            index += int(coords[i] / _compartmentSize[i]);
        }

        try {
            return _compartmentGrid->getCompartment(index);
        }
        catch (const std::exception& e) {
            LOG(ERROR) << "Bad compartment access at index = " << index;
            cout << "Coords = " << coords[0] << " " << coords[1] << " " << coords[2] << endl;
            throw;
        }
    }
    static unsigned int getCompartmentID(const vector<floatingpoint> &coords);
    static Compartment& getCompartment(const int index);
    //@}
    
    // Properties (getters and setters)
    /// Get the center of the grid space
    static const vector<floatingpoint>& getCenter() {return _centerGrid;}
    static const vector<floatingpoint>& getSize() {return _size;}
    static const vector<floatingpoint>& getCompartmentSize() {return _compartmentSize;}
    
    /// Get all compartments within a given range from the specified coordinate
    /// @param ccheck - Compartment to check when initially calling this function
    /// @param compartments - List of compartments that are within range. This will be
    /// populated by the function
    static void findCompartments(const medyan::Vec<3, floatingpoint>& coords,
                                 Compartment* ccheck, floatingpoint dist,
                                 vector<Compartment*>& compartments);
    
    /// Choose a random compartment from the grid (that is activated)
    static Compartment& getRandomCompartment();
    
};

} // namespace medyan

#endif
