
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

#ifndef MEDYAN_CMonomer_h
#define MEDYAN_CMonomer_h

#include "common.h"

#include "Species.h"

/// Represents a container for all Species that could be contained in a
/// particular filament element at a given position.
/*!
 *  CMonomer provides a container to hold all Species that are held at a given
 *  filament position. The species are held in an standard vector.
 */

namespace medyan {

// Forward declarations.
class Compartment;


class CMonomer {
    
friend class ChemManager;
friend class CCylinder;

public:
    struct SpeciesSubarrayIndex {
        medyan::Index start = 0;
        medyan::Index size = 0;
    };

private:
    //@{
    /// Species array
    std::vector< Species* >         _speciesFilament;
    std::vector< SpeciesBound* >    _speciesBound;
    //@}
    
    ///Filament type that this monomer is in
    short _filamentType;
    
    //@{
    /// Species index vectors
    /// These index vectors are used to access the correct species in the actual
    /// species arrays. Each index in the species index vector corresponds to an
    /// offset for that species in the species array.
    // To access the index, use xxxIndex[filament-type][species-xxx-type].
    inline static std::vector<std::vector<SpeciesSubarrayIndex>> speciesFilamentIndex_;
    inline static std::vector<std::vector<SpeciesSubarrayIndex>> speciesBoundIndex_;
    //@}
    
    ///Number of species for each filament type
    static vector<short> _numFSpecies;
    static vector<short> _numBSpecies;
    
public:
    /// Constructor does nothing but memset arrays
    CMonomer(short filamentType) :
        _filamentType(filamentType),
        _speciesFilament(_numFSpecies[filamentType]),
        _speciesBound(_numBSpecies[filamentType])
    {}
    
    
    /// Copy constructor
    /// This constructor will create a new CMonomer, identical to the copied, in a new
    /// compartment. The original species will remain intact, while new identical
    /// species will be initialized.
    CMonomer(const CMonomer& rhs, Compartment* c);
    /// Assignment is not allowed
    CMonomer& operator=(CMonomer &rhs) = delete;
    
    /// Clone a CMonomer. Transfers all [CBounds] (@ref CBound) to the new CMonomer.
    virtual CMonomer* clone(Compartment* c) {
        return new CMonomer(*this, c);
    }
    
    ///Print the Species
    void print();

    //@{
    /// Get Species at a specific index
    /// @note no check on this index. The index value of a species is stored in the
    /// chemical initializer when all reactions are initialized from the chemical input
    /// file.
    Species* speciesFilament (int index);
    Species* speciesPlusEnd  (int index);
    Species* speciesMinusEnd (int index);
    
    SpeciesBound* speciesBound    (int index);
    SpeciesBound* speciesLinker   (int index);
    SpeciesBound* speciesMotor    (int index);
    SpeciesBound* speciesBrancher (int index);
    //@}
    
    //@{
    /// Get the active species of the corresponding type
    /// @note these functions will return -1 if the corresponding type has no marked
    /// species at this point.
    short activeSpeciesFilament();
    short activeSpeciesPlusEnd();
    short activeSpeciesMinusEnd();
    
    short activeSpeciesBrancher();
    //@
    
    /// Check the consistency of the CMonomer for debugging.
    bool isConsistent();

    // Set species index vectors. Should only be called at the beginning of simulation.
    static void setSpeciesFilamentIndex(std::vector<std::vector<SpeciesSubarrayIndex>> speciesFilamentIndex) {
        speciesFilamentIndex_ = std::move(speciesFilamentIndex);
    }
    static void setSpeciesBoundIndex(std::vector<std::vector<SpeciesSubarrayIndex>> speciesBoundIndex) {
        speciesBoundIndex_ = std::move(speciesBoundIndex);
    }
};

} // namespace medyan

#endif
