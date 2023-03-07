
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

#ifndef MEDYAN_Species_h
#define MEDYAN_Species_h

#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <stdexcept>

#include <type_traits>

#include "common.h"

#include "RSpecies.h"

namespace medyan {
//FORWARD DECLARATIONS
class Composite;
class CBound;

///Enumeration for Species types
enum class SpeciesType {
    unspecified,
    general, BULK, DIFFUSING, FILAMENT, BOUND, LINKER, MOTOR, BRANCHER, PLUSEND, MINUSEND,
    singleBinding, pairBinding
};

constexpr const char* speciesSuffix(SpeciesType type) {
    switch(type) {
        case SpeciesType::general:     return "{General}";
        case SpeciesType::BULK:        return "{Bulk}";
        case SpeciesType::DIFFUSING:   return "{Diffusing}";
        case SpeciesType::FILAMENT:    return "{Filament}";
        case SpeciesType::BOUND:       return "{Bound}";
        case SpeciesType::LINKER:      return "{Linker}";
        case SpeciesType::MOTOR:       return "{Motor}";
        case SpeciesType::BRANCHER:    return "{Brancher}";
        case SpeciesType::PLUSEND:     return "{PlusEnd}";
        case SpeciesType::MINUSEND:    return "{MinusEnd}";
        case SpeciesType::singleBinding: return "{SingleBinding}";
        case SpeciesType::pairBinding: return "{PairBinding}";
        default:                       return "{None}";
    }
}


/// Used to associate unique integers with character based names of Species.
/*! Often Species of the same type, let's say "Arp2/3" can be found in different forms, 
 *  for example in cytosol vs bound to a filament. The corresponding specific Species 
 *  will be distinct, however, they will share the same name, as returned by 
 *  Species::getName(). This, in turn, helps to match "Arp2/3" bound to filaments with 
 *  diffusing "Arp2/3". SpeciesNamesDB associates a unique integer with Species of the 
 *  same type (e.g. "Arp2/3", regardless whether it is bound or not), and provides
 *  conversion functions fron integer to  string (SpeciesNamesDB::intToString(int i)) 
 *  and vice versa (SpeciesNamesDB::stringToInt (string name)).
 *
 * SpeciesNamesDB has a function to generate a unique filament name given a string seed.
 * This is particularly useful for when a filament species needs a unique name, say with 
 * the seed "Actin". This unique filament name can be added or removed to match the 
 * corresponding bulk or diffusing species in the system.
 *
 * SpeciesNamesDB also has a function to generate a SpeciesSingleBinding or SpeciesPairBinding
 * name for a Compartment by concatenating two species names, a SpeciesDiffusing/SpeciesBulk with
 * a SpeciesBound, with a "-". This name can then be broken up into its counterparts when needed.
 */  
class SpeciesNamesDB {

private: 
    static unordered_map<string,int> _map_string_int;
    static vector<string> _vec_int_string;
    
    static unsigned long _num; ///<used to generate unique names
public:
    
    /// Given an integer "i", returns the string associated with that integer.
    /// Throws an out of range exception.
    static string intToString(unsigned int i) {
        
        if (i>=_vec_int_string.size())
            throw out_of_range(
                "SpeciesNamesDB::intToString(int i) index error:[" +
                 to_string(i) +"], while the vector size is " +
                 to_string(_vec_int_string.size()));
        
        return _vec_int_string[i];
    }
    
    /// Given a string "name", returns the unique integer associated with that string.
    /// If the string does not exist yet, it is created.
    static int stringToInt (string name) {
        auto mit = _map_string_int.find(name);
        
        if(mit == _map_string_int.end()){
            
            _vec_int_string.push_back(name);
            _map_string_int[name]= _vec_int_string.size()-1;
            
            return _vec_int_string.size()-1;
        }
        else
            return mit->second;
    }

    /// Generate a unique name based on a seed name
    /// (just adds integer value to end of string with a -)
    /// @note - used only for filament species.
    static string genUniqueFilName(string name) {
        string uniqueName = name + to_string(_num);
        if(_map_string_int.find(uniqueName) != _map_string_int.end()) {
        
            return uniqueName;
        }
        else {
            _vec_int_string.push_back(uniqueName);
            _map_string_int[uniqueName] = _vec_int_string.size()-1;
            
            return name + "-" + to_string(++_num);
        }
    }
    
    /// Remove the unique integer value identifier with this filament
    /// species name. @return The name without the integer ending.
    static string removeUniqueFilName(string name) {
        
        //loop through string, get to integer
        for(unsigned long i = 0; i < name.length(); i++) {
            
            if(name.substr(i,1).compare("-") == 0) {
                //return the string, cutting before this value
                return name.substr(0, i);
            }
        }
        //if there was no integer, return the name
        return name;
    }
    
    /// Generate a single or pair binding name based on two strings:
    /// one string being the binding species and the other being
    /// the bound species on the filament.
    static string genBindingName(string binding, string bound) {
        
        string name = binding + "-" + bound;
        
        _vec_int_string.push_back(name);
        _map_string_int[name] = _vec_int_string.size()-1;
        
        return name;
    }
    // Generate binding name, but for pair binding.
    static std::string genBindingName(std::string_view binding1, std::string_view bs1, std::string_view binding2, std::string_view bs2) {
        std::string res;
        res += binding1; res += '-';
        res += bs1;      res += "☣";
        res += binding2; res += '-';
        res += bs2;
        return res;
    }
    
    /// Clear the contents of the database
    static void clear() {
        _map_string_int.clear();
        _vec_int_string.clear();
    }
};
    
/// Represents chemical molecules, tracks their copy number and can be used in
/// [Reactions](@ref Reaction).
/*! This abstract class represents chemical species, such as G-Actin. As a class, 
 *  it provides Species' name, the current copy number of molecules, Species type 
 *  (e.g. SpeciesType::Bulk), and a few other characteristics of a Species.
 *  This class is synergetic with the Reaction, since Species can be added to 
 *  [Reactions](@ref Reaction).
 *  As a type, Species is composed of the following primary fields: type, name, copy 
 *  number. A subclass may added additional primary fields. Two Species are 
 *  mathematically equal, if their primary fields are equal. This means that applying 
 *  the copy constructor will guarantee that primary fields will be equal.
 *
 *  @note Each Species owns an unshared RSpecies object. However, the RSpecies field is
 *  not among the primary fields defining the Species identity (hence, the equality 
 *  operator). In particular, upon copying, the source and target Species will have 
 *  distinct RSpecies fields. However, in a fully constructed program, the source and 
 *  target species (after applying the copy constructor) should eventually be involved 
 *  in the same set of equivalent reactions. This means that their respective reactions 
 *  will return true when the equlaity operator is applied.
 *
 *  @note The Species class allows callbacks (see makeSignaling and related methods).
 */
class Species {
    
protected: //Variables
    int _molecule; ///< unique id identifying the type of the molecule (e.g. the integer id
                   ///< corresponding to "Arp2/3")
    std::unique_ptr<RSpecies> _rspecies; ///< pointer to RSpecies; Species is responsible for creating
                         ///< and destroying RSpecies
    Composite *_parent; ///< pointer to the "container" object holding this Species
                        ///< (could be a nullptr)
     ///< when cloning a reaction, whether this Species
                        ///< should be searched in speciesContainer object of the new Compartment
    

    SpeciesType type_ = SpeciesType::unspecified;

    /// Default Constructor; Should not be used by the end users - only internally
    /// By default, creates a RSpeciesReg.
    Species()  : _parent(nullptr) {
        
        _molecule=SpeciesNamesDB::stringToInt("");
        _rspecies = std::make_unique< RSpeciesReg >(*this);
    }

public:
    /// @param name - a string for the Species name associated with this Species.
    /// For example, "G-Actin" or "Arp2/3"
    /// @param n - copy number
    /// @param ulim - upper limit for this species' copy number
    /// @param sType - the type of the Species
    /// @param type - the type of RSpecies to be created
    /// @param numEvents - the number of events if using averaging
    Species (const string &name, species_copy_t n, species_copy_t ulim, SpeciesType sType, RSpeciesType type)
    
        : _molecule(SpeciesNamesDB::stringToInt(name)), _parent(nullptr), type_(sType) {
        
        //create the appropriate rspecies
        _rspecies = RSpeciesFactory::createRSpecies(*this, n, ulim, type);
    }

    /// Copy constructor
    /// @note The associated RSpecies subpart is not copied, but a new one is created.
    /// This means that the copied destination Species won't be included in any Reaction
    /// interactions of the original source Species. The species copy numbered is copied
    /// to the target. The A Species parent attriute is not copied, but set to nullptr.
    Species (const Species &rhs)
        : _molecule(rhs._molecule), _parent(nullptr), type_(rhs.type_) {
        
        //get type of rhs rspecies
        RSpeciesType t = rhs._rspecies->_type;
            
#ifdef TRACK_UPPER_COPY_N
        _rspecies = RSpeciesFactory::createRSpecies(*this, rhs.getN(), rhs.getUpperLimitForN(), t);
#else
        _rspecies = RSpeciesFactory::createRSpecies(*this, rhs.getN(), max_ulim, t);
#endif

        //transfer signal
        _rspecies->callbacks_ = std::move(rhs._rspecies->callbacks_);
        rhs._rspecies->callbacks_.clear();
            
        //set numevents if averaging
        if(t == RSpeciesType::AVG)
            ((RSpeciesAvg*)_rspecies.get())->_numEvents =
            ((RSpeciesAvg*)rhs._rspecies.get())->_numEvents;
    }
    
    /// Move constructor - makes it possible to easily add Species to STL containers,
    /// such as vector<Species> One has to be careful with vector<Species> as opposed to
    /// vector<Species*>. The latter is "safe" in terms of being internally moved around
    /// by vector.resize, etc., but vector<Species> would copy Species if a move
    /// constructor is not available. This will then destruct the original Species (with
    /// its associated RSpecies), hence, the Reaction network of the Species will be
    /// lost. Since this behavior is most likely not desired or unacceptable, the
    /// internal vector operations should only "move" the Species around, without
    /// copying. Moving transfers the RSpecies pointer from source to target,
    /// stealing resources from the source, leaving it for destruction. The Species
    /// parent attriute is moved it.
    // Jan 21 2020 (Haoran): The _rspecies should not be moved because the
    // Species reference in it cannot be changed, which would be erronous.
    Species (Species &&rhs) = delete;

    /// Assignment operator
    /// An assignment A = B copies the name of B to A. It destroys the Reaction
    /// interactions of A, and resents them to a blank value (i.e. A won't be involved
    /// in any Reactions). B's Reaction interactions are not copied. The copy number of
    /// B is copied to A. The A Species parent attriute is not copied, but set to
    /// nullptr.
    Species& operator=(const Species& rhs)  {
        _molecule = rhs._molecule;
        type_     = rhs.type_;
        
        //get type of rhs rspecies
        RSpeciesType t = rhs._rspecies->_type;
        
#ifdef TRACK_UPPER_COPY_N
        _rspecies = RSpeciesFactory::createRSpecies(*this, rhs.getN(), rhs.getUpperLimitForN(), t);
#else
        _rspecies = RSpeciesFactory::createRSpecies(*this, rhs.getN(), max_ulim, t);
#endif

        //transfer signal
        _rspecies->callbacks_ = std::move(rhs._rspecies->callbacks_);
        rhs._rspecies->callbacks_.clear();
        
        //set numevents if averaging
        if(t == RSpeciesType::AVG)
            ((RSpeciesAvg*)_rspecies.get())->_numEvents =
            ((RSpeciesAvg*)rhs._rspecies.get())->_numEvents;
        
        _parent = nullptr;
        return *this;
    }
    
    /// Move assignment is needed for the same reasons as move constructor.
    /// @see Species (Species &&rhs)
    // Jan 21 2020 (Haoran): Move assignment should be disallowed because the
    // reference of Species in the _rspecies would be wrong.
    Species& operator=(Species&& rhs) = delete;

    virtual Species* clone() {
        return new Species(*this);
    }
    
    Composite* getParent() {return _parent;}
    
    void setParent (Composite *other) {_parent=other;}
    
    bool hasParent() const {return _parent!=nullptr? true : false;}
    
    Composite* getRoot();

    auto getType() const { return type_; }

    /// Return a reference to RSpecies. Notice that value copying won't be allowed 
    /// because RSpecies is not copyable.
    RSpecies& getRSpecies () {return *_rspecies;}
    
    /// Return a reference ptr to RSpecies.
    RSpecies* getRSpeciesPtr () {return _rspecies.get();}
    
    /// Return a constant reference to RSpecies. 
    const RSpecies& getRSpecies () const {return *_rspecies;}
    
    /// Sets the copy number for this Species. 
    /// @param n should be a non-negative number, but no checking is done in run time
    /// @note The operation does not emit any signals about the copy number change.
    void setN(species_copy_t n) {_rspecies->setN(n);}
    
    /// Return the current copy number of this Species
    species_copy_t getN () const {return _rspecies->getN();}
    
#ifdef TRACK_UPPER_COPY_N
    /// Return the upper limit for the copy number of this Species
    species_copy_t getUpperLimitForN() const {return _rspecies->getUpperLimitForN();}
#endif
    
    /// Return this Species' name
    string getName() const {return SpeciesNamesDB::intToString(_molecule);}
    
    /// Return the molecule index associated with this Species' (as int)
    int getMolecule() const {return _molecule;}
    

    // Clear all callbacks associated with this Species.
    void clearSignaling () {_rspecies->clearSignaling();}
    
    /// Connect the callback, rspecies_callback to a signal corresponding to
    /// RSpecies *s.
    /// @param function<void (RSpecies *, int)> callback - a function
    /// object to be called (a slot)
    void connect(std::function<void (RSpecies *, int)> callback);
    
    /// Returns true if two Species objects are equal.
    /// This function would accept derived class of Species, such as SpeciesBulk
    /// Two Species are equal if their SType(s) are equal (i.e. are of the same class),
    /// their names are equal.
    friend bool operator ==(const Species& a, const Species& b)
    {
        if (a.getMolecule() != b.getMolecule() || typeid(a) != typeid(b))
            return false;
        return true;
    }
    
    /// Return true if two Species are not equal.         
    /// @see operator ==(const Species& a, const Species& b) above
    friend bool operator !=(const Species& a, const Species& b){
        return !(a==b);
    }
    
    // virtual methods
    
    /// Virtual destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~Species () = default;

    /// Return the full name of this Species in a string format (e.g. "Arp2/3{Bulk}"
    std::string getFullName() const { return getName() + speciesSuffix(type_); }

    //@{
    /// Change copy number
    void up() {_rspecies->up();}
    void down() {_rspecies->down();}

    
    /// Update the reaction propensities associated with this species.
    /// This will not activate the reactions if currently passivated,
    /// or update any dependents asosciated with these reactions.
    void updateReactantPropensities();
    
    /// Activate the reactions associated with this species. i.e. if
    /// passivated, will activate accordingly and update propensities
    /// and all dependents.
    void activateReactantReactions();
    
    /// Passivate the reactions associated with this species. i.e. if
    /// activated, will passivate accordingly and update propensities
    /// and all dependents.
    void passivateReactantReactions();
};

/// Used for species that can be bound to a Filament.
/// These species can not move cross-compartment.
/// Contains a pointer to a CBound object that this represents.
class SpeciesBound : public Species {
    
protected:
    CBound* _cBound = nullptr; ///< CBound object
    
public:
    /// Default constructor
    SpeciesBound()  : Species() {}
    
    /// The main constructor
    /// @param name - Example, "G-Actin" or "Arp2/3"
    /// @param n - copy number
    /// @param sType - The concrete type of the species
    SpeciesBound (const string &name, species_copy_t n=0, species_copy_t ulim=1, SpeciesType sType=SpeciesType::BOUND)
        :  Species(name, n, ulim, sType, RSpeciesType::REG) {};
    
    /// Copy constructor
    SpeciesBound (const SpeciesBound &rhs)  : Species(rhs){}
    
    /// Move constructor
    SpeciesBound (SpeciesBound &&rhs) = delete;

    /// Regular Assignment
    SpeciesBound& operator=(const SpeciesBound& rhs)  {
        Species::operator=(rhs);
        return *this;
    }
    
    /// Move assignment
    SpeciesBound& operator=(SpeciesBound&& rhs) = delete;

    virtual SpeciesBound* clone() {
        return new SpeciesBound(*this);
    }
    
    /// Default destructor
    ~SpeciesBound () noexcept {};
    
    ///Setter and getter for CBound
    void setCBound(CBound* cBound) {_cBound = cBound;}
    CBound* getCBound() {return _cBound;}
    
    ///remove cBound ptr
    void removeCBound() {_cBound = nullptr;}

};


/// Print self into an iostream
ostream& operator<<(ostream& os, const Species& s);

} // namespace medyan

#endif
