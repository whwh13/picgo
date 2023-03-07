
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

#ifndef MEDYAN_Compartment_h
#define MEDYAN_Compartment_h

#include <algorithm> // count_if
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "common.h"

#include "SpeciesContainer.h"
#include "ReactionContainer.h"
#include "BindingManager.h"
#include "HybridBindingSearchManager.h"
#include "Composite.h"
#include "Chemistry/ChemSim.h"
#ifdef SIMDBINDINGSEARCH
#include "Util/DistModule/dist_driver.h"
#include "Util/DistModule/dist_coords.h"
#endif
#include "MathFunctions.h"
#include "Structure/CellList.hpp"
#include "Util/StableVector.hpp"


namespace medyan {

// Forward declarations.
class BoundaryElement;
class Bead;
class Cylinder;
class SubSystem;
class Triangle;
class Edge;
class Vertex;

/// A container or holding Species and [Reactions](@ref Reactions).

/*! The Compartment class is a container for Species, internal [Reactions](@ref Reactions),
 *  and diffusion [Reactions](@ref Reactions) that can occur. A Compartment object keeps
 *  track of the above while also holding pointers to its neighbors in order to generate
 *  diffusion [Reactions](@ref Reactions) and control other interactions.
 *
 *  The Compartment also keeps Trackable elements in the SubSystem that are in its space, including
 *  [Beads](@ref Bead), [Cylinders](@ref Cylinder), and [BoundaryElements](@ref BoundaryElement).
 *
 *  Lastly, the Compartment holds a container of FilamentBindingManager for updating
 *  binding reactions local to this compartment space only.
 */

class Compartment : public Composite {

public:
    enum class SliceMethod {
        membrane, cylinderBoundary
    };
    enum class ActivateReason {
        Whole, Membrane
    };

    static constexpr int numNeighborDirections = 6;
    static constexpr int getOppositeNeighborDirection(int dir) {
        return dir^1;
    }

protected:
    ///CHEMICAL CONTAINERS
    SpeciesPtrContainerVector _species;  ///< Container with all species
                                         ///< in this compartment
    ReactionPtrContainerVector _internal_reactions; ///< Container with all internal
                                                    ///< reactions in compartment

    // Contains diffusion reactions from this compartment to neighbor compartments. Thus, inward diffusion reactions are stored in neighbor compartments.
    ReactionPtrContainerVector _diffusion_reactions;

    // Diffusion coefficients of Species in compartment. Map from molecule id -> diffusion coefficient.
    std::unordered_map<int, float> diffusionCoefficients_;



    /// All binding managers for this compartment
    vector<unique_ptr<medyan::FilamentBindingManager>> _bindingManagers;
    vector<unique_ptr<medyan::FilamentBindingManager>> _branchingManagers;

#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    //Each compartment has a single instance of hybridbindingmanagers
    medyan::HybridBindingSearchManager* _bindingsearchManagers = NULL;
#endif
    ///ELEMENT CONTAINERS
    unordered_set<BoundaryElement*> _boundaryElements; ///< Set of boundary element
                                                       ///< that are in this compartment

    unordered_set<Bead*> _beads; ///< Set of beads that are in this compartment

    // Indices of neighboring compartments. The neighbors are in the order (x-, x+, y-, y+, z-, z+). Sentinel value -1 indicates that a neighbor does not exist.
    std::array< int, numNeighborDirections > neighborIndices_ { -1, -1, -1, -1, -1, -1 };

    // touch along faces
    vector<Compartment*> _enclosingneighbours; ///< Neighbors that envelop the compartment
    vector<Compartment*> _uniquepermuteneighbours; //Neighbors tht help to parse
    // through unique pairs of compartments to get all necessary binding sites.
    vector<short> _uniquepermuteneighboursstencil;
    vector<short> _enclosingneighboursstencil;///< enclosing compartments  have unique
    // position IDs between 0-26. This ID immediately helps us determine relative
    // position of the Enclosing Compartment with respect to the current compartment.


    ///OTHER COMPARTMENT PROPERTIES
    bool _activated = false; ///< The compartment is activated for diffusion

    floatingpoint _volumeFrac = 1.0; ///< The volume fraction inside the membrane/boundary
    ///< Might be changed to a list or a map when more membranes are involved
    array<floatingpoint, 6> _partialArea {{1.0, 1.0, 1.0, 1.0, 1.0, 1.0}}; ///< The area
    // inside the cell membrane
    ///<In the order of x, y and z, from smaller coordinate value neighbor to larger coordinate value
    ///< Might be changed to a list of arrays or a map of arrays when more membranes are involved

public:
    short _ID = -1;

    // Center coordinate of this compartment.
    medyan::Vec<3, floatingpoint> centerCoord {};

    medyan::CellListHeadUser< Cylinder*, Compartment* > cylinderCell; // Cell of cylinders
    medyan::CellListHeadUser< medyan::StableVector<medyan::Triangle>::Index, Compartment* > triangleCell; // Cell of triangles
    medyan::CellListHeadUser< medyan::StableVector<medyan::Edge>::Index,     Compartment* > edgeCell;     // Cell of edges
    medyan::CellListHeadUser< medyan::StableVector<medyan::Vertex>::Index,   Compartment* > vertexCell;   // Cell of vertices

    /// Default constructor, only takes in number of dimensions
    Compartment() : _species(), _internal_reactions(),
        _diffusion_reactions()
    {}

    /// Constructor which clones another compartment
    Compartment(const Compartment &C) : _species(), _internal_reactions(),
                                        _diffusion_reactions()
    {
        C.cloneSpecies(this);
        C.cloneReactions(this);
        diffusionCoefficients_ = C.diffusionCoefficients_;
        _activated = C._activated;

        // Should eventually clone beads, cylinders, boundary elements.... not clear yet
    }

    /// Assignment operator
    Compartment& operator=(const Compartment &other);

    /// Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~Compartment() noexcept
    {
        clearReactions();
        clearSpecies();

    }

    /// get ID
    int getId() const {return _ID;}
    /// Applies SpeciesVisitor v to every Species* object directly owned by this node.
    /// This method needs to be overriden by descendent classes that contain Species.
    virtual bool apply_impl(SpeciesVisitor &v) override;

    /// Applies ReactionVisitor v to every ReactionBase* object directly owned by this
    /// node. This method needs to be overriden by descendent classes that contain
    /// ReactionBase.
    virtual bool apply_impl(ReactionVisitor &v) override;

    ///Set a compartment as active. Used at initialization.
    void setAsActive() {_activated = true;}
    void setAsInactive() { _activated = false; }


    ///Check if compartment is activated
    bool isActivated() {return _activated;}

    ///Setter and getter for fs
    void setCoordinates(const medyan::Vec<3, floatingpoint>& coords) { centerCoord = coords; }
    const auto& coordinates() const { return centerCoord; }



    /// Removes all reactions from this compartment, diffusing and internal
    void clearReactions() {

        _internal_reactions.reactions().clear();
        _diffusion_reactions.reactions().clear();
    }

    /// Clear all species from this compartment
    void clearSpecies() {
        _species.species().clear();
    }


    /// Returns true if two Compartment objects are equal.
    /// Two Compartment objects are equal if each contains analogous Species and
    /// Reaction objects, in the same order
    friend bool operator==(const Compartment& a, const Compartment& b);

    /// Return true if two Compartment are not equal.
    /// @see operator ==(const Compartment& a, const Compartment& b) above
    friend bool operator !=(const Compartment& a, const Compartment& b)
    {
        return !(a==b);
    }

    /// Returns compartment name
    virtual string getFullName() const override {return "Compartment";};
    /// Returns the number of species in this compartment
    size_t numberOfSpecies() const {return _species.species().size();}
    /// Returns the number of internal reactions in this compartment
    size_t numberOfInternalReactions() const {
        return _internal_reactions.reactions().size();
    }
    /// Returns the total number of reactions in this compartment, diffusing and
    /// internal
    size_t numberOfReactions() const {
        return _internal_reactions.reactions().size() +
               _diffusion_reactions.reactions().size();
    }

    //@{
    /// Species finder functions
    Species* findSpeciesByName(const string &name) {
        return _species.findSpeciesByName(name);
    }
    Species* findSpeciesByIndex (Index index) const {
        return _species.findSpeciesByIndex(index);
    }
    Species* findSpeciesByMolecule (int molecule) const {
        return _species.findSpeciesByMolecule(molecule);
    }
    Species* findSimilarSpecies (const Species &s) {
        return _species.findSimilarSpecies(s);
    }
    //@}

    ///Remove species from this compartment
    size_t removeSpecies(Species* species) {return _species.removeSpecies(species);}

    /// Finds a similar internal reaction, see ReactionBase function
    ReactionBase* findSimilarInternalReaction (const ReactionBase &r) {
        return _internal_reactions.findSimilarReaction(r);
    }

    /// Finds a similar diffusion reaction, see ReactionBase function
    ReactionBase* findSimilarDiffusionReaction (const ReactionBase &r) {
        return _diffusion_reactions.findSimilarReaction(r);
    }

    /// Remove all diffusion reactions that have a given species
    /// @param s - species whose reactions should be removed
    virtual void removeDiffusionReactions (Species* s) {
        _diffusion_reactions.removeReactions(s);
    }

    /// Remove all internal reactions that have a given species
    /// @param s - species whose reactions should be removed
    virtual void removeInternalReactions (Species* s) {
        _internal_reactions.removeReactions(s);
    }

    /// Remove a diffusion reaction
    virtual void removeDiffusionReaction(ReactionBase *r) {
        _diffusion_reactions.removeReaction(r);
    }

    /// Remove an internal reaction
    virtual void removeInternalReaction(ReactionBase *r) {
        _internal_reactions.removeReaction(r);
    }

    /// Add a unique species pointer to this compartment
    Species* addSpeciesUnique (unique_ptr<Species> &&species, float diff_rate = -1.0) {
        Species *sp = _species.addSpeciesUnique(move(species));
        sp->setParent(this);
        diffusionCoefficients_[sp->getMolecule()]=diff_rate;
        return sp;
    }

    /// Add a unique internal reaction pointer to this compartment
    ReactionBase* addInternalReaction(std::unique_ptr<ReactionBase> reaction) {
        ReactionBase *r = _internal_reactions.addReactionUnique(std::move(reaction));
        r->setParent(this);
        return r;
    }

    /// Add an internal reaction pointer to this compartment. Make unique
    ReactionBase* addInternalReaction (ReactionBase* r) {
        _internal_reactions.addReactionUnique(unique_ptr<ReactionBase>(r));
        r->setParent(this);
        return r;
    }

    /// Add a diffusion reaction pointer to this compartment.
    void addDiffusionReaction(std::unique_ptr<ReactionBase> r) {
        r->setParent(this);
        _diffusion_reactions.addReactionUnique(std::move(r));
    }


    /// Add a bound species to this compartment
    /// @param args - any number of SpeciesBound objects
    template<typename ...Args>
    SpeciesBound* addSpeciesBound(Args&& ...args) {
        SpeciesBound *sp =
        (SpeciesBound*)(_species.addSpecies<SpeciesBound>(forward<Args>(args)...));
        sp->setParent(this);
        diffusionCoefficients_[sp->getMolecule()]=-1.0;
        return sp;
    }

    /// Add an internal reaction to this compartment
    template<unsigned short M, unsigned short N, typename ...Args>
    ReactionBase* addInternalReaction (Args&& ...args) {
        ReactionBase *r = _internal_reactions.addReaction<M,N>(forward<Args>(args)...);
        r->setParent(this);
        return r;
    }

    /// Add an internal reaction to this compartment
    /// @param species, rate - specifying the species and rate that should be assigned
    template<template <unsigned short M, unsigned short N> class RXN, unsigned short M, unsigned short N>
    ReactionBase* addInternal(initializer_list<Species*> species, float rate) {
        ReactionBase *r = _internal_reactions.add<RXN,M,N>(species,rate);
        r->setParent(this);
        return r;
    }

    ///Add branching manager to this compartment. Used only in SIMD case
    void addBranchingBindingManager(medyan::FilamentBindingManager* m){
    	_branchingManagers.emplace_back(m);
    }

    /// Add a binding manager to this compartment
    void addFilamentBindingManager(medyan::FilamentBindingManager* m) {
        _bindingManagers.emplace_back(m);
    }
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    void addHybridBindingSearchManager(medyan::HybridBindingSearchManager* bsm){
        if(_bindingsearchManagers == NULL)
            _bindingsearchManagers = bsm;
        else{
            cout<<"Hybrid Binding Search Manager exists. Compartment cannot have multiple"
                    " hybrid binding search managers. Exiting.."<<endl;
            exit(EXIT_FAILURE);
        }
    }

    medyan::HybridBindingSearchManager* getHybridBindingSearchManager(){
        return _bindingsearchManagers;
    }
#endif
#ifdef SIMDBINDINGSEARCH
    dist::Coords bscoords;
    vector<dist::Coords> bscoords_section;
	vector<dist::Coords> bscoords_section_linker;
	vector<dist::Coords> bscoords_section_motor;

    vector<int> Cyldcindexvec;
    vector<int> CylcIDvec;

    template<bool LinkerorMotor>
    dist::Coords& getSIMDcoordsV3(short i, short filamentType){
        if(LinkerorMotor)
            return bscoords_section_linker[filamentType*27 + i];
        else
            return bscoords_section_motor[filamentType*27 + i];
    }
    /*Each compartment is partitioned into 27 sub-volumes. Binding sites are allocated
    to relevant sub-volumes. It is worth noting that the 27 volumes  can be overlapping.
    Binding distance determines if the volumes are overlapping or not.*/
    vector<floatingpoint> partitionedcoordx[27], partitionedcoordy[27], partitionedcoordz[27];
    vector<uint32_t>  cindex_bs_section[27];
    vector<uint32_t> finfo_bs_section[27];

    void SIMDcoordinates_section();
    void SIMDcoordinates4linkersearch_section(bool isvectorizedgather);
    void SIMDcoordinates4motorsearch_section(bool isvectorizedgather);
    void getpartition3Dindex(int (&indices)[3], vector<floatingpoint> coord);
    template<bool rMaxvsCmpsize>
    void getpartitionindex(int (&indices)[3], vector<floatingpoint> coord,
                                floatingpoint (&cmpcornercoords)[6]);

	void addcoord(vector<floatingpoint> coord, uint32_t index, uint32_t cylfinfo, short i);
    bool checkoccupancy(Cylinder* cyl, short it, short _filamentType, short bstatepos);
    bool checkoccupancy(vector<vector<bool>>& boundstate, short bstatepos, int pos);
    void addcoordtopartitons(int (&pindices)[3], vector<floatingpoint> coord, uint32_t
    index, uint32_t cylfinfo);
    void addcoordtopartitons_smallrmax(int (&pindices)[3], vector<floatingpoint> coord,
                                  uint16_t index, uint32_t cylfinfo);
    template<bool rMaxvsCmpsize>
    void addcoordtorMaxbasedpartitons(int (&pindices)[3], vector<floatingpoint> coord,
                                       uint32_t index, uint32_t cylfinfo);

    void deallocateSIMDcoordinates();
#endif
    /// Get binding managers for this compartment
    vector<unique_ptr<medyan::FilamentBindingManager>>& getFilamentBindingManagers() {
        return _bindingManagers;
    }

    //Get BranchingBindingManager
    vector<unique_ptr<medyan::FilamentBindingManager>>& getBranchingManagers() {
        return _branchingManagers;
    }

    /// Get a specific motor binding manager from this compartment
    MotorBindingManager* getMotorBindingManager(int motorType) {

        MotorBindingManager* mp;

        for(auto it = _bindingManagers.begin(); it != _bindingManagers.end(); it++) {

            //find the correct manager and type
            if((mp = dynamic_cast<MotorBindingManager*>((*it).get())) && (*it)->getBoundInt() == motorType)
                return mp;
        }

        return nullptr;
    }

    ///Add a boundary element to this compartment
    void addBoundaryElement(BoundaryElement* be) {_boundaryElements.insert(be);}

    ///Remove a boundary element from this compartment
    ///@note does nothing if boundary element is not in compartment
    void removeBoundaryElement(BoundaryElement* be) {
        auto it = _boundaryElements.find(be);
        if(it != _boundaryElements.end()) _boundaryElements.erase(it);
    }
    ///Check if boundary element is in this container
    bool hasBoundaryElement(BoundaryElement* be) {
        auto it = _boundaryElements.find(be);
        return (it != _boundaryElements.end());
    }

    ///get the boundary elements in this compartment
   unordered_set<BoundaryElement*>& getBoundaryElements() {return _boundaryElements;}

    // Get the cylinders in this compartment.
    auto getCylinders() const { return cylinderCell.manager->getElements(cylinderCell); }

    // Get triangles in this compartment.
    auto getTriangles() const { return triangleCell.manager->getElements(triangleCell); }

    // Get edges in this compartment.
    auto getEdges() const { return edgeCell.manager->getElements(edgeCell); }

    // Get vertices in this compartment.
    auto getVertices() const { return vertexCell.manager->getElements(vertexCell); }
    
    /// Get the diffusion rate of a species
    /// @param - species_name, a string
    auto getDiffusionCoefficient(int molecule) const {
        return diffusionCoefficients_.at(molecule);
    }

    /// Set the diffusion rate of a species in the compartment
    void setDiffusionCoefficient(Species *sp, float diff_rate) {
        int molecule = sp->getMolecule();
        diffusionCoefficients_[molecule]=diff_rate;
    }

    /// Set the diffusion rate of a species in the compartment
    /// @param - molecule, an integer representing the species in the speciesDB
    void setDiffusionCoefficient(int molecule, float diff_rate) {
        diffusionCoefficients_[molecule]=diff_rate;
    }

    /// Set the diffusion rate of a species in the compartment
    /// @param - species_name, a string
    void setDiffusionCoefficient(string species_name, float diff_rate) {
        int molecule = SpeciesNamesDB::stringToInt(species_name);
        diffusionCoefficients_[molecule]=diff_rate;
    }

    /// Add a neighboring compartment to this compartments list of neighbors
    void setNeighbor(int compartmentIndex, int dirIndex) {
        neighborIndices_[dirIndex] = compartmentIndex;
    }

    /// Remove a neighboring compartment
    void resetNeighbor(int dirIndex) {
        neighborIndices_[dirIndex] = -1;
    }

    /// Add a neighboring compartment to this compartments list of neighbors
    void addenclosingNeighbour(Compartment *comp, int stencilpos) {
        auto nit = find(_enclosingneighbours.begin(),_enclosingneighbours.end(), comp);
        if(nit==_enclosingneighbours.end()) {
            _enclosingneighbours.push_back(comp);
            _enclosingneighboursstencil.push_back(stencilpos);
        }
        else
            throw runtime_error(
                    "Compartment::addenclosingNeighbour(): Compartment is already a "
                    "neighbour");
    }


    void adduniquepermuteNeighbour(Compartment *comp, int stencilpos) {
        auto nit = find(_uniquepermuteneighbours.begin(),_uniquepermuteneighbours.end(), comp);
        if(nit==_uniquepermuteneighbours.end()) {
            _uniquepermuteneighbours.push_back(comp);
            _uniquepermuteneighboursstencil.push_back(stencilpos);
        }
        else
            throw runtime_error(
                    "Compartment::addenuniquepermuteNeighbour(): Compartment is already a "
                    "neighbour");
    }

    /// Remove a neighboring compartment
    void removeenclosingNeighbour(Compartment *comp) {
        for(int i = 0; i < _enclosingneighbours.size(); i++){
            if(comp == _enclosingneighbours[i]){
                _enclosingneighbours.erase(_enclosingneighbours.begin() + i);
                _enclosingneighboursstencil.erase(_enclosingneighboursstencil.begin() + i);
                break;
            }
        }
    }

    vector<Compartment*> getenclosingNeighbours(){
        return _enclosingneighbours;
    }

    void removeuniquepermuteNeighbour(Compartment *comp) {
        for(int i = 0; i < _uniquepermuteneighbours.size(); i++){
            if(comp == _uniquepermuteneighbours[i]){
                _uniquepermuteneighbours.erase(_uniquepermuteneighbours.begin() + i);
                _uniquepermuteneighboursstencil.erase(_uniquepermuteneighboursstencil.begin() + i);
                break;
            }
        }
    }
    vector<short> getuniquepermuteneighborsstencil(){
        return _uniquepermuteneighboursstencil;
    }
    vector<Compartment*> getuniquepermuteNeighbours(){
        return _uniquepermuteneighbours;
    }

    /// Clone the species values of another compartment into this one
    void cloneSpecies (Compartment *target) const {
        assert(target->numberOfSpecies()==0);
        for(auto &s : _species.species()){
            Species* sClone = s->clone();
            target->addSpeciesUnique(unique_ptr<Species>(sClone));
        }
    }

    /// Clone the reaction values of another compartment into this one
    void cloneReactions (Compartment *target) const {
        assert(target->numberOfReactions()==0);
        for(auto &r : _internal_reactions.reactions()){

            auto rClone = r->clone(target->_species);
            rClone->setVolumeFrac(target->getVolumeFrac());
            target->addInternalReaction(rClone);
        }
    }

    /// Clone both the species and compartments into this compartment
    void cloneSpeciesReactions(Compartment* C) {
        if(_activated) this->cloneSpecies(C);
        this->cloneReactions(C);
        C->diffusionCoefficients_ = this->diffusionCoefficients_;

    }

    /// Clone a compartment
    /// @note - this does not clone the neighbors, just reactions and species
    virtual Compartment* clone() {
        Compartment *C = new Compartment(*this);
        return C;
    }




    /// Gives the number of neighbors to this compartment
    int numberOfNeighbours() const {
        return std::count_if(neighborIndices_.begin(), neighborIndices_.end(), [](int ni) { return ni != -1; });
    }

    /// Gives the number of enclosing neighbors to this compartment
    size_t numberOfenclosingNeighbours() const {return _enclosingneighbours.size();}


    ///Get the species container vector
    SpeciesPtrContainerVector& getSpeciesContainer() {return _species;}
    const SpeciesPtrContainerVector& getSpeciesContainer() const {return _species;}

    ///Get the internal reaction container vector
    ReactionPtrContainerVector& getInternalReactionContainer() {return _internal_reactions;}
    const ReactionPtrContainerVector& getInternalReactionContainer() const {return _internal_reactions;}

    ///Get the diffusion reaction container vector
    ReactionPtrContainerVector& getDiffusionReactionContainer() {return _diffusion_reactions;}
    const ReactionPtrContainerVector& getDiffusionReactionContainer() const {return _diffusion_reactions;}

    /// Get the vector list of neighbors to this compartment
    auto&       getNeighborIndices()       { return neighborIndices_; }
    const auto& getNeighborIndices() const { return neighborIndices_; }

    /// Print the species in this compartment
    void printSpecies()const {_species.printSpecies();}
    /// Print the reactions in this compartment
    void printReactions()const {
        _internal_reactions.printReactions();
        _diffusion_reactions.printReactions();
    }

    /// Check if all species are unique pointers
    bool areAllSpeciesUnique () {return _species.areAllSpeciesUnique();}

    /// Adds the reactions of this compartment to the ChemSim object
    /// @param - chem, a ChemSim object that runs the reaction-diffusion algorithm
    virtual void addChemSimReactions(medyan::ChemSim* chem) {
        for(auto &r : _internal_reactions.reactions()) chem->addReaction(r.get());
        for(auto &r : _diffusion_reactions.reactions()) chem->addReaction(r.get());
    }

    /// Print properties of this compartment
    virtual void printSelf()const override {
        cout << this->getFullName() << "\n"
        << "Number of neighbors: " << numberOfNeighbours() << endl;
        printSpecies();
        cout << "Reactions:" << endl;
        printReactions();
    }

    //GetType implementation just returns zero (no Compartment types yet)
    virtual int getType() override {return 0;}

    // Helper function for getting the result of geometry from a approximately planar slice
    void computeSlicedVolumeArea(SubSystem& sys, SliceMethod);
    // Helper function that does not scale rates
    void computeNonSlicedVolumeArea();

    // Properties (public variables and getters and setters for private variables)
    bool boundaryInteresting = false; // A marker indicating this compartment is near a certain boundary

    void resetVolumeFrac() { _volumeFrac = 1.0; }
    floatingpoint getVolumeFrac() const { return _volumeFrac; }
    const array<floatingpoint, 6>& getPartialArea()const { return _partialArea; }
    void setPartialArea(const array<floatingpoint, 6>& partialArea) { _partialArea = partialArea; }

};
#ifdef SIMDBINDINGSEARCH
template<>
void Compartment::getpartitionindex<true>(int (&indices)[3], vector<floatingpoint> coord,
                             floatingpoint (&cmpcornercoords)[6]);
template<>
void Compartment::getpartitionindex<false>(int (&indices)[3], vector<floatingpoint> coord,
                              floatingpoint (&cmpcornercoords)[6]);
template<>
void Compartment::addcoordtorMaxbasedpartitons<true>(int (&pindices)[3], vector<floatingpoint>
        coord, uint32_t index, uint32_t cylfinfo);
template<>
void Compartment::addcoordtorMaxbasedpartitons<false>(int (&pindices)[3], vector<floatingpoint>
        coord, uint32_t index, uint32_t cylfinfo);
#endif

} // namespace medyan

#endif
