
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

#ifndef MEDYAN_ChemManager_h
#define MEDYAN_ChemManager_h

#include "common.h"

#include "ReactionTemplate.h"
#include "DissipationTracker.h"

namespace medyan {
//FORWARD DECLARATIONS
class Compartment;
class CompartmentGrid;
class ChemSim;
class Cylinder;
class CCylinder;
class CMonomer;


/// For initailizing chemical reactions based on a specific system
/*!
 *  ChemManager is used for initailizing all chemistry in the system.
 *  Initialized by the CController, ChemManager initializes the chemical components
 *  of the CompartmentGrid. The ChemManager can also update chemical components of
 *  a CCylinder using its stored collection of FilamentReactionTemplate, which is
 *  initialized at startup.
 */
class ChemManager {
    
private:
    chrono::high_resolution_clock::time_point mins, mine;
    //@{
    /// Helper functions. Names are pretty self-explanatory.
    static void setupBindingSites(ChemParams& chemParams, const medyan::SimulConfig& sc);
    
    static void configCMonomer(const medyan::SimulConfig& sc);
    static void initCMonomer(CMonomer* m, short filamentType, Compartment* c, const ChemistryData& chemData);
    
    static void genSpecies(CompartmentGrid& grid, Compartment& protoCompartment, const medyan::SimulConfig& sc);
    
    static void genGeneralReactions(CompartmentGrid& grid, Compartment& protoCompartment, const ChemParams& chemParams, const ChemistryData& chemData);
    static void genBulkReactions(CompartmentGrid& grid, const ChemistryData& chemData);
    
    static void genNucleationReactions(SubSystem& sys, CompartmentGrid& grid, const medyan::SimulConfig& sc);
    static void genFilBindingReactions(SubSystem& sys, const medyan::SimulConfig& sc, DissipationTracker* pdt);
    
    void genFilReactionTemplates(ChemParams& chemParams, const ChemistryData& chemData, DissipationTracker* pdt);
    //@}
    
public:
	static floatingpoint tchemmanager1;
	static floatingpoint tchemmanager2;
	static floatingpoint tchemmanager3;
	static floatingpoint tchemmanager4;

    // An auxiliary function that finds certain species in an array indexed by filament type.
    //
    // Parameters:
    // - species: outer index is filament type.
    // - name:    the name of the species to search for.
    // - catName: the category name for species, used in error message.
    //
    // Returns tuple containing indices to find the species.
    static auto locateSpecies(const vector<vector<string>>& species, const string& name, const string& catName) {
        for(int ft = 0; ft < species.size(); ++ft) {
            auto& speciesVec = species[ft];
            for(int i = 0; i < speciesVec.size(); ++i) {
                if(speciesVec[i] == name) {
                    return tuple{ ft, i };
                }
            }
        }

        LOG(ERROR) << "Cannot find species " << name << " in category " << catName;
        throw runtime_error("Invalid species name.");
    }

    /// Initialize the system, based on the given simulation
    /// Adds all necessary reactions to the ChemSim object
    void initializeSystem(medyan::ChemSim* chemSim, SubSystem& sys, medyan::SimulConfig& sc);
    
    ///Initializer for chem cylinders, based on the given simulation
//    virtual void initializeCCylinder(CCylinder* cc,
//                                     bool extensionFront,
//                                     bool extensionBack,
//                                     bool initialization);

	void initializeCCylinder(CCylinder* cc,
	                                      bool extensionFront,
	                                      bool extensionBack,
	                                      bool initialization,
	                                      int nummonomers = -1,
	                                      int firstmonomer = -1,
	                                      int lastmonomer = -1,
	                                      bool minusendstatus = false,
	                                      bool plusendstatus = false,
	                                      short minusendtype = -1,
	                                      short plusendtype = -1);
    
    /// Update the copy numbers of all species in the chemical network
    /// @note - this only sets the copy number if the simulation time
    ///         tau has passed the release time of the molecule. This
    ///         function is called at every set of chemical steps to check
    ///         if molecules should be released at the current time.
    void updateCopyNumbers(SubSystem& sys, ChemistryData& chemData, const medyan::SimulConfig& sc);

    void restartreleaseandremovaltime(floatingpoint _minimizationTime, ChemistryData& chemData);
    
private:
    ChemistryData chemDataBackup_;

    /// A list of reactions to add to every new CCylinder
    /// @note - is a 2D vector corresponding to different filament types
    vector<vector<unique_ptr<FilamentReactionTemplate>>> _filRxnTemplates =
    vector<vector<unique_ptr<FilamentReactionTemplate>>>(MAX_FILAMENT_TYPES);
};

} // namespace medyan

#endif
