
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

#ifndef MEDYAN_Restart_h
#define MEDYAN_Restart_h

#include "common.h"
#include <random>
#include <chrono>

#include "SubSystem.h"
#include "Boundary.h"
#include "CompartmentGrid.h"

#include "BubbleInitializer.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Linker.h"
#include "MotorGhost.h"
#include "BranchingPoint.h"
#include "Bubble.h"

#include "SysParams.h"
#include "MathFunctions.h"
#include "Cylinder.h"
#include "RestartParams.h"
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include "ChemRNode.h"

namespace medyan {
using namespace mathfunc;

// Forward declarations.
class SubSystem;
class Cylinder;
class FilamentBindingManager;


class Restart {
private:
    SubSystem *_subSystem; ///< A pointer to the subsystem that this controls
    vector<floatingpoint> temp_diffrate_vector; ///vector of diffusion rates
    ChemistryData _chemData;
    vector<floatingpoint> CopyNumbers;
    int  _numChemSteps=0;

	fstream _inputFile; ///< input file being used
	restartSystemData _rsystemdata;
	restartBeadData _rBData;
	vector<restartFilData> _rFDatavec;
    vector<restartCylData> _rCDatavec;
    vector<restartMotorData> _rMDatavec;
    vector<restartLinkerData> _rLDatavec;
    vector<restartBrancherData> _rBDatavec;
    vector<restartCompartmentDiffusingData> _rCmpDDatavec;
	restartBulkData _rbdata;
    restartcopynumberrallyData _rtallydata;

    //Increases copy number of diffusing species corresponding to bound species in each Compartment by number of events.
    void setdiffspeciesnumber(string diffusingspeciesname, Cylinder* c){
        int counter=0;
//        cout<<"Trying to find "<<diffusingspeciesname<<" in Cmp "<<c->getCompartment()
//        ->getId()<<endl;
        for(auto& sd : _chemData.speciesDiffusing) {
            int events=0;
            if(diffusingspeciesname.compare(sd.name)==0){
            	events++;
            }
            CopyNumbers[counter]=CopyNumbers[counter]-events;
            if(CopyNumbers[counter]<0 && SysParams::filamentSetup.USECHEMCOPYNUM == false)
            {cout <<
                "Restart file reaction numbers do not match with diffusing species number."
                << endl;
                exit(EXIT_FAILURE);

            }
/*            if(events > 0) {
                cout << "Cmp " << c->getCompartment()->getId() << " diffusing species "
                     << diffusingspeciesname << "copy number set to " << events +
                        (c->getCompartment()->findSpeciesByName(get<0>(sd)))->getRSpecies().getN()
                     <<" from "<< (c->getCompartment()->findSpeciesByName(get<0>(sd)))
                     ->getRSpecies().getN()<<endl;
            }*/
            events=events+(c->getCompartment()->findSpeciesByName(sd.name))->getRSpecies().getN();

            (c->getCompartment()->findSpeciesByName(sd.name))->getRSpecies().setN(events);
            c->getCompartment()->getDiffusionReactionContainer().updatePropensityComprtment();
            counter++;
        }
    }
public:
    Restart(SubSystem* s, ChemistryData _cd, const string inputFileName)
            : _subSystem(s), _chemData(_cd) {
	    _inputFile.open(inputFileName);
	    if(!_inputFile.is_open()) {
		    cout << "There was an error parsing file " << inputFileName
		         << ". Exiting." << endl;
		    exit(EXIT_FAILURE);
	    }

    }

    ~Restart(){_inputFile.close();}

    void readNetworkSetup();

    void setupInitialNetwork();

    floatingpoint getrestartime(){return _rsystemdata.time;}

    void settorestartphase(){
//STEP #1: Get a copy of diffusion rate and set diffusion rate to 0, reset linker motor and branching managers.
        for(auto& C : _subSystem->getCompartmentGrid()->getCompartments()) {
            for(auto &it: C->getDiffusionReactionContainer().reactions())
            {temp_diffrate_vector.push_back(it->getRate());
                it->setRateMulFactor(0.0f, ReactionBase::RESTARTPHASESWITCH);
            }
            C->getDiffusionReactionContainer().updatePropensityComprtment();
            for(auto &Mgr:C->getFilamentBindingManagers()){
#ifdef NLORIGINAL
                Mgr->clearpossibleBindings();
#else
                Mgr->clearpossibleBindingsstencil();
#endif
            }}
//STEP #1a: Get cylinders, passivate filament reactions.
        for(auto& C : _subSystem->getCompartmentGrid()->getCompartments()) {
            for(auto x : C->getCylinders()) {
                x->getCCylinder()->passivatefilreactions();
                x->getCCylinder()->passivatefilcrossreactions();
            }}
//Step #1b. Passivate general reactions.
	    for(auto& C : _subSystem->getCompartmentGrid()->getCompartments()) {
		    for(auto& rxn : C->getInternalReactionContainer().reactions()) {
		    	if(rxn->getReactionType() == ReactionType::REGULAR)
		    		rxn->passivateReaction();
		    }}
//Step #1c. Get copynumber of diffusing species. This is used later for book keeping
// purposes.
        for(auto& sd : _chemData.speciesDiffusing) {
            string name = sd.name;
	        CopyNumbers.push_back(sd.initialCopyNumber);
            //CopyNumbers.push_back(_subSystem->getCompartmentGrid()
            //->countDiffusingSpecies(name));
        }
        //Set copy number of diffusing species in each compartment to 0.
        for(auto& C : _subSystem->getCompartmentGrid()->getCompartments()) {
            for(auto sd : _chemData.speciesDiffusing) {
                (C->findSpeciesByName(sd.name))->getRSpecies().setN(0);}
            for(auto sd : _chemData.speciesBulk) {
                (C->findSpeciesByName(sd.name))->getRSpecies().setN(0);}
        }
/*//Step #3. Add filament coordinates to be held static during minimization **** NEEEDS TO BE EDITED***
        // coordinates to keep static
        vector<vector<floatingpoint>> staticbeads=get<3>(filaments);
        int nbeads=get<3>(filaments).size();
        for (int ix=0;ix<nbeads;ix++){
            vector<floatingpoint> staticbead=staticbeads[ix];
            vector<floatingpoint> staticbead2;
            staticbead2.insert ( staticbead2.begin(), staticbead.begin()+1, staticbead.begin()+3 );
            for(auto b: Bead::getBeads()) {
                floatingpoint dis=twoPointDistance(b->vcoordinate(),staticbead2);
                if(dis<=0.00001){
                    b->setstaticstate(true);
                }}}*/

//STEP #2 . updating _possbileBindings of Linkers in each compartment.
        //Filter through probable sites in unsortedpairings by making sure the distance between binding sites
        //in them is the same as bound species bond length
        addtoHeapLinkerMotorBrancher();
    }

    int getnumchemsteps(){return _numChemSteps;}

    void addtoHeapLinkerMotorBrancher();

    void CBoundinitializerestart();

    void restartupdateCopyNumbers(){
        // If user requests to use restart file copy numbers after restart.
        if(SysParams::RUNSTATE == false && SysParams::filamentSetup.USECHEMCOPYNUM == false) {
            //Loop through restart data
            // Find corresponding compartment and add species to it.
            for (auto &Cmp:_subSystem->getCompartmentGrid()->getCompartments()) {
                auto cdiffdata = _rCmpDDatavec[Cmp->getId()];
                auto diffspeciesnamevec = cdiffdata.speciesnamevec;
                auto diffspeciescpyvec = cdiffdata.copynumvec;
                if (Cmp->getId() == cdiffdata.id) {
                    for (auto count = 0; count < diffspeciesnamevec.size(); count++) {
                        Species *sp = Cmp->findSpeciesByName(diffspeciesnamevec[count]);
                        if (sp != nullptr) {
                            auto currcopynumber = sp->getRSpecies().getN();
                            auto temptargetcopynumber = diffspeciescpyvec[count];
                            //setN
                            while(currcopynumber > temptargetcopynumber) {
                                sp->down();
                                currcopynumber = sp->getN();
                            }
                            while(currcopynumber < temptargetcopynumber) {
                                sp->up();
                                currcopynumber = sp->getN();
                            }
                            sp->getRSpecies().setN(diffspeciescpyvec[count]);
                        } else {
                            LOG(ERROR) << "Cannot find diffusing species of name "
                                        << diffspeciesnamevec[count]
                                        << " in Compartment ID " << Cmp->getId()
                                        << ". Check chemistry input file. Exiting" << endl;
                            throw std::logic_error("Restart unsuccessful!");
                        }
                    }
                }
            }
            int counter = 0;
            for (auto& C : _subSystem->getCompartmentGrid()->getCompartments()) {
                for (auto &it: C->getDiffusionReactionContainer().reactions()) {
                    it->setRateMulFactor(1.0f, ReactionBase::RESTARTPHASESWITCH);
                    it->updatePropensity();
                }
            }
            //Bulk Species
            for (auto &s : _chemData.speciesBulk) {
                log::error("Bulk species protocols not implemented.");
            }
            //Set diffusing species copy numbers to 0.
            for(auto &s : _chemData.speciesDiffusing) {
                auto name = s.name;
                auto copyNumber = s.initialCopyNumber;
                auto releaseTime = s.releaseTime;
                if (tau() >= releaseTime)
                    s.releaseTime = 0;
            }
    /*            for(auto &s : _chemData.speciesDiffusing) {

                auto name = get<0>(s);
                auto copyNumber = get<1>(s);
                auto releaseTime = get<3>(s);
                cout<<name<<" "<<copyNumber<<endl;
            }*/

        }
        checktallyofcopynumbers();
    }

    void checktallyofcopynumbers(){
        for(int i =0; i < _rtallydata.copynumvec.size(); i++){
            string name = _rtallydata.speciesnamevec.at(i);
            auto copyNum = _subSystem->getCompartmentGrid()->countDiffusingSpecies(name);
            if(_rtallydata.copynumvec[i] != copyNum){
                LOG(ERROR) << "Species copy number does not tally. Species name "
                           << name << " requires total copy number " << _rtallydata.copynumvec[i]
                           << "during restart but current state only has "<<copyNum
                           <<" molecules.. Check chemistry input file. Exiting"
                           << endl;
                throw std::logic_error("Restart unsuccessful!");
            }
        }

    }

    bool crosscheck();
};

} // namespace medyan

#endif
