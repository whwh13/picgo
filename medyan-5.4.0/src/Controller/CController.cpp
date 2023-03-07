
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

#include "Controller/CController.h"

#include "SubSystem.h"

#include "ChemManager.h"

#include "ChemNRMImpl.h"
#include "ChemGillespieImpl.h"
#include "ChemSimpleGillespieImpl.h"

#include "CCylinder.h"
#include "Cylinder.h"

namespace medyan {

void CController::initialize(string& chemAlgorithm, ChemistryData& chem, DissipationTracker* dt, medyan::SimulConfig& sc) {
    
    // Set instance of chemsim algorithm
    if(chemAlgorithm == "NRM") {
        
#if !defined(TRACK_DEPENDENTS)
        cout << "The NRM algorithm relies on tracking dependents. Please set this"
            << " compilation macro and try again. Exiting." << endl;
        exit(EXIT_FAILURE);
#endif
        _subSystem->pChemSim = std::make_unique<ChemNRMImpl>();
    }
    
    else if(chemAlgorithm == "GILLESPIE") {
        
#if !defined(TRACK_DEPENDENTS)
        cout << "The Gillespie algorithm relies on tracking dependents. Please set this"
            << " compilation macro and try again. Exiting." << endl;
        exit(EXIT_FAILURE);
#endif
        _subSystem->pChemSim = std::make_unique<ChemGillespieImpl>();
    }
    
    else if(chemAlgorithm == "SIMPLEGILLESPIE") {
        _subSystem->pChemSim = std::make_unique<ChemSimpleGillespieImpl>();
    }
    else {
        log::error("Chem algorithm not recognized. Exiting.");
        throw std::runtime_error("Chem algorithm not recognized");
    }
    
    //Create manager, intialize
    _chemManager = new ChemManager;
    _subSystem->pChemSim->dt = dt;
    _chemManager->initializeSystem(_subSystem->pChemSim.get(), *_subSystem, sc);
    
    
    // init chemsim
    _subSystem->pChemSim->initialize();
//    extern bool afterchemsiminit;
//    afterchemsiminit = true;
    
    // set some static ptrs
    FilamentReactionTemplate::_ps = _subSystem;
    
    CCylinder::_chemSim = _subSystem->pChemSim.get();
    Cylinder::_chemManager = _chemManager;
}

void CController::initializerestart(floatingpoint time, floatingpoint _minimizationTime, ChemistryData& chemData){
    if(SysParams::RUNSTATE){
        LOG(ERROR) << "initializerestart Function from CController class can "
                      "only be called "
                      "during restart phase. Exiting.";
        throw std::logic_error("Illegal function call pattern");
    }
	_subSystem->pChemSim->initialize();
    _subSystem->pChemSim->initializerestart(time);
    //update copy numbers
    _chemManager->restartreleaseandremovaltime(_minimizationTime, chemData);
}

bool CController::run(FP time, ChemistryData& chemData, const SimulConfig& conf) {
    
    beforeRun();

    //update copy numbers
    _chemManager->updateCopyNumbers(*_subSystem, chemData, conf);
    
    //run the steps
    const auto res = _subSystem->pChemSim->run(time);

    afterRun();

    return res;
}

bool CController::runSteps(int steps, ChemistryData& chemData, const SimulConfig& conf) {
    
    beforeRun();

    //update copy numbers
    _chemManager->updateCopyNumbers(*_subSystem, chemData, conf);
    
    //run the steps
    const auto res = _subSystem->pChemSim->runSteps(steps);

    afterRun();

    return res;
}

bool CController::crosschecktau(){
    return _subSystem->pChemSim->crosschecktau();
}

ofstream CController::_crosscheckdumpFilechem;

} // namespace medyan


