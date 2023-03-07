
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

#include "Controller/DRController.h"

#include "RateChangerImpl.h"

#include "Linker.h"
#include "MotorGhost.h"
#include "Cylinder.h"
#include "BranchingPoint.h"

namespace medyan {
void DRController::initialize(const medyan::SimulConfig& sc) {
    auto& dr = sc.dyRateParams;
    auto& drTypes = dr.dynamicRateType;
    auto& chemParams = sc.chemParams;
    auto& chemData = sc.chemistryData;
    
    //filament polymerization changer
    int filamentIndex = 0;
    
    if(chemParams.numFilaments != 0) {
    
        for(auto &changer : drTypes.dFPolymerizationType) {
        
            if(changer == "BROWRATCHET") {
                //get params
                floatingpoint a = dr.dFilPolymerizationCharLength[filamentIndex];
                Cylinder::_polyChanger.push_back(new BrownianRatchet(a));
            }
            else if(changer == "") {}
            else {
                cout << "Filament polymerization rate changing form not recognized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            filamentIndex++;
        }
    }

    //branching point unbinding changer
    int branchIndex = 0;
    int charLengthIndexbr = 0;
    int charFIndexbr = 0;

    int numBrancherSpecies = 0;
    for(auto& brancherEachFilament : chemData.speciesBrancher) {
        numBrancherSpecies += brancherEachFilament.size();
    }
    if(numBrancherSpecies !=0) {
        
        for(auto &changer : drTypes.dBUnbindingType) {
            
            if(changer == "SLIP") {
                
                //if user did not specify enough parameters, return
                if(charLengthIndexbr >= dr.dBranchUnbindingCharLength.size() )
                    return;
                
                //get the param
                floatingpoint x1 = dr.dBranchUnbindingCharLength[charLengthIndexbr];
                
                //add the rate changer
                BranchingPoint::_unbindingChangers.push_back(new BranchSlip(branchIndex, x1));
                charLengthIndexbr += 1;
            }
            else if(changer == "SLIPF"){
                //if user did not specify enough parameters, return
                if(charFIndexbr >= dr.dBranchUnbindingCharForce
                                                .size() )
                    return;

                //get the param
                floatingpoint x1 = dr
                        .dBranchUnbindingCharForce[charLengthIndexbr];

                //add the rate changer
                BranchingPoint::_unbindingChangers.push_back(new BranchSlipF(branchIndex, x1));
                charFIndexbr += 1;

            }
            else {
                cout << "Branching point unbinding rate changing form not recognized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            
            branchIndex++;
        }
        
    }

    //linker unbinding changer
    int charLengthIndex = 0;
    int ampIndex = 0;
    int linkerIndex = 0;
    
    
    if(sc.chemistryData.numLinkerSpecies() != 0) {
    
        for(auto &changer : drTypes.dLUnbindingType) {
            
            if(changer == "CATCHSLIP") {
                
                //if user did not specify enough parameters, return
                if(ampIndex + 1 >= dr.dLinkerUnbindingAmplitude.size() ||
                   charLengthIndex + 1 >= dr.dLinkerUnbindingCharLength.size())
                    return;
                
                //get two params for each, amp
                floatingpoint a1 = dr.dLinkerUnbindingAmplitude[ampIndex];
                floatingpoint a2 = dr.dLinkerUnbindingAmplitude[ampIndex + 1];
                
                //now char length
                floatingpoint x1 = dr.dLinkerUnbindingCharLength[charLengthIndex];
                floatingpoint x2 = dr.dLinkerUnbindingCharLength[charLengthIndex + 1];
                
                //add the rate changer
                Linker::_unbindingChangers.push_back(new LinkerCatchSlip(linkerIndex, a1, a2, x1, x2));
                
                charLengthIndex += 2;
                ampIndex += 2;
            }
            
            else if(changer == "SLIP") {
                
                //if user did not specify enough parameters, return
                if(charLengthIndex >= dr.dLinkerUnbindingCharLength.size() )
                    return;
                
                //get the param
                floatingpoint x1 = dr.dLinkerUnbindingCharLength[charLengthIndex];
                
                //add the rate changer
                Linker::_unbindingChangers.push_back(new LinkerSlip(linkerIndex, x1));
                charLengthIndex += 1;
            }
            else {
                cout << "Linker unbinding rate changing form not recognized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            
            linkerIndex++;
        }
        
    }
    int forceIndex = 0;
    int motorIndex = 0;
    
    if(sc.chemistryData.numMotorSpecies() != 0) {
    
        //motor unbinding changer
        for(auto &changer : drTypes.dMUnbindingType) {
            
            if(changer == "LOWDUTYCATCH") {
                
                //if user did not specify enough parameters, return
                if(forceIndex >= dr.dMotorUnbindingCharForce.size())
                    return;
                
                //get param
                floatingpoint f = dr.dMotorUnbindingCharForce[forceIndex];
                
                //add the rate changer
                MotorGhost::_unbindingChangers.push_back(new LowDutyMotorCatch(motorIndex, f));
                forceIndex++;
            }
            else if(changer == "HIGHDUTYCATCH") {
                
                //if user did not specify enough parameters, return
                if(forceIndex >= dr.dMotorUnbindingCharForce.size())
                    return;
                
                //get param
                floatingpoint f = dr.dMotorUnbindingCharForce[forceIndex];
                
                //add the rate changer
                MotorGhost::_unbindingChangers.push_back(new HighDutyMotorCatch(motorIndex, f));
                forceIndex++;
            }

            else if(changer == "LOWDUTYCATCHSLIP") {
                cout << "Catch-slip bond implementation of low duty motor not complete. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else {
                cout << "Motor unbinding rate changing form not recognized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            motorIndex++;
        }
        forceIndex = 0;
        motorIndex = 0;
        
        //motor walking 
        for(auto &changer : drTypes.dMWalkingType) {

            if(changer == "LOWDUTYSTALL") {
                
                //if user did not specify enough parameters, return
                if(forceIndex >= dr.dMotorWalkingCharForce.size())
                    return;
                
                //get the param
                floatingpoint f = dr.dMotorWalkingCharForce[forceIndex];
                
                //add the rate changer
                MotorGhost::_walkingChangers.push_back(new LowDutyMotorStall(motorIndex, 0, f));
                forceIndex++;
            }
            else if(changer == "HIGHDUTYSTALL") {
                
                //if user did not specify enough parameters, return
                if(forceIndex >= dr.dMotorWalkingCharForce.size())
                    return;
                
                //get the param
                floatingpoint f = dr.dMotorWalkingCharForce[forceIndex];
                
                //add the rate changer
                MotorGhost::_walkingChangers.push_back(new HighDutyMotorStall(motorIndex, 0, f));
                forceIndex++;
            }
            
            else {
                cout << "Motor walking rate changing form not recognized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            motorIndex++;
        }
    }
}

} // namespace medyan
