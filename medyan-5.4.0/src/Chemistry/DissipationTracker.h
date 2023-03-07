
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

#ifndef MEDYAN_DissipationTracker_h
#define MEDYAN_DissipationTracker_h

#include <fstream>

#include "common.h"

#include "MathFunctions.h"
#include "ReactionBase.h"
#include "Filament.h"
#include "common.h"
#include "SysParams.h"
#include "MotorGhost.h"
#include "Linker.h"
#include "Mechanics/ForceField/Types.hpp"

namespace medyan {
using namespace mathfunc;

/* Dissipation Tracker is used to track the energetics of the system.  It has functions to compute the change in Gibbs free energy of a chemical reaction
 based on the instantaneous compartment copy numbers, using expressions derived in "A Discrete Approximation to Gibbs Free Energy of Chemical Reactions is Needed 
 for Accurately Calculating Entropy Production in Mesoscopic Simulations" by Carlos Floyd, Garegin A. Papoian, and Christopher Jarzynski.  It also tracks the changes in 
 mechanical energy of the system through calls to the ForceFieldManager, which iterates through the force fields and computes the instantaneous mechanical energy of each.  
 It further contains functions to track the spatiotemporal information of motor walking, linker binding, and linker unbinding events.
 */

class DissipationTracker{
    
private:
    
    // reaction counter - not used
    int count;
    
    // cumulative dissipated Gibbs free energy
    floatingpoint cumDissEnergy;
    
    // Mechanical energy before checmial simulation
    floatingpoint G1;
    
    // Mechanical energy after chemical simulation
    floatingpoint GMid;
    
    // Mechanical energy after subsequent mechanical equilibration
    floatingpoint G2;
    
    // Change in chemical free energy after chemical simulation
    floatingpoint GChem;
    
    // cumulative dissiapted chemical energy
    floatingpoint cumDissChemEnergy;
    
    // cumulative dissiapted mechanical energy
    floatingpoint cumDissMechEnergy;

    // cumulative change in chemical energy
    floatingpoint cumGChemEn;
    
    // cumulative change in mechanical energy
    floatingpoint cumGMechEn;
    
    // vector of HRCD elements
    vector<tuple<string, floatingpoint>> HRCDVec;
    
    // vector of HRMD element
    vector<tuple<string, floatingpoint>> HRMDVec1;
    
    // vector of HRMD element
    vector<tuple<string, floatingpoint>> HRMDVecMid;
    
    
    // vector of HRMD element
    vector<tuple<string, floatingpoint>> cumHRMDMechDiss;
    
    // vector of HRMD element
    vector<tuple<string, floatingpoint>> cumHRMDMechEnergy;
    
    // vector of HRMD element
    vector<tuple<string, floatingpoint>> HRMDVec2;
    vector<vector<tuple<string, floatingpoint>>> HRMDMat;
    
    // vector of motor walking data
    //ID, birthtime, walktime, xcoord, ycoord, zcoord,
    vector<tuple<int, floatingpoint, floatingpoint, floatingpoint, floatingpoint,
    floatingpoint>> motorWalkData;

    // vector of motor unbinding data
    //ID, birthtime, deathtime, xcoord, ycoord, zcoord,
    vector<tuple<int, floatingpoint, floatingpoint, floatingpoint, floatingpoint,
            floatingpoint>> motorUnbindingData;
    
    // vector of linker unbinding data
    //ID, birthtime, deathtime, xcoord, ycoord, zcoord,
    vector<tuple<int, floatingpoint, floatingpoint, floatingpoint, floatingpoint,
    floatingpoint>> linkerUnbindingData;
    
    // vector of linker binding data
    //ID, birthtime, xcoord, ycoord, zcoord,
    vector<tuple<int, floatingpoint, floatingpoint, floatingpoint, floatingpoint>>
    linkerBindingData;
    
public:
    
    // Constructor, allow access to objects that information is needed from, set all energy
    // tracking variables to zero
    DissipationTracker(const medyan::SimulConfig& conf):
        // get the stepFrac parameter for myosin walking, assume all fil and motor are type 0
        d_step(conf.chemParams.motorStepSize.empty() ? 0 : conf.chemParams.motorStepSize[0]),
        d_total(conf.geoParams.cylinderSize[0] / conf.chemParams.numBindingSites[0]),
        _stepFrac(d_step / d_total)
    {
        count = 0;
        cumDissEnergy=0;
        cumDissChemEnergy=0;
        cumDissMechEnergy=0;
        cumGChemEn=0;
        cumGMechEn=0;
        G1=0;
        G2=0;
        GMid=0;
        GChem=0;
    };
    
    // increment the reaction counter
    void updateCount(){
        count++;
    };
    
    // return reaction counter
    int getCount(){return count;};
    
    // print reaction counter
    void printCount(ofstream& outputFile){
        outputFile << count << endl;};
    
    // get the number of compartments in each dimension to form the ratio v used later
    int nx = SysParams::Geometry().NX;
    int ny = SysParams::Geometry().NY;
    int nz = SysParams::Geometry().NZ;
    float v = (float)(nx*ny*nz);
    
    // get the stepFrac parameter for myosin walking, assume all fil and motor are type 0
    floatingpoint d_step = 0;
    floatingpoint d_total = 0;
    floatingpoint _stepFrac = 0;
    
    
    // Find the Gibbs free energy change for a reaction using its ReactionBase representation
    floatingpoint getDelGChem(ReactionBase* re){
        using std::log;

        // get the type of reaction
        ReactionType reType = re->getReactionType();
            
        // get the number of reactants
        int M = re->getM();
        
        // get the number of products
        int N = re->getN();
        int sigma = N - M;
        
        float volFrac = re->getVolumeFrac();
        
        // for a vector of stoichiometric coefficients, assumed to be 1 for all
        vector<int> reacNu(M,1);
        vector<int> prodNu(N,1);
        
        // get vector of copy numbers of reactants and products
        vector<species_copy_t> reacN = re->getReactantCopyNumbers();
        vector<species_copy_t> prodN = re->getProductCopyNumbers();
        
        // get vector of names of the reactants and products
        vector<string> reacNames = re->getReactantSpecies();
        vector<string> prodNames = re->getProductSpecies();
        
        // get the HRCDID of this reaction
        string hrcdid;
        hrcdid = re->getHRCDID();
        
        // add the name of the diffusing species to the HRCDID of the diffusion reaction
        if(reType == ReactionType::DIFFUSION){
            hrcdid = "DIF_";
            hrcdid += reacNames[0];
        }
        
        
        float delGZero;
        
        // declare delG and set it to 0
        float delG;
        delG=0;
        
        if(reType == ReactionType::REGULAR) {
            // Regular Reaction
            delGZero =  re->getGNumber();
            delGZero -= sigma*log(volFrac);
            
            delG = delGGenChem(delGZero, reacN, reacNu, prodN, prodNu);
            
        } else if(reType == ReactionType::DIFFUSION){
            // Diffusion Reaction
            
            delG = delGDifChem(reacN[0],prodN[0]);
            
            
        } else if(reType == ReactionType::POLYMERIZATIONPLUSEND){
            // Polymerization Plus End
            
            delGZero = re->getGNumber();
            delGZero += log(volFrac);
            
            species_copy_t nMon = reacN[0];
            delG = delGPolyChem(delGZero,nMon,"P");
            
            
        } else if(reType == ReactionType::POLYMERIZATIONMINUSEND){
            // Polymerization Minus End
            
            delGZero = re->getGNumber();
            delGZero += log(volFrac);
            
            species_copy_t nMon = reacN[0];
            delG = delGPolyChem(delGZero,nMon,"P");
            
            
        } else if(reType == ReactionType::DEPOLYMERIZATIONPLUSEND){
            // Depolymerization Plus End
            
            delGZero = re->getGNumber();
            delGZero -= log(volFrac);
            
            species_copy_t nMon = prodN[0];
            delG = delGPolyChem(delGZero,nMon,"D");
            
        } else if(reType == ReactionType::DEPOLYMERIZATIONMINUSEND){
            // Depolymerization Minus End
            
            delGZero = re->getGNumber();
            delGZero -= log(volFrac);
            
            species_copy_t nMon = prodN[0];
            delG = delGPolyChem(delGZero,nMon,"D");
            
            
        } else if(reType == ReactionType::LINKERBINDING){
            // Linker Binding
            delGZero=re->getGNumber();
            delGZero += log(volFrac);
            
            species_copy_t nMon = reacN[1];
            delG = delGPolyChem(delGZero,nMon,"P");
            

            
        } else if(reType == ReactionType::MOTORBINDING){
            // Motor Binding
            // need to implement volFrac change here
            floatingpoint rn=re->getGNumber();
            
            floatingpoint nh1 = SysParams::Chemistry().motorNumHeadsMin[0];
            floatingpoint nh2 = SysParams::Chemistry().motorNumHeadsMax[0];
            floatingpoint nh = (nh1+nh2)/2.0;
            
            delG = delGMyoChem(nh,rn);
            delG += log(volFrac);
            
            
        } else if(reType == ReactionType::LINKERUNBINDING){
            // Linker Unbinding
            delGZero=re->getGNumber();
            delGZero -= log(volFrac);
            
            species_copy_t nMon = prodN[0];
            delG = delGPolyChem(delGZero,nMon,"D");
            
            CBound* CBound = re->getCBound();
            SpeciesBound* sm1 = CBound->getFirstSpecies();
            Linker* l = ((CLinker*)sm1->getCBound())->getLinker();
            recordLinkerUnbinding(l);
            
            
        } else if(reType == ReactionType::MOTORUNBINDING){
            // Motor Unbinding
            floatingpoint rn=re->getGNumber();
            
            CBound* CBound = re->getCBound();
            SpeciesBound* sm1 = CBound->getFirstSpecies();
            MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();
            int nh = m->getNumHeads();
            
            delG = delGMyoChem(nh,rn);
            delG = -delG;
            delG -= log(volFrac);
            recordMotorUnbinding(m);
            
        } else if(reType == ReactionType::MOTORWALKINGFORWARD){
            // Motor Walking Forward
            delG = re->getGNumber();
            delG = delG*(1/_stepFrac);
            
            
            
        } else if(reType == ReactionType::MOTORWALKINGBACKWARD){
            // Motor Walking Backward
            delG = -(re->getGNumber());
            delG = 0;
            
            
        } else if(reType == ReactionType::AGING){
            // Filament Aging
            
            vector<species_copy_t> numR;
            numR.push_back(Filament::countSpecies(0,reacNames[0]));
            
            vector<species_copy_t> numP;
            numP.push_back(Filament::countSpecies(0,prodNames[0]));
            
            delGZero = (re->getGNumber());
            delGZero -= sigma*log(volFrac);
            
            delG = delGGenChem(delGZero, numR, reacNu, numP, prodNu);
            
        } else if(reType == ReactionType::FILAMENTCREATION){
            // Filament Creation, not currently supported
            
        } else if(reType == ReactionType::FILAMENTDESTRUCTION){
            // Filament Destruction, not currently supported
            
        } else if(reType == ReactionType::BRANCHING){
            // Branching, not currently supported
            
        } else if(reType == ReactionType::BRANCHUNBINDING){
            // Branch Unbinding, not currently supported
            
        } else if(reType == ReactionType::SEVERING){
            // Severing, not currently supported
            
        }
        
        // check for nan
        if(delG<-pow(10,7)){
            delG=0;
        }
        
        if(delG>pow(10,7)){
            delG=0;
        }
        
        if(delG!=delG){
            delG=0;
        }
        
        if(isnan(delG)){
            delG=0;
        }
        
        // record this reaction in the HRCD data
        updateHRCDVec(hrcdid,delG);
        
        if(hrcdid=="DNT"){
            cout<< medyan::underlying(reType) <<endl;
        }


        return delG;
        
        
    }
    
    
    
    // increment the GChem counter when a reaction fires
    void updateDelGChem(ReactionBase* re){
        GChem += getDelGChem(re);
        updateCount();
        
    }
    
    
    // return the cumulative dissipated Gibbs free energy
    floatingpoint getCumDissEnergy(){
        return cumDissEnergy;
    }
    
    // return the cumulative dissipated chemical Gibbs free energy
    floatingpoint getCumDissChemEnergy(){
        return cumDissChemEnergy;
    }
    
    // return the cumulative dissipated mechanical energy
    floatingpoint getCumDissMechEnergy(){
        return cumDissMechEnergy;
    }
    
    // return the cumulative change in chemical Gibbs free energy
    floatingpoint getCumGChemEn(){
        return cumGChemEn;
    }
    
    // return the cumulative change in mechanical energy
    floatingpoint getCumGMechEn(){
        return cumGMechEn;
    }
    
    // return the HRMD cumulative change in mechanical energy
    vector<tuple<string, floatingpoint>> getCumHRMDMechEnergy(){
        return cumHRMDMechEnergy;
    }
    
    // return the HRMD cumulative change in mechanical dissipation
    vector<tuple<string, floatingpoint>> getCumHRMDMechDiss(){
        return cumHRMDMechDiss;
    }

    //get HRMD mat
    vector<vector<tuple<string, floatingpoint>>> getHRMDmat(){
        return HRMDMat;
    }
    //  used to determine if minization should proceed
    floatingpoint getCurrentStress(){
        return GMid-G1;
    };


    //
    
    // set the value of G1 - before chemistry
    void setG1(const EnergyReport& report){
        HRMDVec1.clear();
        
        for(auto i = 0; i < report.individual.size(); i++){
            HRMDVec1.push_back(make_tuple(report.individual[i].name, report.individual[i].energy / kT));
            cumHRMDMechDiss.push_back(make_tuple(report.individual[i].name, 0.0));
            cumHRMDMechEnergy.push_back(make_tuple(report.individual[i].name, 0.0));
        };
        
        G1 = report.total / kT;
      
    }
    
    // set the value of G2 - after minimization
    void setG2(const EnergyReport& report){
        HRMDVec2.clear();
        for(auto i = 0; i < report.individual.size(); i++){
            HRMDVec2.push_back(make_tuple(report.individual[i].name, report.individual[i].energy / kT));
        };
        HRMDMat.push_back(HRMDVec2);
        G2 = report.total / kT;
    }
    
    // set the value of GMid
    void setGMid(const EnergyReport& report){
        HRMDVecMid.clear();
        for(auto i = 0; i < report.individual.size(); i++){
            HRMDVecMid.push_back(make_tuple(report.individual[i].name, report.individual[i].energy / kT));
        };
        HRMDMat.push_back(HRMDVecMid);
        GMid = report.total / kT;
    }
    
    // perform multiple functions to update cumulative energy counters and reset the mechanical energy variables
    void updateAfterMinimization(const EnergyReport& report){
        setG2(report);
        updateCumDissChemEnergy();
        updateCumDissMechEnergy();
        updateCumDissEn();
        updateCumGChemEn();
        updateCumGMechEn();
        updateCumHRMDMechEnergy();
        updateCumHRMDMechDiss();
        resetAfterStep();
    }
    
    // increment cumDissEn after an iteration step has occured
    void updateCumDissEn(){
        cumDissEnergy += GChem + G2 - G1;
    }
    
    // increment cumDissChemEnergy
    void updateCumDissChemEnergy(){
        cumDissChemEnergy += GChem + getCurrentStress();
    }
    
    
    
    // increment cumDissMechEnergy
    void updateCumGChemEn(){
        cumGChemEn += GChem;
    }
    
    // increment cumDissMechEnergy
    void updateCumGMechEn(){
        cumGMechEn += G2-G1;
    }
    
    // increment cumDissMechEnergy
    void updateCumDissMechEnergy(){
        cumDissMechEnergy += G2-GMid;
    }
    
    // increment cumHRMDMechEnergy
    void updateCumHRMDMechEnergy(){
        
        for(auto i = 0; i<cumHRMDMechEnergy.size();i++){
            get<1>(cumHRMDMechEnergy[i]) += get<1>(HRMDVec2[i]) - get<1>(HRMDVec1[i]);
        }
        
    }
    
    // increment cumHRMDMechDiss
    void updateCumHRMDMechDiss(){
        
        for(auto i = 0; i<cumHRMDMechEnergy.size();i++){
            get<1>(cumHRMDMechDiss[i]) += get<1>(HRMDVec2[i]) - get<1>(HRMDVecMid[i]);
        }

    }
    
    
    
    // set new values of energy trackers after an iteration step has occured
    void resetAfterStep(){
        GChem=0;
        GMid=0;
        G1=G2;
        
        for(auto i = 0; i < HRMDVec1.size(); i++){
            HRMDVecMid[i] = make_tuple(get<0>(HRMDVecMid[i]),0.0);
            HRMDVec1[i] = make_tuple(get<0>(HRMDVec1[i]),get<1>(HRMDVec2[i]));
        };
        
    };
    
    // add the changes in the reactions' chemical energy consumptions
    void updateHRCDVec(string hrcdid, floatingpoint delG){
        for(auto i = 0; i<HRCDVec.size(); i++){
            if(get<0>(HRCDVec[i])==hrcdid){
                get<1>(HRCDVec[i]) += delG;
                return;
            }
        };
        HRCDVec.push_back(make_tuple(hrcdid,delG));
    }
    
    vector<tuple<string,floatingpoint>> getHRCDVec(){
        return HRCDVec;
    }
    
    // store the space time information of a motor walking event to motorWalkData
    void recordWalk(MotorGhost* m){
        vector<floatingpoint> mcoords = m->coordinate;
        mcoords.insert(mcoords.begin(), tau());
        //ID birthtime, walktime, coordx coordy coordz
        motorWalkData.push_back(make_tuple(m->getId(), m->getBirthTime(), mcoords[0],
                mcoords[1], mcoords[2], mcoords[3]));

    
    }
    
    vector<tuple<int, floatingpoint, floatingpoint, floatingpoint, floatingpoint,
    floatingpoint>>
    getMotorWalkingData(){
        return motorWalkData;
    }
    
    void clearMotorData(){
        motorWalkData.clear();
    }

    void clearMotorUnbindingData(){
        motorUnbindingData.clear();
    }
    
    // store the space time information of a linker unbinding event to linkerUnbindingData
    void recordLinkerUnbinding(Linker* l){
        vector<floatingpoint> lcoords = l->coordinate;
        lcoords.insert(lcoords.begin(), tau());
        linkerUnbindingData.push_back(make_tuple(l->getId(), l->getBirthTime(), lcoords[0],
                lcoords[1], lcoords[2],lcoords[3]));
    }

    // store the space time information of a linker unbinding event to motorUnbindingData
    void recordMotorUnbinding(MotorGhost* m){
        vector<floatingpoint> mcoords = m->coordinate;
        mcoords.insert(mcoords.begin(), tau());
        //ID birthtime unbindtime coordx coordy coordz
        motorUnbindingData.push_back(make_tuple(m->getId(), m->getBirthTime(), mcoords[0],
                                                 mcoords[1], mcoords[2], mcoords[3]));
    }
    
    vector<tuple<int, floatingpoint, floatingpoint, floatingpoint, floatingpoint,
    floatingpoint>>
    getLinkerUnbindingData(){
        return linkerUnbindingData;
    }

    vector<tuple<int, floatingpoint, floatingpoint, floatingpoint, floatingpoint,
            floatingpoint>>
    getMotorUnbindingData(){
        return motorUnbindingData;
    }
    
    void clearLinkerUnbindingData(){
        linkerUnbindingData.clear();
    }
    
    // store the space time information of a linker binding event to linkerBindingData
    void recordLinkerBinding(Linker* l){
        vector<floatingpoint> lcoords = l->coordinate;
        lcoords.insert(lcoords.begin(), tau());
        linkerBindingData.push_back(make_tuple(l->getId(), lcoords[0], lcoords[1],
                lcoords[2], lcoords[3]));
        

    }
    
    vector<tuple<int, floatingpoint, floatingpoint, floatingpoint, floatingpoint>>
    getLinkerBindingData(){
        return linkerBindingData;
    }
    
    void clearLinkerBindingData(){
        linkerBindingData.clear();
    }

    void clearHRMDMats(){
        HRMDMat.clear();
    }
};

} // namespace medyan

#endif
