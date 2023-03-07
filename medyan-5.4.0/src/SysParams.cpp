
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

#include <sstream>

#include "MathFunctions.h"
#include "SysParams.h"
#include "Util/Io/Log.hpp"

namespace medyan {

void postprocess(SimulConfig& conf) {
    using namespace std;

    auto& geo = conf.geoParams;
    auto& bound = conf.boundParams;
    auto& chem = conf.chemParams;

    //--------------------------------------------------------------------------
    // Geometry.
    //--------------------------------------------------------------------------
    if(geo.cylinderSize.size() != geo.monomerSize.size()) {
        log::error("Number of cylinder sizes and number of monomer sizes must be equal.");
        throw runtime_error("Invalid geometry parameter.");
    }
    for(int i = 0; i < geo.cylinderSize.size(); i++) {

        if(geo.cylinderSize[i] / geo.monomerSize[i] < geo.minCylinderNumMon) {
            log::error("With chemistry, cylinder size specified is too short. Exiting.");
            throw runtime_error("Invalid geometry parameter.");
        }
        geo.cylinderNumMon.push_back(int(geo.cylinderSize[i] / geo.monomerSize[i]));

        geo.minCylinderSize.push_back(geo.minCylinderNumMon * geo.monomerSize[i]);

    }

    //find max compartment side
    geo.largestCompartmentSide = max({
        geo.compartmentSizeX,
        geo.compartmentSizeY,
        geo.compartmentSizeZ
    });
    //find max Cylinder size
    geo.largestCylinderSize = 0;
    for(auto l : geo.cylinderSize)
        geo.largestCylinderSize = max(geo.largestCylinderSize, l);


    //--------------------------------------------------------------------------
    // Boundary parameters.
    //--------------------------------------------------------------------------
    if(const auto& bm = bound.boundaryType.boundaryMove; !bm.empty()) {
        vector<int> leftfrontbottom = {0,0,0};
        vector<int> rightbacktop = {0,0,0};

        for(const auto& eachBM : bm) {
            if(eachBM == "LEFT")
                leftfrontbottom[0] = 1;
            else if(eachBM == "BOTTOM")
                leftfrontbottom[1] = 1;
            else if(eachBM == "FRONT")
                leftfrontbottom[2] = 1;
            else if(eachBM == "RIGHT")
                rightbacktop[0] = 1;
            else if(eachBM == "TOP")
                rightbacktop[1] = 1;
            else if(eachBM == "BACK")
                rightbacktop[2] = 1;
        }

        for(int i = 0; i < 3; i++){
            int addthemup = leftfrontbottom[i] + rightbacktop[i];
            if(addthemup > 0)
                bound.transfershareaxis = i;
            if(addthemup == 2)
                bound.planestomove = 2;
            else if(leftfrontbottom[i] == 1)
                bound.planestomove = 1;
            else if(rightbacktop[i] == 1)
                bound.planestomove = 0;
        }
    }


    //--------------------------------------------------------------------------
    // Chemistry parameters.
    //--------------------------------------------------------------------------
    // Figure out the binding sites.
    chem.bindingSites.clear();
    for(int i = 0; i < chem.numBindingSites.size(); i++) {

        chem.maxbindingsitespercylinder = max(chem.maxbindingsitespercylinder, chem.numBindingSites[i]);

        vector<short> tempBindingSites;

        int deltaBinding = geo.cylinderNumMon[i] /
                           chem.numBindingSites[i];

        int firstBindingSite = deltaBinding / 2 + 1;
        int bindingCount = firstBindingSite;

        //add all other binding sites
        while(bindingCount < geo.cylinderNumMon[i]) {
            tempBindingSites.push_back(bindingCount);
            bindingCount += deltaBinding;
        }

        // Push to binding sites.
        chem.bindingSites.push_back(move(tempBindingSites));
    }

    // Find the maximum allowed Cindex and shift operator.
    {
        auto np2 = mathfunc::nextPowerOf2(uint32_t(chem.maxbindingsitespercylinder));

        if(np2 == chem.maxbindingsitespercylinder)
            np2 *= 2;

        cout<<"np2 "<<np2<<" shift "<< chem.shiftbybits<<endl;
        chem.shiftbybits = log2(np2);
        chem.maxStableIndex = numeric_limits<uint32_t>::max()/chem.shiftbybits -1;
    }

    //--------------------------------------------------------------------------
    // Chemistry data.
    //--------------------------------------------------------------------------
    conf.membraneMeshChemistryInfo = MembraneMeshChemistryInfo::fromChemistryData(conf.chemistryData);

} // End of postprocess.


bool SysParams::RUNSTATE=true;
bool SysParams::INITIALIZEDSTATUS=false;
bool SysParams::DURINGCHEMISTRY=false;

int SysParams::exvolcounter[3] = {0,0,0};
long long SysParams::exvolcounterz[3] = {0,0,0};
#ifdef NLSTENCILLIST
short SysParams::numcylcylNL = 0;
#endif
vector<float> SysParams::BUBBareRate ={};
ChemParams   SysParams::CParams;
GeoParams    SysParams::GParams;
BoundParams  SysParams::BParams;
DyRateParams SysParams::DRParams;
#ifdef TRACKDIDNOTMINIMIZE
MinimizationParams SysParams::MinParams;
#endif

} // namespace medyan
