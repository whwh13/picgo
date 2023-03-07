
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2017-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Mechanics_ForceField_Volume_TriangleBeadVolumeFF_hpp
#define MEDYAN_Mechanics_ForceField_Volume_TriangleBeadVolumeFF_hpp

#include <stdexcept>
#include <string>
#include <vector>

#include "common.h"
#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Volume/TriangleBeadExclVolume.hpp"

namespace medyan {

struct TriangleBeadVolumeFFFactory {
    auto operator()(
        SubSystem&         sys,
        const std::string& type
    ) const {
        using namespace std;

        vector< unique_ptr< ForceField > > res;

        if (type == "REPULSION")
            res.emplace_back(
                new TriangleBeadExclVolume(sys));
        else if(type == "") {}
        else {
            LOG(ERROR) << "Volume FF " << type << " is not recognized.";
            throw runtime_error("Volume FF type not recognized");
        }

        return res;
    }

};

} // namespace medyan

#endif
