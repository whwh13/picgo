
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------
#ifndef cross_check_h
#define cross_check_h

#include "Cylinder.h"
#include "Bead.h"
#include "Filament.h"
#include <stdio.h>
#include <limits>

namespace medyan {
typedef std::numeric_limits< floatingpoint > dbl;
class cross_checkclass{
    public:
    inline static bool Aux = false;
};

namespace cross_check{
inline bool crosscheckforces(floatingpoint* force){
    bool state=false;
    for(auto b: Bead::getBeads()) {
        
        //set bead index
        auto idx=3*b->getStableIndex();
        cout.precision(dbl::max_digits10);
//        if(force[idx]!=b->force[0])
//            std::cout<<"0"<<endl;
//        if(force[idx+1]!=b->force[1])
//            std::cout<<"1"<<endl;
//        if(force[idx+2]!=b->force[2])
//            std::cout<<"2"<<endl;
        
        /* Haoran 03/18/2019
        if(force[idx]==b->force[0] && force[idx+1]==b->force[1] && force[idx+2]==b->force[2])
            state=true;
        else{
            state=false;
            std::cout<<"vectorized "<<force[idx]<<" "<<force[idx+1]<<" "<<force[idx+2]<<endl;
            std::cout<<"old way "<< b->force() <<endl;
            exit(EXIT_FAILURE);
        }
        */
    }
    return state;
}
    inline bool crosscheckAuxforces(floatingpoint* force){
        bool state=false;
        for(auto b: Bead::getBeads()) {
            
            //set bead index
            auto idx=3*b->getStableIndex();
            cout.precision(dbl::max_digits10);
            //        if(force[idx]!=b->force[0])
            //            std::cout<<"0"<<endl;
            //        if(force[idx+1]!=b->force[1])
            //            std::cout<<"1"<<endl;
            //        if(force[idx+2]!=b->force[2])
            //            std::cout<<"2"<<endl;
            
            if(force[idx]==b->force[0] && force[idx+1]==b->force[1] && force[idx+2]==b->force[2])
                state=true;
            else{
                state=false;
                std::cout<<"vectorized "<<force[idx]<<" "<<force[idx+1]<<" "<<force[idx+2]<<endl;
                std::cout<<"old way "<< b->force <<endl;
                exit(EXIT_FAILURE);
            }
        }
        return state;
    }
}

} // namespace medyan

#endif /* cross_check_h */
