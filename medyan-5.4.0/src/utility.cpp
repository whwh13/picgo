
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

#include "utility.h"

#include "Util/Environment.hpp"

#ifdef COMPILER_MSVC
    #include <intrin.h>
    #pragma intrinsic(__rdtsc)
#endif

unsigned long long rdtsc(){

#ifdef COMPILER_MSVC
    return __rdtsc();

#else
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;

#endif
}
