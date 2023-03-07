
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

#ifndef MEDYAN_common_h
#define MEDYAN_common_h

#include <array>
#include <cstddef>
#include <iostream>
#include <vector>

#include <fmt/core.h>

#include "utility.h"
#include "Util/Io/Log.hpp"

#ifdef CUDAACCL
#include <cuda.h>
#include <cuda_runtime.h>
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

#endif

namespace medyan {

using Index = std::ptrdiff_t;
using Size  = std::ptrdiff_t;

using fmt::format;


// Constants.

constexpr double snan = std::numeric_limits<double>::signaling_NaN();
constexpr float  snanf = std::numeric_limits<float>::signaling_NaN();
constexpr FP     snanfp = std::numeric_limits<FP>::signaling_NaN();

constexpr double inf = std::numeric_limits<double>::infinity();
constexpr float  inff = std::numeric_limits<float>::infinity();
constexpr FP     inffp = std::numeric_limits<FP>::infinity();

} // namespace medyan

///Species constants
typedef unsigned int species_copy_t;
const species_copy_t max_ulim = 1000000;

///Global time
inline floatingpoint global_time = 0;

inline floatingpoint tau() {return global_time;}
inline void resetglobaltime() {global_time=0.0;}
///Some constants
const floatingpoint kT = 4.1; //in pN * nm

const int cylinder_cache = 500;
const int bead_cache = 1000;//total number of beads that can be appended before
// revectorization

///To use STL containers, libraries, etc
using namespace std;

///Num filament types maximum
#define MAX_FILAMENT_TYPES 10

//@{
/// Constant Species index identifiers
/// @note - DO NOT CHANGE!!!

#define SPECIESFILAMENT       0
#define SPECIESPLUSEND        1
#define SPECIESMINUSEND       2

#define SPECIESBOUND          0
#define SPECIESLINKER         1
#define SPECIESMOTOR          2
#define SPECIESBRANCHER       3
//@}

//@{
/// Constant Reaction reactant and product numbers
/// @note - DO NOT CHANGE!!!

/// Polymerization
#define POLYREACTANTS         2
#define POLYPRODUCTS          2

/// Depolymerization
#define DEPOLYREACTANTS       2
#define DEPOLYPRODUCTS        2

/// Linker and motor binding
#define LMBINDINGREACTANTS    3
#define LMBINDINGPRODUCTS     2

/// Linker and motor unbinding
#define LMUNBINDINGREACTANTS  2
#define LMUNBINDINGPRODUCTS   3

/// Motor walking
#define MWALKINGREACTANTS     2
#define MWALKINGPRODUCTS      2

/// Nucleation
#define NUCLEATIONREACTANTS   2
#define NUCLEATIONPRODUCTS    3

/// Destruction
#define DESTRUCTIONREACTANTS  2
#define DESTRUCTIONPRODUCTS   2

/// Aging
#define AGINGREACTANTS        1
#define AGINGPRODUCTS         1

/// Severing
#define SEVERINGREACTANTS     1
#define SEVERINGPRODUCTS      0

/// Branching
#define BRANCHINGREACTANTS    3
#define BRANCHINGPRODUCTS     2

/// Branch unbinding
#define BUNBINDINGREACTANTS   1
#define BUNBINDINGPRODUCTS    2
//@}

#endif
