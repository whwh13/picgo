//
// Created by aravind on 9/19/17.
//

#include "CUDAcommon.h"

namespace medyan {

Callbacktime CUDAcommon::ctime;
Callbackcount CUDAcommon::ccount;
PolyPlusEndTemplatetime CUDAcommon::ppendtime;
timeminimization CUDAcommon::tmin;
motorwalking CUDAcommon::mwalk;
chemdetails CUDAcommon::cdetails;
#if defined(CUDAACCL)
CUDAvars  CUDAcommon::cudavars;
CylCylNLvars CUDAcommon::cylcylnlvars;
SERLtime CUDAcommon::serltime;
CUDAtime CUDAcommon::cudatime;
#endif

} // namespace medyan
