#ifndef MEDYAN_Side_SideProcedures_hpp
#define MEDYAN_Side_SideProcedures_hpp

#include <string>

#include "Side/SurfaceMeshDiffusion.hpp"
#include "Util/Io/Log.hpp"

namespace medyan::side {


inline void runSideProcedure(std::string procName) {
    if(procName == "mesh-diffusion") {
        surfaceMeshDiffusion();
    }
    else {
        log::error("Unrecognized side procedure {}", procName);
    }
}

} // namespace medyan::side

#endif
