#ifndef MEDYAN_Structure_MemFilConnOp_hpp
#define MEDYAN_Structure_MemFilConnOp_hpp

#include "Structure/Cylinder.h"
#include "Structure/MemFilConn.hpp"

namespace medyan {
inline auto getCoordinateOnFilamentEnd(const MemFilConn& conn) {
    return (1 - conn.coordCylinder) * conn.cylinder->getFirstBead()->coordinate()
              + conn.coordCylinder  * conn.cylidner->getSecondBead()->coordinate();
}

inline auto getCoordinateOnMembraneEnd(const MemFilConn& conn) {

}

} // namespace medyan

#endif
