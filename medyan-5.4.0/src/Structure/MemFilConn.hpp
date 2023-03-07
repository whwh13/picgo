#ifndef MEDYAN_Structure_MemFilConn_hpp
#define MEDYAN_Structure_MemFilConn_hpp

#include "Structure/Database.h"

namespace medyan {
// Forward declarations
class Cylinder;
class Membrane;
class Triangle;

struct MemFilConn : Database< MemFilConn, false > {
    // The connection data
    //---------------------------------

    // On filament
    Cylinder* cylinder = nullptr;
    double coordCylinder = 0.0;

    // On membrane triangle
    Membrane* membrane = nullptr;
    Triangle* triangle = nullptr;
    double coordMembrane[2] {};

    // Mechanical properties
    //---------------------------------
    double springConstant = 0.0;

};

} // namespace medyan

#endif
