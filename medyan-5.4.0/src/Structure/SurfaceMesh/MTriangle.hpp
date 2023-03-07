#ifndef MEDYAN_Structure_SurfaceMesh_MTriangle_Hpp
#define MEDYAN_Structure_SurfaceMesh_MTriangle_Hpp

/******************************************************************************

Storing some mechanical properties of the triangle patches

******************************************************************************/

namespace medyan {

struct MTriangle {

    // Local area elasticity, applicable only in material coordinates
    double kArea = 0.0;
    double eqArea = 1.0;
};

} // namespace medyan

#endif
