#ifndef MEDYAN_Structure_SurfaceMesh_MMembrane_hpp
#define MEDYAN_Structure_SurfaceMesh_MMembrane_hpp

namespace medyan {

/******************************************************************************
Stores mechanical properties of the membrane.
******************************************************************************/

struct MMembrane {

    double kVolume  = 0; // The elastic constant of volume
    double eqVolume = 1; // The equilibrium volume
    // The final volume will be the computed volume plus this offset.
    // This offset should only be used in open, fixed-border membranes, where the enclosing volume is ill-defined, but the tetrahedra-tiled volume is still defined.
    // Example:
    //   In a box of size 1x1x1, there is an open membrane at z=0.5, with border vertices fixed.
    //   The tetrahedra-tiled volume is 1/6, but if we need the entire volume "below" the membrane, we need to offset the volume by 2/6.
    double volumeOffset = 0;

    double eqArea   = 1; // The equilibrium area, used in quadratic global stretching energy
    double kArea    = 0; // Area elasticity, used in quadratic global stretching energy

    double tension  = 0; // Membrane tension, used in linear stretching energy.
                         // Note that the tension can also be defined in quadratic energy case,
                         // but is not stored in this variable.

    // Bending constant and spontaneous curvature in the Helfrich model.
    double kBending = 0;
    double eqCurv   = 0;

    // Area per lipid (for both layers).
    // Ref: Area per Lipid and Cholesterol Interactions in Membranes from Separated Local-Field 13C NMR Spectroscopy (2014) Leftin et al. Biophys J.
    double areaPerLipid = 0.6;
};

// This structure will be used for every curvature inducing protein on the membrane.
struct ProteinCurvatureMismatchParams {
    // The index of surface diffusing species (same for all membranes). -1 is sentinel value for invalid index.
    int    speciesIndex = -1;
    double kBending = 0;
    double eqCurv = 0;
    double area = 0;
};

} // namespace medyan

#endif
