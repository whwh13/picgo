#ifndef MEDYAN_Parameter_Membrane_hpp
#define MEDYAN_Parameter_Membrane_hpp

#include "common.h"
#include "Structure/SurfaceMesh/Types.hpp"

namespace medyan {

/// Struct to hold membrane setup information
struct MembraneSetup {
    // meta info
    //---------------------------------
    std::string name;
    int type = 0;

    // mechanical parameters
    //---------------------------------
    // mesh mode
    MembraneMeshVertexSystem vertexSystem = MembraneMeshVertexSystem::general;

    // Vertex pinning.
    MembraneVertexPinning vertexPinning = MembraneVertexPinning::border1;

    // Currently, the vertex system requires certain reservoir property.
    // - material: no reservoir.
    // - normal: has reservoir.
    // - general: has reservoir.
    bool hasLipidReservoir = true;

    // elasticity
    double areaElasticity = 400;
    double bendingElasticity = 100;
    double eqMeanCurv = 0;

    // tension
    double tension = 0;

    // enclosed volume conservation
    double volumeElasticity = 0.8;
};

// Each membrane object is initialized using the following parameters.
struct MembraneInit {
    struct SpeciesInit {
        std::string name;
        int         copyNumber = 0;
    };

    // Membrane profile name.
    std::string                name;

    // Mesh information, initialized by shape or from file.
    std::vector< std::string > meshParams;

    // initial membrane area factor compared to equilibrium area
    double                     eqAreaFactor = 1.0;

    // Volume offset for open, fixed-border membranes. See MMembrane::volumeOffset.
    double                     volumeOffset = 0.0;

    // Chemistry information.
    std::vector< SpeciesInit > speciesInitVec;
};


// Parameters for initializing a FixedVertexAttachment.
struct FixedVertexAttachmentInitSearch {
    // Initial coordinate.
    Vec<3, FP> coord {};

    // If no vertices exist within the range, the attachment is not created.
    FP range = 0;

    // Mechanical properties.
    FP kStretch = 0;
};

// Parameters for initializing many FixedVertexAttachments.
struct FixedVertexAttachmentInitSearchMultiple {
    // Initial coordinates.
    std::vector<FixedVertexAttachmentInitSearch> attachments;

    // If no vertices exist within the range, warn or throw.
    bool throwOnNotFound = true;
};

struct MembraneSettings {
    std::vector< MembraneSetup > setupVec;
    std::vector< MembraneInit >  initVec;

    FixedVertexAttachmentInitSearchMultiple fixedVertexAttachmentInits;
};

} // namespace medyan

#endif
