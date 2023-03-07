#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneFF_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneFF_hpp

#include <memory> // unique_ptr
#include <stdexcept>
#include <string>
#include <vector>

#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Membrane/FixedVertexAttachmentStretching.hpp"
#include "Mechanics/ForceField/Membrane/FixedVertexCoordinates.hpp"
#include "Mechanics/ForceField/Membrane/MembraneBending.hpp"
#include "Mechanics/ForceField/Membrane/MembraneBendingHelfrich.hpp"
#include "Mechanics/ForceField/Membrane/MembraneMeshless.hpp"
#include "Mechanics/ForceField/Membrane/MembraneStretchingLocal.hpp"
#include "Mechanics/ForceField/Membrane/MembraneStretchingImpl.hpp"
#include "Mechanics/ForceField/Membrane/MembraneStretchingGlobal.hpp"
#include "Mechanics/ForceField/Membrane/MembraneTension.hpp"
#include "Mechanics/ForceField/Membrane/MembraneTriangleProtect.hpp"
#include "Mechanics/ForceField/Types.hpp"
#include "Util/Io/Log.hpp"

namespace medyan {

inline auto createMembraneForceFields(
    std::string_view               stretchingType,
    std::string_view               tensionType,
    std::string_view               bendingType,
    const SurfaceGeometrySettings& gs
) {
    using namespace std;

    std::vector< std::unique_ptr< ForceField > > forceFields;

    // Fixed vertex coordinates are always present.
    forceFields.push_back(std::make_unique<FixedVertexCoordinates>());

    // Whether triangle protection is enabled.
    bool enableTriangleProtect = false;

    if (stretchingType == "LOCAL_HARMONIC") {
        // In material coordinates, it is mandatory, and is applicable to
        //   non-reservior-touching membrane.
        // In normal coordinates, it cannot be modeled without more degrees of freedom.
        forceFields.push_back(
            std::make_unique< MembraneStretchingLocal >()
        );
    }
    else if(stretchingType == "GLOBAL_HARMONIC") {
        // Assumption: surface tension is uniform on the membrane
        // Only applicable to non-reservior-touching membrane in general coordinates.
        forceFields.push_back(
            std::make_unique< MembraneStretchingGlobal >()
        );

        enableTriangleProtect = true;
    }
    else if(stretchingType == "") {}
    else {
        log::error("Membrane stretching FF type {} is not recognized.", stretchingType);
        throw std::runtime_error("Membrane stretching FF type not recognized");
    }

    if(tensionType == "CONSTANT") {
        // In material coordinates, it is applicable to reservior-touching
        // border triangles.
        //
        // In normal or general corodinates, it is applicable to the whole
        // reservior-touching membrane, assuming that surface tension is
        // constant globally.

        forceFields.push_back(
            std::make_unique< MembraneTension >()
        );

        enableTriangleProtect = true;
    }
    else if(tensionType == "") {}
    else {
        log::error("Membrane tension FF type {} is not recognized.", tensionType);
        throw std::runtime_error("Membrane tension FF type not recognized");
    }
    
    if (bendingType == "HELFRICH") {
        if(gs.curvPol == SurfaceCurvaturePolicy::withSign || gs.curvPol == SurfaceCurvaturePolicy::mem3dg) {
            forceFields.emplace_back(
                new MembraneBending<MembraneBendingHelfrich>()
            );
        }
        else {
            forceFields.emplace_back(
                new MembraneBending<MembraneBendingHelfrichQuadratic>()
            );
        }
    }
    else if(bendingType == "") {}
    else {
        log::error("Membrane bending FF type {} is not recognized.", bendingType);
        throw std::runtime_error("Membrane bending FF type not recognized");
    }

    // Enable the protective force field
    if(enableTriangleProtect) {
        forceFields.emplace_back(
            new MembraneTriangleProtect< MembraneTriangleProtectFene, true >()
        );
    }

    // Add vertex attachment stretching force field.
    forceFields.emplace_back(
        new FixedVertexAttachmentStretching()
    );


    return forceFields;
}

} // namespace medyan

#endif
