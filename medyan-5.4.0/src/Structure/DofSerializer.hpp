#ifndef MEDYAN_Structure_DofSerializer_hpp
#define MEDYAN_Structure_DofSerializer_hpp

// The file provides functions to do data copy between elements in the system
// and the vectorized data in the mechanics energy/force calculations.
//
// Functions:
//
//   - serializeDof(...)
//
//     Copy the degree-of-freedom data from all system elements (such as the
//     bead coordinates) to the coordinate array.
//
//     Returns the starting indices of different types of elements, which is
//     useful when building interactions in force fields.
//
//   - deserializeDof(...)
//
//     Copy the vectorized coordinate and force data to all the element
//     instances in the system.

#include <algorithm>

#include "Mechanics/ForceField/Types.hpp"
#include "Structure/Bead.h"
#include "Structure/Bubble.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Vertex.hpp"
#include "Structure/SubSystem.h"

namespace medyan {

// Update pinned state of all vertices.
// This is an auxiliary function before DOF serialization.
inline void updateVertexPinning(SubSystem& sys) {
    // Reset pinning.
    for(auto& v : sys.vertices) {
        v.pinned = false;
    }

    // Update pinning based on border conditions.
    for(auto& m : sys.membranes) {
        auto& mesh = m.getMesh();
        if(m.getSetup().vertexPinning == MembraneVertexPinning::border1 || m.getSetup().vertexPinning == MembraneVertexPinning::border2) {
            // Pin border vertices.
            for(auto& v : mesh.getVertices()) {
                if(mesh.isVertexOnBorder(v)) {
                    v.attr.vertex(sys).pinned = true;

                    // Further pin all neighbors if border2 is selected.
                    if(m.getSetup().vertexPinning == MembraneVertexPinning::border2) {
                        mesh.forEachHalfEdgeTargetingVertex(v, [&](auto hei) {
                            auto vin = mesh.target(mesh.opposite(hei));
                            mesh.attribute(vin).vertex(sys).pinned = true;
                        });
                    }
                }
            }
        }
    }
}

// Copies all the system data to the CGMethod data vector
inline FFCoordinateStartingIndex serializeDof(
    SubSystem&         sys,
    std::vector< FP >& coord
) {
    FFCoordinateStartingIndex si {};
    si.ps = &sys;
    Index curIdx = 0;
    coord.clear();

    updateVertexPinning(sys);

    //---------------------------------
    // Copy all the coordinate information here
    // Also initializes the force tolerance vector

    // Bead coord
    si.bead = curIdx;
    coord.reserve(coord.size() + 3 * Bead::getBeads().size());
    for(auto pb : Bead::getBeads()) {
        coord.insert(coord.end(), pb->coord.begin(), pb->coord.end());
        curIdx += 3;
    }

    // (Moveble) bubbles.
    si.movableBubble = curIdx;
    {
        Index loopIndex = 0;
        for(auto& b : sys.bubbles) {
            if(!b.fixed) {
                coord.insert(coord.end(), b.coord.begin(), b.coord.end());
                // Update contiguous looping index inside the bubble as well.
                b.loopIndex = loopIndex++;
                curIdx += 3;
            }
        }
    }

    // (Movable) vertex coords.
    si.movableVertex = curIdx;
    coord.reserve(coord.size() + 3 * sys.vertices.size());
    {
        Index loopIndex = 0;
        for(auto& m : sys.membranes) {
            auto& mesh = m.getMesh();
            for(Membrane::MeshType::VertexIndex vi{0}; vi < mesh.numVertices(); ++vi) {
                auto& v = mesh.attribute(vi).vertex(sys);
                if(!v.pinned) {
                    // Also caches coordinate index in the mesh structure.
                    mesh.attribute(vi).cachedCoordIndex = curIdx;

                    // Copy coordinates.
                    coord.insert(coord.end(), v.coord.begin(), v.coord.end());

                    curIdx += 3;
                }

                // Update contiguous looping index inside the vertex as well.
                v.loopIndex = loopIndex++;
            }
        }
    }

    // Meshless spin vertex coord.
    si.meshlessSpinVertex = curIdx;
    coord.reserve(coord.size() + 5 * sys.meshlessSpinVertices.size());
    {
        Index loopIndex = 0;
        for(auto& v : sys.meshlessSpinVertices) {
            coord.insert(coord.end(), v.coord.data(), v.coord.data() + v.coord.size());
            coord.push_back(v.theta);
            coord.push_back(v.phi);
            // Update contiguous looping index inside the vertex as well.
            v.loopIndex = loopIndex++;
            curIdx += 5;
        }
    }

    // Membrane 2d coord
    si.mem2d = curIdx;
    // Add things for membrane 2d coordinates here

    //---------------------------------
    // The independent variables have all been assigned.
    si.ndof = curIdx;

    //---------------------------------
    // Now assign space for dependent variables.

    // Fixed bubbles.
    si.fixedBubble = curIdx;
    {
        Index loopIndex = 0;
        for(auto& b : sys.bubbles) {
            if(b.fixed) {
                coord.insert(coord.end(), b.coord.begin(), b.coord.end());
                // Update contiguous looping index inside the bubble as well.
                b.loopIndex = loopIndex++;
                curIdx += 3;
            }
        }
    }

    // Fixed vertex coords.
    si.fixedVertex = curIdx;
    for(auto& m : sys.membranes) {
        auto& mesh = m.getMesh();
        for(Membrane::MeshType::VertexIndex vi{0}; vi < mesh.numVertices(); ++vi) {
            auto& v = mesh.attribute(vi).vertex(sys);
            if(v.pinned) {
                // Also caches coordinate index in the mesh structure.
                mesh.attribute(vi).cachedCoordIndex = curIdx;

                // Copy coordinates.
                coord.insert(coord.end(), v.coord.begin(), v.coord.end());

                curIdx += 3;
            }
        }
    }

    //---------------------------------
    // Return the starting index information for vectorizing the force fields
    return si;
}

// Copies all the CGMethod data back to the system
//
// Note:
//   - The copying must be in the same order with the initCGMethodData
//     function.
inline void deserializeDof(
    SubSystem&               sys,
    const std::vector< FP >& coord,
    const std::vector< FP >& force
) {
    std::size_t curIdx = 0;

    // Copy coord and force data to beads
    for(auto pb : Bead::getBeads()) {
        std::copy(coord.begin() + curIdx, coord.begin() + curIdx + 3, pb->coord.begin());
        std::copy(force.begin() + curIdx, force.begin() + curIdx + 3, pb->force.begin());
        curIdx += 3;
    }

    // Copy coord and force data to movable bubbles
    for(auto& b : sys.bubbles) {
        if(!b.fixed) {
            std::copy(coord.begin() + curIdx, coord.begin() + curIdx + 3, b.coord.begin());
            std::copy(force.begin() + curIdx, force.begin() + curIdx + 3, b.force.begin());
            curIdx += 3;
        }
    }

    // Movable vertex coords.
    for(auto& m : sys.membranes) {
        auto& mesh = m.getMesh();
        for(Membrane::MeshType::VertexIndex vi {0}; vi < mesh.numVertices(); ++vi) {
            if(!mesh.attribute(vi).vertex(sys).pinned) {
                auto& v = sys.vertices[mesh.attribute(vi).vertexSysIndex];
                std::copy(coord.begin() + curIdx, coord.begin() + curIdx + 3, v.coord.begin());
                std::copy(force.begin() + curIdx, force.begin() + curIdx + 3, v.force.begin());
                curIdx += 3;
            }
        }
    }

    // Copy coord data back to meshless spin vertices.
    for(auto& v : sys.meshlessSpinVertices) {
        std::copy(coord.begin() + curIdx, coord.begin() + curIdx + 3, v.coord.data());
        v.theta = coord[curIdx + 3];
        v.phi   = coord[curIdx + 4];
        curIdx += 5;
    }

    // Copy coord and force data to Membrane 2d points


    //---------------------------------
    // Now copy back dependent variables.

    // Fixed bubbles.
    for(auto& b : sys.bubbles) {
        if(b.fixed) {
            // Do not copy any data since the bubble is fixed anyway.
            curIdx += 3;
        }
    }

    // Fixed vertex coords.
    for(auto& m : sys.membranes) {
        auto& mesh = m.getMesh();
        for(Membrane::MeshType::VertexIndex vi{0}; vi < mesh.numVertices(); ++vi) {
            if(mesh.attribute(vi).vertex(sys).pinned) {
                // Do not copy any data since the vertex is fixed anyway.
                curIdx += 3;
            }
        }
    }
}


// Helper function to find the starting index of bead coordinate in the minimizer.
inline int findBeadCoordIndex(
    Bead&                            b,
    const FFCoordinateStartingIndex& si
) {
    // if(SysParams::Mechanics().globalFilamentModel == FilamentModel::bezier) {
    //     const auto& fil = *static_cast< Filament* >(b.getParent());
    //     return 
    //         si.bezierFilPos[fil.getIndex()].start
    //         + fil.getCoordIndexInKvecrWithBeadPosition(b.getPosition());
    // }
    // else if(SysParams::Mechanics().globalFilamentModel == FilamentModel::geodesic) {
    //     const auto& fil = *static_cast< Filament* >(b.getParent());
    //     const auto& gcInfo = si.gcFil[fil.getIndex()];
    //     const int   relPos = fil.getRelativePositionFromMinusEnd(b.getPosition());
    //     return relPos == 0
    //         // bead is the minus end bead
    //         ? gcInfo.rMinusStart
    //         // bead is not the minus end
    //         : gcInfo.rNotMinusStart + 3 * (relPos - 1);
    // }
    // else {
    {
        return b.getIndex() * 3 + si.bead;
    }
}

// Helper function to find the starting index of bubble coordinate in the minimizer.
// Requires that the loopIndex is up-to-date.
inline Index findBubbleCoordIndex(Bubble& b, const FFCoordinateStartingIndex& si) {
    return (b.fixed ? si.fixedBubble : si.movableBubble) + b.loopIndex * 3;
}
// Not providing a function to find the starting index of vertex coordinate in the minimizer.
// Because the vertex coordinate index is already cached in mesh structures, and is not associated with loopIndex.

// Helper function to find the starting index of meshless vertex coordinate in the minimizer.
// Requires that loopIndex is up-to-date.
inline Index findMeshlessSpinVertexCoordIndex(MeshlessSpinVertex& v, const FFCoordinateStartingIndex& si) {
    return si.meshlessSpinVertex + v.loopIndex * 5;
}

} // namespace medyan

#endif
