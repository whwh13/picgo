#ifndef MEDYAN_Mechanics_ForceField_Types_hpp
#define MEDYAN_Mechanics_ForceField_Types_hpp

// Notes
//   - This file is used to define various types about the force field, used by various files. Do not include project files here unless necessary, to avoid circular includes.

#include <ostream>
#include <string>
#include <tuple>
#include <vector>

#include <Eigen/Core>

#include "common.h"

namespace medyan {
// Forward declaration.
class SubSystem;

struct ForceFieldTypes {
    enum class LoadForceEnd { Plus, Minus };
};

struct EnergyReport {
    struct EachEnergy {
        string name;
        floatingpoint energy;
    };

    floatingpoint total;
    std::vector< EachEnergy > individual;
};

// Format and print EnergyReport with std::ostream.
// Note: This function should be replaced by specializing std::formatter in C++20.
inline std::ostream& operator<<(std::ostream& os, const EnergyReport& er) {
    const int spaceBeforeName = 4;
    const int spaceAfterName = 1;
    const int nameMaxLength = 17;
    os << "Total energy = " << er.total << '\n';
    for(auto& indi : er.individual) {
        os << std::string(spaceBeforeName, ' ') << indi.name;
        if(indi.name.size() > nameMaxLength) {
            os << '\n' << std::string(spaceBeforeName + nameMaxLength + spaceAfterName, ' ');
        }
        else {
            os << std::string(nameMaxLength - indi.name.size() + spaceAfterName, ' ');
        }
        os << indi.energy << '\n';
    }
    return os;
}

// This information is used in the force field vectorization.
struct FFCoordinateStartingIndex {
    // The Bezier model of filament requires a vector (matrix) indicating bead
    // positions, angles and control nodes.
    struct FilVec {
        // Starting index in the vectorized array
        std::size_t start = 0;
        // Length of data in the vectorized array
        std::size_t length = 0;
    };

    // The independent variable indices for the GC filament model.
    //
    // Number of segments: N
    //
    // Independent variables:
    // Length of a: 3(N+1)
    // Length of zeta: 3N
    // Length of dof: 3 + 3(N+1) + 3N
    //
    // Dependent variables:
    // Length of dependent knot coordinates: 3N
    struct FilGCVec {
        int rMinusStart = 0;  // length is 3
        int numSegments = 0;

        // The starting index of dependent bead coordinates
        int rNotMinusStart = 0;

        // Auxiliary functions for data lengths
        int rMinusLength()    const { return 3; }
        int aLength()         const { return 3 * (numSegments + 1); }
        int zetaLength()      const { return 3 * numSegments; }
        int rNotMinusLength() const { return 3 * numSegments; }

        // Auxiliary functions for data starting indices
        int aStart()    const { return rMinusStart + rMinusLength(); }
        int zetaStart() const { return aStart() + aLength(); }

        template< typename Float >
        auto mapData(Float* data) const {
            return std::tuple {
                // rMinus
                Eigen::Vector3d::Map(data + rMinusStart),
                // a
                Eigen::Matrix3Xd::Map(data + aStart(),    3, numSegments + 1),
                // zeta
                Eigen::Matrix3Xd::Map(data + zetaStart(), 3, numSegments),
            };
        }
        template< typename Float >
        auto mapDepData(Float* data) const {
            // rNotMinus
            return Eigen::Matrix3Xd::Map(data + rNotMinusStart, 3, numSegments);
        }
    };

    // Pointer to the system.
    SubSystem* ps = nullptr;

    // If Bezier filament is used, beads will NOT be stored according to their
    // indices. Instead, the information of filamental beads will be stored
    // with individual filaments.
    std::size_t bead = 0;
    Index movableBubble = 0;
    Index fixedBubble = 0;
    Index movableVertex = 0;
    Index fixedVertex = 0;
    Index meshlessSpinVertex = 0;
    std::size_t mem2d = 0;

    // The following will be used only when the Bezier filament model is used.
    std::vector< FilVec > bezierFilPos;
    std::vector< FilVec > bezierFilAng;

    // The GC filament model data indices.
    std::vector< FilGCVec > gcFil;

    // Total number of degrees of freedom. It is not necessarily the size of vectorized data, as the vectorized data may contain dependent variables.
    int ndof = 0;

};

} // namespace medyan

#endif
