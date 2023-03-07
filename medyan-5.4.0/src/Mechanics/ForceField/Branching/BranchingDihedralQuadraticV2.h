#ifndef MEDYAN_Mechanics_ForceField_Branching_BranchingDihedralQuadraticV2_H
#define MEDYAN_Mechanics_ForceField_Branching_BranchingDihedralQuadraticV2_H

#include "common.h" // floatingpoint

namespace medyan {
struct BranchingDihedralQuadraticV2 {
    floatingpoint energy(
        const floatingpoint *coord, size_t nint,
        const unsigned int *beadSet, const floatingpoint *kdih, const floatingpoint *pos
    ) const;

    void forces(
        const floatingpoint *coord, floatingpoint *f, size_t nint,
        const unsigned int *beadSet, const floatingpoint *kdih, const floatingpoint *pos,
        floatingpoint *stretchforce) const;
};

} // namespace medyan

#endif
