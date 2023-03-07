#ifndef MEDYAN_Structure_Output_OMembrane_hpp
#define MEDYAN_Structure_Output_OMembrane_hpp

#include "Structure/OutputStruct.hpp"
#include "Structure/SubSystem.h"
#include "Util/Io/H5.hpp"

namespace medyan {

// Extract membrane meta data.
inline void extract(OutputStructMembraneMeta& outMemMeta, const MembraneMeshChemistryInfo& meshChemInfo) {

    outMemMeta.vertexColumnNamesFloat64 = { "coord.x", "coord.y", "coord.z", "area.", "curvature.curv" };
    for(int i = 0; i < meshChemInfo.diffusingSpeciesNames.size(); ++i) {
        auto& name = meshChemInfo.diffusingSpeciesNames[i];
        outMemMeta.vertexColumnNamesFloat64.push_back("conc." + name);
    }

    outMemMeta.vertexColumnNamesInt64.clear();
    for(int i = 0; i < meshChemInfo.diffusingSpeciesNames.size(); ++i) {
        auto& name = meshChemInfo.diffusingSpeciesNames[i];
        outMemMeta.vertexColumnNamesInt64.push_back("copyNumber." + name);
    }
}
inline void extract(OutputStructMembrane& outMem, const MembraneMeshChemistryInfo& meshChemInfo, const Membrane& mem, const SubSystem& sys) {
    using HI = Membrane::MeshType::HalfEdgeIndex;
    using TI = Membrane::MeshType::TriangleIndex;
    using VI = Membrane::MeshType::VertexIndex;
    auto& mesh = mem.getMesh();

    outMem.type = mem.getSetup().type;

    outMem.numVertices = mesh.numVertices();
    outMem.numTriangles = mesh.numTriangles();
    outMem.numBorders = mesh.numBorders();

    // Vertex data.
    const Size numVertAttrF64 = 5 + meshChemInfo.diffusingSpeciesNames.size();
    outMem.vertexDataFloat64 = Eigen::MatrixXd(numVertAttrF64, outMem.numVertices);
    for(int i = 0; i < outMem.numVertices; ++i) {
        auto&      gv = mesh.attribute(VI{i}).gVertex;
        auto&      v = mesh.attribute(VI{i}).vertex(sys);
        auto&      c = v.coord;
        auto&      cv = v.cVertex;
        const auto area = gv.astar / 3;

        for(int dim = 0; dim < 3; ++dim) {
            outMem.vertexDataFloat64(dim, i) = c[dim];
        }
        outMem.vertexDataFloat64(3, i) = area;
        outMem.vertexDataFloat64(4, i) = gv.curv;
        for(int si = 0; si < meshChemInfo.diffusingSpeciesNames.size(); ++si) {
            outMem.vertexDataFloat64(5 + si, i) = cv.species.findSpeciesByIndex(si)->getRSpecies().getTrueN() / area;
        }
    }

    const auto vertexDataInt64NumRows = meshChemInfo.diffusingSpeciesNames.size();
    outMem.vertexDataInt64 = Eigen::MatrixXi(vertexDataInt64NumRows, outMem.numVertices);
    for(int i = 0; i < outMem.numVertices; ++i) {
        auto& cv = mesh.attribute(VI{i}).vertex(sys).cVertex;
        for(int si = 0; si < meshChemInfo.diffusingSpeciesNames.size(); ++si) {
            outMem.vertexDataInt64(si, i) = cv.species.findSpeciesByIndex(si)->getRSpecies().getTrueN();
        }
    }

    // Triangle data.
    outMem.triangleDataInt64 = Eigen::MatrixXi(3, outMem.numTriangles);
    for(TI ti {0}; ti < outMem.numTriangles; ++ti) {
        Index i = 0;
        mesh.forEachHalfEdgeInTriangle(ti, [&](HI hei) {
            outMem.triangleDataInt64(i++, ti.index) = mesh.target(hei).index;
        });
    }

    // Packed border vertices.
    outMem.packedBorderVertices.clear();
    for(Index bi = 0; bi < outMem.numBorders; ++bi) {
        auto& border = mesh.getBorders()[bi];
        Index currentBorderCountIndex = outMem.packedBorderVertices.size();

        // Use the next element to store number of vertices in this border.
        // To be counted later.
        outMem.packedBorderVertices.push_back(0);

        Size nv = 0;
        mesh.forEachHalfEdgeInPolygon(border, [&](HI hei) {
            outMem.packedBorderVertices.push_back(mesh.target(hei).index);
            ++nv;
        });
        outMem.packedBorderVertices[currentBorderCountIndex] = nv;
    }
}
inline void extract(std::vector<OutputStructMembrane>& outMems, const MembraneMeshChemistryInfo& meshChemInfo, const SubSystem& sys) {
    outMems.clear();
    outMems.reserve(sys.membranes.size());
    for(auto& m : sys.membranes) {
        OutputStructMembrane outMem;
        extract(outMem, meshChemInfo, m, sys);
        outMems.push_back(std::move(outMem));
    }
}

// Write membrane meta data.
inline void write(h5::Group& grpMemMeta, const OutputStructMembraneMeta& memMeta) {
    h5::writeDataSet(grpMemMeta, "vertexColumnNamesFloat64", memMeta.vertexColumnNamesFloat64);
    h5::writeDataSet(grpMemMeta, "vertexColumnNamesInt64", memMeta.vertexColumnNamesInt64);
}
inline void write(h5::Group& grpMems, const OutputStructMembrane& outMem, Index memIndex) {
    auto grpMem = grpMems.createGroup(to_string(memIndex));
    h5::writeDataSet<std::int64_t>(grpMem, "type",         outMem.type);
    h5::writeDataSet<std::int64_t>(grpMem, "numVertices",  outMem.numVertices);
    h5::writeDataSet<std::int64_t>(grpMem, "numTriangles", outMem.numTriangles);
    h5::writeDataSet<std::int64_t>(grpMem, "numBorders",   outMem.numBorders);
    h5::writeDataSet(grpMem, "vertexDataFloat64", outMem.vertexDataFloat64);
    h5::writeDataSet(grpMem, "vertexDataInt64", outMem.vertexDataInt64);
    h5::writeDataSet(grpMem, "triangleDataInt64", outMem.triangleDataInt64);
    h5::writeDataSet(grpMem, "packedBorderVertices", outMem.packedBorderVertices);
}
inline void write(h5::Group& grpSnapshot, const std::vector<OutputStructMembrane>& outMems) {
    using namespace h5;
    auto grpMem = grpSnapshot.createGroup("membranes");
    h5::writeDataSet<std::int64_t>(grpMem, "count", outMems.size());
    for(Index i = 0; i < outMems.size(); ++i) {
        write(grpMem, outMems[i], i);
    }
}

// Read membrane meta data.
inline void read(OutputStructMembraneMeta& outMemMeta, const h5::Group& grpMemMeta) {
    h5::readDataSet(outMemMeta.vertexColumnNamesFloat64, grpMemMeta, "vertexColumnNamesFloat64");
    h5::readDataSet(outMemMeta.vertexColumnNamesInt64, grpMemMeta, "vertexColumnNamesInt64");
}
inline void read(OutputStructMembrane& outMem, const h5::Group& grpMems, Index memIndex) {
    auto grpMem = grpMems.getGroup(toString(memIndex));
    outMem.type         = h5::readDataSet<std::int64_t>(grpMem, "type");
    outMem.numVertices  = h5::readDataSet<std::int64_t>(grpMem, "numVertices");
    outMem.numTriangles = h5::readDataSet<std::int64_t>(grpMem, "numTriangles");
    outMem.numBorders   = h5::readDataSet<std::int64_t>(grpMem, "numBorders");
    h5::readDataSet(outMem.vertexDataFloat64, grpMem, "vertexDataFloat64");
    h5::readDataSet(outMem.vertexDataInt64, grpMem, "vertexDataInt64");
    h5::readDataSet(outMem.triangleDataInt64, grpMem, "triangleDataInt64");
    h5::readDataSet(outMem.packedBorderVertices, grpMem, "packedBorderVertices");
}
inline void read(std::vector<OutputStructMembrane>& outMems, const h5::Group& grpSnapshot) {
    using namespace h5;

    Group grpMem = grpSnapshot.getGroup("membranes");
    auto count = h5::readDataSet<std::int64_t>(grpMem, "count");
    outMems.resize(count);
    for(Index i = 0; i < count; ++i) {
        read(outMems[i], grpMem, i);
    }
}

} // namespace medyan

#endif
