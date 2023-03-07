#ifndef MEDYAN_Structure_SurfaceMesh_AdaptiveMeshGeometryManager_hpp
#define MEDYAN_Structure_SurfaceMesh_AdaptiveMeshGeometryManager_hpp

namespace medyan::adaptive_mesh {

template< typename Mesh > struct GeometryManager {
    template< typename Context >
    static void computeAllTriangleNormals(Context& sys, Mesh& mesh) {
        const size_t numTriangles = mesh.numTriangles();

        for(typename Mesh::TriangleIndex ti { 0 }; ti < numTriangles; ++ti) {
            medyan::adaptiveComputeTriangleNormal(sys, mesh, ti);
        }
    }

    template< typename Context >
    static void computeAllAngles(Context& sys, Mesh& mesh) {
        const size_t numHalfEdges = mesh.numHalfEdges();

        for(typename Mesh::HalfEdgeIndex hei { 0 }; hei < numHalfEdges; ++hei) {
            medyan::adaptiveComputeAngle(sys, mesh, hei);
        }
    }

    // Requires
    //   - Unit normals in triangles (geometric)
    //   - Angles in halfedges (geometric)
    static void computeAllVertexNormals(Mesh& mesh) {
        const size_t numVertices = mesh.numVertices();

        for(typename Mesh::VertexIndex vi {0}; vi < numVertices; ++vi) {
            medyan::adaptiveComputeVertexNormal(mesh, vi);
        }
    }
}; // End GeometryManager

} // namespace medyan::adaptive_mesh

#endif
