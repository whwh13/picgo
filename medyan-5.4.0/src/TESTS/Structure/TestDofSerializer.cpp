#include <catch2/catch.hpp>

#include "Structure/DofSerializer.hpp"
#include "Structure/SubSystemFunc.hpp"

namespace medyan {

TEST_CASE("DOF serialization", "[DofSerializer]") {
    // Initialize a system to be serialized.
    SubSystem sys;

    const auto addMembraneSheet = [&]() {
        MembraneSetup memSetup;

        // Initialize 6x6 grid.
        const int nx = 7;
        const int ny = 7;
        std::vector<Membrane::CoordinateType> coords;
        coords.reserve(nx * ny);
        for(int j = 0; j < ny; ++j) {
            for(int i = 0; i < nx; ++i) {
                coords.push_back({ 1.0 * i, 1.0 * j, 0.0 });
            }
        }
        std::vector<std::array<int, 3>> triangles;
        triangles.reserve(2 * (nx - 1) * (ny - 1));
        for(int j = 0; j < ny-1; ++j) {
            for(int i = 0; i < nx-1; ++i) {
                triangles.push_back({ j * nx + i, j * nx + i + 1, (j + 1) * nx + i + 1 });
                triangles.push_back({ j * nx + i, (j + 1) * nx + i + 1, (j + 1) * nx + i });
            }
        }

        // Initialize membrane.
        return SubSystemFunc{}.emplaceTrackable<Membrane>(sys, memSetup, coords, triangles);
    };

    const auto clearMembranes = [&]() {
        for(auto it = sys.membranes.begin(); it != sys.membranes.end(); ++it) {
            SubSystemFunc{}.removeTrackable<Membrane>(sys, sys.membranes.indexat(it));
        }
    };

    ScopeGuard clearMembranesGuard{ clearMembranes };


    SECTION("Vertex pinning") {
        // Add a membrane.
        auto memIndex = addMembraneSheet();
        auto& mem = sys.membranes.at(memIndex);

        const auto countPinned = [&]() {
            int count = 0;
            for(auto& v : mem.getMesh().getVertices()) {
                if(v.attr.vertex(sys).pinned) {
                    ++count;
                }
            }
            return count;
        };

        REQUIRE(mem.getMesh().getVertices().size() == 49);

        // Do not pin.
        mem.setup.vertexPinning = MembraneVertexPinning::none;
        updateVertexPinning(sys);
        REQUIRE(countPinned() == 0);

        // Pin border vertices.
        mem.setup.vertexPinning = MembraneVertexPinning::border1;
        updateVertexPinning(sys);
        REQUIRE(countPinned() == 24);

        // Pin border vertices and their neighbors.
        mem.setup.vertexPinning = MembraneVertexPinning::border2;
        updateVertexPinning(sys);
        REQUIRE(countPinned() == 40);
    }
}

} // namespace medyan
