#ifndef MEDYAN_Side_SurfaceMeshDiffusion_hpp
#define MEDYAN_Side_SurfaceMeshDiffusion_hpp

#include <fstream>

#include "Chemistry/ChemNRMImpl.h"
#include "Controller/GController.h"
#include "Structure/SubSystem.h"
#include "Structure/SubSystemFunc.hpp"
#include "Structure/SurfaceMesh/AdaptiveMesh.hpp"
#include "Structure/SurfaceMesh/SurfaceMeshGeneratorPreset.hpp"
#include "SysParams.h"

namespace medyan::side {

inline void surfaceMeshDiffusion() {
    log::info("Running side procedure: surface mesh diffusion test.");


    // Background and other stuff
    //---------------------------------
    log::info("Initializing.");
    SubSystem subSystem;

    GeoParams geoParams;
    geoParams.compartmentSizeX = 500;
    geoParams.compartmentSizeY = 500;
    geoParams.compartmentSizeZ = 500;
    geoParams.NX = 1;
    geoParams.NY = 1;
    geoParams.NZ = 1;

    ChemistryData chemData;
    chemData.speciesMembraneDiffusing = {
        { "mem-diffu-test-a", 100.0 },
        { "mem-diffu-test-b", 0.5 },
        { "diffu-potential-test", 100.0 },
        { "adsorption-potential-test", 0.0 },
    };
    chemData.speciesBulk = {
        {
            "adsorption-test-bulk",
            20000,
            0.0,                    // Release time.
            0.0,                    // Removal time.
            "REG",
        },
    };
    chemData.reactionsSurface = {
        {
            { "mem-diffu-test-a", "mem-diffu-test-b", "mem-diffu-test-b" },
            { "mem-diffu-test-b", "mem-diffu-test-b", "mem-diffu-test-b" },
            1.0,
        },
    };
    chemData.reactionsAdsorptionDesorption = {
        { "adsorption-test-bulk", "adsorption-potential-test", 100, 0.1 },
    };

    GController gController(&subSystem);
    gController.initializeGrid(geoParams);
    auto& grid = *subSystem.getCompartmentGrid();

    // Initialize membrane with surface proteins
    //---------------------------------
    MembraneSetup memSetup;
    MembraneInit  memInit;
    memInit.meshParams = {
        "PLANE",
        "0", "0", "50",         // a point on the plane
        "0", "0", "1",          // facing z direction
        "10", "10", "0",        // box origin
        "480", "480", "100"     // box size
    };

    const auto meshData = mesh_gen::generateMeshViaParams< floatingpoint >(memInit.meshParams);
    auto flatMembraneIndex = SubSystemFunc{}.emplaceTrackable<Membrane>(
        subSystem,
        memSetup,
        meshData.vertexCoordinateList,
        meshData.triangleList
    );
    auto& flatMembrane = subSystem.membranes[flatMembraneIndex];

    ScopeGuard guard { [&] {
        clearChemistry(subSystem, flatMembrane.getMesh());
        SubSystemFunc{}.removeTrackable<Membrane>(subSystem, flatMembraneIndex);
    } };


    adaptive_mesh::MembraneMeshAdapter meshAdapter(MeshAdapterSettings {});
    meshAdapter.adapt(subSystem, flatMembrane.getMesh());

    updateGeometryValueForSystem(subSystem, flatMembrane.getMesh());
    // The geometry of the membrane should be fixed beyond this point.

    log::info("Adding surface chemistry...");
    {
        const auto memChemInfo = MembraneMeshChemistryInfo::fromChemistryData(chemData);
        setChemistry(subSystem, flatMembrane.getMesh(), memChemInfo);
    }

    log::info("Adding adsorption/desorption chemistry...");
    {
        for(auto& b : chemData.speciesBulk) {
            grid.addSpeciesBulk(
                b.name,
                b.initialCopyNumber,
                max_ulim,
                SpeciesType::BULK,
                b.rspeciesType == "REG" ? RSpeciesType::REG : RSpeciesType::CONST
            );
        }

        setAdsorptionDesorptionReactions(subSystem, chemData.reactionsAdsorptionDesorption);
    }


    // Setup vertex energies, by abusing surface mesh curvature mismatch energies.
    // "test-a" and "test-b" do not use curvature mismatch energy, while "diffu-potential-test" does use the 0th curvature mismatch setup.
    subSystem.indexMapMembraneDiffusingSpeciesToCurvatureMismatchSetup = { -1, -1, 0, 0 };
    subSystem.indexMapDesorptionReactionToCurvatureMismatchSetup = { 0 };
    subSystem.allVertexCurvatureMismatchParams.resize(1, flatMembrane.getMesh().numVertices());
    for(Index vi = 0; vi < flatMembrane.getMesh().numVertices(); ++vi) {
        auto& v = flatMembrane.getMesh().attribute(Membrane::MeshType::VertexIndex{vi}).vertex(subSystem);
        v.loopIndex = vi;

        subSystem.allVertexCurvatureMismatchParams(0, vi).energy = -5 * std::sin(2 * M_PI / 1000 * v.coord[0]);
    }

    // Update and fix reaction rates
    setReactionRates(subSystem, flatMembrane.getMesh());
    setAdsorptionDesorptionReactionRates(subSystem, chemData);

    // Set initial species. 10000 copies of species 0 and 2 on vertex 0.
    flatMembrane.getMesh().getVertices()[0].attr.vertex(subSystem).cVertex.species.findSpeciesByIndex(0)->setN(10000);
    flatMembrane.getMesh().getVertices()[0].attr.vertex(subSystem).cVertex.species.findSpeciesByIndex(2)->setN(10000);


    // Initialize chem sim
    //---------------------------------
    ChemNRMImpl nrm;

    // Activate all reactions
    medyan::forEachReactionInMesh(
        subSystem,
        flatMembrane.getMesh(),
        [&](ReactionDy& r) {
            nrm.addReaction(&r);
            r.activateReaction();
        }
    );

    // Setup output
    //---------------------------------
    std::ofstream outChem("side-surface-mesh-diffusion.txt");

    const auto printChemistry = [&]() {
        outChem << "FRAME " << nrm.getTime() << '\n';
        for(Membrane::MeshType::VertexIndex vi {0}; vi < flatMembrane.getMesh().numVertices(); ++vi) {
            auto& vertex = flatMembrane.getMesh().attribute(vi).vertex(subSystem);
            outChem << vi.index << ' '
                << flatMembrane.getMesh().attribute(vi).gVertex.astar / 3 << ' '
                << vertex.coord[0] << ' '
                << vertex.coord[1] << ' '
                << vertex.coord[2] << ' ';
            for(auto& sp : vertex.cVertex.species.species()) {
                outChem << sp->getN() << ' ';
            }
            outChem << '\n';
        }

        outChem << std::flush;
    };


    std::ofstream outTriangle("side-surface-mesh-diffusion-triangle.txt");

    const auto printTriangle = [&]() {
        auto triangleInfo = flatMembrane.getMesh().extract< Membrane::MeshType::VertexTriangleInitializer >(subSystem).triangleVertexIndexList;
        outTriangle << triangleInfo.size() << '\n';
        for(auto& v : triangleInfo) {
            outTriangle << v[0] << ' ' << v[1] << ' ' << v[2] << ' ';
        }

        outTriangle << std::flush;
    };



    // Start simulation
    //---------------------------------
    log::info("Starting simulation.");

    printTriangle();

    printChemistry();
    for (int i = 0; i < 100; ++i) {
        nrm.run(100);
        printChemistry();
        log::info("Current time: {}", nrm.getTime());
    }

}

} // namespace medyan::side

#endif
