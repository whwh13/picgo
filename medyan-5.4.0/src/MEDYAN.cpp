
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

/*! \mainpage MEDYAN software package
 
 \section intro_sec Introduction
 
The cell cytoskeleton plays a key role in human biology and disease, contributing ubiquitously
 to such important processes as embryonic development, wound repair and cancer 
 metastasis. Papoian laboratory is interested in gaining deeper understanding of the
 physical chemistry behind these complex, far-from-equilibrium mechanochemical 
 processes. Their approach and model, named Mechanochemical Dynamics of Active Networks
 (**MEDYAN**), is based on combining stochastic reaction-diffusion treatment
 of cellular biochemical processes with polymer physics of cytoskeletal filament network 
 growth, while explicitly coupling chemistry and mechanics.
 
 Papoian laboratory has developed a third-generation software package based on
 the **MEDYAN** model, to simulate growth dynamics of actin based filamentous networks *in vitro* 
 and *in vivo*. Recent papers where **MEDYAN** or its second-generation predecessor, **StochTools**,
 were used can be found on the publication section of [the Papoian group's main web page]
 (http://papoian.chem.umd.edu/ ) or on [the MEDYAN web page] (http://www.medyan.org ).
 The **MEDYAN** package can also be extended to simulate the dynamics of any active
 matter network.
 
 \section install_sec Installation
 
    See docs/manual/installation.md for details.
 
 */

#include <thread>

#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"

#include "common.h"
#include "MedyanArgs.hpp"
#include "MedyanConfig.hpp"
#include "MedyanMeta.hpp"
#include "Analysis/Io/ReadSnapshot.hpp"
#include "Controller/Controller.h"
#include "Side/SideProcedures.hpp"
#include "Visualization.hpp"


int main(int argc, char **argv) {
    using namespace medyan;

    std::cout.precision(8);
    printHeader();

    auto cmdRes = medyanInitFromCommandLine(argc, argv);
    auto& cmdConfig = cmdRes.cmdConfig;

    log::debug("MEDYAN program started.");

    int returnCode = 0;

    switch(cmdConfig.runMode) {

        case MedyanRunMode::simulation:
            {

                const auto runController = [&]{
                    Controller c;
                    auto conf = c.initialize(cmdConfig);
                    c.run(cmdConfig, conf);
                };

                if(cmdConfig.guiEnabled) {
                    // Run controller on another thread.
                    std::thread simul(runController);

                    guiRunRealtime();
                    simul.join();
                }
                else {
                    // Run controller in this main thread.
                    runController();
                }
            }
            break;

        case MedyanRunMode::gui:
            guiRunTrajectory();
            break;

        case MedyanRunMode::analyze:
            {
                auto inputFilePath = cmdConfig.inputDirectory / "snapshot.traj";
                auto pdbFilePath = cmdConfig.outputDirectory / "snapshot.pdb";
                analysis::SnapshotReader sr(inputFilePath, pdbFilePath, cmdConfig.outputDirectory, "snapshot");
                sr.readAndConvertToVmd(0, cmdConfig.runAnalyzeParams);
            }
            break;

        case MedyanRunMode::config:
            if(cmdConfig.inputFile.empty()) {
                interactiveConfig();
            } else {
                normalizeConfig(cmdConfig.inputFile, cmdConfig.inputDirectory, cmdConfig.outputDirectory);
            }
            break;

        case MedyanRunMode::test:
            returnCode = Catch::Session().run(argc - (cmdRes.argpNext - 1), argv + (cmdRes.argpNext - 1));
            break;

        case MedyanRunMode::side:
            medyan::side::runSideProcedure(cmdConfig.sideProcName);
            break;
    }

    log::debug("MEDYAN program finished.");

    return returnCode;
}