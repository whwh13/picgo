#ifndef MEDYAN_MedyanArgs_hpp
#define MEDYAN_MedyanArgs_hpp

#include <string>
#include <thread>

#include "Rand.h"
#include "SysParams.h"
#include "Util/Io/CmdParse.hpp"
#include "Util/Io/Log.hpp"
#include "utility.h" // rdtsc

namespace medyan {

struct MedyanCmdInitResult {

    // Remaining arguments
    int argpNext = 0; // If the parsing is complete, this points to invalid location

    medyan::CommandLineConfig cmdConfig;

};

inline auto medyanInitFromCommandLine(int argc, char** argv) {
    using namespace cmdparse;

    // The initialization result output
    MedyanCmdInitResult res;

    // Parsing the command line
    {
        Command cmdMain("MEDYAN", "");

        cmdMain.addOption(Option('s', "", "file", "System input file", false,
            [&](const Command&, const std::string& arg) {
                res.cmdConfig.inputFile = arg;
                // Also sets the input directory to the directory of this file if not set.
                if(res.cmdConfig.inputDirectory.empty()) {
                    res.cmdConfig.inputDirectory = res.cmdConfig.inputFile.parent_path();
                    // If still empty, set to current directory.
                    if(res.cmdConfig.inputDirectory.empty()) {
                        res.cmdConfig.inputDirectory = ".";
                    }
                }
            }
        ));
        cmdMain.addOption(makeOptionWithVar('i', "", "path", "Input directory", false, res.cmdConfig.inputDirectory));
        cmdMain.addOption(makeOptionWithVar('o', "", "path", "Output directory", false, res.cmdConfig.outputDirectory));
        cmdMain.addOption(Option(0, "seed-fixed", "seed", "Fixed random generator seed", false,
            [&](const Command&, const std::string& arg) {
                res.cmdConfig.rngSeedFixed = true;
                VariableWrite<unsigned long long>{std::string("seed")}(res.cmdConfig.rngSeed, arg);
            }
        ));
        cmdMain.addOption(Option(0, "gui", "Enable GUI", false,
            [&](const Command&) {
                res.cmdConfig.guiEnabled = true;
            }
        ));
        cmdMain.addOption(makeOptionWithVar('t', "", "int", "Thread count (0 for auto)", false, res.cmdConfig.numThreads));
        cmdMain.addOption(Option(0, "trap-fp-invalid", "Enable trapping of floating point invalid operations", false,
            [&](const Command&) {
                res.cmdConfig.trapInvalidFP = true;
            }
        ));
        cmdMain.addHelp();

        // Add gui command
        Command& cmdGui = cmdMain.addCommand(
            "gui", "GUI mode",
            [&] { res.cmdConfig.runMode = medyan::MedyanRunMode::gui; }
        );
        cmdGui.addHelp();

        // Add analyze command
        Command& cmdAnalyze = cmdMain.addCommand(
            "analyze", "Analyze simulation output",
            [&] { res.cmdConfig.runMode = medyan::MedyanRunMode::analyze; }
        );
        cmdAnalyze.addOption(Option(0, "bond-frame", "frame", "Frame of membrane topology information", false,
            [&](const Command&, const std::string& arg) {
                if(arg == "all") {
                    res.cmdConfig.runAnalyzeParams.analyzeMembraneBondAllFrames = true;
                } else {
                    VariableWrite< size_t >{"bond-frame"}(res.cmdConfig.runAnalyzeParams.analyzeMembraneBondFrame, arg);
                }
            }
        ));
        cmdAnalyze.addOption(makeOptionWithVar(0, "frame-interval", "int", "Interval of frames", false, res.cmdConfig.runAnalyzeParams.analyzeFrameInterval));
        cmdAnalyze.addHelp();

        // Add interactive configuration command
        auto& cmdConfig = cmdMain.addCommand(
            "config", "Interactive system configuration, or to normalize an existing config file.",
            [&] { res.cmdConfig.runMode = medyan::MedyanRunMode::config; }
        );

        // Add test command
        auto& cmdTest = cmdMain.addCommand(
            "test", "Run MEDYAN tests.",
            [&] { res.cmdConfig.runMode = medyan::MedyanRunMode::test; }
        );
        cmdTest.setTerminating(true);

        // Add side routines.
        {
            auto& cmdSide = cmdMain.addCommand(
                "side", "Run MEDYAN side procedures.",
                [&] { res.cmdConfig.runMode = medyan::MedyanRunMode::side; }
            );
            cmdSide.addPosArgForVar(
                "name", "Name of side procedure.", true,
                res.cmdConfig.sideProcName
            );
            cmdSide.addHelp();
        }

        // Add validation
        cmdMain.setValidation([&] {
            if(res.cmdConfig.runMode == medyan::MedyanRunMode::simulation) {
                if(res.cmdConfig.inputFile.empty()) {
                    throw ValidationError("Must specify the system input file (-s)");
                }
                if(res.cmdConfig.inputDirectory.empty()) {
                    throw ValidationError("Must specify the input directory (-i)");
                }
                if(res.cmdConfig.outputDirectory.empty()) {
                    throw ValidationError("Must specify the output directory (-o)");
                }
            }
        });

        try {

            cmdMain.ruleCheck();
            res.argpNext = cmdMain.parse(argc, argv);

        } catch (const CommandLogicError& e) {
            std::cerr << e.what() << std::endl;
            // Internal error, no help message generated.
            throw;
        } catch (const ParsingError& e) {
            std::cerr << e.what() << std::endl;
            cmdMain.printUsage();
            throw;
        } catch (const ValidationError& e) {
            std::cerr << e.what() << std::endl;
            cmdMain.printUsage();
            throw;
        }
    }

    //-------------------------------------------------------------------------
    // Post processing
    //-------------------------------------------------------------------------

    // Initialize loggers.
    //---------------------------------
    log::configureLoggers(res.cmdConfig.outputDirectory);

    // Number of threads
    //---------------------------------
    // The actual number of threads will also be affected by the result of
    // std::thread::hardware_concurrency()
    //
    // Let N1 be the user input, N2 be the result of calling hardward_concurrency,
    // and N be the final number of threads. Then
    //
    //   - N1 < 0: (Treat as unspecified) N = 1
    //   - N1 == 0
    //       - N2 == 0: N = 1; Show warning
    //       - N2 > 0: N = N2
    //   - N1 > 0
    //       - N2 == 0: N = N1
    //       - 0 < N2 < N1: N = N1; Show warning
    //       - N2 >= N1: N = N1
    if(res.cmdConfig.numThreads < 0) {
        res.cmdConfig.numThreads = 1;
    } else {
        unsigned hc = std::thread::hardware_concurrency();

        if(res.cmdConfig.numThreads == 0) {
            if(hc == 0) {
                res.cmdConfig.numThreads = 1;
                log::warn("Cannot auto determine the number of threads. Using 1 thread.");
            } else {
                res.cmdConfig.numThreads = hc;
            }
        } else {
            if(hc > 0 && hc < res.cmdConfig.numThreads) {
                log::warn("The number of threads specified is greater than hardware concurrency ({}).", hc);
            }
        }
    }
    log::info("The program will utilize {} thread(s) for computing.", res.cmdConfig.numThreads);

    // Seed global random generator
    //---------------------------------
    if(!res.cmdConfig.rngSeedFixed) {
        res.cmdConfig.rngSeed = rdtsc();
    }
    log::debug("Global RNG seed: {}", res.cmdConfig.rngSeed);
    Rand::eng.seed(res.cmdConfig.rngSeed);


    return res;

}

} // namespace medyan

#endif
