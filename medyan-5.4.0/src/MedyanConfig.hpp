#ifndef MEDYAN_MedyanConfig_hpp
#define MEDYAN_MedyanConfig_hpp

/*
This file includes all one needs for interacting with MEDYAN user inputs.

Features:
  - Simulation config storage (SysParams.h)
  - Simulation config validation (SysParamsValidate.hpp)
  - Reading config from input files (Parser.h)
  - Reading config from command line guide
  - Generating config input files
  - GUI config display (future)
  - GUI config interactive modification (far future)

Documentation:

  Data structures
  ---------------
    medyan::SimulConfig

      All the simulation configuration is stored in an object of this type.

*/

#include <filesystem>
#include <iostream>
#include <stdexcept>

#include "Parser.h"
#include "SysParams.h"
#include "Util/Io/Log.hpp"

namespace medyan {


//---------------------------------
// Auxiliary functions
//---------------------------------

// read the whole file into string
inline std::string readFileToString(std::filesystem::path file) {
    using namespace std;

    ifstream ifs(file);
    if(!ifs.is_open()) {
        log::error("Cannot open file {}", file);
        throw std::runtime_error("Cannot open input file.");
    }
    log::info("Loading file {}", file);

    stringstream ss;
    ss << ifs.rdbuf(); // Read the file
    return ss.str();
}

// Read from stdin (with or without default values)
inline int cinIntDef(int def) {
    using namespace std;
    string line; getline(cin, line);
    if(line.empty()) return def;
    while(true) {
        try {
            return stoi(line);
        }
        catch(const exception& e) {
            cout << "Invalid input: " << line << ". Please correct it: ";
            getline(cin, line);
            if(line.empty()) return def;
        }
    }
}
inline double cinDoubleDef(double def) {
    using namespace std;
    string line; getline(cin, line);
    if(line.empty()) return def;
    while(true) {
        try {
            return stod(line);
        }
        catch(const exception& e) {
            cout << "Invalid input: " << line << ". Please correct it: ";
            getline(cin, line);
            if(line.empty()) return def;
        }
    }
}
inline std::string cinStringDef(std::string def) {
    using namespace std;
    string line; getline(cin, line);
    if(line.empty()) return def;
    return line;
}
inline std::string cinChooseString(std::vector<std::string> choices) {
    using namespace std;
    if(choices.empty()) return "";
    // The first one is default
    auto s = cinStringDef(choices[0]);
    while(find(choices.begin(), choices.end(), s) == choices.end()) {
        cout << "Invalid choice: " << s << ". Please correct the input: ";
        s = cinStringDef(choices[0]);
    }
    return s;
}
inline bool cinYesNo() {
    using namespace std;
    string line;
    getline(cin, line);
    while(line != "y" && line != "n") {
        cout << "Please indicate yes or no using (y/n): ";
        getline(cin, line);
    }
    return line == "y";
}
inline bool cinYesNoDef(bool def) {
    using namespace std;
    string line;
    getline(cin, line);
    if(line.empty()) return def;
    while(line != "y" && line != "n") {
        cout << "Please indicate yes or no using (y/n): ";
        getline(cin, line);
        if(line.empty()) return def;
    }
    return line == "y";
}


// Read bubble input
inline void readBubbleConfig(SimulConfig& sc, std::istream& is) {
    sc.bubbleData = BubbleParser::readBubbles(is);
}

// Read filament input
inline void readFilamentConfig(SimulConfig& sc, std::istream& is) {
    sc.filamentData = FilamentParser::readFilaments(is);
}

enum class InputGenOverwriteAction {
    overwrite, // warn and overwrite
    ignore,    // warn and ignore
    confirm    // confirm with stdin input. should not be used in non-interactive mode
};

// Read simulation configuration from file
inline SimulConfig getSimulConfigFromInput(
    const std::filesystem::path& systemInputFile,
    const std::filesystem::path& inputDirectory
) {
    using namespace std;
    using namespace std::filesystem;

    // Auxiliary struct handling file streams
    struct ReadFile {
        ifstream ifs;
        ReadFile(const path& p) : ifs(p) {
            if(!ifs.is_open()) {
                log::error("There was an error parsing file {}", p);
                throw std::runtime_error("Cannot open input file.");
            }
            log::info("Loading input file {}", p);
        }
    };

    SimulConfig conf;
    conf.metaParams.systemInputFile = systemInputFile;
    conf.metaParams.inputDirectory  = inputDirectory;

    // Read system input
    parseKeyValueList(
        conf,
        lexTokenList(sExprTokenize(readFileToString(systemInputFile))),
        systemParser(),
        KeyValueParserUnknownKeyAction::warn
    );

    // Read chemistry input
    if(conf.chemParams.chemistrySetup.inputFile.empty()) {
        log::error("Need to specify a chemical input file. Exiting.");
        throw std::runtime_error("No chemistry input file specified.");
    }
    else {
        parseKeyValueList(
            conf,
            lexTokenList(sExprTokenize(readFileToString(inputDirectory / conf.chemParams.chemistrySetup.inputFile))),
            chemDataParser(),
            KeyValueParserUnknownKeyAction::warn
        );
    }

    // Read bubble input
    if(!conf.bubbleSetup.inputFile.empty()) {
        ReadFile file(inputDirectory / conf.bubbleSetup.inputFile);
        readBubbleConfig(conf, file.ifs);
    }

    // Read filament input
    if(!conf.filamentSetup.inputFile.empty()) {
        ReadFile file(inputDirectory / conf.filamentSetup.inputFile);
        readFilamentConfig(conf, file.ifs);
    }

    // Do post-processing.
    postprocess(conf);

    return conf;
}

// Output configuration to files
//
// altOstream is used when the file name is not specified.
inline void generateInput(
    const SimulConfig&      conf,
    std::ostream&           altOstream = std::cout,
    InputGenOverwriteAction overwriteAction = InputGenOverwriteAction::overwrite
) {
    using namespace std;

    // Generate system input
    {
        const auto& p = conf.metaParams.systemInputFile;
        bool useAlt = false;
        if(p.empty()) {
            log::info("The system input file is not specified in the config.");
            useAlt = true;
        }
        else if(filesystem::exists(p)) {
            switch(overwriteAction) {

            case InputGenOverwriteAction::overwrite:
                log::warn("The file {} already exists and will be overwritten.", p.string());
                break;

            case InputGenOverwriteAction::ignore:
                log::warn("The file {} already exists.", p.string());
                useAlt = true;
                break;

            case InputGenOverwriteAction::confirm:
                cout << "The file " << p.string() << " already exists. Overwrite (y/n)? ";
                useAlt = !cinYesNo();
                break;
            }
        }

        if(useAlt) {
            log::info("Using alternative output target.");
            outputTokenList(altOstream, buildTokens(conf, systemParser()));
            log::info("----- End of system input -----");
        } else {
            ofstream ofs(p);
            outputTokenList(ofs, buildTokens(conf, systemParser()));
            log::info("System input file written to {}.", p.string());
        }
    }

    // Generate chemistry input
    {
        const auto& name = conf.chemParams.chemistrySetup.inputFile;
        const auto p = conf.metaParams.inputDirectory / name;
        bool useAlt = false;
        if(name.empty()) {
            log::info("The chemistry input file is not specified in the config.");
            useAlt = true;
        }
        else {
            if(filesystem::exists(p)) {
                switch(overwriteAction) {

                case InputGenOverwriteAction::overwrite:
                    log::warn("The file {} already exists and will be overwritten.", p.string());
                    break;

                case InputGenOverwriteAction::ignore:
                    log::warn("The file {} already exists.", p.string());
                    useAlt = true;
                    break;

                case InputGenOverwriteAction::confirm:
                    cout << "The file " << p.string() << " already exists. Overwrite (y/n)? ";
                    useAlt = !cinYesNo();
                    break;
                }
            }
        }

        if(useAlt) {
            log::info("Using alternative output target.");
            outputTokenList(altOstream, buildTokens(conf, chemDataParser()));
            log::info("----- End of chemistry input -----");
        } else {
            ofstream ofs(p);
            outputTokenList(ofs, buildTokens(conf, chemDataParser()));
            log::info("Chemistry input file written to {}.", p.string());
        }
    }

    // Generate bubble input
    if(!conf.bubbleSetup.inputFile.empty()) {
        log::warn("Additional bubble input is specified, but an input file is not generated.");
        log::info("It will be supported in future versions.");
    }

    // Generate filament input
    if(!conf.filamentSetup.inputFile.empty()) {
        log::warn("Additional filament input is specified, but an input file is not generate.");
        log::info("It will be supported in future versions.");
    }


}


// Starts an interactive session of the command line configuration generation
//
// altOstream is used when the file name is not specified.
inline void interactiveConfig(std::ostream& altOstream = std::cout) {
    using namespace std;


    // States
    //-------------------------------------------------------------------------
    bool shouldAddMyosin = true;
    bool shouldAddLinker = true;
    bool shouldAddBrancher = true;

    filesystem::path fileDirectory;

    // Initialize a default configuration
    //-------------------------------------------------------------------------
    SimulConfig conf;

    conf.geoParams.compartmentSizeX = 500.0;
    conf.geoParams.compartmentSizeY = 500.0;
    conf.geoParams.compartmentSizeZ = 500.0;
    conf.geoParams.monomerSize = {2.7};
    conf.geoParams.cylinderSize = {108.0};

    conf.chemParams.motorNumHeadsMin = {15};
    conf.chemParams.motorNumHeadsMax = {30};
    conf.chemParams.motorStepSize = {6.0};
    conf.chemParams.numBindingSites = {4};
    conf.chemParams.numFilaments = 1;
    conf.chemParams.chemistryAlgorithm.algorithm = "NRM";

    conf.mechParams.mechanicsFFType.BoundaryFFType = "REPULSIONEXP";
    conf.boundParams.BoundaryCutoff = 300.0;
    conf.boundParams.BoundaryK = 41.0;
    conf.boundParams.BScreenLength = 2.7;

    conf.mechParams.mechanicsFFType.FStretchingType = "HARMONIC";
    conf.mechParams.FStretchingK = {100.0};
    conf.mechParams.mechanicsFFType.FBendingType = "COSINE";
    conf.mechParams.FBendingK = {672.0};
    conf.mechParams.FBendingTheta = {0.0};

    conf.mechParams.mechanicsAlgorithm.ConjugateGradient = "POLAKRIBIERE";
    conf.mechParams.mechanicsAlgorithm.gradientTolerance = 5.0;
    conf.mechParams.mechanicsAlgorithm.maxDistance = 1.0;
    conf.mechParams.mechanicsAlgorithm.lambdaMax = 0.01;

    conf.filamentSetup.projectionType = "STRAIGHT";

    conf.chemistryData.speciesFilament[0].push_back("AF");
    conf.chemistryData.speciesPlusEnd[0].push_back("PA");
    conf.chemistryData.speciesMinusEnd[0].push_back("MA");

    const auto addMyosin = [&] {
        conf.mechParams.mechanicsFFType.MStretchingType = "HARMONIC";
        conf.mechParams.MStretchingK = {2.5};

        conf.chemistryData.speciesMotor[0].push_back("MOA");
        conf.chemistryData.speciesBound[0].push_back("MEA");
        conf.chemistryData.M_BINDING_INDEX[0] = "MEA";

        conf.chemistryData.motorReactions.push_back({
            { "MEA:BOUND:1", "MEA:BOUND:2", "MD:DIFFUSING" },
            { "MOA:MOTOR:1", "MOA:MOTOR:2" },
            0.2, 1.7, 175.0, 225.0,
            0.0, "NA"
        });
        conf.chemistryData.motorWalkingReactions[0].push_back({
            "MOA:MOTOR:N", "MEA:BOUND:N+1",
            "MOA:MOTOR:N+1", "MEA:BOUND:N",
            0.2,
            0.0, "NA"
        });
    };
    const auto addLinker = [&] {
        conf.mechParams.mechanicsFFType.LStretchingType = "HARMONIC";
        conf.mechParams.LStretchingK = {8.0};

        conf.chemistryData.speciesLinker[0].push_back("LA");
        conf.chemistryData.speciesBound[0].push_back("LEA");
        conf.chemistryData.L_BINDING_INDEX[0] = "LEA";

        conf.chemistryData.linkerReactions.push_back({
            { "LEA:BOUND:1", "LEA:BOUND:2", "LD:DIFFUSING" },
            { "LA:LINKER:1", "LA:LINKER:2" },
            0.01, 0.3, 30.0, 40.0,
            0.0, "NA"
        });
    };
    const auto addBrancher = [&] {
        conf.mechParams.mechanicsFFType.BrStretchingType = "HARMONIC";
        conf.mechParams.BrStretchingK = {100.0};
        conf.mechParams.BrStretchingL = {6.0};
        conf.mechParams.mechanicsFFType.BrBendingType = "COSINE";
        conf.mechParams.BrBendingK = {10.0};
        conf.mechParams.BrBendingTheta = {1.22};
        conf.mechParams.mechanicsFFType.BrDihedralType = "COSINE";
        conf.mechParams.BrDihedralK = {10.0};
        conf.mechParams.mechanicsFFType.BrPositionType = "COSINE";
        conf.mechParams.BrPositionK = {20.0};

        conf.chemistryData.speciesBrancher[0].push_back("BA");
        conf.chemistryData.speciesBound[0].push_back("BEA");
        conf.chemistryData.B_BINDING_INDEX[0] = "BEA";

        conf.chemistryData.branchingReactions[0].push_back({
            { "BD:DIFFUSING", "A:DIFFUSING", "BEA:BOUND" },
            { "BA:BRANCHER", "PA:PLUSEND" },
            0.001, 0.01, "ALL", 100.0
        });
    };

    // The interactive session
    //-------------------------------------------------------------------------

    log::info("-----------------------------------------------------");
    log::info("You are now in the interactive configuration session.");
    log::info("Note that this only provides the very basic configuration of the system."
        " For more detailed configuration, read the documents and manually edit the input files.");

    log::info("##################################################");
    log::info("  Units in MEDYAN are nm, second, pN, and pN*nm");
    log::info("##################################################");

    //---------- basics ----------
    cout << endl;
    log::info("Basic simulation parameters");
    log::info("The geometry information");
    log::info("Note: Network size = compartment size (500nm by default) .* (NX, NY, NZ)");
    cout << "Num compartments in x direction (default 2): "; conf.geoParams.NX = cinIntDef(2);
    cout << "Num compartments in y direction (default 2): "; conf.geoParams.NY = cinIntDef(2);
    cout << "Num compartments in z direction (default 2): "; conf.geoParams.NZ = cinIntDef(2);

    cout << "Boundary shape {CUBIC (default), SPHERICAL, CYLINDER}: ";
    conf.boundParams.boundaryType.boundaryShape = cinChooseString({ "CUBIC", "SPHERICAL", "CYLINDER" });
    if(
        const auto& s = conf.boundParams.boundaryType.boundaryShape;
        s == "SPHERICAL" || s == "CYLINDER"
    ) {
        cout << "Diameter of boundary in nm (default 1000): "; conf.boundParams.diameter = cinDoubleDef(1000.0);
    }

    cout << endl;
    log::info("Run time specification");
    cout << "Total running time in second (default 1000): "; conf.chemParams.chemistryAlgorithm.runTime = cinDoubleDef(1000.0);
    cout << "Snapshot output time interval in second (default 1.0): "; conf.chemParams.chemistryAlgorithm.snapshotTime = cinDoubleDef(1.0);

    cout << endl;
    log::info("Initial filaments (applicable to type 0 filament only)");
    cout << "Number of filaments (default 30): "; conf.filamentSetup.numFilaments = cinIntDef(30);
    log::info("The cylinder length is {} nm.", conf.geoParams.cylinderSize[0]);
    cout << "Number of cylinders of each filament (default 1): "; conf.filamentSetup.filamentLength = cinIntDef(1);

    //---------- chemistry ----------
    cout << endl;
    log::info("Basic chemistry parameters");
    cout << "Initial diffusing actin copy number (default 5000): ";
    conf.chemistryData.speciesDiffusing.push_back({
        "A", cinIntDef(5000), 80.0, 0.0, 0.0, "REG",
        0, "NONE", 0.0
    });
    cout << "Enable filament polymerization/depolymerization (y/n) (default y): ";
    if(cinYesNoDef(true)) {
        conf.chemistryData.polymerizationReactions[0].push_back({
            "A:DIFFUSING", "PA:PLUSEND",
            "AF:FILAMENT", "PA:PLUSEND",
            0.154, 0.0, "NA"
        });
        conf.chemistryData.polymerizationReactions[0].push_back({
            "A:DIFFUSING", "MA:MINUSEND",
            "AF:FILAMENT", "MA:MINUSEND",
            0.017, 0.0, "NA"
        });
        conf.chemistryData.depolymerizationReactions[0].push_back({
            "AF:FILAMENT", "PA:PLUSEND",
            "A:DIFFUSING", "PA:PLUSEND",
            1.4, 0.0, "NA"
        });
        conf.chemistryData.depolymerizationReactions[0].push_back({
            "AF:FILAMENT", "MA:MINUSEND",
            "A:DIFFUSING", "MA:MINUSEND",
            0.8, 0.0, "NA"
        });
    }

    cout << endl;
    log::info("Myosin - Non-muscle mysoin IIA");
    log::info("  Includes force fields, diffusing/bound species, and binding/unbinding/walking reactions.");
    cout << "Add myosin (y/n) (default y): "; shouldAddMyosin = cinYesNoDef(true);
    if(shouldAddMyosin) {
        cout << "Initial diffusing myosin copy number (default 50): ";
        conf.chemistryData.speciesDiffusing.push_back({
            "MD", cinIntDef(50), 0.8, 0.0, 0.0, "REG",
            0, "NONE", 0.0
        });
    }

    log::info("Linker - alpha-actinin crosslinker");
    log::info("  Includes force fields, diffusing/bound species, and binding/unbinding reactions.");
    cout << "Add linker (y/n) (default y): "; shouldAddLinker = cinYesNoDef(true);
    if(shouldAddLinker) {
        cout << "Initial diffusing linker copy number (default 500): ";
        conf.chemistryData.speciesDiffusing.push_back({
            "LD", cinIntDef(500), 8.0, 0.0, 0.0, "REG",
            0, "NONE", 0.0
        });
    }

    log::info("Bracher - Arp2/3 brancher");
    log::info("  Includes force fields, diffusing/bound species, and binding/unbinding reactions.");
    cout << "Add brancher (y/n) (default y): "; shouldAddBrancher = cinYesNoDef(true);
    if(shouldAddBrancher) {
        cout << "Initial diffusing brancher copy number (default 20): ";
        conf.chemistryData.speciesDiffusing.push_back({
            "BD", cinIntDef(20), 80.0, 1.0, 0.0, "REG",
            0, "NONE", 0.0
        });
    }

    //---------- mechanics ----------
    cout << endl;
    log::info("Mechanics parameters");
    log::info("Energy minimization time interval:");
    log::info("  Use a lower value if: 1. simulation fails or generates warnings");
    log::info("                        2. has very fast chemical reactions");
    log::info("  Recommend value: 0.001 - 0.05");
    cout << "Energy minimization interval in second (default 0.01): ";
    conf.chemParams.chemistryAlgorithm.neighborListTime =
        conf.chemParams.chemistryAlgorithm.minimizationTime = cinDoubleDef(0.01);

    cout << endl;
    log::info("Actin force fields (stretching and bending are already enabled)");
    cout << "Enable volume exclusion (y/n) (default y): ";
    if(cinYesNoDef(true)) {
        conf.mechParams.mechanicsFFType.cylinderVolumeExclusionFFType = CylinderVolumeExclusionFFType::integral;
        conf.mechParams.VolumeCutoff = 108.0;
        conf.mechParams.VolumeK = {1e5};
    }

    //---------- dynamic rates ----------
    cout << endl;
    log::info("Dynamic rate parameters");
    cout << "Enable Brownian Ratchet (y/n) (default y): ";
    if(cinYesNoDef(true)) {
        conf.dyRateParams.dynamicRateType.dFPolymerizationType = { "BROWRATCHET" };
        conf.dyRateParams.dFilPolymerizationCharLength = {2.7};
    }

    if(shouldAddMyosin) {
        cout << "Enable motor catch bond (y/n) (default y): ";
        if(cinYesNoDef(true)) {
            conf.dyRateParams.dynamicRateType.dMUnbindingType = { "LOWDUTYCATCH" };
            conf.dyRateParams.dMotorUnbindingCharForce = {12.62};
        }
        cout << "Enable motor walking stall (y/n) (default y): ";
        if(cinYesNoDef(true)) {
            conf.dyRateParams.dynamicRateType.dMWalkingType = { "LOWDUTYSTALL" };
            conf.dyRateParams.dMotorWalkingCharForce = {90.0};
        }
    }
    if(shouldAddLinker) {
        cout << "Enable linker slip bond (y/n) (default y): ";
        if(cinYesNoDef(true)) {
            conf.dyRateParams.dynamicRateType.dLUnbindingType = { "SLIP" };
            conf.dyRateParams.dLinkerUnbindingCharLength = {0.24};
        }
    }

    //---------- post processing ----------
    if(shouldAddMyosin) addMyosin();
    if(shouldAddLinker) addLinker();
    if(shouldAddBrancher) addBrancher();

    cout << endl;
    log::info("Output files settings");
    cout << "The directory for generated files: "; fileDirectory = cinStringDef("");
    while(fileDirectory.empty()) {
        fileDirectory = filesystem::current_path();
        cout << "Will use " << fileDirectory.string() << " to generate files. Proceed (y/n)? ";
        if(!cinYesNo()) {
            cout << "The directory for generated files: "; fileDirectory = cinStringDef("");
        }
    }
    conf.metaParams.inputDirectory = fileDirectory;
    cout << "The name of the system input file (default system.txt): ";
    conf.metaParams.systemInputFile = fileDirectory / cinStringDef("system.txt");
    cout << "The name of the chemistry input file (default chemistry.txt): ";
    conf.chemParams.chemistrySetup.inputFile = cinStringDef("chemistry.txt");

    // File generation
    //-------------------------------------------------------------------------
    cout << endl;
    log::info("Configuration complete! Generating files...");
    generateInput(
        conf,
        altOstream,
        InputGenOverwriteAction::confirm
    );

    cout << endl;
    log::info("Input file generation complete.");
    log::info("Thanks for using MEDYAN interactive configuration!");

}

// Given input files, generate the normalized input files in the output directory.
inline void normalizeConfig(
    const std::filesystem::path& systemInputFile,
    const std::filesystem::path& inputDirectory,
    const std::filesystem::path& outputDirectory
) {
    auto simulConfig = getSimulConfigFromInput(systemInputFile, inputDirectory);

    auto sysFileName = systemInputFile.filename();
    auto chemFileName = simulConfig.chemParams.chemistrySetup.inputFile.filename();

    // Set new sys/chem file names.
    auto newSysFileName = std::filesystem::path("normalized." + sysFileName.string());
    auto newChemFileName = std::filesystem::path("normalized." + chemFileName.string());

    // Update simulConfig.
    simulConfig.metaParams.systemInputFile = outputDirectory / newSysFileName;
    simulConfig.metaParams.inputDirectory = outputDirectory;
    simulConfig.chemParams.chemistrySetup.inputFile = newChemFileName;

    // Generate new files.
    generateInput(simulConfig);

    log::info("Normalized files generated in {}", outputDirectory);
}

} // namespace medyan

#endif
