
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

#include <stdexcept> // runtime_error
#include <utility> // move

#include "Parser.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

#include "SysParams.h"
#include "Util/Io/Log.hpp"
#include "Util/Parser/BubbleParser.hpp"
#include "Util/Parser/MembraneParser.hpp"
#include "Util/Parser/OutputParser.hpp"

namespace medyan {

KeyValueParser<SimulConfig> buildSystemParser() {
    using namespace std;

    KeyValueParser<SimulConfig> sysParser;

    //--------------------------------------------------------------------------
    // Headers.
    //--------------------------------------------------------------------------
    {
        sysParser.addComment("==================================================");
        sysParser.addComment(" Important notes:");
        sysParser.addComment(" 1. Units in MEDYAN are nm, second, pN, and pN*nm");
        sysParser.addComment("==================================================");
        sysParser.addEmptyLine();
    }

    //--------------------------------------------------------------------------
    // System element definition.
    //--------------------------------------------------------------------------
    {
        sysParser.addComment("==================================================");
        sysParser.addComment(" System elements.");
        sysParser.addComment("==================================================");
        sysParser.addEmptyLine();

        sysParser.addComment(" Membrane setup.");
        sysParser.addArgs(
            "membrane",
            [] (SimulConfig& sc, const SExpr::ListType& sel) {
                if(sel.size() < 2) {
                    log::error("Membrane setup specification is incorrect.");
                    log::info("Membrane setup usage: (membrane <setup-name> [<item1> ...]");
                    throw runtime_error("Membrane setup specification is incorrect.");
                }
                else {
                    auto setupName = get<SExpr::StringType>(sel[1].data);
                    auto& setupVec = sc.membraneSettings.setupVec;
                    // Add this new setup if the name is not already in the list.
                    if(find_if(setupVec.begin(), setupVec.end(), [&](const MembraneSetup& eachSetup) { return eachSetup.name == setupName; }) == setupVec.end()) {
                        MembraneSetup ms;
                        ms.name = setupName;
                        for(int i = 2; i < sel.size(); ++i) {
                            parseKeyValue<MembraneSetup>(ms, sel[i], membraneSetupParser().dict, KeyValueParserUnknownKeyAction::error);
                        }
                        sc.membraneSettings.setupVec.push_back(move(ms));
                    }
                    else {
                        log::error("Membrane setup {} already exists.", setupName);
                        throw runtime_error("Membrane setup already exists.");
                    }
                }
            },
            [] (const SimulConfig& sc) {
                vector<list<SExprToken>> res;

                for(auto& eachSetup : sc.membraneSettings.setupVec) {
                    list<SExprToken> tokenList;
                    // Add setup name.
                    tokenList.push_back(SExprToken::makeString(eachSetup.name));
                    // All all setup items.
                    tokenList.splice(tokenList.end(), buildTokens<MembraneSetup>(eachSetup, membraneSetupParser()));
                    res.push_back(move(tokenList));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();
    }

    //--------------------------------------------------------------------------
    // System geometries.
    //--------------------------------------------------------------------------
    {
        sysParser.addComment("==================================================");
        sysParser.addComment(" Geometric parameters");
        sysParser.addComment("==================================================");
        sysParser.addEmptyLine();

        sysParser.addComment(" Set network sizes and shape");
        sysParser.addComment(" - Set the number of compartments in x, y and z directions");
        sysParser.addComment(" - Network size = compartment size (500nm by default) * (NX, NY, NZ)");

        sysParser.addStringArgsWithAliases(
            "NX", { "NX:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "There was an error parsing input file at grid dimensions. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if(lineVector.size() == 2)
                    sc.geoParams.NX = stoi(lineVector[1]);
                else {}
            },
            [] (const SimulConfig& sc) {
                vector< string > res;
                if(sc.geoParams.NX) res.push_back(toString(sc.geoParams.NX));
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "NY", { "NY:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "There was an error parsing input file at grid dimensions. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if(lineVector.size() == 2)
                    sc.geoParams.NY = stoi(lineVector[1]);
                else {}
            },
            [] (const SimulConfig& sc) {
                vector< string > res;
                if(sc.geoParams.NY) res.push_back(toString(sc.geoParams.NY));
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "NZ", { "NZ:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "There was an error parsing input file at grid dimensions. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if(lineVector.size() == 2)
                    sc.geoParams.NZ = stoi(lineVector[1]);
                else {}
            },
            [] (const SimulConfig& sc) {
                vector< string > res;
                if(sc.geoParams.NZ) res.push_back(toString(sc.geoParams.NZ));
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" The compartment size");
        sysParser.addComment(" - Based on Kuramoto length, see Popov et al., PLoS Comp Biol, 2016 ");
        sysParser.addComment(" - Some chemical reaction rates are scaled based on compartment size");

        sysParser.addStringArgsWithAliases(
            "COMPARTMENTSIZEX", { "COMPARTMENTSIZEX:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "There was an error parsing input file at compartment size. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if(lineVector.size() == 2)
                    sc.geoParams.compartmentSizeX = stod(lineVector[1]);
                else {}
            },
            [] (const SimulConfig& sc) {
                vector< string > res;
                if(sc.geoParams.compartmentSizeX) res.push_back(toString(sc.geoParams.compartmentSizeX));
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "COMPARTMENTSIZEY", { "COMPARTMENTSIZEY:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "There was an error parsing input file at compartment size. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if(lineVector.size() == 2)
                    sc.geoParams.compartmentSizeY = stod(lineVector[1]);
                else {}
            },
            [] (const SimulConfig& sc) {
                vector< string > res;
                if(sc.geoParams.compartmentSizeY) res.push_back(toString(sc.geoParams.compartmentSizeY));
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "COMPARTMENTSIZEZ", { "COMPARTMENTSIZEZ:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "There was an error parsing input file at compartment size. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if(lineVector.size() == 2)
                    sc.geoParams.compartmentSizeZ = stod(lineVector[1]);
                else {}
            },
            [] (const SimulConfig& sc) {
                vector< string > res;
                if(sc.geoParams.compartmentSizeZ) res.push_back(toString(sc.geoParams.compartmentSizeZ));
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "FILCREATIONBOUNDS", { "FILCREATIONBOUNDS:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size()!=7) {
                    log::error("FILCREATIONBOUNDS should have 6 elements. Exiting.");
                    throw std::runtime_error("Invalid arguments.");
                }
                else {
                    auto& s = sc.geoParams.fracGridSpan;
                    parse(s[0][0], lineVector[1]);
                    parse(s[0][1], lineVector[2]);
                    parse(s[0][2], lineVector[3]);
                    parse(s[1][0], lineVector[4]);
                    parse(s[1][1], lineVector[5]);
                    parse(s[1][2], lineVector[6]);
                }
            },
            [] (const SimulConfig& sc) {
                auto& s = sc.geoParams.fracGridSpan;
                return vector<string> {
                    toString(s[0][0]),
                    toString(s[0][1]),
                    toString(s[0][2]),
                    toString(s[1][0]),
                    toString(s[1][1]),
                    toString(s[1][2]),
                };
            }
        );

        sysParser.addEmptyLine();

        sysParser.addComment(" Cylinder setup");
        sysParser.addComment(" - Changes not recommended");

        sysParser.addStringArgsWithAliases(
            "MONOMERSIZE", { "MONOMERSIZE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.geoParams.monomerSize.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.geoParams.monomerSize.push_back(atof((lineVector[i].c_str())));
                }
                else {}
            },
            [] (const SimulConfig& sc) {
                vector< string > res;
                for(const auto& s : sc.geoParams.monomerSize)
                    res.push_back(toString(s));
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "CYLINDERSIZE", { "CYLINDERSIZE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.geoParams.cylinderSize.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.geoParams.cylinderSize.push_back(atof((lineVector[i].c_str())));
                }
                else {}
            },
            [] (const SimulConfig& sc) {
                vector< string > res;
                for(const auto& s : sc.geoParams.cylinderSize)
                    res.push_back(toString(s));
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Surface curvature policy.");
        sysParser.addComment(" - Choose in with-sign, squared or mem3dg");
        sysParser.addSingleArg(
            "surface-curvature-policy",
            [](auto&& conf) -> auto& { return conf.geoParams.surfaceGeometrySettings.curvPol; }
        );

        sysParser.addComment(" Configure mesh adapter.");
        sysParser.addArgs(
            "mesh-adapter",
            [](SimulConfig& conf, const SExpr::ListType& sel) {
                auto& s = conf.geoParams.meshAdapterSettings;
                for(int i = 1; i < sel.size(); ++i) {
                    parseKeyValue<MeshAdapterSettings>(s, sel[i], meshAdapterParser().dict, KeyValueParserUnknownKeyAction::error);
                }
            },
            [](const SimulConfig& conf) {
                auto& s = conf.geoParams.meshAdapterSettings;
                return buildTokens<MeshAdapterSettings>(s, meshAdapterParser());
            }
        );

        sysParser.addEmptyLine();
    }

    //--------------------------------------------------------------------------
    // Boundary parameters.
    //--------------------------------------------------------------------------
    {
        sysParser.addComment("==================================================");
        sysParser.addComment(" Boundary parameters");
        sysParser.addComment("==================================================");
        sysParser.addEmptyLine();

        sysParser.addComment(" Define network boundary geometry (CUBIC, SPHERICAL, CYLINDER)");

        sysParser.addStringArgsWithAliases(
            "BOUNDARYSHAPE", { "BOUNDARYSHAPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout << "A boundary shape needs to be specified. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    sc.boundParams.boundaryType.boundaryShape = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { sc.boundParams.boundaryType.boundaryShape };
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Define how network boundary moves");
        sysParser.addComment(" - Usage: BOUNDARYMOVE LEFT/RIGHT/BOTTOM/TOP/FRONT/BACK");
        sysParser.addComment(" - Changes are not recommended.");

        sysParser.addStringArgsWithAliases(
            "BOUNDARYMOVE", { "BOUNDARYMOVE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout << "A boundary move type needs to be specified. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    sc.boundParams.boundaryType.boundaryMove.push_back(lineVector[1]);
                }
            },
            [] (const SimulConfig& sc) {
                vector<vector<string>> res;
                for(const auto& s : sc.boundParams.boundaryType.boundaryMove) {
                    res.push_back({s});
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Set diameter for SPHERICAL or CYLINDER type");
        sysParser.addComment(" - CUBIC: No need to set, boundary is the same as network size");

        sysParser.addStringArgsWithAliases(
            "BOUNDARYDIAMETER", { "BOUNDARYDIAMETER:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() == 2) {
                    sc.boundParams.diameter = atof((lineVector[1].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                return sc.boundParams.diameter ?
                    vector<string> { toString(sc.boundParams.diameter) } :
                    vector<string> {};
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Boundary interactions ");
        sysParser.addComment(" - Repulsion: Popov et al, 2016, PLoS Comp Biol");

        sysParser.addStringArgsWithAliases(
            "BOUNDARYCUTOFF", { "BOUNDARYCUTOFF:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                        "There was an error parsing input file at Boundary parameters. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    sc.boundParams.BoundaryCutoff = atof((lineVector[1].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.boundParams.BoundaryCutoff) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "BOUNDARYINTERACTIONK", { "BOUNDARYINTERACTIONK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                        "There was an error parsing input file at Boundary parameters. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    sc.boundParams.BoundaryK = atof((lineVector[1].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.boundParams.BoundaryK) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "BOUNDARYSCREENLENGTH", { "BOUNDARYSCREENLENGTH:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                        "There was an error parsing input file at Boundary parameters. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    sc.boundParams.BScreenLength = atof((lineVector[1].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.boundParams.BScreenLength) };
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Set boundary move parameters");
        sysParser.addComment(" - Changes are not recommended.");

        sysParser.addStringArgsWithAliases(
            "BMOVESPEED", { "BMOVESPEED:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() == 2) {
                    sc.boundParams.moveSpeed = atof((lineVector[1].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                return sc.boundParams.moveSpeed ?
                    vector<string> { toString(sc.boundParams.moveSpeed) } :
                    vector<string> {};
            }
        );
        sysParser.addStringArgsWithAliases(
            "BMOVESTARTTIME", { "BMOVESTARTTIME:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() == 2) {
                    sc.boundParams.moveStartTime = atof((lineVector[1].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                return sc.boundParams.moveStartTime ?
                    vector<string> { toString(sc.boundParams.moveStartTime) } :
                    vector<string> {};
            }
        );
        sysParser.addStringArgsWithAliases(
            "BMOVEENDTIME", { "BMOVEENDTIME:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() == 2) {
                    sc.boundParams.moveEndTime = atof((lineVector[1].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                return sc.boundParams.moveEndTime ?
                    vector<string> { toString(sc.boundParams.moveEndTime) } :
                    vector<string> {};
            }
        );
        sysParser.addStringArgsWithAliases(
            "TRANSFERSHAREAXIS", { "TRANSFERSHAREAXIS:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 3) {
                    cout <<
                        "There was an error parsing input file at Chemistry parameters. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }

                else{
                    cout<<"TRANSFERSHARE AXIS "<<lineVector[2]<<endl;
                    if(lineVector[2]=="X")
                        sc.boundParams.transfershareaxis=0;
                    else if(lineVector[2]=="Y")
                        sc.boundParams.transfershareaxis=1;
                    else if(lineVector[2]=="Z")
                        sc.boundParams.transfershareaxis=2;
                    else if(lineVector[2]=="RADIAL") {
                        sc.boundParams.transfershareaxis = 3;
                        cout<<"RADIAL transfer not implemented. Change paramters. Exiting"
                                "."<<endl;
                        exit(EXIT_FAILURE);
                    }
                    else{
                        cout << "There was an error parsing input file at Chemistry parameters. Exiting."
                            << endl;
                        exit(EXIT_FAILURE);
                    }
                }
            },
            [] (const SimulConfig& sc) {
                vector<vector<string>> res;
                switch(sc.boundParams.transfershareaxis) {
                    // TODO: figure out the 1st element
                case 0:
                    res.push_back({ "???", "X" });
                    break;
                case 1:
                    res.push_back({ "???", "Y" });
                    break;
                case 2:
                    res.push_back({ "???", "Z" });
                    break;
                case 3:
                    res.push_back({ "???", "RADIAL" });
                    break;
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

    }

    //--------------------------------------------------------------------------
    // Mechanical parameters.
    //--------------------------------------------------------------------------
    {
        sysParser.addComment("==================================================");
        sysParser.addComment(" Mechanical parameters");
        sysParser.addComment("==================================================");
        sysParser.addEmptyLine();

        sysParser.addComment("====== Mechanical algorithm ======");
        sysParser.addEmptyLine();

        sysParser.addStringArgsWithAliases(
            "CONJUGATEGRADIENT", { "CONJUGATEGRADIENT:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                        "A conjugate gradient method must be specified. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsAlgorithm.ConjugateGradient = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<vector<string>> res;
                if(const auto& s = sc.mechParams.mechanicsAlgorithm.ConjugateGradient; !s.empty()) {
                    res.push_back({ s });
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "MD", {},
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout << "A Mechanics algorithm must be specified. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsAlgorithm.MD = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<vector<string>> res;
                if(const auto& s = sc.mechParams.mechanicsAlgorithm.MD; !s.empty()) {
                    res.push_back({ s });
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "GRADIENTTOLERANCE", { "GRADIENTTOLERANCE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsAlgorithm.gradientTolerance = atof(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& x = sc.mechParams.mechanicsAlgorithm.gradientTolerance; x) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addSingleArg(
            "energy-change-relative-tolerance",
            [](auto&& conf) -> auto& { return conf.mechParams.mechanicsAlgorithm.energyChangeRelativeTolerance; }
        );
        sysParser.addStringArgsWithAliases(
            "MAXDISTANCE", { "MAXDISTANCE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsAlgorithm.maxDistance = atof(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& x = sc.mechParams.mechanicsAlgorithm.maxDistance; x) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "LAMBDAMAX", { "LAMBDAMAX:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsAlgorithm.lambdaMax = atof(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& x = sc.mechParams.mechanicsAlgorithm.lambdaMax; x) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "LAMBDARUNNINGAVERAGEPROBABILITY", { "LAMBDARUNNINGAVERAGEPROBABILITY:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsAlgorithm.lambdarunningaverageprobability = atof(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& x = sc.mechParams.mechanicsAlgorithm.lambdarunningaverageprobability; x) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "LINESEARCHALGORITHM", { "LINESEARCHALGORITHM:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsAlgorithm.linesearchalgorithm = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<vector<string>> res;
                if(const auto& s = sc.mechParams.mechanicsAlgorithm.linesearchalgorithm; !s.empty()) {
                    res.push_back({ s });
                }
                return res;
            }
        );
        sysParser.addComment(" Whether to try to restore when line search error happens. Default to true.");
        sysParser.addStringArgs(
            "try-to-recover-in-line-search-error",
            [] (SimulConfig& conf, const vector<string>& lineVector) {
                if(lineVector.size() == 2) {
                    if(lineVector[1] == "true") {
                        conf.mechParams.mechanicsAlgorithm.tryToRecoverInLineSearchError = true;
                    }
                    else if(lineVector[1] == "false") {
                        conf.mechParams.mechanicsAlgorithm.tryToRecoverInLineSearchError = false;
                    }
                    else {
                        log::error("Must specify \"true\" or \"false\" for try-to-recover-in-line-search-error.");
                        throw runtime_error("Invalid input variable.");
                    }
                }
            },
            [] (const SimulConfig& conf) {
                return vector<string> { conf.mechParams.mechanicsAlgorithm.tryToRecoverInLineSearchError ? "true" : "false" };
            }
        );
        sysParser.addSingleArg(
            "adapt-mesh-during-minimization",
            [](auto&& conf) -> auto& { return conf.mechParams.mechanicsAlgorithm.adaptMeshDuringMinimization; }
        );
        sysParser.addEmptyLine();

        sysParser.addComment("====== Force fields ======");
        sysParser.addEmptyLine();

        sysParser.addComment(" Actin filaments ");
        sysParser.addEmptyLine();

        sysParser.addComment(" Stretching: Popov et al, 2016, PLoS Comp Biol");
        sysParser.addStringArgsWithAliases(
            "FSTRETCHINGFFTYPE", { "FSTRETCHINGFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Filament stretching FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.FStretchingType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.FStretchingType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "FSTRETCHINGK", { "FSTRETCHINGK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.FStretchingK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.FStretchingK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.FStretchingK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Bending: Ott et al, 1993, Phys Rev E");
        sysParser.addStringArgsWithAliases(
            "FBENDINGFFTYPE", { "FBENDINGFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Filament bending FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.FBendingType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.FBendingType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "FBENDINGK", { "FBENDINGK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.FBendingK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.FBendingK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.FBendingK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "FBENDINGTHETA", { "FBENDINGTHETA:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.FBendingTheta.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.FBendingTheta.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.FBendingTheta) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Twisting: Currently not implemented");
        sysParser.addStringArgsWithAliases(
            "FTWISTINGFFTYPE", { "FTWISTINGFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Filament twisting FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.FTwistingType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.FTwistingType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "FTWISTINGK", { "FTWISTINGK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.FTwistingK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.FTwistingK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.FTwistingK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "FTWISTINGPHI", { "FTWISTINGPHI:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.FTwistingPhi.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.FTwistingPhi.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.FTwistingPhi) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Linkers");
        sysParser.addEmptyLine();
        sysParser.addComment(" Stretching: Didonna et al, Phys Rev E, 2007");
        sysParser.addStringArgsWithAliases(
            "LSTRETCHINGFFTYPE", { "LSTRETCHINGFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Linker stretching FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.LStretchingType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.LStretchingType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "LSTRETCHINGK", { "LSTRETCHINGK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.LStretchingK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.LStretchingK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.LStretchingK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();
        sysParser.addComment(" Bending/Twisting: not implemented");
        sysParser.addStringArgsWithAliases(
            "LBENDINGFFTYPE", { "LBENDINGFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Linker bending FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.LBendingType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.LBendingType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "LBENDINGK", { "LBENDINGK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.LBendingK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.LBendingK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.LBendingK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "LBENDINGTHETA", { "LBENDINGTHETA:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.LBendingTheta.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.LBendingTheta.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.LBendingTheta) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "LTWISTINGFFTYPE", { "LTWISTINGFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Linker twisting FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.LTwistingType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.LTwistingType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "LTWISTINGK", { "LTWISTINGK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.LTwistingK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.LTwistingK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.LTwistingK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "LTWISTINGPHI", { "LTWISTINGPHI:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.LTwistingPhi.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.LTwistingPhi.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.LTwistingPhi) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Myosin motors");
        sysParser.addEmptyLine();
        sysParser.addComment(" Stretching: Vilfan, Biophys J, 2010");
        sysParser.addStringArgsWithAliases(
            "MSTRETCHINGFFTYPE", { "MSTRETCHINGFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Motor stretching FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.MStretchingType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.MStretchingType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "MSTRETCHINGK", { "MSTRETCHINGK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.MStretchingK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.MStretchingK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.MStretchingK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();
        sysParser.addComment(" Bending/Twisting: not implemented");
        sysParser.addStringArgsWithAliases(
            "MBENDINGFFTYPE", { "MBENDINGFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Motor bending FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.MBendingType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.MBendingType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "MBENDINGK", { "MBENDINGK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.MBendingK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.MBendingK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.MBendingK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "MBENDINGTHETA", { "MBENDINGTHETA:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.MBendingTheta.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.MBendingTheta.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.MBendingTheta) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "MTWISTINGFFTYPE", { "MTWISTINGFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Motor twisting FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.MTwistingType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.MTwistingType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "MTWISTINGK", { "MTWISTINGK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.MTwistingK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.MTwistingK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.MTwistingK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "MTWISTINGPHI", { "MTWISTINGPHI:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.MTwistingPhi.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.MTwistingPhi.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.MTwistingPhi) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Arp2/3 brancher");
        sysParser.addComment(" - 4 force fields: Popov et al, Plos Comp Biol, 2016");
        sysParser.addComment(" - No reliable literature values for this FF");
        sysParser.addEmptyLine();
        sysParser.addStringArgsWithAliases(
            "BRSTRETCHINGFFTYPE", { "BRSTRETCHINGFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Brancher stretching FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.BrStretchingType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.BrStretchingType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "BRSTRETCHINGK", { "BRSTRETCHINGK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.BrStretchingK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.BrStretchingK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.BrStretchingK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "BRSTRETCHINGL", { "BRSTRETCHINGL:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.BrStretchingL.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.BrStretchingL.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.BrStretchingL) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();
        sysParser.addStringArgsWithAliases(
            "BRBENDINGFFTYPE", { "BRBENDINGFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Brancher bending FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.BrBendingType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.BrBendingType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "BRBENDINGK", { "BRBENDINGK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.BrBendingK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.BrBendingK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.BrBendingK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "BRBENDINGTHETA", { "BRBENDINGTHETA:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.BrBendingTheta.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.BrBendingTheta.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.BrBendingTheta) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();
        sysParser.addStringArgsWithAliases(
            "BRDIHEDRALFFTYPE", { "BRDIHEDRALFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Brancher dihedral FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.BrDihedralType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.BrDihedralType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "BRDIHEDRALK", { "BRDIHEDRALK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.BrDihedralK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.BrDihedralK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.BrDihedralK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();
        sysParser.addStringArgsWithAliases(
            "BRPOSITIONFFTYPE", { "BRPOSITIONFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Brancher position FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.BrPositionType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.BrPositionType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "BRPOSITIONK", { "BRPOSITIONK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.BrPositionK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.BrPositionK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.BrPositionK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Boundary");
        sysParser.addComment(" Boundary repulsion force field. Parameters are set in the boundary section.");
        sysParser.addStringArgsWithAliases(
            "BOUNDARYFFTYPE", { "BOUNDARYFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Boundary FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.BoundaryFFType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.BoundaryFFType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Volume exclusion: Popov et al, 2016, PLoS Comp Biol");
        sysParser.addStringArgsWithAliases(
            "VOLUMEFFTYPE", { "VOLUMEFFTYPE:" },
            [] (SimulConfig& conf, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    log::error("There was an error parsing input file at Volume FF type. Exiting.");
                    throw runtime_error("Error parsing volume FF type.");
                }
                else if (lineVector.size() == 2) {
                    parse(conf.mechParams.mechanicsFFType.cylinderVolumeExclusionFFType, lineVector[1]);
                }
            },
            [] (const SimulConfig& conf) {
                return vector<string> { toString(conf.mechParams.mechanicsFFType.cylinderVolumeExclusionFFType) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "VOLUMEK", { "VOLUMEK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.VolumeK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.VolumeK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.VolumeK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "VOLUMECUTOFF", { "VOLUMECUTOFF:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                        "Error reading Volume cutoff. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    sc.mechParams.VolumeCutoff = stod(lineVector[1]);
                }
            },
            [] (const SimulConfig& sc) {
                vector<vector<string>> res;
                if(sc.mechParams.VolumeCutoff) {
                    res.push_back({ toString(sc.mechParams.VolumeCutoff) });
                }
                return res;
            }
        );
        sysParser.addStringArgs(
            "volume-exclusion-monomer-interval",
            [](SimulConfig& conf, const vector<string>& args) {
                conf.mechParams.volumeExclusionMonomerInterval.clear();
                for(int i = 1; i < args.size(); ++i) {
                    conf.mechParams.volumeExclusionMonomerInterval.push_back(parse<Size>(args[i]));
                }
            },
            [](const SimulConfig& conf) {
                vector<string> res;
                for(const auto& i : conf.mechParams.volumeExclusionMonomerInterval) {
                    res.push_back(toString(i));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Bubble interactions");
        sysParser.addStringArgsWithAliases(
            "BUBBLEFFTYPE", { "BUBBLEFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Bubble FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.BubbleFFType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.BubbleFFType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "BUBBLEINTERACTIONK", { "BUBBLEINTERACTIONK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.BubbleK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.BubbleK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.BubbleK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "BUBBLESCREENLENGTH", { "BUBBLESCREENLENGTH:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.BubbleScreenLength.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.BubbleScreenLength.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.BubbleScreenLength) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "BUBBLECUTOFF", { "BUBBLECUTOFF:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                        "Error reading Bubble cutoff. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    sc.mechParams.BubbleCutoff = stod(lineVector[1]);
                }
            },
            [] (const SimulConfig& sc) {
                vector<vector<string>> res;
                if(sc.mechParams.BubbleCutoff) {
                    res.push_back({ toString(sc.mechParams.BubbleCutoff) });
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "BUBBLERADIUS", { "BUBBLERADIUS:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.BubbleRadius.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.BubbleRadius.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.BubbleRadius) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "NUMBUBBLETYPES", { "NUMBUBBLETYPES:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                        "Error reading number of Bubble types. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    sc.mechParams.numBubbleTypes = atoi(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                vector<vector<string>> res;
                if(sc.mechParams.numBubbleTypes) {
                    res.push_back({ toString(sc.mechParams.numBubbleTypes) });
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" MTOC");
        sysParser.addStringArgsWithAliases(
            "MTOCFFTYPE", { "MTOCFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at MTOC FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.MTOCFFType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.MTOCFFType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "MTOCBENDINGK", { "MTOCBENDINGK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.MTOCBendingK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.MTOCBendingK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.MTOCBendingK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" AFM bead");
        sysParser.addStringArgsWithAliases(
            "AFMFFTYPE", { "AFMFFTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at AFM FF type. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.mechParams.mechanicsFFType.AFMFFType = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(const auto& s = sc.mechParams.mechanicsFFType.AFMFFType; !s.empty()) {
                    res.push_back(s);
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "AFMBENDINGK", { "AFMBENDINGK:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.mechParams.AFMBendingK.clear();
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.mechParams.AFMBendingK.push_back(atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& k : sc.mechParams.AFMBendingK) {
                    res.push_back(toString(k));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Special bead constraint force field only used in GC filament model.");
        sysParser.addStringArgs(
            "bead-constraint-k",
            [](SimulConfig& conf, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    LOG(ERROR) << "Error reading bead-constraint-k.";
                    throw runtime_error("Error reading bead-constraint-k.");
                }
                conf.mechParams.beadConstraintK = stod(lineVector[1]);
            },
            [] (const SimulConfig& conf) {
                return vector<string> { toString(conf.mechParams.beadConstraintK) };
            }
        );
        sysParser.addStringArgs(
            "bead-constraint-rmax",
            [](SimulConfig& conf, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    LOG(ERROR) << "Error reading bead-constraint-rmax.";
                    throw runtime_error("Error reading bead-constraint-rmax.");
                }
                conf.mechParams.beadConstraintRmax = stod(lineVector[1]);
            },
            [] (const SimulConfig& conf) {
                return vector<string> { toString(conf.mechParams.beadConstraintRmax) };
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Membrane related force fields");
        sysParser.addSingleArg("membrane-stretching-ff-type", [](auto&& sc) -> auto& { return sc.mechParams.mechanicsFFType.memStretchingFFType; });
        sysParser.addSingleArg("membrane-tension-ff-type", [](auto&& sc) -> auto& { return sc.mechParams.mechanicsFFType.memTensionFFType; });
        sysParser.addSingleArg("membrane-bending-ff-type", [](auto&& sc) -> auto& { return sc.mechParams.mechanicsFFType.memBendingFFType; });
        sysParser.addSingleArg("volume-conservation-ff-type", [](auto&& sc) -> auto& { return sc.mechParams.mechanicsFFType.volumeConservationFFType; });
        sysParser.addSingleArg("triangle-bead-volume-ff-type", [](auto&& sc) -> auto& { return sc.mechParams.mechanicsFFType.triangleBeadVolumeFFType; });

        sysParser.addSingleArg("triangle-bead-volume-k", [](auto&& sc) -> auto& { return sc.mechParams.triangleBeadVolume.k; });
        sysParser.addSingleArg("triangle-bead-volume-cutoff", [](auto&& sc) -> auto& { return sc.mechParams.triangleBeadVolume.cutoff; });
        sysParser.addSingleArg("triangle-bead-volume-cutoff-mech", [](auto&& sc) -> auto& { return sc.mechParams.triangleBeadVolume.cutoffMech; });
        sysParser.addEmptyLine();

        sysParser.addComment(" Membrane protein curvature mismatch");
        sysParser.addComment(" - Usage: (curvature-mismatch <species-name> <k-bending> <eq-curvature>)");
        sysParser.addStringArgs(
            "curvature-mismatch",
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 4) {
                    LOG(ERROR) << "Curvature mismatch must have 3 arguments.";
                    throw runtime_error("Curvature mismatch must have 3 arguments.");
                }
                else {
                    sc.mechParams.proteinCurvatureMismatchSetups.push_back({
                        lineVector[1],
                        parse<double>(lineVector[2]),
                        parse<double>(lineVector[3]),
                    });
                }
            },
            [] (const SimulConfig& sc) {
                vector<vector<string>> res;
                for(const auto& setup : sc.mechParams.proteinCurvatureMismatchSetups) {
                    res.push_back({
                        setup.speciesName,
                        toString(setup.kBending),
                        toString(setup.eqCurv),
                    });
                }
                return res;
            }
        );
        sysParser.addComment(" The following options can turn on/off curvature sensing/generation using curvature mismatch model. Typically, they should both be set to true.");
        sysParser.addSingleArg("curvature-mismatch-enable-curvature-sensing",    [](auto&& conf) -> auto& { return conf.mechParams.curvatureMismatchEnableCurvatureSensing; });
        sysParser.addSingleArg("curvature-mismatch-enable-curvature-generation", [](auto&& conf) -> auto& { return conf.mechParams.curvatureMismatchEnableCurvatureGeneration; });
        sysParser.addEmptyLine();

        sysParser.addComment("====== Protocols ======");
        sysParser.addEmptyLine();

        sysParser.addComment(" Hessian tracking");
        sysParser.addStringArgsWithAliases(
            "HESSIANTRACKING", { "HESSIANTRACKING:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 5) {
                    cout <<
                    "There was an error parsing input file at Hessian tracking. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 5) {
                    sc.mechParams.hessTracking = true;
                    //sc.mechParams.hessDelta = atof(lineVector[1].c_str());
                    sc.mechParams.hessSkip = atof(lineVector[1].c_str());
                    int dense = atoi(lineVector[2].c_str());
                    if(dense == 1){
                        sc.mechParams.denseEstimationBool = true;
                    }else{
                        sc.mechParams.denseEstimationBool = false;
                    }
                    int rocksnapbool = atoi(lineVector[3].c_str());
                    if(rocksnapbool == 1){
                        sc.mechParams.rockSnapBool = true;
                    }else{
                        sc.mechParams.rockSnapBool = false;
                    }
                    int hessmatprintbool = atoi(lineVector[4].c_str());
                    if(hessmatprintbool == 1){
                        sc.mechParams.hessMatrixPrintBool = true;
                    }else{
                        sc.mechParams.hessMatrixPrintBool = false;
                    }
                }
            },
            [] (const SimulConfig& sc) {
                vector<vector<string>> res;
                if(sc.mechParams.hessTracking) {
                    res.push_back({
                        toString(sc.mechParams.hessSkip),
                        toString(sc.mechParams.denseEstimationBool ? 0 : 1),
                        toString(sc.mechParams.rockSnapBool ? 0 : 1),
                        toString(sc.mechParams.hessMatrixPrintBool ? 0 : 1),
                    });
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "EIGENTRACKING", { "EIGENTRACKING:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                        "There was an error parsing input file at Chemistry algorithm. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    const char * testStr1 = "OFF";
                    const char * testStr2 = lineVector[1].c_str();
                    if(strcmp(testStr1, testStr2) == 0)
                        sc.mechParams.eigenTracking = false;
                    else
                        sc.mechParams.eigenTracking = true;
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                res.push_back(sc.mechParams.eigenTracking ? "ON" : "OFF");
                return res;
            }
        );

        sysParser.addEmptyLine();

        sysParser.addComment(" Special protocols");
        sysParser.addStringArgsWithAliases(
            "SPECIALPROTOCOL", { "SPECIALPROTOCOL:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() <= 1) {
                    // Ignore empty protocol.
                }
                else if(lineVector[1] == "PINBOUNDARYFILAMENTS") {
                    if(lineVector.size() != 5) {
                        log::error("There was an error parsing input file at pinning boundary filaments. Exiting.");
                        log::info("Usage: (SPECIALPROTOCOL PINBOUNDARYFILAMENTS <k> <distance> <time>)");
                        throw runtime_error("Invalid boundary filament pinning.");
                    }
                    else{
                        sc.mechParams.pinBoundaryFilaments = true;
                        sc.mechParams.pinK = atof(lineVector[2].c_str());
                        sc.mechParams.pinDistance = atof(lineVector[3].c_str());
                        sc.mechParams.pinTime = atof(lineVector[4].c_str());
                    }

                }
                else if (lineVector[1] == "PINLOWERBOUNDARYFILAMENTS") {

                    if(lineVector.size() != 5) {
                        log::error("There was an error parsing input file at pinning lower boundary filaments. Exiting.");
                        log::info("Usage: (SPECIALPROTOCOL PINLOWERBOUNDARYFILAMENTS <k> <time> <fraction>)");
                        throw runtime_error("Invalid lower boundary filament pinning.");
                    }
                    else {
                        sc.mechParams.pinLowerBoundaryFilaments = true;
                        sc.mechParams.pinK = atof(lineVector[2].c_str());
                        sc.mechParams.pinTime = atof(lineVector[3].c_str());
                        sc.mechParams.pinFraction = atof(lineVector[4].c_str());
                    }
                }
                else if(lineVector[1] == "pin-initial-filament-below-z") {
                    if(lineVector.size() != 4) {
                        log::error("Invalid pin-initial-filament-below-z arguments.");
                        log::info("Usage: (SPECIALPROTOCOL pin-initial-filament-below-z <k> <z>)");
                        throw runtime_error("Invalid pin-initial-filament-below-z arguments.");
                    }
                    else {
                        sc.mechParams.pinInitialFilamentBelowZ = true;
                        sc.mechParams.pinK = parse<double>(lineVector[2]);
                        sc.mechParams.pinInitialFilamentBelowZValue = parse<double>(lineVector[3]);
                    }
                }
                else if(lineVector[1] == "MAKELINKERSSTATIC") {
                    if (lineVector.size() != 3) {
                        log::error("Invalid MAKELINKERSTATIC arguments.");
                        log::info("Usage: (SPECIALPROTOCOL MAKELINKERSSTATIC <time>)");
                        throw runtime_error("Invalid MAKELINKERSTATIC arguments.");
                    } else {
                        sc.chemParams.makeLinkersStatic = true;
                        sc.chemParams.makeLinkersStaticTime = atof(lineVector[2].c_str());
                    }
                }
                else if(lineVector[1] == "MAKEFILAMENTSSTATIC") {
                    if (lineVector.size() != 3) {
                        log::error("Invalid MAKEFILAMENTSSTATIC arguments.");
                        log::info("Usage: (SPECIALPROTOCOL MAKEFILAMENTSSTATIC <time>)");
                        throw runtime_error("Invalid MAKEFILAMENTSSTATIC arguments.");
                    } else {
                        sc.chemParams.makeFilamentsStatic = true;
                        sc.chemParams.makeFilamentsStaticTime = atof(lineVector[2].c_str());
                    }
                    
                }
                else if(lineVector[1]  == "RATEDEPEND") {
                    if (lineVector.size() != 4) {
                        log::error("Invalid RATEDEPEND arguments.");
                        log::info("Usage: (SPECIALPROTOCOL RATEDEPEND <time> <force>)");
                        throw runtime_error("Invalid RATEDEPEND arguments.");
                    } else {
                        sc.chemParams.makeRateDepend = true;
                        sc.chemParams.makeRateDependTime = atof(lineVector[2].c_str());
                        sc.chemParams.makeRateDependForce = atof(lineVector[3].c_str());
                    }
                }
                else if(lineVector[1]  == "AFM") {
                    if (lineVector.size() != 7) {
                        log::error("Invalid AFM arguments.");
                        log::info("Usage: (SPECIALPROTOCOL AFM <step1> <step2> <iter-change> <step-total> <step-time>)");
                        throw runtime_error("Invalid AFM arguments.");
                    } else {
                        sc.chemParams.makeAFM = true;
                        //displacement of each pull
                        sc.chemParams.AFMStep1 = atof(lineVector[2].c_str());
                        sc.chemParams.AFMStep2 = atof(lineVector[3].c_str());
                        //change dispalcement from 1 to 2
                        sc.chemParams.IterChange = atof(lineVector[4].c_str());
                        //total step of each AFM pull
                        sc.chemParams.StepTotal = atof(lineVector[5].c_str());
                        //time between each pull
                        sc.chemParams.StepTime = atof(lineVector[6].c_str());
                    }
                }
                else if(lineVector[1] == "scale-membrane-eq-area") {
                    if(lineVector.size() == 4) {
                        auto& change = sc.mechParams.membraneEqAreaChange.emplace();
                        parse(change.rate,      lineVector[2]);
                        parse(change.minEqArea, lineVector[3]);
                    }
                    else {
                        log::error("Invalid scale-membrane-eq-area arguments.");
                        log::info("Usage: (SPECIALPROTOCOL scale-membrane-eq-area <rate> <min-eq-area>)");
                        throw runtime_error("Invalid scale-membrane-eq-area arguments.");
                    }
                }
                else if(lineVector[1] == "scale-membrane-eq-volume") {
                    if(lineVector.size() == 4) {
                        auto& change = sc.mechParams.membraneEqVolumeChange.emplace();
                        parse(change.rate,        lineVector[2]);
                        parse(change.minEqVolume, lineVector[3]);
                    }
                    else {
                        log::error("Invalid scale-membrane-eq-volume arguments.");
                        log::info("Usage: (SPECIALPROTOCOL scale-membrane-eq-volume <rate> <min-eq-volume>)");
                        throw runtime_error("Invalid scale-membrane-eq-volume arguments.");
                    }
                }
                else {
                    // Invalid special protocol.
                    log::error("Invalid special protocol: {}", lineVector[1]);
                    throw runtime_error("Invalid special protocol.");
                }
            },
            [] (const SimulConfig& sc) {
                vector<vector<string>> res;
                if(sc.mechParams.pinBoundaryFilaments) {
                    res.push_back({
                        "PINBOUNDARYFILAMENTS",
                        toString(sc.mechParams.pinK),
                        toString(sc.mechParams.pinDistance),
                        toString(sc.mechParams.pinTime)
                    });
                }
                if(sc.mechParams.pinLowerBoundaryFilaments) {
                    res.push_back({
                        "PINLOWERBOUNDARYFILAMENTS",
                        toString(sc.mechParams.pinK),
                        toString(sc.mechParams.pinTime),
                        toString(sc.mechParams.pinFraction)
                    });
                }
                if(sc.mechParams.pinInitialFilamentBelowZ) {
                    res.push_back({
                        "pin-initial-filament-below-z",
                        toString(sc.mechParams.pinK),
                        toString(sc.mechParams.pinInitialFilamentBelowZValue),
                    });
                }
                if(sc.chemParams.makeLinkersStatic) {
                    res.push_back({ "MAKELINKERSSTATIC", toString(sc.chemParams.makeLinkersStaticTime) });
                }
                if(sc.chemParams.makeFilamentsStatic) {
                    res.push_back({ "MAKEFILAMENTSSTATIC", toString(sc.chemParams.makeFilamentsStaticTime) });
                }
                if(sc.chemParams.makeRateDepend) {
                    res.push_back({
                        "RATEDEPEND",
                        toString(sc.chemParams.makeRateDependTime),
                        toString(sc.chemParams.makeRateDependForce)
                    });
                }
                if(sc.chemParams.makeAFM) {
                    res.push_back({
                        "AFM",
                        toString(sc.chemParams.AFMStep1),
                        toString(sc.chemParams.AFMStep2),
                        toString(sc.chemParams.IterChange),
                        toString(sc.chemParams.StepTotal),
                        toString(sc.chemParams.StepTime)
                    });
                }
                if(sc.mechParams.membraneEqAreaChange.has_value()) {
                    auto& change = sc.mechParams.membraneEqAreaChange.value();
                    res.push_back({
                        "scale-membrane-eq-area",
                        toString(change.rate),
                        toString(change.minEqArea),
                    });
                }
                if(sc.mechParams.membraneEqVolumeChange.has_value()) {
                    auto& change = sc.mechParams.membraneEqVolumeChange.value();
                    res.push_back({
                        "scale-membrane-eq-volume",
                        toString(change.rate),
                        toString(change.minEqVolume),
                    });
                }
                return res;
            }
        );
        sysParser.addEmptyLine();
    }

    //--------------------------------------------------------------------------
    // Chemistry parameters.
    //--------------------------------------------------------------------------
    {
        sysParser.addComment("==================================================");
        sysParser.addComment(" Chemistry parameters");
        sysParser.addComment("==================================================");
        sysParser.addEmptyLine();

        sysParser.addComment("====== Chemistry algorithm ======");
        sysParser.addEmptyLine();

        sysParser.addStringArgsWithAliases(
            "CALGORITHM", { "CALGORITHM:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                        "There was an error parsing input file at Chemistry algorithm. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.chemistryAlgorithm.algorithm = lineVector[1];
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { sc.chemParams.chemistryAlgorithm.algorithm };
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Use either time mode or step mode");
        sysParser.addStringArgsWithAliases(
            "RUNTIME", { "RUNTIME:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Chemistry algorithm. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.chemistryAlgorithm.runTime = atof(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.chemParams.chemistryAlgorithm.runTime) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "SNAPSHOTTIME", { "SNAPSHOTTIME:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Chemistry algorithm. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.chemistryAlgorithm.snapshotTime = atof(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.chemParams.chemistryAlgorithm.snapshotTime) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "MINIMIZATIONTIME", { "MINIMIZATIONTIME:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                        "There was an error parsing input file at Chemistry algorithm. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.chemistryAlgorithm.minimizationTime = atof(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.chemParams.chemistryAlgorithm.minimizationTime) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "NEIGHBORLISTTIME", { "NEIGHBORLISTTIME:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                        "There was an error parsing input file at Chemistry algorithm. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.chemistryAlgorithm.neighborListTime = atof(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.chemParams.chemistryAlgorithm.neighborListTime) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "INITIALSLOWDOWNTIME", { "INITIALSLOWDOWNTIME:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                        "There was an error parsing input file at Chemistry algorithm. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.chemistryAlgorithm.initialSlowDownTime = atof(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.chemParams.chemistryAlgorithm.initialSlowDownTime) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "DATADUMPTIME", { "DATADUMPTIME:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Chemistry algorithm. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.chemistryAlgorithm.datadumpTime = atof(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.chemParams.chemistryAlgorithm.datadumpTime) };
            }
        );
        sysParser.addEmptyLine();
        sysParser.addStringArgsWithAliases(
            "RUNSTEPS", { "RUNSTEPS:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Chemistry algorithm. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.chemistryAlgorithm.runSteps = atoi(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.chemParams.chemistryAlgorithm.runSteps) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "SNAPSHOTSTEPS", { "SNAPSHOTSTEPS:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Chemistry algorithm. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.chemistryAlgorithm.snapshotSteps = atoi(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.chemParams.chemistryAlgorithm.snapshotSteps) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "MINIMIZATIONSTEPS", { "MINIMIZATIONSTEPS:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                        "There was an error parsing input file at Chemistry algorithm. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.chemistryAlgorithm.minimizationSteps = atoi(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.chemParams.chemistryAlgorithm.minimizationSteps) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "NEIGHBORLISTSTEPS", { "NEIGHBORLISTSTEPS:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                        "There was an error parsing input file at Chemistry algorithm. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.chemistryAlgorithm.neighborListSteps = atoi(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.chemParams.chemistryAlgorithm.neighborListSteps) };
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment("====== Chemistry setup ======");
        sysParser.addEmptyLine();
        sysParser.addStringArgsWithAliases(
            "CHEMISTRYFILE", { "CHEMISTRYFILE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "Error reading chemistry input file. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }

                else if (lineVector.size() == 2)
                    sc.chemParams.chemistrySetup.inputFile = lineVector[1];

                else if(lineVector.size() < 2) {
                    cout << "Must specify a chemistry input file. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { sc.chemParams.chemistrySetup.inputFile.string() };
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment("====== Numbers and protocols ======");
        sysParser.addEmptyLine();

        sysParser.addComment(" Numbers");
        sysParser.addStringArgsWithAliases(
            "NUMFILAMENTTYPES", { "NUMFILAMENTTYPES:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout <<
                        "There was an error parsing input file at Chemistry parameters. Exiting."
                        << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.numFilaments = atoi(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.chemParams.numFilaments) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "NUMBINDINGSITES", { "NUMBINDINGSITES:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.chemParams.numBindingSites.push_back(atoi(lineVector[i].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(auto x : sc.chemParams.numBindingSites) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "NUMMOTORHEADSMIN", { "NUMMOTORHEADSMIN:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.chemParams.motorNumHeadsMin.push_back(atoi(lineVector[i].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(auto x : sc.chemParams.motorNumHeadsMin) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "NUMMOTORHEADSMAX", { "NUMMOTORHEADSMAX:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.chemParams.motorNumHeadsMax.push_back(atoi(lineVector[i].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(auto x : sc.chemParams.motorNumHeadsMax) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "MOTORSTEPSIZE", { "MOTORSTEPSIZE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.chemParams.motorStepSize.push_back(atof(lineVector[i].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(auto x : sc.chemParams.motorStepSize) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "LINKERBINDINGSKIP", { "LINKERBINDINGSKIP:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                    "There was an error parsing input file at Chemistry algorithm. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.linkerbindingskip = atoi(lineVector[1].c_str());
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.chemParams.linkerbindingskip) };
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Whether to allow all pairwise bindings on the same filament.");
        sysParser.addComment(" - Even if bindings on the same filament is allowed, pairwise bindings will not happen on the same or neighboring cylinders.");
        sysParser.addSingleArg(
            "allow-same-filament-pair-binding",
            [](auto&& conf) -> auto& { return conf.chemParams.allowPairBindingOnSameFilament; }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Dissipation tracking");
        sysParser.addStringArgsWithAliases(
            "DISSIPATIONTRACKING", { "DISSIPATIONTRACKING:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                    "There was an error parsing input file at Chemistry algorithm. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.dissTracking = (lineVector[1] == "ON");
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { sc.chemParams.dissTracking ? "ON" : "OFF" };
            }
        );
        sysParser.addStringArgsWithAliases(
            "EVENTTRACKING", { "EVENTTRACKING:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() != 2) {
                    cout <<
                    "There was an error parsing input file at Chemistry algorithm. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2) {
                    sc.chemParams.eventTracking = (lineVector[1] == "ON");
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { sc.chemParams.eventTracking ? "ON" : "OFF" };
            }
        );
        sysParser.addEmptyLine();
    }

    //--------------------------------------------------------------------------
    // Dynamic rate parameters.
    //--------------------------------------------------------------------------
    {
        sysParser.addComment("==================================================");
        sysParser.addComment(" Dynamic rate parameters");
        sysParser.addComment("==================================================");
        sysParser.addEmptyLine();

        sysParser.addComment(" Brownian ratchet model");
        sysParser.addComment(" - Footer et al, PNAS, 2007");
        sysParser.addStringArgsWithAliases(
            "DFPOLYMERIZATIONTYPE", { "DFPOLYMERIZATIONTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.dyRateParams.dynamicRateType.dFPolymerizationType.push_back(lineVector[i]);
                }
            },
            [] (const SimulConfig& sc) {
                return sc.dyRateParams.dynamicRateType.dFPolymerizationType;
            }
        );
        sysParser.addStringArgsWithAliases(
            "DFPOLYMERIZATIONLEN", { "DFPOLYMERIZATIONLEN:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.dyRateParams.dFilPolymerizationCharLength.push_back(
                                atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& x : sc.dyRateParams.dFilPolymerizationCharLength) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Motor catch-bond model");
        sysParser.addComment(" - Erdmann et al, JCP, 2013");
        sysParser.addStringArgsWithAliases(
            "DMUNBINDINGTYPE", { "DMUNBINDINGTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.dyRateParams.dynamicRateType.dMUnbindingType.push_back(lineVector[i]);
                }
            },
            [] (const SimulConfig& sc) {
                return sc.dyRateParams.dynamicRateType.dMUnbindingType;
            }
        );
        sysParser.addStringArgsWithAliases(
            "DMUNBINDINGFORCE", { "DMUNBINDINGFORCE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.dyRateParams.dMotorUnbindingCharForce.push_back(
                                atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& x : sc.dyRateParams.dMotorUnbindingCharForce) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Motor walking");
        sysParser.addComment(" - Komianos & Papoian, PRX, 2018");
        sysParser.addComment(" - Tunable parameters based on different studies");
        sysParser.addComment(" - recommended values 24pN - 100 pN");
        sysParser.addStringArgsWithAliases(
            "DMWALKINGTYPE", { "DMWALKINGTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.dyRateParams.dynamicRateType.dMWalkingType.push_back(lineVector[i]);
                }
            },
            [] (const SimulConfig& sc) {
                return sc.dyRateParams.dynamicRateType.dMWalkingType;
            }
        );
        sysParser.addStringArgsWithAliases(
            "DMWALKINGFORCE", { "DMWALKINGFORCE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.dyRateParams.dMotorWalkingCharForce.push_back(
                                atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& x : sc.dyRateParams.dMotorWalkingCharForce) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Linker slip-bond model");
        sysParser.addComment(" - Ferrer et al, PNAS, 2008");
        sysParser.addStringArgsWithAliases(
            "DLUNBINDINGTYPE", { "DLUNBINDINGTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.dyRateParams.dynamicRateType.dLUnbindingType.push_back(lineVector[i]);
                }
            },
            [] (const SimulConfig& sc) {
                return sc.dyRateParams.dynamicRateType.dLUnbindingType;
            }
        );
        sysParser.addStringArgsWithAliases(
            "DLUNBINDINGLEN", { "DLUNBINDINGLEN:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.dyRateParams.dLinkerUnbindingCharLength.push_back(
                                atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& x : sc.dyRateParams.dLinkerUnbindingCharLength) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "DLUNBINDINGAMP", { "DLUNBINDINGAMP:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.dyRateParams.dLinkerUnbindingAmplitude.push_back(
                                atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& x : sc.dyRateParams.dLinkerUnbindingAmplitude) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Brancher slip-bond model");
        sysParser.addComment(" - Source unclear.");
        sysParser.addStringArgsWithAliases(
            "DBUNBINDINGTYPE", { "DBUNBINDINGTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.dyRateParams.dynamicRateType.dBUnbindingType.push_back(lineVector[i]);
                }
            },
            [] (const SimulConfig& sc) {
                return sc.dyRateParams.dynamicRateType.dBUnbindingType;
            }
        );
        sysParser.addStringArgsWithAliases(
            "DBUNBINDINGLEN", { "DBUNBINDINGLEN:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.dyRateParams.dBranchUnbindingCharLength.push_back(
                                atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& x : sc.dyRateParams.dBranchUnbindingCharLength) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "DBUNBINDINGF", { "DBUNBINDINGF:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    for(int i = 1; i < lineVector.size(); i++)
                        sc.dyRateParams.dBranchUnbindingCharForce.push_back(
                                atof((lineVector[i].c_str())));
                }
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                for(const auto& x : sc.dyRateParams.dBranchUnbindingCharForce) {
                    res.push_back(toString(x));
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        /// Manual Rate Changer
        // It takes 5 inputs as start_time, plusend_poly, plusend_depoly, minusend_poly, minusend_depoly
        // Currently it applies type 0 to all filament types
        sysParser.addComment(" Manual rate changer.");
        sysParser.addComment(" - Do not change unless you know what you are doing.");
        sysParser.addEmptyLine();
        sysParser.addStringArgsWithAliases(
            "MANUALSTARTTIME", { "MANUALSTARTTIME:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    sc.dyRateParams.manualCharStartTime = atof((lineVector[1].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.dyRateParams.manualCharStartTime) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "MANUALPLUSPOLYRATIO", { "MANUALPLUSPOLYRATIO:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    sc.dyRateParams.manualPlusPolyRate = atof((lineVector[1].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.dyRateParams.manualPlusPolyRate) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "MANUALPLUSDEPOLYRATIO", { "MANUALPLUSDEPOLYRATIO:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    sc.dyRateParams.manualPlusDepolyRate = atof((lineVector[1].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.dyRateParams.manualPlusDepolyRate) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "MANUALMINUSPOLYRATIO", { "MANUALMINUSPOLYRATIO:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    sc.dyRateParams.manualMinusPolyRate = atof((lineVector[1].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.dyRateParams.manualMinusPolyRate) };
            }
        );
        sysParser.addStringArgsWithAliases(
            "MANUALMINUSDEPOLYRATIO", { "MANUALMINUSDEPOLYRATIO:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if (lineVector.size() >= 2) {
                    sc.dyRateParams.manualMinusDepolyRate = atof((lineVector[1].c_str()));
                }
            },
            [] (const SimulConfig& sc) {
                return vector<string> { toString(sc.dyRateParams.manualMinusDepolyRate) };
            }
        );
        sysParser.addEmptyLine();
    }

    //--------------------------------------------------------------------------
    // System initializations.
    //--------------------------------------------------------------------------
    {
        sysParser.addComment("==================================================");
        sysParser.addComment(" Initial network setup");
        sysParser.addComment("==================================================");
        sysParser.addEmptyLine();

        sysParser.addComment(" Initial bubble setup from external input");
        sysParser.addStringArgsWithAliases(
            "BUBBLEFILE", { "BUBBLEFILE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "Error reading bubble input file. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2)
                    sc.bubbleSetup.inputFile = lineVector[1];
                else {}
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(!sc.bubbleSetup.inputFile.empty()) {
                    res.push_back( sc.bubbleSetup.inputFile.string() );
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Generate initial bubbles");
        sysParser.addStringArgsWithAliases(
            "NUMBUBBLES", { "NUMBUBBLES:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "Error reading number of bubbles. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2)
                    sc.bubbleSetup.numBubbles = atoi(lineVector[1].c_str());
                else {}
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(sc.bubbleSetup.numBubbles) {
                    res.push_back( toString(sc.bubbleSetup.numBubbles) );
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "BUBBLETYPE", { "BUBBLETYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "Error reading bubble type. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2)
                    sc.bubbleSetup.bubbleType = atoi(lineVector[1].c_str());
                else {}
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(sc.bubbleSetup.bubbleType) {
                    res.push_back( toString(sc.bubbleSetup.bubbleType) );
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" MTOC/AFM initialization.");
        sysParser.addArgs(
            "init-mtoc",
            [](SimulConfig& conf, const SExpr::ListType& sel) {
                MTOCInit bi;
                for(int i = 1; i < sel.size(); ++i) {
                    parseKeyValue<MTOCInit>(bi, sel[i], mtocInitParser().dict, KeyValueParserUnknownKeyAction::error);
                }
                conf.mtocSettings.initVec.push_back(move(bi));
            },
            [](const SimulConfig& conf) {
                vector<list<SExprToken>> res;
                for(auto& eachInit : conf.mtocSettings.initVec) {
                    list<SExprToken> tokenList;
                    // Add all setup items.
                    tokenList.splice(tokenList.end(), buildTokens<MTOCInit>(eachInit, mtocInitParser()));
                    res.push_back(move(tokenList));
                }
                return res;
            }
        );
        sysParser.addArgs(
            "init-afm",
            [](SimulConfig& conf, const SExpr::ListType& sel) {
                AFMInit bi;
                for(int i = 1; i < sel.size(); ++i) {
                    parseKeyValue<AFMInit>(bi, sel[i], afmInitParser().dict, KeyValueParserUnknownKeyAction::error);
                }
                conf.afmSettings.initVec.push_back(move(bi));
            },
            [](const SimulConfig& conf) {
                vector<list<SExprToken>> res;
                for(auto& eachInit : conf.afmSettings.initVec) {
                    list<SExprToken> tokenList;
                    // Add all setup items.
                    tokenList.splice(tokenList.end(), buildTokens<AFMInit>(eachInit, afmInitParser()));
                    res.push_back(move(tokenList));
                }
                return res;
            }
        );

        sysParser.addComment(" Membrane initialization.");
        sysParser.addArgs(
            "init-membrane",
            [](SimulConfig& conf, const SExpr::ListType& sel) {
                if(sel.size() < 2) {
                    log::error("Membrane initialization is incorrect.");
                    log::info("Membrane initialization usage: (init-membrane <setup-name> [<item1> ...]");
                    throw runtime_error("Membrane initialization is incorrect.");
                }
                else {
                    MembraneInit mi;
                    mi.name = get<SExpr::StringType>(sel[1].data);
                    for(int i = 2; i < sel.size(); ++i) {
                        parseKeyValue<MembraneInit>(mi, sel[i], membraneInitParser().dict, KeyValueParserUnknownKeyAction::error);
                    }
                    conf.membraneSettings.initVec.push_back(move(mi));
                }
            },
            [](const SimulConfig& conf) {
                vector<list<SExprToken>> res;
                for(auto& eachInit : conf.membraneSettings.initVec) {
                    list<SExprToken> tokenList;
                    // Add setup name.
                    tokenList.push_back(SExprToken::makeString(eachInit.name));
                    // Add all setup items.
                    tokenList.splice(tokenList.end(), buildTokens<MembraneInit>(eachInit, membraneInitParser()));
                    res.push_back(move(tokenList));
                }
                return res;
            }
        );

        sysParser.addComment(" Fixed vertex attachment initialization.");
        sysParser.addComment(" - Each attachment is initialized using \"coord_xyz..., search_range, k_stretch\".");
        sysParser.addStringArgs(
            "init-fixed-vertex-attachment",
            [](SimulConfig& conf, const vector<string>& lineVector) {
                const int partsize = 5;
                if(lineVector.size() % partsize == 1) {
                    // Clear original.
                    auto& att = conf.membraneSettings.fixedVertexAttachmentInits.attachments;
                    att.clear();
                    for(int i = 1; i < lineVector.size(); i += 5) {
                        att.push_back({
                            {
                                parse<FP>(lineVector[i + 0]),
                                parse<FP>(lineVector[i + 1]),
                                parse<FP>(lineVector[i + 2]),
                            },
                            parse<FP>(lineVector[i + 3]),
                            parse<FP>(lineVector[i + 4]),
                        });
                    }
                }
                else {
                    log::error("Fixed vertex attachment initialization is incorrect.");
                    throw runtime_error("Fixed vertex attachment initialization is incorrect.");
                }
            },
            [](const SimulConfig& conf) {
                vector<string> res;
                for(auto& eachAtt : conf.membraneSettings.fixedVertexAttachmentInits.attachments) {
                    res.push_back(toString(eachAtt.coord[0]));
                    res.push_back(toString(eachAtt.coord[1]));
                    res.push_back(toString(eachAtt.coord[2]));
                    res.push_back(toString(eachAtt.range));
                    res.push_back(toString(eachAtt.kStretch));
                }
                return res;
            }
        );

        sysParser.addComment(" Initial filament setup from external input");
        sysParser.addStringArgsWithAliases(
            "FILAMENTFILE", { "FILAMENTFILE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "Error reading filament input file. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2)
                    sc.filamentSetup.inputFile = lineVector[1];
                else {}
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(!sc.filamentSetup.inputFile.empty()) {
                    res.push_back( sc.filamentSetup.inputFile.string() );
                }
                return res;
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Generate initial filaments");
        sysParser.addStringArgsWithAliases(
            "NUMFILAMENTS", { "NUMFILAMENTS:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "Error reading number of filaments. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2)
                    sc.filamentSetup.numFilaments = atoi(lineVector[1].c_str());
                else {}
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(sc.filamentSetup.numFilaments) {
                    res.push_back( toString(sc.filamentSetup.numFilaments) );
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "FILAMENTLENGTH", { "FILAMENTLENGTH:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "Error reading filament length. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2)
                    sc.filamentSetup.filamentLength = atoi(lineVector[1].c_str());
                else {}
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(sc.filamentSetup.filamentLength) {
                    res.push_back( toString(sc.filamentSetup.filamentLength) );
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "FILAMENTTYPE", { "FILAMENTTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "Error reading filament type. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2)
                    sc.filamentSetup.filamentType = atoi(lineVector[1].c_str());
                else {}
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(sc.filamentSetup.filamentType) {
                    res.push_back( toString(sc.filamentSetup.filamentType) );
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "PROJECTIONTYPE", { "PROJECTIONTYPE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "Error reading filament projection type. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2)
                    sc.filamentSetup.projectionType = lineVector[1];
                else {}
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(!sc.filamentSetup.projectionType.empty()) {
                    res.push_back( sc.filamentSetup.projectionType );
                }
                return res;
            }
        );
        sysParser.addComment(" In cylindrical boundary quasi 2D systems, allow angles between cylinder and cylindrical boundary to be within -20 and 20 degrees.");
        sysParser.addStringArgs(
            "restrict-cylinder-boundary-angle",
            [] (SimulConfig& conf, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    log::error("Error specifying restrict-cylinder-boundary-angle.");
                    log::info("Usage: (restrict-cylinder-boundary-angle <true/false>)");
                    throw runtime_error("Error specifying restrict-cylinder-boundary-angle.");
                }
                else if (lineVector.size() == 2)
                    conf.filamentSetup.restrictCylinderBoundaryAngle = lineVector[1] == "true";
                else {}
            },
            [](const SimulConfig& conf) {
                return vector<string> { conf.filamentSetup.restrictCylinderBoundaryAngle ? "true" : "false" };
            }
        );
        sysParser.addComment(" In cylindrical boundary quasi 2D systems, make filaments only initialize in a ring region.");
        sysParser.addStringArgs(
            "form-rings",
            [](SimulConfig& conf, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    log::error("Error specifying form-rings.");
                    log::info("Usage: (form-rings <true/false>)");
                    throw runtime_error("Error specifying form-rings.");
                }
                else if (lineVector.size() == 2)
                    conf.filamentSetup.formRings = lineVector[1] == "true";
                else {}
            },
            [](const SimulConfig& conf) {
                return vector<string> { conf.filamentSetup.formRings ? "true" : "false" };
            }
        );
        sysParser.addEmptyLine();

        sysParser.addComment(" Restart settings");
        sysParser.addStringArgsWithAliases(
            "RESTARTPHASE", { "RESTARTPHASE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                sc.metaParams.isRestart = true;
                if(lineVector.size() > 2) {
                    log::error("Error reading restart params. Exiting.");
                    throw runtime_error("Error reading restart params.");
                }
                else if (lineVector.size() == 2){
                    if(lineVector[1].find("USECHEMCOPYNUM"))
                        sc.filamentSetup.USECHEMCOPYNUM = true;
                }

                // Also sets the system run state. Future: this should be moved to post processing.
                SysParams::RUNSTATE = false;
            },
            [] (const SimulConfig& sc) {
                vector<vector<string>> res;
                if(sc.metaParams.isRestart) {
                    vector<string> line;
                    if(sc.filamentSetup.USECHEMCOPYNUM) {
                        line.push_back("USECHEMCOPYNUM");
                    }
                    res.push_back(move(line));
                }
                return res;
            }
        );
        sysParser.addStringArgsWithAliases(
            "PINRESTARTFILE", { "PINRESTARTFILE:" },
            [] (SimulConfig& sc, const vector<string>& lineVector) {
                if(lineVector.size() > 2) {
                    cout << "Error reading filament projection type. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 2)
                    sc.filamentSetup.pinRestartFile = lineVector[1];
                else {}
            },
            [] (const SimulConfig& sc) {
                vector<string> res;
                if(!sc.filamentSetup.pinRestartFile.empty()) {
                    res.push_back(sc.filamentSetup.pinRestartFile);
                }
                return res;
            }
        );
        sysParser.addEmptyLine();
    }

    //--------------------------------------------------------------------------
    // Various runtime configuration options.
    //--------------------------------------------------------------------------
    {
        sysParser.addComment("==================================================");
        sysParser.addComment(" Runtime configuration options.");
        sysParser.addComment("==================================================");
        sysParser.addEmptyLine();

        sysParser.addComment(" Output settings.");
        sysParser.addArgs(
            "output",
            [](SimulConfig& conf, const SExpr::ListType& sel) {
                for(int i = 1; i < sel.size(); ++i) {
                    parseKeyValue<OutputParams>(conf.outputParams, sel[i], outputParamsParser().dict, KeyValueParserUnknownKeyAction::error);
                }
            },
            [](const SimulConfig& conf) {
                return buildTokens<OutputParams>(conf.outputParams, outputParamsParser());
            }
        );
        sysParser.addEmptyLine();
    }

    return sysParser;
}


FilamentData FilamentParser::readFilaments(std::istream& is) {
    is.clear();
    is.seekg(0);

    FilamentData fd;
    std::string line;
    
    while(getline(is, line)) {

        if(line.find("#") != string::npos) { continue; }

        vector<string> lineVector = split<string>(line);
        if(lineVector.size() >= 8) {
            if(lineVector[0]=="FILAMENT"){
                auto type = parse<int>(lineVector[1]);
                std::vector<Vec<3, FP>> coords;
                for(int i = 2; i < lineVector.size(); i += 3) {
                    if(i + 3 > lineVector.size()) {
                        log::error("In reading filament coordinates, the number of coordinates is not divisible by 3.");
                        throw std::runtime_error("Error reading filament coordinates.");
                    }
                    coords.push_back({
                        parse<FP>(lineVector[i]),
                        parse<FP>(lineVector[i+1]),
                        parse<FP>(lineVector[i+2]),
                    });
                }

                fd.filaments.push_back({ type, std::move(coords) });
            }
        }
    }

    return fd;
}

vector<MembraneParser::MembraneInfo> MembraneParser::readMembranes(std::istream& is) {

    is.clear();
    is.seekg(0);

    vector<MembraneInfo> res;
    
    bool wasEmpty = true;
    int stage; // 0: init, number of vertices; 1: vertex coordinate; 2: triangle vertex index
    size_t numVertices;
    
    string line;
    while(getline(is, line)) {
        
        bool isEmpty = (line.empty() || line.find("#") != string::npos); // empty line or line with '#'

        if(wasEmpty && !isEmpty) { // New membrane
            res.emplace_back();
            stage = 0;
        }

        wasEmpty = isEmpty;
        if(isEmpty) continue;

        auto& activeMem = res.back();

        vector<string> lineVector = split<string>(line);
        size_t lineVectorSize = lineVector.size();

        switch(stage) {
        case 0: // init, number of vertices
            if(lineVectorSize != 1) cout << "First line of membrane should be number of vertices." << endl;
            numVertices = atoi(lineVector[0].c_str());
            activeMem.vertexCoordinateList.reserve(numVertices);
            stage = 1;
            break;
        case 1: // vertex coordinate
            {
                if(lineVectorSize != 3) {
                    cout << "Each vertex should have 3 coordinates" << endl;
                }

                activeMem.vertexCoordinateList.emplace_back();
                auto& activeCoord = activeMem.vertexCoordinateList.back();

                // Parse coordinate information
                for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
                    activeCoord[coordIdx] = atof(lineVector[coordIdx].c_str());
                }
            }

            if(activeMem.vertexCoordinateList.size() == numVertices) stage = 2;
            break;
        case 2: // triangle vertex indices
            {
                if(lineVectorSize != 3) {
                    cout << "Each triangle should have 3 indices" << endl;
                }

                activeMem.triangleVertexIndexList.emplace_back();
                auto& activeTriangle = activeMem.triangleVertexIndexList.back();

                for(size_t i = 0; i < 3; ++i)
                    activeTriangle[i] = atoi(lineVector[i].c_str());
            }
            break;
        }
    }

    return res;
}

BubbleData BubbleParser::readBubbles(std::istream& is) {
    
    is.clear();
    is.seekg(0);
    
    BubbleData ret;
    string line;

    while(getline(is, line)) {

        if(line.find("#") != string::npos) { continue; }

        vector<string> lineVector = split<string>(line);
        if(lineVector.size() == 5) {
            auto& newBubble = ret.bubbles.emplace_back();
            parse(newBubble.type, lineVector[1]);
            parse(newBubble.coord[0], lineVector[2]);
            parse(newBubble.coord[1], lineVector[3]);
            parse(newBubble.coord[2], lineVector[4]);
        }
    }
    return ret;
}


KeyValueParser<SimulConfig> buildChemDataParser() {
    using namespace std;

    KeyValueParser<SimulConfig> chemDataParser;

    chemDataParser.addComment("============================ SPECIES ============================");
    chemDataParser.addEmptyLine();

    chemDataParser.addStringArgs(
        "species-general",
        [](SimulConfig& conf, const vector<string>& lineVector) {
            if(lineVector.size() != 3) {
                log::error("Error reading species-general. Exiting.");
                log::info("Usage: (species-general <name> <rspecies-type>)");
                throw std::runtime_error("Error reading species-general.");
            }
            conf.chemistryData.speciesGeneral.push_back({
                lineVector[1],
                lineVector[2],
            });
        },
        [] (const SimulConfig& conf) {
            vector<vector<string>> res;
            for(auto& sinfo : conf.chemistryData.speciesGeneral) {
                res.push_back({
                    sinfo.name,
                    sinfo.rspeciesType,
                });
            }
            return res;
        }
    );

    chemDataParser.addStringArgsWithAliases(
        "SPECIESBULK", { "SPECIESBULK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() !=  6 && lineVector.size() !=  8) {
                log::error("Error reading a bulk species. Exiting.");
                throw std::runtime_error("Error reading a bulk species.");
            }
            else if (lineVector.size() == 6) {

                if(lineVector[5] != "CONST" && lineVector[5] != "REG") {

                    log::error("Rspecies type \"{}\" for bulk species not valid. Exiting.", lineVector[5]);
                    throw std::runtime_error("Rspecies type for bulk species not valid.");
                }

                sc.chemistryData.speciesBulk.push_back({
                    lineVector[1], atoi(lineVector[2].c_str()),
                                atof(lineVector[3].c_str()), atof(lineVector[4].c_str()),
                                lineVector[5], "NONE", 0.0
                });
            }
            else if (lineVector.size() == 8) {

                if(lineVector[5] != "CONST" && lineVector[5] != "REG") {

                    cout << "Option for bulk species not valid. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }

                sc.chemistryData.speciesBulk.push_back({
                    lineVector[1], atoi(lineVector[2].c_str()),
                                atof(lineVector[3].c_str()), atof(lineVector[4].c_str()),
                                lineVector[5],lineVector[6], atof(lineVector[7].c_str())
                });
            }
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            for(const auto& sb : sc.chemistryData.speciesBulk) {
                // WTF is 5? Just use named members!
                if(sb.copyNumberManipulationType == "NONE") {
                    res.push_back({
                        sb.name,
                        toString(sb.initialCopyNumber),
                        toString(sb.releaseTime),
                        toString(sb.removalTime),
                        sb.rspeciesType
                    });
                } else {
                    res.push_back({
                        sb.name,
                        toString(sb.initialCopyNumber),
                        toString(sb.releaseTime),
                        toString(sb.removalTime),
                        sb.rspeciesType,
                        sb.copyNumberManipulationType,
                        toString(sb.holdMolarity)
                    });
                }
            }
            return res;
        }
    );
    chemDataParser.addEmptyLine();

    chemDataParser.addStringArgsWithAliases(
        "SPECIESDIFFUSING", { "SPECIESDIFFUSING:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() >  9 || lineVector.size() < 7) {
                cout << "Error reading a diffusing species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 8) {

                if(lineVector[6] != "AVG") {

                    cout << "Too many arguments for a non AVG-qualified diffusing "
                            "species. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }

                if(lineVector[6] == "AVG")
                    sc.chemistryData.speciesDiffusing.push_back({
                        lineVector[1], atoi(lineVector[2].c_str()),
                     atof(lineVector[3].c_str()), atof(lineVector[4].c_str()),
                     atof(lineVector[5].c_str()), lineVector[6], atoi(lineVector[7].c_str
                             ()),"NONE", 0.0
                    });
            }
            else if (lineVector.size() == 7) {

                if(lineVector[6] != "REG") {

                    cout << "Not enough arguments for a non REG-qualified diffusing species. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }

                sc.chemistryData.speciesDiffusing.push_back({
                    lineVector[1], atoi(lineVector[2].c_str()),
                     atof(lineVector[3].c_str()), atof(lineVector[4].c_str()),
                     atof(lineVector[5].c_str()), lineVector[6], 0, "NONE", 0.0
                });
            }
            else if (lineVector.size() == 9) {

                if(lineVector[6] != "REG") {

                    cout << "Not enough arguments for a non REG-qualified diffusing species. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }

                sc.chemistryData.speciesDiffusing.push_back({
                    lineVector[1], atoi(lineVector[2].c_str()),
                                                         atof(lineVector[3].c_str()), atof(lineVector[4].c_str()),
                                                         atof(lineVector[5].c_str()),
                                                         lineVector[6], 0, lineVector[7],
                                                         atof(lineVector[8].c_str())
                });
            }
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            for(const auto& sd : sc.chemistryData.speciesDiffusing) {
                if(sd.rspeciesType == "AVG") {
                    res.push_back({
                        sd.name,
                        toString(sd.initialCopyNumber),
                        toString(sd.diffusionCoefficient),
                        toString(sd.releaseTime),
                        toString(sd.removalTime),
                        sd.rspeciesType,
                        toString(sd.numEvents)
                    });
                } else if(sd.copyNumberManipulationType == "NONE") {
                    res.push_back({
                        sd.name,
                        toString(sd.initialCopyNumber),
                        toString(sd.diffusionCoefficient),
                        toString(sd.releaseTime),
                        toString(sd.removalTime),
                        sd.rspeciesType
                    });
                } else {
                    res.push_back({
                        sd.name,
                        toString(sd.initialCopyNumber),
                        toString(sd.diffusionCoefficient),
                        toString(sd.releaseTime),
                        toString(sd.removalTime),
                        sd.rspeciesType,
                        sd.copyNumberManipulationType,
                        toString(sd.holdMolarity)
                    });
                }
            }
            return res;
        }
    );
    chemDataParser.addComment("Membrane surface diffusing species.");
    chemDataParser.addComment("- Usage: (species-membrane-diffusing <name> <diffusion-coeff> <area>)");
    chemDataParser.addComment("  - diffusion-coeff: diffusion coefficient of the species, in nm^2/s.");
    chemDataParser.addComment("  - area: projected area of this protein on the membrane, in nm^2.");
    chemDataParser.addStringArgs(
        "species-membrane-diffusing",
        [](SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() == 4) {
                ChemistryData::SpeciesMembraneDiffusingInfo sinfo;
                parse(sinfo.name, lineVector[1]);
                parse(sinfo.diffusionCoefficient, lineVector[2]);
                parse(sinfo.area, lineVector[3]);
                sc.chemistryData.speciesMembraneDiffusing.push_back(move(sinfo));
            }
            else {
                log::error("Usage: (species-membrane-diffusing <name> <diffusion-coeff> <area>)");
                throw runtime_error("Invalid paramters for membrane-diffusing.");
            }
        },
        [](const SimulConfig& sc) {
            vector<vector<string>> res;
            for(auto& sinfo : sc.chemistryData.speciesMembraneDiffusing) {
                res.push_back({
                    sinfo.name,
                    toString(sinfo.diffusionCoefficient),
                    toString(sinfo.area),
                });
            }
            return res;
        }
    );
    chemDataParser.addEmptyLine();

    chemDataParser.addStringArgsWithAliases(
        "SPECIESFILAMENT", { "SPECIESFILAMENT:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() !=  3) {
                cout << "Error reading a filament species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {

                sc.chemistryData.speciesFilament[atoi(lineVector[2].c_str())].push_back(lineVector[1]);

            }
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& sf = sc.chemistryData.speciesFilament;
            for(int i = 0; i < sf.size(); ++i) {
                for(const auto& eachS : sf[i]) {
                    res.push_back({ eachS, toString(i) });
                }
            }
            return res;
        }
    );
    chemDataParser.addEmptyLine();

    chemDataParser.addStringArgsWithAliases(
        "SPECIESBOUND", { "SPECIESBOUND:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() !=  3) {
                cout << "Error reading a filament bound species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {

                sc.chemistryData.speciesBound[atoi(lineVector[2].c_str())].push_back(lineVector[1]);

            }
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& sb = sc.chemistryData.speciesBound;
            for(int i = 0; i < sb.size(); ++i) {
                for(const auto& eachS : sb[i]) {
                    res.push_back({ eachS, toString(i) });
                }
            }
            return res;
        }
    );
    chemDataParser.addEmptyLine();

    chemDataParser.addStringArgsWithAliases(
        "SPECIESLINKER", { "SPECIESLINKER:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() !=  3) {
                cout << "Error reading a filament linker species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {

                sc.chemistryData.speciesLinker[atoi(lineVector[2].c_str())].push_back(lineVector[1]);

            }
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& ss = sc.chemistryData.speciesLinker;
            for(int i = 0; i < ss.size(); ++i) {
                for(const auto& eachS : ss[i]) {
                    res.push_back({ eachS, toString(i) });
                }
            }
            return res;
        }
    );
    chemDataParser.addStringArgsWithAliases(
        "SPECIESMOTOR", { "SPECIESMOTOR:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() !=  3) {
                cout << "Error reading a filament motor species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {

                sc.chemistryData.speciesMotor[atoi(lineVector[2].c_str())].push_back(lineVector[1]);

            }
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& ss = sc.chemistryData.speciesMotor;
            for(int i = 0; i < ss.size(); ++i) {
                for(const auto& eachS : ss[i]) {
                    res.push_back({ eachS, toString(i) });
                }
            }
            return res;
        }
    );
    chemDataParser.addStringArgsWithAliases(
        "SPECIESBRANCHER", { "SPECIESBRANCHER:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() !=  3) {
                cout << "Error reading a filament brancher species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {

                sc.chemistryData.speciesBrancher[atoi(lineVector[2].c_str())].push_back(lineVector[1]);

            }
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& ss = sc.chemistryData.speciesBrancher;
            for(int i = 0; i < ss.size(); ++i) {
                for(const auto& eachS : ss[i]) {
                    res.push_back({ eachS, toString(i) });
                }
            }
            return res;
        }
    );
    chemDataParser.addEmptyLine();

    chemDataParser.addStringArgsWithAliases(
        "SPECIESPLUSEND", { "SPECIESPLUSEND:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() !=  3) {
                cout << "Error reading a filament plus end species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {

                sc.chemistryData.speciesPlusEnd[atoi(lineVector[2].c_str())].push_back(lineVector[1]);

            }
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& ss = sc.chemistryData.speciesPlusEnd;
            for(int i = 0; i < ss.size(); ++i) {
                for(const auto& eachS : ss[i]) {
                    res.push_back({ eachS, toString(i) });
                }
            }
            return res;
        }
    );
    chemDataParser.addStringArgsWithAliases(
        "SPECIESMINUSEND", { "SPECIESMINUSEND:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() !=  3) {
                cout << "Error reading a filament minus end species. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {

                sc.chemistryData.speciesMinusEnd[atoi(lineVector[2].c_str())].push_back(lineVector[1]);

            }
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& ss = sc.chemistryData.speciesMinusEnd;
            for(int i = 0; i < ss.size(); ++i) {
                for(const auto& eachS : ss[i]) {
                    res.push_back({ eachS, toString(i) });
                }
            }
            return res;
        }
    );
    chemDataParser.addEmptyLine();

    chemDataParser.addStringArgsWithAliases(
        "BRANCHERBINDINGSITE", { "BRANCHERBINDINGSITE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() !=  3) {
                cout << "Error reading a brancher binding site. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3)
                sc.chemistryData.B_BINDING_INDEX[atoi(lineVector[2].c_str())] = lineVector[1];
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& ss = sc.chemistryData.B_BINDING_INDEX;
            for(int i = 0; i < ss.size(); ++i) {
                if (!ss[i].empty()) {
                    res.push_back({ ss[i], toString(i) });
                }
            }
            return res;
        }
    );
    chemDataParser.addStringArgsWithAliases(
        "LINKERBINDINGSITE", { "LINKERBINDINGSITE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() !=  3) {
                cout << "Error reading a linker binding site. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3)
                sc.chemistryData.L_BINDING_INDEX[atoi(lineVector[2].c_str())] = lineVector[1];
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& ss = sc.chemistryData.L_BINDING_INDEX;
            for(int i = 0; i < ss.size(); ++i) {
                if (!ss[i].empty()) {
                    res.push_back({ ss[i], toString(i) });
                }
            }
            return res;
        }
    );
    chemDataParser.addStringArgsWithAliases(
        "MOTORBINDINGSITE", { "MOTORBINDINGSITE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() !=  3) {
                cout << "Error reading a motor binding site. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3)
                sc.chemistryData.M_BINDING_INDEX[atoi(lineVector[2].c_str())] = lineVector[1];
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& ss = sc.chemistryData.M_BINDING_INDEX;
            for(int i = 0; i < ss.size(); ++i) {
                if (!ss[i].empty()) {
                    res.push_back({ ss[i], toString(i) });
                }
            }
            return res;
        }
    );
    chemDataParser.addEmptyLine();

    chemDataParser.addComment("============================ REACTIONS ============================");
    chemDataParser.addEmptyLine();

    chemDataParser.addComment(" General reactions");
    chemDataParser.addStringArgsWithAliases(
        "GENREACTION", { "GENREACTION:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            vector<string> reactants;
            vector<string> products;

            // check parameters related to dissipation tracking if it is enabled
            float gnum = 0.0;
            int dissOffSet = 0;
            string HRCDID = "NA";
            if(sc.chemParams.dissTracking){
                string dissString = lineVector[1];
                istringstream iss(dissString);
                string token;
                vector<string> dissTokens;
                dissOffSet = 1;

                while (std::getline(iss, token, ':')) {
                    if (!token.empty())
                        dissTokens.push_back(token);
                }

                gnum = atof(dissTokens[0].c_str());

                if(dissTokens.size()!=1){
                    HRCDID = dissTokens[1];
                } else {
                    HRCDID = "NA";
                }
            }


            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {
                for(auto it  = lineVector.begin() + 1 + dissOffSet; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }
                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }

                sc.chemistryData.genReactions.push_back(
                        tuple<vector<string>, vector<string>, floatingpoint, floatingpoint,
                        string>(reactants, products, atof(lineVector[lineVector.size() - 1].c_str()),gnum,HRCDID));
                
            }
            else {
                cout << "Error reading a general reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            for(const auto& r : sc.chemistryData.genReactions) {
                vector<string> line;
                if(sc.chemParams.dissTracking) {
                    line.push_back(toString(get<3>(r)) + ":" + get<4>(r));
                }

                const auto& reactants = get<0>(r);
                const auto& products = get<1>(r);
                for(int i = 0; i < reactants.size(); ++i) {
                    if(i) { line.push_back("+"); line.push_back(reactants[i]); }
                    else  { line.push_back(reactants[i]); }
                }
                line.push_back("->");
                for(int i = 0; i < products.size(); ++i) {
                    if(i) { line.push_back("+"); line.push_back(products[i]); }
                    else  { line.push_back(products[i]); }
                }

                line.push_back(toString(get<2>(r)));

                res.push_back(move(line));
            }
            return res;
        }
    );
    chemDataParser.addEmptyLine();

    chemDataParser.addComment(" General reactions on membrane surfaces.");
    chemDataParser.addStringArgs(
        "reaction-surface",
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            ChemistryData::ReactionSurfaceInfo rinfo;

            if(lineVector.size() < 3) {
                log::error("Usage: (reaction-surface <rate> [s1 [s2 [...]]] -> [s3 [s4 [...]]])");
                throw runtime_error("Invalid reaction on surface.");
            }
            const auto arrowIndex = find(lineVector.begin(), lineVector.end(), "->") - lineVector.begin();
            if(arrowIndex == lineVector.size()) {
                log::error("Cannot find \"->\" in reaction-surface.");
                throw runtime_error("Invalid reaction on surface.");
            }

            parse(rinfo.rate, lineVector[1]);
            for(int i = 2; i < arrowIndex; ++i) {
                rinfo.reactants.push_back(lineVector[i]);
            }
            for(int i = arrowIndex + 1; i < lineVector.size(); ++i) {
                rinfo.products.push_back(lineVector[i]);
            }

            sc.chemistryData.reactionsSurface.push_back(move(rinfo));                
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            for(auto& rinfo : sc.chemistryData.reactionsSurface) {
                vector<string> line;

                line.push_back(toString(rinfo.rate));
                for(auto& s : rinfo.reactants) line.push_back(s);
                line.push_back("->");
                for(auto& s : rinfo.products)  line.push_back(s);

                res.push_back(move(line));
            }
            return res;
        }
    );
    chemDataParser.addEmptyLine();

    chemDataParser.addComment(" Adsorption/desorption on membrane surface.");
    chemDataParser.addStringArgs(
        "reaction-adsorption-desorption",
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            ChemistryData::ReactionAdsorptionDesorptionInfo rinfo;

            if(lineVector.size() != 6 || lineVector[4] != "<->") {
                log::error("Usage: (reaction-adsorption-desorption <on rate> <off rate> <3D diffusing> <-> <2D diffusing>)");
                throw runtime_error("Invalid adsorption/desorption reaction.");
            }
            parse(rinfo.onRate, lineVector[1]);
            parse(rinfo.offRate, lineVector[2]);
            rinfo.speciesName3D = lineVector[3];
            rinfo.speciesName2D = lineVector[5];

            sc.chemistryData.reactionsAdsorptionDesorption.push_back(move(rinfo));
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            for(auto& rinfo : sc.chemistryData.reactionsAdsorptionDesorption) {
                vector<string> line;

                line.push_back(toString(rinfo.onRate));
                line.push_back(toString(rinfo.offRate));
                line.push_back(rinfo.speciesName3D);
                line.push_back("<->");
                line.push_back(rinfo.speciesName2D);

                res.push_back(move(line));
            }
            return res;
        }
    );
    chemDataParser.addEmptyLine();

    chemDataParser.addComment(" Bulk reactions");
    chemDataParser.addStringArgsWithAliases(
        "BULKREACTION", { "BULKREACTION:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            vector<string> reactants;
            vector<string> products;

            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {

                for(auto it  = lineVector.begin() + 1; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }

                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }

                sc.chemistryData.bulkReactions.push_back(
                tuple<vector<string>, vector<string>, floatingpoint>
                (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
                
            }
            else {
                cout << "Error reading a bulk reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            for(const auto& r : sc.chemistryData.bulkReactions) {
                vector<string> line;

                const auto& reactants = get<0>(r);
                const auto& products = get<1>(r);
                for(int i = 0; i < reactants.size(); ++i) {
                    if(i) { line.push_back("+"); line.push_back(reactants[i]); }
                    else  { line.push_back(reactants[i]); }
                }
                line.push_back("->");
                for(int i = 0; i < products.size(); ++i) {
                    if(i) { line.push_back("+"); line.push_back(products[i]); }
                    else  { line.push_back(products[i]); }
                }

                line.push_back(toString(get<2>(r)));

                res.push_back(move(line));
            }
            return res;
        }
    );
    chemDataParser.addEmptyLine();

    chemDataParser.addComment(" Filament reactions");
    chemDataParser.addStringArgsWithAliases(
        "NUCLEATIONREACTION", { "NUCLEATIONREACTION:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            vector<string> reactants;
            vector<string> products;

            int filType = atoi(lineVector[1].c_str());

            auto arrowIt = find(lineVector.begin(), lineVector.end(), "->");
            if(arrowIt != lineVector.end()) {

                for(auto it  = lineVector.begin() + 2; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }

                for(auto it = arrowIt + 1; it != lineVector.end() - 1; it++) {
                    if(*it != "+")  products.push_back((*it));
                }

                sc.chemistryData.nucleationReactions[filType].push_back(
                tuple<vector<string>, vector<string>, floatingpoint>
                (reactants, products, atof(lineVector[lineVector.size() - 1].c_str())));
            }
            else {
                cout << "Error reading a nucleation reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& rxns = sc.chemistryData.nucleationReactions;
            for(int k = 0; k < rxns.size(); ++k) {
                for(const auto& r : rxns[k]) {
                    vector<string> line;

                    line.push_back(toString(k));

                    const auto& reactants = get<0>(r);
                    const auto& products = get<1>(r);
                    for(int i = 0; i < reactants.size(); ++i) {
                        if(i) { line.push_back("+"); line.push_back(reactants[i]); }
                        else  { line.push_back(reactants[i]); }
                    }
                    line.push_back("->");
                    for(int i = 0; i < products.size(); ++i) {
                        if(i) { line.push_back("+"); line.push_back(products[i]); }
                        else  { line.push_back(products[i]); }
                    }

                    line.push_back(toString(get<2>(r)));

                    res.push_back(move(line));
                }
            }
            return res;
        }
    );
    chemDataParser.addStringArgsWithAliases(
        "DEPOLYMERIZATIONREACTION", { "DEPOLYMERIZATIONREACTION:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            vector<string> reactants;
            vector<string> products;

            // check parameters related to dissipation tracking if it is enabled
            float gnum = 0.0;
            int dissOffset = 0;
            string HRCDID = "NA";
            if(sc.chemParams.dissTracking){
                string dissString = lineVector[1].c_str();
                istringstream iss(dissString);
                string token;
                vector<string> dissTokens;
                dissOffset = 1;

                while (std::getline(iss, token, ':')) {
                    if (!token.empty())
                        dissTokens.push_back(token);
                }

                gnum = atof(dissTokens[0].c_str());

                if(dissTokens.size()!=1){
                    HRCDID = dissTokens[1];
                } else {
                    HRCDID = "NA";
                }
            }

            const int filType = parse<int>(lineVector[1 + dissOffset]);

            if(lineVector.size() == 10 + dissOffset) {

                sc.chemistryData.depolymerizationReactions[filType].push_back({
                    lineVector[2 + dissOffset],
                    lineVector[4 + dissOffset],
                    lineVector[6 + dissOffset],
                    lineVector[8 + dissOffset],
                    parse<floatingpoint>(lineVector[9 + dissOffset]),
                    gnum,
                    HRCDID,
                });
            }
            else {
                log::error("Error reading a depolymerization reaction.");
                throw std::runtime_error("Error reading depolymerization reaction.");
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& rxns = sc.chemistryData.depolymerizationReactions;
            for(int fi = 0; fi < rxns.size(); ++fi) {
                for(const auto& r : rxns[fi]) {
                    vector<string> line;

                    if(sc.chemParams.dissTracking) {
                        line.push_back(toString(r.gnum) + ":" + r.hrcdid);
                    }

                    line.push_back(toString(fi));

                    line.push_back(r.speciesReactantFilament);
                    line.push_back("+");
                    line.push_back(r.speciesReactantPlusEndMinusEnd);
                    line.push_back("->");
                    line.push_back(r.speciesProductDiffusingBulk);
                    line.push_back("+");
                    line.push_back(r.speciesProductPlusEndMinusEnd);

                    line.push_back(toString(r.rate));

                    res.push_back(move(line));
                }
            }
            return res;
        }
    );
    chemDataParser.addStringArgsWithAliases(
        "POLYMERIZATIONREACTION", { "POLYMERIZATIONREACTION:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            vector<string> reactants;
            vector<string> products;

            // check parameters related to dissipation tracking if it is enabled
            float gnum = 0.0;
            int dissOffset = 0;
            string HRCDID = "NA";
            if(sc.chemParams.dissTracking){
                string dissString = lineVector[1].c_str();
                istringstream iss(dissString);
                string token;
                vector<string> dissTokens;
                dissOffset = 1;

                while (std::getline(iss, token, ':')) {
                    if (!token.empty())
                        dissTokens.push_back(token);
                }

                gnum = atof(dissTokens[0].c_str());

                if(dissTokens.size()!=1){
                    HRCDID = dissTokens[1];
                } else {
                    HRCDID = "NA";
                }
            }

            const int filType = parse<int>(lineVector[1 + dissOffset]);

            if(lineVector.size() == 10 + dissOffset) {

                sc.chemistryData.polymerizationReactions[filType].push_back({
                    lineVector[2 + dissOffset],
                    lineVector[4 + dissOffset],
                    lineVector[6 + dissOffset],
                    lineVector[8 + dissOffset],
                    parse<floatingpoint>(lineVector[9 + dissOffset]),
                    gnum,
                    HRCDID,
                });
            }
            else {
                log::error("Error reading a polymerization reaction.");
                throw std::runtime_error("Error reading polymerization reaction.");
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& rxns = sc.chemistryData.polymerizationReactions;
            for(int fi = 0; fi < rxns.size(); ++fi) {
                for(const auto& r : rxns[fi]) {
                    vector<string> line;

                    if(sc.chemParams.dissTracking) {
                        line.push_back(toString(r.gnum) + ":" + r.hrcdid);
                    }

                    line.push_back(toString(fi));

                    line.push_back(r.speciesReactantDiffusingBulk);
                    line.push_back("+");
                    line.push_back(r.speciesReactantPlusEndMinusEnd);
                    line.push_back("->");
                    line.push_back(r.speciesProductFilament);
                    line.push_back("+");
                    line.push_back(r.speciesProductPlusEndMinusEnd);

                    line.push_back(toString(r.rate));

                    res.push_back(move(line));
                }
            }
            return res;
        }
    );
    chemDataParser.addEmptyLine();

    chemDataParser.addComment(" Linker/motor/branching reactions");
    chemDataParser.addStringArgsWithAliases(
        "LINKERREACTION", { "LINKERREACTION:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {

            // lineVector[dissOffset + 1] was for filament type but it is no longer required. It is only kept for backward compatibility.

            // check parameters related to dissipation tracking if it is enabled
            floatingpoint gnum = 0.0;
            int dissOffset = 0;
            string HRCDID = "NA";
            if(sc.chemParams.dissTracking){
                string dissString = lineVector[1].c_str();
                istringstream iss(dissString);
                string token;
                vector<string> dissTokens;
                dissOffset = 1;

                while (std::getline(iss, token, ':')) {
                    if (!token.empty())
                        dissTokens.push_back(token);
                }

                gnum = atof(dissTokens[0].c_str());

                if(dissTokens.size()!=1){
                    HRCDID = dissTokens[1];
                } else {
                    HRCDID = "NA";
                }
            }

            const int expectedSize = dissOffset + 15;
            if(expectedSize == lineVector.size()) {

                sc.chemistryData.linkerReactions.push_back({
                    // Reactants.
                    {
                        lineVector[dissOffset + 2],
                        lineVector[dissOffset + 4],
                        lineVector[dissOffset + 6],
                    },
                    // Products.
                    {
                        lineVector[dissOffset + 8],
                        lineVector[dissOffset + 10],
                    },
                    atof(lineVector[dissOffset + 11].c_str()),
                    atof(lineVector[dissOffset + 12].c_str()),
                    atof(lineVector[dissOffset + 13].c_str()),
                    atof(lineVector[dissOffset + 14].c_str()),
                    gnum, HRCDID
                });
                
            }
            else {
                cout << "Error reading a linker reaction. Exiting." << endl;
                throw runtime_error("Error parsing a linker reaction.");
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& rxns = sc.chemistryData.linkerReactions;
            for(auto& r : rxns) {
                vector<string> line;

                if(sc.chemParams.dissTracking) {
                    line.push_back(toString(r.gnum) + ":" + r.hrcdid);
                }

                // Dummy variable for filament type, for backward compatibility.
                line.push_back("0");

                line.push_back(r.reactantInfo.speciesBound1);
                line.push_back("+");
                line.push_back(r.reactantInfo.speciesBound2);
                line.push_back("+");
                line.push_back(r.reactantInfo.linkerDiffusing);
                line.push_back("<->");
                line.push_back(r.productInfo.linkerBound1);
                line.push_back("+");
                line.push_back(r.productInfo.linkerBound2);

                line.push_back(toString(r.onRate));
                line.push_back(toString(r.offRate));
                line.push_back(toString(r.rMin));
                line.push_back(toString(r.rMax));

                res.push_back(move(line));
            }
            return res;
        }
    );
    chemDataParser.addStringArgsWithAliases(
        "MOTORREACTION", { "MOTORREACTION:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {

            // lineVector[dissOffset + 1] was for filament type but it is no longer required. It is only kept for backward compatibility.

            // check parameters related to dissipation tracking if it is enabled
            floatingpoint gnum = 0.0;
            int dissOffset = 0;
            string HRCDID = "NA";
            if(sc.chemParams.dissTracking){
                string dissString = lineVector[1].c_str();
                istringstream iss(dissString);
                string token;
                vector<string> dissTokens;
                dissOffset = 1;

                while (std::getline(iss, token, ':')) {
                    if (!token.empty())
                        dissTokens.push_back(token);
                }

                gnum = atof(dissTokens[0].c_str());

                if(dissTokens.size()!=1){
                    HRCDID = dissTokens[1];
                } else {
                    HRCDID = "NA";
                }
            }

            const int expectedSize = dissOffset + 15;
            if(expectedSize == lineVector.size()) {

                sc.chemistryData.motorReactions.push_back({
                    // Reactants.
                    {
                        lineVector[dissOffset + 2],
                        lineVector[dissOffset + 4],
                        lineVector[dissOffset + 6],
                    },
                    // Products.
                    {
                        lineVector[dissOffset + 8],
                        lineVector[dissOffset + 10],
                    },
                    atof(lineVector[dissOffset + 11].c_str()),
                    atof(lineVector[dissOffset + 12].c_str()),
                    atof(lineVector[dissOffset + 13].c_str()),
                    atof(lineVector[dissOffset + 14].c_str()),
                    gnum, HRCDID
                });
                
            }
            else {
                cout << "Error reading a motor reaction. Exiting." << endl;
                throw runtime_error("Error parsing a motor reaction.");
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& rxns = sc.chemistryData.motorReactions;
            for(auto& r : rxns) {
                vector<string> line;

                if(sc.chemParams.dissTracking) {
                    line.push_back(toString(r.gnum) + ":" + r.hrcdid);
                }

                // Dummy variable for filament type, for backward compatibility.
                line.push_back("0");

                line.push_back(r.reactantInfo.speciesBound1);
                line.push_back("+");
                line.push_back(r.reactantInfo.speciesBound2);
                line.push_back("+");
                line.push_back(r.reactantInfo.linkerDiffusing);
                line.push_back("<->");
                line.push_back(r.productInfo.linkerBound1);
                line.push_back("+");
                line.push_back(r.productInfo.linkerBound2);

                line.push_back(toString(r.onRate));
                line.push_back(toString(r.offRate));
                line.push_back(toString(r.rMin));
                line.push_back(toString(r.rMax));

                res.push_back(move(line));
            }
            return res;
        }
    );
    chemDataParser.addStringArgsWithAliases(
        "MOTORWALKINGREACTION", { "MOTORWALKINGREACTION:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            vector<string> reactants;
            vector<string> products;

            // check parameters related to dissipation tracking if it is enabled
            floatingpoint gnum = 0.0;
            int dissOffset = 0;
            string HRCDID = "NA";
            if(sc.chemParams.dissTracking){
                string dissString = lineVector[1].c_str();
                istringstream iss(dissString);
                string token;
                vector<string> dissTokens;
                dissOffset = 1;

                while (std::getline(iss, token, ':')) {
                    if (!token.empty())
                        dissTokens.push_back(token);
                }

                gnum = atof(dissTokens[0].c_str());

                if(dissTokens.size()!=1){
                    HRCDID = dissTokens[1];
                } else {
                    HRCDID = "NA";
                }
            }
            const int filType = parse<int>(lineVector[1 + dissOffset]);

            if(lineVector.size() == 10 + dissOffset) {

                sc.chemistryData.motorWalkingReactions[filType].push_back({
                    lineVector[2 + dissOffset],
                    lineVector[4 + dissOffset],
                    lineVector[6 + dissOffset],
                    lineVector[8 + dissOffset],
                    parse<floatingpoint>(lineVector[9 + dissOffset]),
                    gnum,
                    HRCDID,
                });
            }
            else {
                LOG(ERROR) << "Error reading a motor walking reaction.";
                throw std::runtime_error("Error reading motor walking reaction.");
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& rxns = sc.chemistryData.motorWalkingReactions;
            for(int fi = 0; fi < rxns.size(); ++fi) {
                for(const auto& r : rxns[fi]) {
                    vector<string> line;

                    if(sc.chemParams.dissTracking) {
                        line.push_back(toString(r.gnum) + ":" + r.hrcdid);
                    }

                    line.push_back(toString(fi));

                    line.push_back(r.speciesReactantMotor);
                    line.push_back("+");
                    line.push_back(r.speciesReactantEmptySite);
                    line.push_back("->");
                    line.push_back(r.speciesProductMotor);
                    line.push_back("+");
                    line.push_back(r.speciesProductEmptySite);

                    line.push_back(toString(r.rate));

                    res.push_back(move(line));
                }
            }
            return res;
        }
    );
    chemDataParser.addStringArgsWithAliases(
        "BRANCHINGREACTION", { "BRANCHINGREACTION:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            vector<string> reactants;
            vector<string> products;

            int filType = atoi(lineVector[1].c_str());

            auto arrowIt = find(lineVector.begin(), lineVector.end(), "<->");
            if(arrowIt != lineVector.end()) {

                for(auto it  = lineVector.begin() + 2; it != arrowIt; it++) {
                    if(*it != "+") reactants.push_back((*it));
                }

                for(auto it = arrowIt + 1; it != lineVector.end() - 4; it++) {
                    if(*it != "+")  products.push_back((*it));
                }

                sc.chemistryData.branchingReactions[filType].push_back(
                tuple<vector<string>, vector<string>, floatingpoint, floatingpoint, string, floatingpoint>
                (reactants, products, atof(lineVector[lineVector.size() - 4].c_str()),
                                      atof(lineVector[lineVector.size() - 3].c_str()),
                                           lineVector[lineVector.size() - 2].c_str(),
                                      atof(lineVector[lineVector.size() - 1].c_str())));
            }
            else {
                cout << "Error reading a branching reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& rxns = sc.chemistryData.branchingReactions;
            for(int k = 0; k < rxns.size(); ++k) {
                for(const auto& r : rxns[k]) {
                    vector<string> line;

                    line.push_back(toString(k));

                    const auto& reactants = get<0>(r);
                    const auto& products = get<1>(r);
                    for(int i = 0; i < reactants.size(); ++i) {
                        if(i) { line.push_back("+"); line.push_back(reactants[i]); }
                        else  { line.push_back(reactants[i]); }
                    }
                    line.push_back("<->");
                    for(int i = 0; i < products.size(); ++i) {
                        if(i) { line.push_back("+"); line.push_back(products[i]); }
                        else  { line.push_back(products[i]); }
                    }

                    line.push_back(toString(get<2>(r)));
                    line.push_back(toString(get<3>(r)));
                    line.push_back(          get<4>(r) );
                    line.push_back(toString(get<5>(r)));

                    res.push_back(move(line));
                }
            }
            return res;
        }
    );
    chemDataParser.addEmptyLine();

    chemDataParser.addComment(" Filament aging/destruction/severing");
    chemDataParser.addStringArgsWithAliases(
        "AGINGREACTION", { "AGINGREACTION:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            vector<string> reactants;
            vector<string> products;

            // check parameters related to dissipation tracking if it is enabled 
            floatingpoint gnum = 0.0;
            int dissOffset = 0;
            string HRCDID = "NA";
            if(sc.chemParams.dissTracking){
                string dissString = lineVector[1].c_str();
                istringstream iss(dissString);
                string token;
                vector<string> dissTokens;
                dissOffset = 1;

                while (std::getline(iss, token, ':')) {
                    if (!token.empty())
                        dissTokens.push_back(token);
                }

                gnum = atof(dissTokens[0].c_str());

                if(dissTokens.size()!=1){
                    HRCDID = dissTokens[1];
                } else {
                    HRCDID = "NA";
                }
            }

            const int filType = parse<int>(lineVector[1 + dissOffset]);

            if(lineVector.size() == 6 + dissOffset) {

                sc.chemistryData.agingReactions[filType].push_back({
                    lineVector[2 + dissOffset],
                    lineVector[4 + dissOffset],
                    parse<floatingpoint>(lineVector[5 + dissOffset]),
                    gnum,
                    HRCDID,
                });
            }
            else {
                log::error("Error reading an aging reaction.");
                throw std::runtime_error("Error reading aging reaction.");
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& rxns = sc.chemistryData.agingReactions;
            for(int fi = 0; fi < rxns.size(); ++fi) {
                for(const auto& r : rxns[fi]) {
                    vector<string> line;

                    if(sc.chemParams.dissTracking) {
                        line.push_back(toString(r.gnum) + ":" + r.hrcdid);
                    }

                    line.push_back(toString(fi));

                    line.push_back(r.speciesReactant);
                    line.push_back("->");
                    line.push_back(r.speciesProduct);

                    line.push_back(toString(r.rate));

                    res.push_back(move(line));
                }
            }
            return res;
        }
    );
    chemDataParser.addStringArgsWithAliases(
        "DESTRUCTIONREACTION", { "DESTRUCTIONREACTION:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            int filType = atoi(lineVector[1].c_str());

            if(lineVector.size() == 10) {

                sc.chemistryData.destructionReactions[filType].push_back({
                    lineVector[2],
                    lineVector[4],
                    lineVector[6],
                    lineVector[8],
                    parse<floatingpoint>(lineVector[9]),
                });
            }
            else {
                log::error("Error reading a destruction reaction.");
                throw std::runtime_error("Error reading destruction reaction.");
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& rxns = sc.chemistryData.destructionReactions;
            for(int fi = 0; fi < rxns.size(); ++fi) {
                for(const auto& r : rxns[fi]) {
                    vector<string> line;

                    line.push_back(toString(fi));

                    line.push_back(r.speciesReactantPlusEnd);
                    line.push_back("+");
                    line.push_back(r.speciesReactantMinusEnd);
                    line.push_back("->");
                    line.push_back(r.speciesProduct1);
                    line.push_back("+");
                    line.push_back(r.speciesProduct2);

                    line.push_back(toString(r.rate));

                    res.push_back(move(line));
                }
            }
            return res;
        }
    );
    chemDataParser.addStringArgsWithAliases(
        "SEVERINGREACTION", { "SEVERINGREACTION:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            int filType = atoi(lineVector[1].c_str());

            auto arrowIt = find(lineVector.begin(), lineVector.end(), "AT");
            if(arrowIt != lineVector.end()) {

                auto it = arrowIt + 1;
                
                sc.chemistryData.severingReactions[filType].push_back(tuple<string, floatingpoint>
                ((*it), atof(lineVector[lineVector.size() - 1].c_str())));
            }
            else {
                cout << "Error reading a severing reaction. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            const auto& rxns = sc.chemistryData.severingReactions;
            for(int k = 0; k < rxns.size(); ++k) {
                for(const auto& r : rxns[k]) {
                    vector<string> line;

                    line.push_back(toString(k));
                    line.push_back("AT");
                    line.push_back(get<0>(r));
                    line.push_back(toString(get<1>(r)));

                    res.push_back(move(line));
                }
            }
            return res;
        }
    );
    chemDataParser.addEmptyLine();

    return chemDataParser;
}



void PinRestartParser::resetPins() {

    //loop through filaments
    for(auto &f: Filament::getFilaments()) {

        _inputFile.clear();
        _inputFile.seekg(0);

        // Get minus end bead
        auto b1 = f->getMinusEndCylinder()->getFirstBead();
        auto b2 = f->getPlusEndCylinder()->getSecondBead();

        int filID = f->getId();
        string searchID = "FILAMENT " + toString(filID) + ":";

        string line;

        while(getline(_inputFile, line)) {

            if(line.find("#") != string::npos) { continue; }

            else if(line.find(searchID) != string::npos) {

                vector<string> lineVector = split<string>(line);
                if(lineVector.size() !=  8) {
                    cout << "Error reading a restart pin position. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                else if (lineVector.size() == 8) {
                    
                    b1->pinnedPosition = vector<floatingpoint>{stof(lineVector[2].c_str()
                                                               ), stof(lineVector[3].c_str()), stof(lineVector[4].c_str())};
                    b2->pinnedPosition = vector<floatingpoint>{stof(lineVector[5].c_str()
                                                               ), stof(lineVector[6].c_str()), stof(lineVector[7].c_str())};
                    
                    if(!areEqual(b1->pinnedPosition[0],0.0) && !areEqual(b1->pinnedPosition[1],0.0) && !areEqual(b1->pinnedPosition[2],0.0)) {
                        b1->addAsPinned();

//                        cout << "Pinned filament! coordinates = " << b1->coordinate[0] << " " << b1->coordinate[1] << " " << b1->coordinate[2] << endl;
//                        cout << "Pin position = " << b1->pinnedPosition[0] << " " << b1->pinnedPosition[1] << " " << b1->pinnedPosition[2] << endl;
                    }

                    if(!areEqual(b2->pinnedPosition[0],0.0) && !areEqual(b2->pinnedPosition[1],0.0) && !areEqual(b2->pinnedPosition[2],0.0)) {
                        b2->addAsPinned();

//                        cout << "Pinned filament! coordinates = " << b2->coordinate[0] << " " << b2->coordinate[1] << " " << b2->coordinate[2] << endl;
//                        cout << "Pin position = " << b2->pinnedPosition[0] << " " << b2->pinnedPosition[1] << " " << b2->pinnedPosition[2] << endl;
                    }
                }
            }
        }
    }
}

} // namespace medyan
