#ifndef MEDYAN_Util_Parser_BubbleParser_hpp
#define MEDYAN_Util_Parser_BubbleParser_hpp

#include "Parser.h"
#include "Parameter/Bubble.hpp"
#include "Util/Parser/StringParser.hpp"

namespace medyan {

inline const auto& mtocInitParser() {
    using namespace std;

    static auto res = [] {
        KeyValueParser< MTOCInit > p;

        p.addEmptyLine();
        p.addComment(" Bubble parameters.");
        p.addSingleArg("bubble-type", [](auto&& bi) -> auto& { return bi.bubbleType; });
        p.addStringArgs(
            "bubble-coord",
            [](MTOCInit& bi, const vector<string>& input) {
                if(input.size() != 4) {
                    log::error("Expected 3 coordinates for bubble-coord.");
                    throw runtime_error("Expected 3 coordinates for bubble-coord");
                }
                else {
                    parse(bi.bubbleCoord[0], input[1]);
                    parse(bi.bubbleCoord[1], input[2]);
                    parse(bi.bubbleCoord[2], input[3]);
                }
            },
            [](const MTOCInit& bi) {
                return vector<string> {
                    toString(bi.bubbleCoord[0]),
                    toString(bi.bubbleCoord[1]),
                    toString(bi.bubbleCoord[2]),
                };
            }
        );
        p.addStringArgs(
            "bubble-fixed",
            [](MTOCInit& bi, const vector<string>& input) {
                if(input.size() != 2 || (input[1] != "true" && input[1] != "false")) {
                    log::error("Expected true or false for bubble-fixed.");
                    throw runtime_error("Expected true or false for bubble-fixed");
                }
                else {
                    bi.bubbleFixed = input[1] == "true";
                }
            },
            [](const MTOCInit& bi) {
                return vector<string> { bi.bubbleFixed ? "true" : "false" };
            }
        );
        p.addEmptyLine();

        p.addComment(" Attached filaments.");
        p.addSingleArg("filament-type", [](auto&& bi) -> auto& { return bi.filamentType; });
        p.addSingleArg("num-filaments", [](auto&& bi) -> auto& { return bi.numFilaments; });
        p.addSingleArg("num-cylinders-per-filament", [](auto&& bi) -> auto& { return bi.numCylindersPerFilament; });
        p.addComment(" Parameters for initial filament orientation.");
        p.addSingleArg("theta1", [](auto&& bi) -> auto& { return bi.theta1; });
        p.addSingleArg("theta2", [](auto&& bi) -> auto& { return bi.theta2; });
        p.addSingleArg("phi1", [](auto&& bi) -> auto& { return bi.phi1; });
        p.addSingleArg("phi2", [](auto&& bi) -> auto& { return bi.phi2; });
        p.addEmptyLine();

        p.addComment(" Mechanical parameters.");
        p.addSingleArg("attachment-stretching-k", [](auto&& bi) -> auto& { return bi.attachmentStretchingK; });

        p.addComment(" Chemical parameters.");
        p.addComment(" Emission-absorption reactions. Usage: (reaction-emi-abs <init-copy-number> <species-name-general> <species-name-diffusing-bulk> <diffusing/bulk> <emission-rate> <absorption-rate> <dynamic-rate-type>)");
        p.addComment(" - Emission rate is of unit 1/s.");
        p.addComment(" - Absorption rate is of unit nm^3/s.");
        p.addComment(" Example: (reaction-emi-abs A B diffusing 0.1 0.2 force)");
        p.addStringArgs(
            "reaction-emi-abs",
            [](MTOCInit& bi, const vector<string>& input) {
                if(input.size() != 8) {
                    log::error("Expected 7 arguments for reaction-emi-abs.");
                    throw runtime_error("Invalid arguments for reaction-emi-abs");
                }
                else {
                    ReactionEmissionAbsorptionSetupInit ri {};
                    parse(ri.initNumSpecies1, input[1]);
                    auto& r = ri.setup;
                    r.speciesName1 = input[2];
                    r.speciesName2 = input[3];
                    r.speciesType2 = input[4];
                    parse(r.emiRateConst, input[5]);
                    parse(r.absRateConst, input[6]);
                    parse(r.dyRateType, input[7]);
                    bi.vecRxnEmiAbs.push_back(std::move(ri));
                }
            },
            [](const MTOCInit& bi) {
                vector<vector<string>> res;
                for(auto& ri : bi.vecRxnEmiAbs) {
                    auto& r = ri.setup;
                    res.push_back({
                        toString(ri.initNumSpecies1),
                        r.speciesName1,
                        r.speciesName2,
                        r.speciesType2,
                        toString(r.emiRateConst),
                        toString(r.absRateConst),
                        toString(r.dyRateType),
                    });
                }
                return res;
            }
        );
        p.addSingleArg("emi-force-1", [](auto&& bi) -> auto& { return bi.emiForce1; });
        p.addSingleArg("emi-force-2", [](auto&& bi) -> auto& { return bi.emiForce2; });

        return p;
    }();

    return res;
}

inline const auto& afmInitParser() {
    using namespace std;

    static auto res = [] {
        KeyValueParser< AFMInit > p;

        p.addEmptyLine();
        p.addComment(" Bubble parameters.");
        p.addSingleArg("bubble-type", [](auto&& bi) -> auto& { return bi.bubbleType; });
        p.addStringArgs(
            "bubble-coord",
            [](AFMInit& bi, const vector<string>& input) {
                if(input.size() != 4) {
                    log::error("Expected 3 coordinates for bubble-coord.");
                    throw runtime_error("Expected 3 coordinates for bubble-coord");
                }
                else {
                    parse(bi.bubbleCoord[0], input[1]);
                    parse(bi.bubbleCoord[1], input[2]);
                    parse(bi.bubbleCoord[2], input[3]);
                }
            },
            [](const AFMInit& bi) {
                return vector<string> {
                    toString(bi.bubbleCoord[0]),
                    toString(bi.bubbleCoord[1]),
                    toString(bi.bubbleCoord[2]),
                };
            }
        );
        p.addStringArgs(
            "bubble-fixed",
            [](AFMInit& bi, const vector<string>& input) {
                if(input.size() != 2 || (input[1] != "true" && input[1] != "false")) {
                    log::error("Expected true or false for bubble-fixed.");
                    throw runtime_error("Expected true or false for bubble-fixed");
                }
                else {
                    bi.bubbleFixed = input[1] == "true";
                }
            },
            [](const AFMInit& bi) {
                return vector<string> { bi.bubbleFixed ? "true" : "false" };
            }
        );
        p.addEmptyLine();

        p.addComment(" Attached filaments.");
        p.addSingleArg("filament-type", [](auto&& bi) -> auto& { return bi.filamentType; });
        p.addSingleArg("num-filaments", [](auto&& bi) -> auto& { return bi.numFilaments; });
        p.addSingleArg("num-cylinders-per-filament", [](auto&& bi) -> auto& { return bi.numCylindersPerFilament; });
        p.addComment(" Parameters for initial filament orientation.");
        p.addSingleArg("theta1", [](auto&& bi) -> auto& { return bi.theta1; });
        p.addSingleArg("theta2", [](auto&& bi) -> auto& { return bi.theta2; });
        p.addSingleArg("phi1", [](auto&& bi) -> auto& { return bi.phi1; });
        p.addSingleArg("phi2", [](auto&& bi) -> auto& { return bi.phi2; });
        p.addEmptyLine();

        p.addComment(" Mechanical parameters.");
        p.addSingleArg("attachment-stretching-k", [](auto&& bi) -> auto& { return bi.attachmentStretchingK; });

        p.addComment(" Chemical parameters.");
        p.addComment(" Emission-absorption reactions. Usage: (reaction-emi-abs <init-copy-number> <species-name-general> <species-name-diffusing-bulk> <diffusing/bulk> <emission-rate> <absorption-rate> <dynamic-rate-type>)");
        p.addComment(" - Emission rate is of unit 1/s.");
        p.addComment(" - Absorption rate is of unit nm^3/s.");
        p.addComment(" Example: (reaction-emi-abs A B diffusing 0.1 0.2 force)");
        p.addStringArgs(
            "reaction-emi-abs",
            [](AFMInit& bi, const vector<string>& input) {
                if(input.size() != 8) {
                    log::error("Expected 7 arguments for reaction-emi-abs.");
                    throw runtime_error("Invalid arguments for reaction-emi-abs");
                }
                else {
                    ReactionEmissionAbsorptionSetupInit ri {};
                    parse(ri.initNumSpecies1, input[1]);
                    auto& r = ri.setup;
                    r.speciesName1 = input[2];
                    r.speciesName2 = input[3];
                    r.speciesType2 = input[4];
                    parse(r.emiRateConst, input[5]);
                    parse(r.absRateConst, input[6]);
                    parse(r.dyRateType, input[7]);
                    bi.vecRxnEmiAbs.push_back(std::move(ri));
                }
            },
            [](const AFMInit& bi) {
                vector<vector<string>> res;
                for(auto& ri : bi.vecRxnEmiAbs) {
                    auto& r = ri.setup;
                    res.push_back({
                        toString(ri.initNumSpecies1),
                        r.speciesName1,
                        r.speciesName2,
                        r.speciesType2,
                        toString(r.emiRateConst),
                        toString(r.absRateConst),
                        toString(r.dyRateType),
                    });
                }
                return res;
            }
        );
        p.addSingleArg("emi-force-1", [](auto&& bi) -> auto& { return bi.emiForce1; });
        p.addSingleArg("emi-force-2", [](auto&& bi) -> auto& { return bi.emiForce2; });

        return p;
    }();

    return res;
}

} // namespace medyan

#endif
