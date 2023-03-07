#ifndef MEDYAN_Util_Parser_MembraneParser_hpp
#define MEDYAN_Util_Parser_MembraneParser_hpp

#include "Parser.h"
#include "SysParams.h"

namespace medyan {

inline const auto& membraneSetupParser() {
    using namespace std;

    static auto res = [] {
        KeyValueParser< MembraneSetup > p;

        p.addEmptyLine();
        p.addSingleArg("vertex-system", [](auto&& ms) -> auto& { return ms.vertexSystem; });
        p.addSingleArg("has-lipid-reservoir", [](auto&& ms) -> auto& { return ms.hasLipidReservoir; });
        p.addSingleArg("area-k",        [](auto&& ms) -> auto& { return ms.areaElasticity; });
        p.addSingleArg("bending-k",     [](auto&& ms) -> auto& { return ms.bendingElasticity; });
        p.addSingleArg("eq-curv",       [](auto&& ms) -> auto& { return ms.eqMeanCurv; });
        p.addSingleArg("tension",       [](auto&& ms) -> auto& { return ms.tension; });
        p.addSingleArg("volume-k",      [](auto&& ms) -> auto& { return ms.volumeElasticity; });

        return p;
    }();

    return res;
}

inline const auto& membraneInitParser() {
    using namespace std;

    static auto res = [] {
        KeyValueParser< MembraneInit > p;

        p.addEmptyLine();
        p.addStringArgs(
            "mesh",
            [](MembraneInit& mi, const vector<string>& input) {
                mi.meshParams.assign(input.begin() + 1, input.end());
            },
            [](const MembraneInit& mi) {
                return mi.meshParams;
            }
        );
        p.addSingleArg("eq-area-factor", [](auto&& mi) -> auto& { return mi.eqAreaFactor; });
        p.addSingleArg("volume-offset",  [](auto&& mi) -> auto& { return mi.volumeOffset; });
        p.addStringArgs(
            "species",
            [](MembraneInit& mi, const vector<string>& input) {
                if(input.size() != 3) {
                    log::error("Membrane species init specification is incorrect.");
                    log::info("Membrane init species usage: (species <name> <copy-number>)");
                    throw runtime_error("Membrane species init specification is incorrect.");
                }
                else {
                    mi.speciesInitVec.push_back({ input[1], parse<int>(input[2]) });
                }
            },
            [](const MembraneInit& mi) {
                vector<vector<string>> res;
                for(const auto& s : mi.speciesInitVec) {
                    res.push_back({ s.name, toString(s.copyNumber) });
                }
                return res;
            }
        );

        return p;
    }();

    return res;
}

inline const auto& meshAdapterParser() {
    using namespace std;

    static auto res = [] {
        KeyValueParser< MeshAdapterSettings > p;

        p.addEmptyLine();
        p.addSingleArg("max-size", [](auto&& s) -> auto& { return s.maxSize; });

        return p;
    }();

    return res;
}

} // namespace medyan

#endif
