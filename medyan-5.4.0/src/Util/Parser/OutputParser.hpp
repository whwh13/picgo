#ifndef MEDYAN_Util_Parser_OutputParser_hpp
#define MEDYAN_Util_Parser_OutputParser_hpp

#include "Parser.h"
#include "Parameter/Output.hpp"
#include "Util/Parser/StringParser.hpp"

namespace medyan {

inline const auto& outputParamsParser() {
    using namespace std;

    static auto res = [] {
        KeyValueParser< OutputParams > p;

        p.addEmptyLine();
        p.addComment(" Logging behavior adjustments.");
        p.addSingleArg(
            "log-brief-profiling-each-snapshot",
            [](auto&& params) -> auto& { return params.logBriefProfilingEachSnapshot; }
        );
        p.addSingleArg(
            "log-simulation-time-each-cycle",
            [](auto&& params) -> auto& { return params.logSimulationTimeEachCycle; }
        );

        return p;
    }();

    return res;
}

} // namespace medyan

#endif
