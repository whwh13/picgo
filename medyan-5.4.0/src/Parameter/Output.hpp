#ifndef MEDYAN_Parameter_Output_hpp
#define MEDYAN_Parameter_Output_hpp

namespace medyan {

struct OutputParams {
    // Controls some of the logging behavior.
    bool logBriefProfilingEachSnapshot = false;
    bool logSimulationTimeEachCycle = false;
};

} // namespace medyan

#endif
