#ifndef MEDYAN_Controller_SpecialProtocols_hpp
#define MEDYAN_Controller_SpecialProtocols_hpp

#include <functional>
#include <stdexcept>

#include "Structure/SubSystem.h"

namespace medyan {

struct SpecialProtocol {
    enum class Schedule {
        // The protocol will not run.
        none,
        // The protocol will run when current simulation time >= set time for the first time.
        once,
        // The protocol will run at every iteration for current simulation time < set time.
        before,
        // The protocol will run at every iteration for current simulation time >= set time.
        after,
    };

    // Settings.
    //----------------------------------
    Schedule schedule = Schedule::none;
    double   setTime = 0;
    std::function<void(SubSystem&, const SimulConfig&)> run;

    // States.
    //----------------------------------

    // Will be incremented after every run.
    int      runCount = 0;
    // Will be updated after every run.
    int      lastRunIteration = -1;
    double   lastRunTime = 0;
};

struct SpecialProtocolManager {
    std::vector<SpecialProtocol> protocols;

    // Will be updated before every run.
    int      currentIteration = -1;

    // Auxiliary function to execute a protocol unconditionally.
    void executeUnconditionally(SubSystem& sys, const SimulConfig& conf, SpecialProtocol& protocol) const {
        protocol.run(sys, conf);
        protocol.runCount++;
        protocol.lastRunIteration = currentIteration;
        protocol.lastRunTime = sys.tau();
    }

    // Execute all protocols checking conditions.
    void execute(SubSystem& sys, const SimulConfig& conf) {
        ++currentIteration;

        for(auto& protocol : protocols) {
            switch(protocol.schedule) {
                case SpecialProtocol::Schedule::none:
                    break;
                case SpecialProtocol::Schedule::once:
                    if(protocol.runCount == 0 && sys.tau() >= protocol.setTime) {
                        executeUnconditionally(sys, conf, protocol);
                    }
                    break;
                case SpecialProtocol::Schedule::before:
                    if(sys.tau() < protocol.setTime) {
                        executeUnconditionally(sys, conf, protocol);
                    }
                    break;
                case SpecialProtocol::Schedule::after:
                    if(sys.tau() >= protocol.setTime) {
                        executeUnconditionally(sys, conf, protocol);
                    }
                    break;
                default:
                    log::error("Unknown protocol schedule ID: {}. Ignored.", underlying(protocol.schedule));
                    throw runtime_error("Invalid special protocol schedule.");
            }
        }
    }

    // Auxiliary function to add a special protocol.
    //
    // run: (SubSystem&, const SimulConfig&) -> void.
    template< typename RunFunc >
    void addProtocol(SpecialProtocol::Schedule schedule, double setTime, RunFunc&& run) {
        SpecialProtocol protocol;
        protocol.schedule = schedule;
        protocol.setTime = setTime;
        protocol.run = std::move(run);
        protocols.push_back(std::move(protocol));
    }
};

} // namespace medyan

#endif
