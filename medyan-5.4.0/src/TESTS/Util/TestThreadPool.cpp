#ifdef THIS_THREADPOOL_IS_CURRENTLY_NOT_USED_IN_MEDYAN
#include <chrono>

#include "catch2/catch.hpp"

#include "Util/ThreadPool.hpp"

TEST_CASE("Thread Pool test", "[ThreadPool]") {
    using namespace std;

    SECTION("Job submission and job counter") {
        size_t numThreads = 4;

        // Create stuff needed for the test
        atomic_bool hold {true};
        atomic_int tasksEngaged {0};
        ThreadPool tp(numThreads);
        REQUIRE(tp.numThreads() == numThreads);
        CHECK(tp.numWorkingThreads() == 0);
        CHECK(tp.numIdleThreads() == numThreads);

        auto f = [&] {
            volatile unsigned int a = 0;
            ++tasksEngaged;
            // Loop forever until hold is false
            while(hold) {
                ++a;
                this_thread::yield();
            }
        };
        vector< future< void > > futs;

        // Push a job to the queue
        futs.push_back(tp.submit(f));
        while(tasksEngaged < 1) this_thread::yield();
        REQUIRE(tp.numThreads() == numThreads);
        CHECK(tp.numWorkingThreads() == 1);
        CHECK(tp.numIdleThreads() == numThreads - 1);
        CHECK(tp.numTasks() == 0);

        // Push 10 jobs to the queue
        for(int i = 0; i < 10; ++i) {
            futs.push_back(tp.submit(f));
        }
        while(tasksEngaged < numThreads) this_thread::yield();
        REQUIRE(tp.numThreads() == numThreads);
        CHECK(tp.numWorkingThreads() == numThreads);
        CHECK(tp.numIdleThreads() == 0);
        CHECK(tp.numTasks() == 11 - numThreads);

        // Finish the job
        hold = false;
    }

    SECTION("Pressure test") {
        size_t numThreads = 100;

        ThreadPool tp(numThreads);

        int sum = 0;
        vector< future< int > > futs;
        for(int i = 0; i < 10000; ++i) {
            futs.push_back(tp.submit(
                [](int x) { return x; },
                i + 1
            ));
        }
        for(auto& fut : futs) sum += fut.get();
        CHECK(sum == 50005000);

    }
}
#endif
