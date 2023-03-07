#include <algorithm> // transform
#include <string>
#include <vector>

#include "catch2/catch.hpp"

#include "Util/Io/CmdParse.hpp"

namespace medyan {
namespace {

void parseForCommand(cmdparse::Command& cmd, std::vector< std::string > args) {
    using namespace std;
    size_t s = args.size();
    vector< char* > arr(s);
    for(size_t i = 0; i < s; ++i) arr[i] = &args[i][0];
    cmd.parse(s, arr.data());
}

} // namespace

TEST_CASE("Command line argument parsing", "[CmdParse]") {
    using namespace std;
    using namespace cmdparse;

    Command cmd("cmd", "desc");

    SECTION("Variable parsing check") {
        int argInt = 0;
        long argLong = 0;
        float argFloat = 0;
        double argDouble = 0;
        string argString;
        vector<int> argIntVector;

        cmd.addPosArgForVar("int", "desc", true, argInt);
        cmd.addPosArgForVar("long", "desc", true, argLong);
        cmd.addPosArgForVar("float", "desc", true, argFloat);
        cmd.addPosArgForVar("double", "desc", true, argDouble);
        cmd.addPosArgForVar("string", "desc", true, argString);
        cmd.addPosArgForVector("long", "desc", true, argIntVector);
        parseForCommand(cmd, {"cmd", "--", "-1", "2000", "-3e10", "4e-100", " 5 test string ", "6", "70", "-800"});

        CHECK(argInt == -1);
        CHECK(argLong == 2000);
        CHECK(argFloat == Approx(-3e10));
        CHECK(argDouble == Approx(4e-100));
        CHECK(argString == " 5 test string ");
        REQUIRE(argIntVector.size() == 3);
        CHECK(argIntVector[0] == 6);
        CHECK(argIntVector[1] == 70);
        CHECK(argIntVector[2] == -800);
    }

    SECTION("Positional arguments") {
        SECTION("No positional argument") {
            CHECK_THROWS_AS(parseForCommand(cmd, {"cmd", "arg"}), ParsingError);
        }

        SECTION("1 positional argument") {
            double arg = 1.0;

            SECTION("Required argument") {
                cmd.addPosArgForVar("arg", "desc", true, arg);
                CHECK_THROWS_AS(parseForCommand(cmd, {"cmd"}), ValidationError);
            }
            SECTION("Optional argument") {
                cmd.addPosArgForVar("arg", "desc", false, arg);
                parseForCommand(cmd, {"cmd", "--", "-2.0"});
                CHECK(arg == Approx(-2.0));

                SECTION("No required argument after optional argument.") {
                    cmd.addPosArgForVar("arg2", "desc2", true, arg);
                    CHECK_THROWS_AS(cmd.ruleCheck(), CommandLogicError);
                }
            }
        }

        SECTION("Positional argument list") {
            vector< float > args;

            SECTION("Required argument") {
                cmd.addPosArgForVector("args", "desc", true, args);
                CHECK_THROWS_AS(parseForCommand(cmd, {"cmd"}), ValidationError);
            }
            SECTION("Optional argument") {
                cmd.addPosArgForVector("args", "desc", false, args);

                SECTION("Supply none") {
                    parseForCommand(cmd, {"cmd"});
                    CHECK(args.size() == 0);
                }
                SECTION("Supply multiple") {
                    parseForCommand(cmd, {"cmd", "1.0", "2"});
                    REQUIRE(args.size() == 2);
                    CHECK(args[1] == Approx(2.0f));
                }
            }
        }

        SECTION("Mixing positional arguments") {
            // TODO
        }

    }
}

} // namespace medyan
