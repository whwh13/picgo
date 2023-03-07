#ifndef MEDYAN_Util_Io_CmdParse_Hpp
#define MEDYAN_Util_Io_CmdParse_Hpp

#include <cstdlib> // exit
#include <functional>
#include <iostream> // ostream, cout
#include <memory> // unique_ptr
#include <sstream> // string format
#include <stdexcept> // custom exception
#include <string>
#include <utility> // forward, move
#include <vector>

#include "Util/Parser/StringParser.hpp"

namespace medyan::cmdparse {

// Forward decl
class Command;
template< typename T > struct VariableWrite;
template< typename T > struct VectorAppend;

// Exceptions
struct CommandLogicError : std::logic_error {
    CommandLogicError( const std::string& what_arg ) : std::logic_error("Command logic error: " + what_arg) {}
};
struct ParsingError : std::runtime_error {
    ParsingError( const std::string& what_arg ) : std::runtime_error("Parsing error: " + what_arg) {}
};
struct ValidationError : std::runtime_error {
    ValidationError( const std::string& what_arg ) : std::runtime_error("Validation error: " + what_arg) {}
};

// Positional argument
class PosArg {
private:
    const char* _name;
    std::string _description;
    bool _required;
    bool _list;

    struct State {
        size_t _occurenceCount = 0;
    } state_;

    std::function<void(const std::string&)> activate_;

public:

    PosArg(const char* name, const std::string& description, bool required, bool list, const std::function<void(const std::string&)>& activate) :
        _name(name), _description(description), _required(required), _list(list), activate_(activate) {}

    PosArg(PosArg&&) = default;

    bool isList() const { return _list; }
    bool isRequired()const { return _required; }

    const char* getName()const { return _name; }
    const std::string& getDescription()const { return _description; }

    void occur() { ++state_._occurenceCount; }
    size_t getOccurenceCount()const { return state_._occurenceCount; }

    const auto& activate() const { return activate_; }

};

class Option {
private:

    char _short = 0; // without "-". 0 for no short name
    std::string _long; // without "--". "" for no long name
    std::string _description;

    bool _hasVariable;
    std::string _variableName; // Useless if _hasVariable is false

    bool _required;

    struct State {
        size_t _occurenceCount = 0;
    } state_;

    std::function<void(const Command &, const std::string &)> activateFunc_;

public:

    template< typename ActivateFuncWithoutVar > // void(const Command&)
    Option(
        char                     shortName,
        const std::string&       longName,
        const std::string&       description,
        bool                     required,
        ActivateFuncWithoutVar&& act
    ) :
        _short(shortName),
        _long(longName),
        _description(description),
        _hasVariable(false),
        _variableName(),
        _required(required),
        activateFunc_([act](const Command& cmd, const std::string&) { act(cmd); })
    {}
    template< typename ActivateFuncWithVar > // void(const Command&, const std::string&)
    Option(
        char                  shortName,
        const std::string&    longName,
        const std::string&    variableName,
        const std::string&    description,
        bool                  required,
        ActivateFuncWithVar&& act
    ):
        _short(shortName),
        _long(longName),
        _description(description),
        _hasVariable(true),
        _variableName(variableName),
        _required(required),
        activateFunc_(std::forward<ActivateFuncWithVar>(act))
    {}

    Option(Option&&) = default;

    char getShortName()const { return _short; }
    const std::string& getLongName()const { return _long; }
    std::string getReadableName()const;
    std::string getUsageName()const;
    const std::string& getDescription()const { return _description; }
    const std::string& getVariableName()const { return _variableName; }

    bool hasVariable()const { return _hasVariable; }
    bool isRequired()const { return _required; }

    void occur() { ++state_._occurenceCount; }
    size_t getOccurenceCount()const { return state_._occurenceCount; }

    void activate(const Command& cmd, const std::string& var) const { return activateFunc_(cmd, var); }
};

// Option factory
inline Option makeOptionAsFlag(char shortName, const std::string& longName, const std::string& description, bool required) {
    // Add an option without argument
    return Option(shortName, longName, description, required, [](const Command&){});
}

template< typename T >
inline Option makeOptionWithVar(char shortName, const std::string& longName, const std::string& variableName, const std::string& description, bool required, T& var) {
    return Option(shortName, longName, variableName, description, required, [variableName, &var](const Command&, const std::string& arg) {
        VariableWrite<T>{variableName}(var, arg);
        // Caveat: if () is used instead of {}, then T(a)(b, c) will be considered as the definition of variable a with type T initialized with (b, c).
    });
}


class Command {
private:

    std::string _name;
    std::string _inheritedName;
    std::string _description;

    bool terminating_ = false; // When this is true, the command will
                               // immediately end parsing, ignoring the rest of
                               // the arguments.

    std::vector<std::unique_ptr<PosArg>>  _posArgs;
    std::vector<Option>                   options_;
    std::vector<std::unique_ptr<Command>> subcommands_;
    Command* parent_ = nullptr;

    std::function<void()> userValidation_;

    // Real parsing function
    int parse_(std::vector<std::string>& feed, size_t argp);
    void _parsePosArg(const std::string& arg); // helper
    // Internal validation. Recursive check for invoked subcommands.
    void validate_() const;

    // State variable
    struct State {
        bool _parseAsPosArg = false; // After the delimiter, this will be set to true.
        size_t _posArgCount = 0; // Number of positional argument encountered.
        size_t _posArgIndex = 0; // The index for the next positional argument.
                                 // Normally same with _posArgCount except on arg list.
        size_t _occurenceCount = 0;
    } state_;

    std::function< void() > activateFunc_;

    // Private constructor, allowing parent assignment and name inheriting
    template< typename ActivateFunc > // void()
    Command(
        Command*           parent,
        const std::string& inheritedName,
        const std::string& name,
        const std::string& description,
        ActivateFunc&&     activateFunc
    ) :
        _name(name),
        _inheritedName(inheritedName),
        _description(description),
        parent_(parent),
        activateFunc_(activateFunc)
    {}

public:


    // Public constructor, where parent and inheritedName are default
    template< typename ActivateFunc > // void()
    Command(const std::string& name, const std::string& description, ActivateFunc&& activateFunc) :
        _name(name), _description(description), activateFunc_(activateFunc) {}
    Command(const std::string& name, const std::string& description) :
        _name(name), _description(description) {}

    Command(Command&&) = default;

    // Check validity of specified rules. Recursively check for subcommands.
    void ruleCheck() const;

    const auto& getName()        const { return _name; }
    auto        getFullName()    const { return _inheritedName + ' ' + _name; }
    const auto& getDescription() const { return _description; }

    void setTerminating(bool b) { terminating_ = b; }

    void occur() { ++state_._occurenceCount; }

    void activate() const { if(activateFunc_) activateFunc_(); }

    // Specification
    template< typename... Args >
    PosArg* addPosArg(Args&&... args) {
        _posArgs.emplace_back(new PosArg(std::forward<Args>(args)...));
        return _posArgs.back().get();
    }
    template< typename T >
    PosArg* addPosArgForVar(const char* name, const std::string& description, bool required, T& var) {
        // Implies non-list
        return addPosArg(name, description, required, false, [name, &var](const std::string& arg) {
            VariableWrite<T>(std::string(name))(var, arg);
        });
    }
    template< typename T >
    PosArg* addPosArgForVector(const char* name, const std::string& description, bool required, std::vector<T>& var) {
        // Implies list
        return addPosArg(name, description, required, true, [name, &var](const std::string& arg) {
            VectorAppend<T>(std::string(name))(var, arg);
        });
    }

    Option& addOption(Option&& op) {
        options_.push_back(std::move(op));
        return options_.back();
    }
    Option& addHelp() {
        // Auxilliary for
        return addOption(Option('h', "help", "Print usage and exit", false, [](const Command& cmd) {
            cmd.printUsage();
            std::exit(EXIT_SUCCESS);
        }));
    }
    template< typename... Args >
    Command& addCommand(Args&&... args) {
        subcommands_.emplace_back(new Command(
            this,
            (_inheritedName.length() ? _inheritedName + ' ' + _name : _name),
            std::forward<Args>(args)...
        ));
        return *subcommands_.back();
    }
    Command& addCommandSimple(const std::string& name, const std::string& description) {
        return addCommand(name, description, []{});
    }

    template< typename Func > // Must be void()
    void setValidation(Func&& validation) { userValidation_ = std::forward<Func>(validation); }

    // The parsing function used only by the root command.
    // States will be changed in this function
    int parse(int argc, char** argv) {
        occur();
        activate();

        std::vector<std::string> inputFeed(argc);
        for(size_t i = 0; i < argc; ++i) inputFeed[i] = argv[i];
        int argpNext = parse_(inputFeed, 1);

        validate_();

        return argpNext;
    }

    // Auxilliary function that prints usage help message
    void printUsage(std::ostream& os = std::cout) const;
};

// Helper functions

// Helper function to write to a variable. Throws on error
template< typename T >
struct VariableWrite {
    std::string argName;
    VariableWrite() = default;
    VariableWrite(const std::string& argName): argName(argName) {}

    void operator()(T& var, const std::string& s) const {
        try {
            medyan::parse(var, s);
        }
        catch(const std::exception&) {
            throw ParsingError("Cannot understand argument value " + s + (argName.length() ? " for " + argName : std::string{}));
        }
    }
};

// Helper function to append to a vector. Throws on error
template< typename T >
struct VectorAppend {
    std::string argName;
    VectorAppend() = default;
    VectorAppend(const std::string& argName): argName(argName) {}

    void operator()(std::vector<T>& var, const std::string& s) const {
        T tmp;
        try {
            medyan::parse(tmp, s);
        }
        catch(const std::exception&) {
            throw ParsingError("Cannot append the argument value " + s + (argName.length() ? " for " + argName : std::string{}));
        }
        var.push_back(std::move(tmp));
    }
};

} // namespace medyan::cmdparse

#endif // include guard
