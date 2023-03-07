#include "Util/Io/CmdParse.hpp"

#include <algorithm>

namespace medyan::cmdparse {

namespace {

void usagePairFormatter(const std::string& key, const std::string& description, std::ostream& os) {
    static const size_t leftMargin = 4;
    static const size_t maxKeySize = 21;
    static const size_t midMargin = 1;

    os << std::string(leftMargin, ' ');

    size_t len = key.length();
    os << key;
    if(len > maxKeySize) os << '\n' << std::string(leftMargin + maxKeySize + midMargin, ' ');
    else os << std::string(maxKeySize + midMargin - len, ' ');

    os << description << '\n';
}

} // namespace


std::string Option::getReadableName()const {
    std::ostringstream oss;
    if(_short) {
        oss << '-' << _short;
        if(_long.length())
            oss << " (--" << _long << ')';
    } else if(_long.length())
        oss << "--" << _long;
    return oss.str();
}

std::string Option::getUsageName()const {
    std::ostringstream oss;
    if(_short) {
        oss << '-' << _short;
        if(_long.length()) {
            oss << ", --" << _long;
            if(_hasVariable) {
                oss << "=<" << _variableName << '>';
            }
        } else {
            if(_hasVariable) {
                oss << " <" << _variableName << '>';
            }
        }
    } else if(_long.length()) {
        oss << "--" << _long;
        if(_hasVariable) {
            oss << "=<" << _variableName << '>';
        }
    }
    return oss.str();
}


void Command::ruleCheck()const {

    // Check positional argument positions
    {
        size_t num = _posArgs.size();
        size_t numOpt = 0;
        size_t numList = 0;
        for(size_t i = 0; i < num; ++i) {
            if(_posArgs[i]->isRequired()) {
                if(numOpt)
                    throw CommandLogicError("Required positional argument should not appear after optional positional argument.");
            } else {
                ++numOpt;
            }
            if(_posArgs[i]->isList()) {
                ++numList;
                if(numList > 1)
                    throw CommandLogicError("Command should not contain more than one positional argument list.");
            }
        }
    }

    // Check option names
    for(const auto& op : options_) {
        if(op.getShortName() == 0 && op.getLongName() == "")
            throw CommandLogicError("Command should not contain options without name.");
    }

    // Recursively check all subcommands
    for(const auto& sc : subcommands_) {
        sc->ruleCheck();
    }
}

void Command::_parsePosArg(const std::string& arg) {
    if(state_._posArgIndex >= _posArgs.size())
        throw ParsingError("More positional argument specified than required: " + arg);

    ++state_._posArgCount;

    _posArgs[state_._posArgIndex]->occur();
    _posArgs[state_._posArgIndex]->activate()(arg);
    if(!_posArgs[state_._posArgIndex]->isList())
        ++state_._posArgIndex;

}

int Command::parse_(std::vector<std::string>& feed, size_t argp) {
    if(terminating_) {
        return argp;
    }

    for(; argp < feed.size(); ++argp) {
        std::string& thisArg = feed[argp];

        // Deduce the argument type
        if(state_._parseAsPosArg) {
            _parsePosArg(thisArg);

            continue;
        }
        
        if (thisArg == "--") {
            state_._parseAsPosArg = true;
            continue;
        }
        
        {
            // if no positional argument shows up, check if the argument is a subcommand
            auto cmdMatch = std::find_if(
                subcommands_.begin(), subcommands_.end(),
                [&thisArg](const std::unique_ptr<Command>& sc) { return sc->getName() == thisArg; }
            );
            if(state_._posArgCount == 0 && cmdMatch != subcommands_.end()) {
                // hand over the parsing to the subcommand
                (*cmdMatch)->occur();
                (*cmdMatch)->activate();
                return (*cmdMatch)->parse_(feed, argp + 1);
            }
        }

        // Check if the argument is a short option.
        if(thisArg.length() >= 2 && thisArg[0] == '-' && thisArg[1] != '-') {
            char shortName = thisArg[1];
            Command *p = this;
            auto shortNameMatch = [shortName](const Option& op) { return op.getShortName() == shortName; };
            auto shortOpMatch = std::find_if(
                p->options_.begin(), p->options_.end(),
                shortNameMatch
            );
            while(shortOpMatch == p->options_.end() && p->parent_) {
                p = p->parent_;
                shortOpMatch = std::find_if(
                    p->options_.begin(), p->options_.end(),
                    shortNameMatch
                );
            }

            if(shortOpMatch == p->options_.end())
                throw ParsingError(std::string("Unrecognized short option -") + shortName);

            // Short Option found
            shortOpMatch->occur();

            if(shortOpMatch->hasVariable()) {
                if(thisArg.length() > 2) {
                    // Use the rest as argument
                    std::string content = thisArg.substr(2);
                    shortOpMatch->activate(*this, content);
                } else {
                    // Use the next as argument
                    if(argp + 1 == feed.size()) {
                        // No next argument
                        throw ParsingError(std::string("Argument required for option -") + shortName);
                    } else {
                        ++argp;
                        shortOpMatch->activate(*this, feed[argp]);
                    }
                }
            } else {
                shortOpMatch->activate(*this, "");
                if(thisArg.length() > 2) {
                    // Prepare the rest for next round
                    thisArg = '-' + thisArg.substr(2);
                    --argp; // To keep argp unchanged in the next iteration
                } // else do nothing
            }

            continue;
        }
        
        // Check if the argument is a long option.
        if(thisArg.length() > 2 && thisArg[0] == '-' && thisArg[1] == '-') {
            size_t eqIdx = thisArg.find('=');
            std::string longName = (
                eqIdx == std::string::npos ?
                thisArg.substr(2) :
                thisArg.substr(2, eqIdx - 2)
            );

            if(longName.length() == 0)
                throw ParsingError("Invalid option " + thisArg);
            
            Command *p = this;
            auto longNameMatch = [&longName](const Option& op) { return op.getLongName() == longName; };
            auto longOpMatch = std::find_if(
                p->options_.begin(), p->options_.end(),
                longNameMatch
            );
            while(longOpMatch == p->options_.end() && p->parent_) {
                p = p->parent_;
                longOpMatch = std::find_if(
                    p->options_.begin(), p->options_.end(),
                    longNameMatch
                );
            }

            if(longOpMatch == p->options_.end())
                throw ParsingError(std::string("Unrecognized long option --") + longName);

            // Long Option found
            longOpMatch->occur();

            if(longOpMatch->hasVariable()) {
                if(eqIdx != std::string::npos) {
                    // Use everything after '=' as the argument
                    std::string content = thisArg.substr(eqIdx + 1);
                    longOpMatch->activate(*this, content);
                } else {
                    // Use the next as argument
                    if(argp + 1 == feed.size()) {
                        // No next argument
                        throw ParsingError(std::string("Argument required for option --") + longName);
                    } else {
                        ++argp;
                        longOpMatch->activate(*this, feed[argp]);
                    }
                }
            } else {
                longOpMatch->activate(*this, "");
                if(eqIdx != std::string::npos)
                    throw ParsingError(std::string("Option --") + longName + std::string(" should not take variable"));
            }

            continue;
        }

        // Default as positional argument
        _parsePosArg(thisArg);

    } // End argp for loop

    return argp;
}


void Command::validate_() const {
    // Check for required positional argument
    for(auto& pa : _posArgs) {
        if(pa->isRequired() && pa->getOccurenceCount() == 0)
            throw ValidationError(std::string("Must specify positional argument ") + pa->getName());
    }

    // Check for required option
    for(const auto& op : options_) {
        if(op.isRequired() && op.getOccurenceCount() == 0)
            throw ValidationError(std::string("Must specify option ") + op.getReadableName());
    }

    // Run user validation
    if(userValidation_) userValidation_();

    // Recursively validate subcommands
    for(const auto& sc : subcommands_) {
        if(sc->state_._occurenceCount)
            sc->validate_();
    }
}

void Command::printUsage(std::ostream& os)const {
    os << "Usage:\n";

    std::ostringstream ossFullName;
    if(_inheritedName.length()) ossFullName << _inheritedName << ' ';
    ossFullName << _name;

    size_t numReqOp = std::count_if(
        options_.begin(), options_.end(),
        [](const Option& op) { return op.isRequired(); }
    );
    size_t numOptOp = options_.size() - numReqOp;

    // Usage of command
    os << "    " << ossFullName.str();
    if(numReqOp)
        os << " <options>";
    if(numOptOp)
        os << " [<options>]";
    if(_posArgs.size()) {
        os << " [--]";
        for(auto& pa : _posArgs) {
            os << ' ';
            if(!pa->isRequired()) os << '[';
            os << '<' << pa->getName() << '>';
            if(pa->isList()) os << "...";
            if(!pa->isRequired()) os << ']';
        }
    }
    os << '\n';

    // Usage with subcommand
    if(subcommands_.size())
        os << "    " << ossFullName.str() << " <command>\n";

    // Section of subcommand
    if(subcommands_.size()) {
        os << "\nCommands:\n";
        for(const auto& sc : subcommands_)
            usagePairFormatter(sc->getName(), sc->getDescription(), os);
    }

    // Section of option
    if(options_.size()) {
        os << "\nOptions:\n";
        for(const auto& op : options_)
            if(op.isRequired())
                usagePairFormatter(op.getUsageName(), "[Required] " + op.getDescription(), os);
        for(const auto& op : options_)
            if(!op.isRequired())
                usagePairFormatter(op.getUsageName(), op.getDescription(), os);
    }

    os << std::endl;
}

} // namespace medyan::cmdparse
