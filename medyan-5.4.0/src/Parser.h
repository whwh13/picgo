
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

#ifndef MEDYAN_Parser_h
#define MEDYAN_Parser_h

#include <algorithm>
#include <cctype> // isspace
#include <charconv>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <sstream>
#include <type_traits>
#include <variant>
#include <vector>

#include "common.h"
#include "SysParams.h"
#include "utility.h"
#include "Util/Math/Vec.hpp"
#include "Util/Parser/StringParser.hpp"
#include "Util/SExprParser.hpp"

namespace medyan {



inline std::vector< std::string > getStringVector(const SExpr::ListType& sel) {
    // Precondition:
    //   - Every element in sel is of string type
    using namespace std;

    vector< string > res;
    res.reserve(sel.size());
    for(const auto& eachSE : sel) {
        // Check whether eachSE.data is a string.
        if(holds_alternative< SExpr::StringType >(eachSE.data)) {
            res.push_back(get< SExpr::StringType >(eachSE.data));
        } else {
            throw std::runtime_error("getStringVector: sel contains non-string element");
        }
    }
    return res;
}



// Key-value parser
//
// Treats an input s-expression as a list of key-value pairs and parse by
// matching keywords
//
// Features:
//   - Convert system input (s-expression) to params
//   - Build formatted tokens (for output) from params
template< typename Params >
struct KeyValueParser {

    using ParserFunc = std::function< void(Params&, const SExpr::ListType&) >;
    using TokenBuildFunc = std::function< std::list< SExprToken >(const Params&) >;

    using ParserFuncDict = std::map< std::string, ParserFunc >;
    using TokenBuildList = std::vector< TokenBuildFunc >;

    // Two data structures are required for the parser to work.
    //   - A dictionary to match keywords to the appropriate function for
    //     parsing the arguments
    //   - A list which instructs how an input file should be generated from
    //     system params
    ParserFuncDict dict;
    TokenBuildList tokenBuildList;

    // Handy functions
    //---------------------------------

    // Parse or print a specific single-valued parameter.
    //
    // Template parameters
    //   - LocateParam:
    //     (Params&)       -> T&,       AND
    //     (const Params&) -> const T&
    template< typename LocateParam >
    void addSingleArg(
        std::string name,
        LocateParam&& locate
    ) {
        using namespace std;

        addStringArgs(
            name,
            [name, locate](Params& params, const vector<string>& lineVector) {
                auto& param = locate(params);

                if (lineVector.size() != 2) {
                    log::error("{} must have exactly one value.", name);
                    throw runtime_error("Invalid argument.");
                }

                auto& arg = lineVector[1];

                parse(param, arg);
            },
            [name, locate](const Params& params) {
                auto& param = locate(params);

                return vector<string> { toString(param) };
            }
        );
    }

    // Template parameters
    //   - FuncParse: void(Params&, const vector<string>&), including key
    //   - FuncBuild:
    //     vector<string>(const Params&), excluding key, OR
    //     vector<vector<string>>(const Params&), excluding key (for repeating items)
    template<
        typename FuncParse,
        typename FuncBuild
    >
    void addStringArgs(
        std::string name,
        FuncParse&& funcParse,
        FuncBuild&& funcBuild
    ) {
        using namespace std;

        return addStringArgsWithAliases(
            move(name),
            {},
            forward<FuncParse>(funcParse),
            forward<FuncBuild>(funcBuild)
        );
    }

    // With alias
    // Template parameters
    //   - FuncParse: void(Params&, const vector<string>&), including key
    //   - FuncBuild:
    //     vector<string>(const Param&), excluding key, OR
    //     vector<vector<string>>(const Param&), excluding key (for repeating items)
    template<
        typename FuncParse,
        typename FuncBuild
    >
    void addStringArgsWithAliases(
        std::string name,
        std::vector< std::string > aliases,
        FuncParse&& funcParse,
        FuncBuild&& funcBuild
    ) {
        using namespace std;

        return addArgsWithAliases(
            name,
            move(aliases),
            [funcParse] (Params& params, const SExpr::ListType& keyAndArgs) {
                funcParse(params, getStringVector(keyAndArgs));
            },
            [funcBuild, name] (const Params& params) {

                if constexpr(is_same_v< invoke_result_t< FuncBuild, const Params& >, vector<string> >) {

                    vector<string> tempResult = funcBuild(params);

                    list< SExprToken > res;
                    for(auto& s : tempResult) {
                        res.push_back(SExprToken::makeString(move(s)));
                    }

                    return res;

                } else {
                    vector<vector<string>> tempResult = funcBuild(params);

                    vector<list< SExprToken >> res;
                    for(auto& eachVS : tempResult) {
                        res.emplace_back();
                        for(auto& s : eachVS) {
                            res.back().push_back(SExprToken::makeString(move(s)));
                        }
                    }

                    return res;
                }
            }
        );
    }

    // Template parameters
    //   - FuncParse: (Params&, const SExpr::ListType&) -> void, including key
    //   - FuncBuild: (const Params&) -> list< ConfigFileToken >, excluding key
    template<
        typename FuncParse,
        typename FuncBuild
    >
    void addArgs(
        std::string name,
        FuncParse&& funcParse,
        FuncBuild&& funcBuild
    ) {
        using namespace std;
        return addArgsWithAliases(
            move(name),
            {},
            forward<FuncParse>(funcParse),
            forward<FuncBuild>(funcBuild)
        );
    }

    // With alias
    // Template parameters
    //   - FuncParse: (Params&, const SExpr::ListType&) -> void, including key
    //   - FuncBuild: (const Params&) -> list< SExprToken >, excluding key
    template<
        typename FuncParse,
        typename FuncBuild
    >
    void addArgsWithAliases(
        std::string name,
        std::vector< std::string > aliases,
        FuncParse&& funcParse,
        FuncBuild&& funcBuild
    ) {
        using namespace std;

        const auto insertResult = dict.insert({ name, forward<FuncParse>(funcParse) });
        if(!insertResult.second) {
            log::error("Duplicate keyword: {} in parser.", name);
            throw runtime_error("Duplicate keyword.");
        }
        auto& refFuncParse = insertResult.first->second;
        // Insert for aliases
        for(auto& alias : aliases) {
            const auto aliasInsertResult = dict.insert({ move(alias), refFuncParse });
            if(!aliasInsertResult.second) {
                log::error("Duplicate keyword: {} in parser.", alias);
                throw runtime_error("Duplicate keyword.");
            }
        }

        tokenBuildList.push_back(
            [funcBuild, name] (const Params& params) {

                list< SExprToken > res;

                if constexpr(
                    is_same_v< invoke_result_t< FuncBuild, const Params& >, list<SExprToken> >
                ) {
                    // single entry
                    res.push_back({ SExprToken::Type::parenthesisLeft });
                    res.push_back(SExprToken::makeString(name));

                    res.splice(res.end(), funcBuild(params));

                    res.push_back({ SExprToken::Type::parenthesisRight });
                    res.push_back({ SExprToken::Type::lineBreak });

                } else {
                    // multiple entries
                    vector<list<SExprToken>> tempResult = funcBuild(params);
                    for(auto& eachList : tempResult) {
                        res.push_back({ SExprToken::Type::parenthesisLeft });
                        res.push_back(SExprToken::makeString(name));

                        res.splice(res.end(), move(eachList));

                        res.push_back({ SExprToken::Type::parenthesisRight });
                        res.push_back({ SExprToken::Type::lineBreak });
                    }
                }

                return res;
            }
        );
    }

    void addComment(std::string comment) {
        // The comment must start with valid comment specifiers such as ';'
        tokenBuildList.push_back(
            [comment{ std::move(comment) }] (const Params&) {
                std::list< SExprToken > res {
                    SExprToken {
                        SExprToken::Type::comment,
                        comment
                    }
                };
                res.push_back({ SExprToken::Type::lineBreak });
                return res;
            }
        );
    }
    void addEmptyLine() {
        tokenBuildList.push_back(
            [] (const Params&) {
                return std::list< SExprToken > {
                    SExprToken { SExprToken::Type::lineBreak }
                };
            }
        );
    }

};

// What to do when the key is not in the parser function dictionary.
enum class KeyValueParserUnknownKeyAction {
    ignore, warn, error
};

// Parse an s-expr as a key-value pair, using the parser func dictionary.
template< typename Params >
inline void parseKeyValue(
    Params&                                                params,
    const SExpr&                                           se,
    const typename KeyValueParser<Params>::ParserFuncDict& dict,
    KeyValueParserUnknownKeyAction                         unknownKeyAction = KeyValueParserUnknownKeyAction::ignore
) {
    using namespace std;
    using SES = SExpr::StringType;
    using SEL = SExpr::ListType;

    // se.data must be a list type, or an exception will be thrown
    // se contains the (key arg1 arg2 ...) data
    // The key must be a string type, or an exception will be thrown

    // Here the assignments are by value because of possible temporary values on the right hand side.
    const SES key = get<SES>(car(se).data);
    SEL keyAndArgs = get<SEL>(se.data);

    // Check against the parsing dictionary
    if(auto it = dict.find(key); it == dict.end()) {
        switch(unknownKeyAction) {
            case KeyValueParserUnknownKeyAction::ignore:
                break;
            case KeyValueParserUnknownKeyAction::warn:
                log::warn("In the input file, {} cannot be recognized.", key);
                break;
            case KeyValueParserUnknownKeyAction::error:
                log::error("In the input file, {} cannot be recognized.", key);
                throw runtime_error("Unknown key in parser");
        }
    }
    else {
        // Execute the settings
        (it->second)(params, move(keyAndArgs));
    }
}

// Parse a list of s-expr as key-value pairs, using the parser func dictionary.
template< typename Params >
inline void parseKeyValueList(
    Params&                                                params,
    const SExpr::ListType&                                 sel,
    const typename KeyValueParser<Params>::ParserFuncDict& dict,
    KeyValueParserUnknownKeyAction                         unknownKeyAction = KeyValueParserUnknownKeyAction::warn
) {
    // se.data must be a list type, or an exception will be thrown
    for(const SExpr& eachList : sel) {
        parseKeyValue(params, eachList, dict, unknownKeyAction);
    }
}

// Parsing an s-expr as key-value pairs.
template< typename Params >
inline void parseKeyValueList(
    Params&                        params,
    const SExpr::ListType&         sel,
    const KeyValueParser<Params>&  parser,
    KeyValueParserUnknownKeyAction unknownKeyAction = KeyValueParserUnknownKeyAction::warn
) {
    return parseKeyValueList(params, sel, parser.dict, unknownKeyAction);
}

// Build the token list for key-value inputs.
template< typename Params >
inline auto buildTokens(
    const Params&                                          params,
    const typename KeyValueParser<Params>::TokenBuildList& tokenBuildList
) {
    std::list< SExprToken > res;

    for(const auto& eachTokenBuild : tokenBuildList) {
        res.splice(res.end(), eachTokenBuild(params));
    }
    return res;
}
template< typename Params >
inline auto buildTokens(
    const Params&                 params,
    const KeyValueParser<Params>& parser
) {
    return buildTokens(params, parser.tokenBuildList);
}


//-----------------------------------------------------------------------------
// The actual parsers dealing with all medyan system parameters.
//-----------------------------------------------------------------------------

KeyValueParser<SimulConfig> buildSystemParser();
inline const auto& systemParser() {
    static const auto parser = buildSystemParser();
    return parser;
}

KeyValueParser<SimulConfig> buildChemDataParser();
inline const auto& chemDataParser() {
    static const auto parser = buildChemDataParser();
    return parser;
}


/// A general parser
/*!
 *  A parser object, when initialized, opens an input file. Upon destruction, it 
 *  closes the file.
 */
class Parser {
protected:
    fstream _inputFile; ///< input file being used
    
public:
    Parser(string inputFileName) {
        _inputFile.open(inputFileName);
        if(!_inputFile.is_open()) {
            cout << "There was an error parsing file " << inputFileName
                 << ". Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        cout << "Loading file " << inputFileName << endl;
    }
    ~Parser() {_inputFile.close();}
};

/// Used to parse initial Filament information, initialized by the Controller.
struct FilamentParser {
  
    /// Reads filament input file. Returns a vector of tuples containing
    /// filament type and positions (start and end points).
    /// @note - Does not check for coordinate correctness.
    static medyan::FilamentData readFilaments(std::istream&);
};

/// Used to parse initial membrane vertex and neighbor information, initialized by the Controller.
struct MembraneParser {

    struct MembraneInfo {
        using coordinate_type = medyan::Vec< 3, floatingpoint >;
        std::vector< coordinate_type > vertexCoordinateList;
        std::vector< std::array< int, 3 > > triangleVertexIndexList;
    };

    /// Reads membrane vertex input file.
    /// @note - Does not check for coordinate correctness.
    static std::vector<MembraneInfo> readMembranes(std::istream&);
};

/// Used to parse initial Bubble information, initialized by the Controller.
struct BubbleParser {    
    /// Reads bubble input file. Returns a vector of tuples containing
    /// bubble type and position.
    /// @note - Does not check for coordinate correctness.
    static BubbleData readBubbles(std::istream&);
};



/// Used to parse pin positions if needed upon restart
class PinRestartParser: public Parser {
    
public:
    PinRestartParser(string inputFileName) : Parser(inputFileName) {}
    ~PinRestartParser() {}
    
    /// Reads pin positions from file, and sets filaments
    void resetPins();
};

} // namespace medyan

#endif
