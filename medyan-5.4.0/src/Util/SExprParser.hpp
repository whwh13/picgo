#ifndef MEDYAN_Util_SExprParser_hpp
#define MEDYAN_Util_SExprParser_hpp

#include <algorithm>
#include <ostream>
#include <sstream>
#include <list>

#include "Util/SExpr.hpp"

namespace medyan {

// Tokenizer and lexer for s-expressions.
//
// This extends the input config to s-expressions, which allows more extensibility, but it does not include the lisp syntax.
//
// The input file will be treated as a quoted list of s-expressions. For backwards compatibility and simplicity, several special specifications exist:
//   - All symbols are parsed as strings.
//   - '#' and ';' both start a line comment.
//   - If a top level (not within any parentheses) token is not a parenthesis, an implicit pair of parentheses will be added before the token and before the next closest top level line break or end of input.
//   - Double quotation mark is parsed as string, but nested quotations, escapings are not allowed. Line comment marker in a quoted string has no effect.
//   - cons that is not a list is currently not supported.
//   - Most syntactic sugar is currently not supported.

struct SExprToken {
    // Spaces are not recorded as tokens. However, inline comments and line breaks will be recorded because
    // - Comments are used when generating input files.
    // - Line breaks are used to determine implicit parentheses.
    enum class Type {
        string,            // Surrounding quotes are removed.
        comment,           // Needed for generating input files. Leading marker is removed.
        parenthesisLeft,
        parenthesisRight,
        lineBreak,         // Needed for implicit parenthesis.
        unknown
    };

    Type        type = Type::string;
    std::string content;

    // Use default in C++20.
    bool operator==(const SExprToken& other) const {
        return type == other.type && content == other.content;
    }
    bool operator!=(const SExprToken& other) const {
        return !(*this == other);
    }

    static const char* defaultString(Type type) {
        switch(type) {
            case Type::string:           return "";
            case Type::comment:          return ";";
            case Type::parenthesisLeft:  return "(";
            case Type::parenthesisRight: return ")";
            case Type::lineBreak:        return "\n";
            case Type::unknown:          return "";
            default:                     return "";
        }
    }

    static auto makeString(std::string_view sv) {
        return SExprToken { Type::string, std::string(sv) };
    }
};

inline bool isTokenSpecialChar(char x) {
    return std::isspace(static_cast<unsigned char>(x)) ||
        x == '#' || x == ';' ||
        x == '(' || x == ')' ||
        x == '"';
}


//-----------------------------------------------------------------------------
// From string to token list.
//-----------------------------------------------------------------------------

// Tokenize and the input file.
inline std::list< SExprToken > sExprTokenize(std::string_view sv) {
    const auto len = sv.length();

    std::list< SExprToken > tokens;

    for(int i = 0; i < len; ++i) {
        const char a = sv[i];
        if(a == '\n') {
            tokens.push_back({ SExprToken::Type::lineBreak });
        }
        else if(std::isspace(static_cast<unsigned char>(a))) {
            // Do nothing
        }
        else if(a == '#' || a == ';') {
            int j = i + 1;
            while(j < len && sv[j] != '\n') ++j;

            // Now j points to either the end of string, or the next line break.
            // The leading token and the line break will not be stored.
            tokens.push_back({SExprToken::Type::comment, std::string(sv.substr(i+1, j-i-1))});
            i = j - 1;
        }
        else if(a == '(') {
            tokens.push_back({ SExprToken::Type::parenthesisLeft });
        }
        else if(a == ')') {
            tokens.push_back({ SExprToken::Type::parenthesisRight });
        }
        else if(a == '"') {
            int j = i + 1;
            while(j < len && sv[j] != '"') ++j;

            // Now j points to either end of string, or the next double quotation
            if(j < len) {
                tokens.push_back(SExprToken::makeString(sv.substr(i+1, j-i-1)));
                i = j;
            } else {
                log::error("Quotation marks do not match");
                throw std::runtime_error("Quotation marks do not match.");
            }
        }
        else {
            int j = i + 1;
            while(j < len && !isTokenSpecialChar(sv[j])) ++j;

            // Now j points to either the end of string, or the next special char
            tokens.push_back(SExprToken::makeString(sv.substr(i, j-i)));
            i = j - 1;
        }
    }

    return tokens;
}

//-----------------------------------------------------------------------------
// Serialize token list to stream.
//-----------------------------------------------------------------------------


inline void outputTokenList(std::ostream& os, const std::list< SExprToken >& tokens) {
    using TT = SExprToken::Type;
    const int indent = 4;

    auto lastType = TT::unknown;
    int indentLevel = 0;
    for(const auto& tok : tokens) {
        // Add indentation.
        if(lastType == TT::lineBreak) {
            os << std::string(std::max(0, indentLevel * indent), ' ');
        }

        // Print token.
        if(tok.type == TT::string) {

            auto content = tok.content;
            // Strip double quotations
            content.erase(
                std::remove(content.begin(), content.end(), '"'),
                content.end()
            );

            // Check special characters exist
            const bool hasSpecial = (
                std::find_if(content.begin(), content.end(), isTokenSpecialChar)
                    != content.end());

            if(lastType == TT::string || lastType == TT::parenthesisRight) {
                os << ' ';
            }

            if(hasSpecial || content.empty()) {
                os << '"' << content << '"';
            } else {
                os << content;
            }
        }
        else if(tok.type == TT::comment) {
            if(lastType != TT::lineBreak && lastType != TT::unknown) {
                os << ' ';
            }
            os << SExprToken::defaultString(tok.type) << tok.content;
        }
        else if(tok.type == TT::parenthesisLeft) {
            ++indentLevel;
            if(lastType == TT::string || lastType == TT::parenthesisRight) {
                os << ' ';
            }

            os << SExprToken::defaultString(tok.type);
        }
        else if(tok.type == TT::parenthesisRight) {
            --indentLevel;

            os << SExprToken::defaultString(tok.type);
        }
        else {
            os << SExprToken::defaultString(tok.type);
        }

        lastType = tok.type;
    }
}

inline std::string toString(const std::list< SExprToken >& tokens) {
    std::ostringstream oss;
    outputTokenList(oss, tokens);
    return oss.str();
}


//-----------------------------------------------------------------------------
// From token list to a list of s-expressions.
//-----------------------------------------------------------------------------

inline SExpr::ListType lexTokenList(
    const std::list< SExprToken >&           tokens,
    std::list< SExprToken >::const_iterator& tokenIter,
    const int                                depth,
    const bool                               implicitParentheses
) {
    using namespace std;
    using SEL = SExpr::ListType;

    SEL res;

    while(tokenIter != tokens.cend()) {
        if(tokenIter->type == SExprToken::Type::string) {
            // Add to the current list
            res.push_back(
                SExpr { tokenIter->content }
            );
            ++tokenIter;
        }
        else if(tokenIter->type == SExprToken::Type::parenthesisLeft) {
            res.push_back(
                SExpr { lexTokenList(tokens, ++tokenIter, depth + 1, false) }
            );
        }
        else if(tokenIter->type == SExprToken::Type::parenthesisRight) {
            if(implicitParentheses) {
                LOG(ERROR) << "Unexpected ')' in implicit parentheses.";
                throw runtime_error("Unmatched parentheses");
            } else {
                // End current list
                ++tokenIter;
                return res;
            }
        }
        else if(tokenIter->type == SExprToken::Type::lineBreak) {
            if(implicitParentheses) {
                // End current list
                ++tokenIter;
                return res;
            }
            else {
                ++tokenIter;
            }
        }
        else {
            ++tokenIter;
        }
    }

    if(implicitParentheses) {
        // Reaching the end
        return res;
    }
    else {
        LOG(ERROR) << "')' is expected, but not found.";
        throw runtime_error("Unmatched parentheses");
    }
}

inline SExpr::ListType lexTokenList(const list< SExprToken >& tokens) {
    using namespace std;
    using SEL = SExpr::ListType;

    SEL res;

    for(auto tokenIter = tokens.cbegin(); tokenIter != tokens.cend(); ) {

        if (tokenIter->type == SExprToken::Type::string) {
            // Implicit parentheses assumed
            res.push_back(
                SExpr { lexTokenList(tokens, tokenIter, 1, true) }
            );
        }
        else if (tokenIter->type == SExprToken::Type::parenthesisLeft) {
            res.push_back(
                SExpr { lexTokenList(tokens, ++tokenIter, 1, false) }
            );
        }
        else if (tokenIter->type == SExprToken::Type::parenthesisRight) {
            throw runtime_error("Unmatched parentheses: unexpected ')'.");
        }
        else {
            ++tokenIter;
        }

    }

    return res;
}


//-----------------------------------------------------------------------------
// From s-expression to token list.
//-----------------------------------------------------------------------------

// Build tokens from s-expression.
// No comment or line break will be created. No formatting will be applied.
inline std::list< SExprToken > buildTokenList(const SExpr& se) {
    using namespace std;
    using TokenList = list< SExprToken >;

    TokenList tokens;

    struct TokenBuildVisitor {
        TokenList& theList;
        void operator()(const SExpr::StringType& str) const {
            theList.push_back(SExprToken::makeString( str ));
        }
        void operator()(const SExpr::ListType& seList) const {
            theList.push_back({ SExprToken::Type::parenthesisLeft });
            for(const auto& eachSE : seList) {
                visit(TokenBuildVisitor{theList}, eachSE.data);
            }
            theList.push_back({ SExprToken::Type::parenthesisRight });
        }
    };

    visit(TokenBuildVisitor{ tokens }, se.data);
    return tokens;
}

inline std::list< SExprToken > buildTokenList(const SExpr::ListType& sel) {
    using namespace std;
    using TokenList = list< SExprToken >;

    TokenList tokens;

    for(auto& eachSE : sel) {
        tokens.splice(tokens.end(), buildTokenList(eachSE));
    }
    return tokens;
}


} // namespace medyan

#endif
