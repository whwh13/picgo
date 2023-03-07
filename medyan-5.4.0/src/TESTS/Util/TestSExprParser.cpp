#include <vector>

#include <catch2/catch.hpp>

#include "Util/SExprParser.hpp"

TEST_CASE("s-expression parser", "[SExpr]") {
    using namespace std;
    using namespace medyan;

    const string strGood =
        "prop1         val1              \n"
        "                                \n"
        "; This is a comment             \n"
        "(prop2                          \n"
        "  (subprop1   \"val2 1\")       \n"
        "  (subprop2   \"val2 2\"))      \n"
        "(prop3        val3) ; comment3  \n";

    const string strBad =
        "))( (\n"
        "; (prop whatever) (да\n"
        "（） 草";

    SECTION("String tokenization") {
        {
            auto tl = sExprTokenize(strGood);

            REQUIRE(tl.size() == 26);
            auto tlv = vector<SExprToken>(tl.begin(), tl.end());

            // Check some strings.
            CHECK(tlv[0].type == SExprToken::Type::string);
            CHECK(tlv[0].content == "prop1");
            CHECK(tlv[7].type == SExprToken::Type::string);
            CHECK(tlv[7].content == "prop2");
            CHECK(tlv[11].type == SExprToken::Type::string);
            CHECK(tlv[11].content == "val2 1");

            // Check some comments.
            CHECK(tlv[4].type == SExprToken::Type::comment);
            CHECK(tlv[4].content == " This is a comment             ");
            CHECK(tlv[24].type == SExprToken::Type::comment);
            CHECK(tlv[24].content == " comment3  ");

            // Check some parentheses.
            CHECK(tlv[6].type == SExprToken::Type::parenthesisLeft);
            CHECK(tlv[9].type == SExprToken::Type::parenthesisLeft);
            CHECK(tlv[12].type == SExprToken::Type::parenthesisRight);
            CHECK(tlv[18].type == SExprToken::Type::parenthesisRight);

            // Check some line breaks.
            CHECK(tlv[2].type == SExprToken::Type::lineBreak);
            CHECK(tlv[3].type == SExprToken::Type::lineBreak);
            CHECK(tlv[25].type == SExprToken::Type::lineBreak);


            // Check string rebuild.
            auto str2 = toString(tl);
            auto tl2 = sExprTokenize(str2);
            CHECK(tl == tl2);

            auto str3 = toString(tl2);
            CHECK(str2 == str3);
        }
        {
            auto tl = sExprTokenize(strBad);

            REQUIRE(tl.size() == 9);
            auto tlv = vector<SExprToken>(tl.begin(), tl.end());

            CHECK(tlv[0].type == SExprToken::Type::parenthesisRight);
            CHECK(tlv[2].type == SExprToken::Type::parenthesisLeft);
            CHECK(tlv[5].type == SExprToken::Type::comment);
            CHECK(tlv[5].content == " (prop whatever) (да");
            CHECK(tlv[6].type == SExprToken::Type::lineBreak);
            CHECK(tlv[7].type == SExprToken::Type::string);
            CHECK(tlv[7].content == "（）");

            // Check string rebuild.
            auto str2 = toString(tl);
            auto tl2 = sExprTokenize(str2);
            CHECK(tl == tl2);

            auto str3 = toString(tl2);
            CHECK(str2 == str3);
        }
    }

    SECTION("s-expression build") {
        {
            auto tl = sExprTokenize(strGood);
            auto sel = lexTokenList(tl);

            REQUIRE(sel.size() == 3);

            REQUIRE(holds_alternative<SExpr::ListType>(sel[0].data));
            auto& l0 = get<SExpr::ListType>(sel[0].data);
            REQUIRE(l0.size() == 2);
            CHECK(holds_alternative<SExpr::StringType>(l0[0].data));
            CHECK(holds_alternative<SExpr::StringType>(l0[1].data));

            REQUIRE(holds_alternative<SExpr::ListType>(sel[1].data));
            auto& l1 = get<SExpr::ListType>(sel[1].data);
            REQUIRE(l1.size() == 3);
            CHECK(holds_alternative<SExpr::StringType>(l1[0].data));
            CHECK(holds_alternative<SExpr::ListType>(l1[1].data));
            CHECK(holds_alternative<SExpr::ListType>(l1[2].data));

            // Check token rebuild.
            auto tl2 = buildTokenList(sel);
            auto sel2 = lexTokenList(tl2);
            CHECK(sel == sel2);
        }
        {
            auto tl = sExprTokenize(strBad);
            CHECK_THROWS(lexTokenList(tl));
        }
    }
}
