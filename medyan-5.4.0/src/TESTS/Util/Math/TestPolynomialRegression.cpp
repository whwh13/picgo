#include <vector>

#include <catch2/catch.hpp>

#include "Util/Math/PolynomialRegression.hpp"

TEST_CASE("Polynomial regression test", "[PolynomialRegression]") {
    using namespace std;
    using namespace medyan;

    // General case 1 (exact match).
    {
        INFO("General case 1");
        // y = 3x^2 + 2x + 1
        double x[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        double y[] { 6, 17, 34, 57, 86, 121, 162, 209, 262, 321 };
        auto reg2 = polynomialRegression<2>(x, y);
        REQUIRE(reg2.size() == 3);
        CHECK(reg2[0] == Approx(1.0f));
        CHECK(reg2[1] == Approx(2.0f));
        CHECK(reg2[2] == Approx(3.0f));

        auto reg5 = polynomialRegression<5>(x, y);
        REQUIRE(reg5.size() == 6);
        CHECK(reg5[0] == Approx(1.0f));
        CHECK(reg5[1] == Approx(2.0f));
        CHECK(reg5[2] == Approx(3.0f));
        CHECK(reg5[3] == Approx(0.0f).margin(1e-8));
        CHECK(reg5[4] == Approx(0.0f).margin(1e-8));
        CHECK(reg5[5] == Approx(0.0f).margin(1e-8));
    }

    // General case 2 (exact match).
    {
        INFO("General case 2");
        array<double, 5> coeff { -0.5, 2.0, 0, 1, 4 };
        vector<double> x(50);
        vector<double> y(50);
        for(int i = 0; i < 50; ++i) {
            x[i] = 0.1 * i;
            y[i] = 0;
            for(int j = 0; j < coeff.size(); ++j) {
                y[i] += coeff[j] * pow(x[i], j);
            }
        }

        auto reg4 = polynomialRegression<4>(x, y);
        REQUIRE(reg4.size() == 5);
        CHECK(reg4[0] == Approx(coeff[0]));
        CHECK(reg4[1] == Approx(coeff[1]));
        CHECK(reg4[2] == Approx(coeff[2]).margin(1e-5));
        CHECK(reg4[3] == Approx(coeff[3]));
        CHECK(reg4[4] == Approx(coeff[4]));
    }

    // Single number case
    {
        INFO("Single number case");
        float x[] { 10 };
        float y[] { 20 };
        auto reg0 = polynomialRegression<0>(x, y);
        REQUIRE(reg0.size() == 1);
        CHECK(reg0[0] == Approx(20.0f));

        CHECK_THROWS(polynomialRegression<1>(x, y));
    }

    // Linear case
    {
        INFO("Linear case");
        double x[] { -2, 0, 2, 4, 6 };
        float  y[] { 100, 200, 300, 200, 100 };
        auto reg1 = polynomialRegression<1>(x, y);
        REQUIRE(reg1.size() == 2);
        CHECK(reg1[1] == Approx(0.0f).margin(1e-8));

        CHECK_THROWS(polynomialRegression<5>(x, y));
    }
}
