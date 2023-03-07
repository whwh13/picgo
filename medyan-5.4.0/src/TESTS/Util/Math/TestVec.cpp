#include <cstddef> // size_t
#include <vector>

#include <catch2/catch.hpp>

#include "Util/Math/Vec.hpp"

using namespace std;

TEST_CASE("Vec and RefVec tests", "[Vec]") {
    using namespace medyan;

    Vec< 3, float > v3f_1 = {0.0f, 1.0f, 2.0f};
    vector< float > vector3x2f = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
    auto v3f_2 = makeRefVec< 3 >(&vector3x2f[3]);  // 3 4 5
    auto v3f_3 = makeRefVec< 3 >(&vector3x2f[0]);  // 0 1 2

    Vec< 4, double > v4d_1 = {-4.0, -3.0, -2.0, -1.0};
    VecArray< 4, double > va4d = { {-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0} };
    auto v4d_2 = va4d[0];                            // -1 -2 -3 -4
    auto v4d_3 = makeRefVec< 4 >(va4d.data() + 4);   // -5 -6 -7 -8

    // Check element access
    REQUIRE(v3f_1[1] == Approx(1.0f));
    REQUIRE(v3f_2[1] == Approx(4.0f));
    REQUIRE(v3f_3[1] == Approx(1.0f));

    REQUIRE(v4d_1[2] == Approx(-2.0));
    REQUIRE(v4d_2[2] == Approx(-3.0));
    REQUIRE(v4d_3[2] == Approx(-7.0));

    SECTION("Conversion and copy assignments") {
        // copy assignment between Vec's
        {
            Vec3f v3f { 3.0f, 4.0f, 5.0f };
            v3f = v3f_1;
            CHECK(v3f[1] == Approx(1.0f));
        }
        // copy assignment between VecMap's
        {
            const double pool_1_restore[] { -1.0, -2.0, -3.0, -4.0 };
            double pool_1[] { -1.0, -2.0, -3.0, -4.0 };
            double pool_2[] { -5.0, -6.0, -7.0, -8.0 };
            const auto vcm4d_1 = makeRefVec< 4 >(pool_1_restore);
            const auto vm4d_1 = makeRefVec< 4 >(pool_1);
            const auto vm4d_2 = makeRefVec< 4 >(pool_2);

            // lvalue = lvalue
            vm4d_1 = vm4d_2;
            CHECK(pool_1[2] == Approx(-7.0));
            vm4d_1 = vcm4d_1;                      // restore
            REQUIRE(pool_1[2] == Approx(-3.0));    // check restore

            // rvalue = lvalue
            makeRefVec<4>(pool_1) = vm4d_2;
            CHECK(pool_1[2] == Approx(-7.0));
            vm4d_1 = vcm4d_1;                      // restore

            // lvalue = rvalue
            vm4d_1 = makeRefVec<4>(pool_2);
            CHECK(pool_1[2] == Approx(-7.0));
            vm4d_1 = vcm4d_1;                      // restore

            // rvalue = rvalue
            makeRefVec<4>(pool_1) = makeRefVec<4>(pool_2);
            CHECK(pool_1[2] == Approx(-7.0));
            vm4d_1 = vcm4d_1;                      // restore
        }
        // assignment between different Vec's
        {
            Vec3f v3f   { 0.0f, 1.0f, 2.0f };
            Vec3d v3d_1 {  6.0,  7.0,  8.0 };
            Vec3d v3d_2 { -1.0, -2.0, -3.0 };

            v3f = v3d_1;
            CHECK(v3f[1] == Approx(7.0f));

            v3d_2 = v3f;
            CHECK(v3d_2[1] == Approx(7.0));
        }
        // assignment between Vec and RefVec
        {
            Vec3f v3f { 0.0f, 1.0f, 2.0f };
            std::vector<double> vv3d { 6.0, 7.0, 8.0, -1.0, -2.0, -3.0 };
            auto v3d_1 = VecMap< 3, double >{ vv3d.data() };     //  6  7  8
            auto v3d_2 = makeRefVec< 3 >(vv3d.data() + 3);       // -1 -2 -3

            v3f = v3d_1;
            CHECK(v3f[1] == Approx(7.0f));

            v3d_2 = v3f;
            CHECK(v3d_2[1] == Approx(7.0));
        }
        // assignment between different RefVec's
        {
            std::vector<float> vv3f { 0.0f, 1.0f, 2.0f };
            std::vector<double> vv3d { 6.0, 7.0, 8.0, -1.0, -2.0, -3.0 };
            auto v3f = makeRefVec< 3 >(vv3f.data());                //  0  1  2
            auto v3d_1 = VecMap< 3, double >{ vv3d.data() };        //  6  7  8
            auto v3d_2 = makeRefVec< 3 >(vv3d.data() + 3);          // -1 -2 -3

            v3f = v3d_1;
            CHECK(v3f[1] == Approx(7.0f));

            v3d_2 = v3f;
            CHECK(v3d_2[1] == Approx(7.0));
        }
        // construction of Vec using conversions
        {
            // This is to ensure no ambiguity of construction from the same type
            Vec3f v3f (v3f_1);
            CHECK(v3f[1] == Approx(1.0f));

            Vec3d v3d (v3f_1);
            CHECK(v3d[1] == Approx(1.0));

            Vec< 4, float > v4f (v4d_3);
            CHECK(v4f[2] == Approx(-7.0));
        }
    }

    SECTION("Compound operators") {
        v3f_1 += v3f_2;
        CHECK(v3f_1[0] == Approx(3.0f));
        CHECK(v3f_1[1] == Approx(5.0f));
        CHECK(v3f_1[2] == Approx(7.0f));

        v3f_2 -= v3f_1;
        CHECK(v3f_2[0] == Approx(-0.0f));
        CHECK(v3f_2[1] == Approx(-1.0f));
        CHECK(v3f_2[2] == Approx(-2.0f));

        v3f_3 += v3f_2;
        CHECK(v3f_3[0] == Approx(0.0f));
        CHECK(v3f_3[1] == Approx(0.0f));
        CHECK(v3f_3[2] == Approx(0.0f));

        v4d_1 *= 2.0;
        CHECK(v4d_1[0] == Approx(-8.0));
        CHECK(v4d_1[1] == Approx(-6.0));
        CHECK(v4d_1[2] == Approx(-4.0));
        CHECK(v4d_1[3] == Approx(-2.0));

        v4d_2 /= 2.0;
        CHECK(v4d_2[0] == Approx(-0.5));
        CHECK(v4d_2[1] == Approx(-1.0));
        CHECK(v4d_2[2] == Approx(-1.5));
        CHECK(v4d_2[3] == Approx(-2.0));

        v4d_3 *= -1.0;
        CHECK(v4d_3[0] == Approx(5.0));
        CHECK(v4d_3[1] == Approx(6.0));
        CHECK(v4d_3[2] == Approx(7.0));
        CHECK(v4d_3[3] == Approx(8.0));
    }

    SECTION("Arithmetic operators") {
        auto res1 = v3f_1 + v3f_2;
        CHECK(res1[0] == Approx(3.0f));
        CHECK(res1[1] == Approx(5.0f));
        CHECK(res1[2] == Approx(7.0f));

        auto res2 = v3f_3 - v3f_2;
        CHECK(res2[0] == Approx(-3.0f));
        CHECK(res2[1] == Approx(-3.0f));
        CHECK(res2[2] == Approx(-3.0f));

        auto res3 = v4d_1 * 2.0;
        CHECK(res3[0] == Approx(-8.0));
        CHECK(res3[1] == Approx(-6.0));
        CHECK(res3[2] == Approx(-4.0));
        CHECK(res3[3] == Approx(-2.0));

        auto res4 = v4d_3 / 2.0;
        CHECK(res4[0] == Approx(-2.5));
        CHECK(res4[1] == Approx(-3.0));
        CHECK(res4[2] == Approx(-3.5));
        CHECK(res4[3] == Approx(-4.0));
    }

    SECTION("Vector products") {
        auto res1 = cross(v3f_1, v3f_2);
        CHECK(res1[0] == Approx(-3.0f));
        CHECK(res1[1] == Approx( 6.0f));
        CHECK(res1[2] == Approx(-3.0f));

        auto res2 = dot(v4d_1, v4d_2);
        CHECK(res2 == Approx(20.0));
    }

    SECTION("Factory functions") {
        std::vector< double > source {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
        auto res = makeVec<3>(&source[1]);
        CHECK(res[0] == Approx(1.0));
        CHECK(res[1] == Approx(2.0));
        CHECK(res[2] == Approx(3.0));
    }
}

TEST_CASE("VecArray tests", "[VecArray]") {
    using namespace medyan;

    VecArray< 3, float > v3f;
    VecArray< 4, double > v4d;

    // Appending elements
    const Vec< 3, float > push_v3f[] {
        { 0.0f, 0.0f, 0.0f },
        { 1.0f, 1.0f, 1.0f },
        { 2.0f, 2.0f, 2.0f },
        { 3.0f, 3.0f, 3.0f },
        { 4.0f, 4.0f, 4.0f }
    };
    const Vec< 4, double > push_v4d[] {
        { -4.0, -4.0, -4.0, -4.0 },
        { -3.0, -3.0, -3.0, -3.0 },
        { -2.0, -2.0, -2.0, -2.0 },
        { -1.0, -1.0, -1.0, -1.0 }
    };

    for(auto& x : push_v3f) v3f.push_back(x);
    for(auto& x : push_v4d) v4d.push_back(x);

    REQUIRE(v3f.size() == 5);
    REQUIRE(v3f.size_raw() == 15);
    REQUIRE(v4d.size() == 4);
    REQUIRE(v4d.size_raw() == 16);

    SECTION("RefVec accessing and iterating") {
        // Accessing
        auto v3f_2 = v3f[2];
        CHECK(v3f_2[0] == Approx(2.0f));
        CHECK(v3f_2[1] == Approx(2.0f));
        CHECK(v3f_2[2] == Approx(2.0f));

        Vec< 3, float > v3f_2_new { 2.1f, 2.1f, 2.1f };
        v3f[2] = v3f_2_new;
        CHECK(v3f_2[0] == Approx(2.1f));
        CHECK(v3f_2[1] == Approx(2.1f));
        CHECK(v3f_2[2] == Approx(2.1f));

        const auto& v4d_cref = v4d;
        auto v4d_1c = v4d_cref[1];
        CHECK(v4d_1c[0] == Approx(-3.0));
        CHECK(v4d_1c[1] == Approx(-3.0));
        CHECK(v4d_1c[2] == Approx(-3.0));
        CHECK(v4d_1c[3] == Approx(-3.0));

        // Iterating
        auto v4d_2 = v4d[2];
        size_t v4d_2_cnt = 0;
        double v4d_2_sum = 0.0;
        for(auto x : v4d_2) {
            ++v4d_2_cnt;
            v4d_2_sum += x;
        }
        REQUIRE(v4d_2_cnt == 4);
        CHECK(v4d_2_sum == Approx(-8.0));

        auto v4d_3 = v4d[3];
        double v4d_3_accu = 0.0;
        for(auto& x : v4d_3) {
            v4d_3_accu += 0.1;
            x += v4d_3_accu;
        }
        CHECK(v4d_3[0] == Approx(-0.9));
        CHECK(v4d_3[1] == Approx(-0.8));
        CHECK(v4d_3[2] == Approx(-0.7));
        CHECK(v4d_3[3] == Approx(-0.6));
    }

    SECTION("Conversion and assignment of VecMap") {
        // Using conversion from RefVec to Vec
        Vec< 3, float > v3f_2_copy (v3f[2]);
        CHECK(v3f_2_copy[0] == 2.0f); // Check copy is successful

        v3f_2_copy[0] = 2.1f;
        CHECK(v3f[2][0] == Approx(2.0f)); // Should not change original value

        // Using copy assignment operator from Vec to RefVec
        Vec< 3, float > v3f_2_another_copy;
        v3f_2_another_copy = v3f[2];
        CHECK(v3f_2_another_copy[0] == Approx(2.0f)); // Check copy is successful

        v3f_2_another_copy[0] = 2.1f;
        CHECK(v3f[2][0] == Approx(2.0f)); // Should not change original value

        // Copy assignment operator between RefVecs
        v3f[3] = v3f[2];
        CHECK(v3f[3][0] == Approx(2.0f)); // Check copy is successful

    }

    SECTION("Iterator usage") {
        // Iterator pointer, dereference, op[]
        auto v3f_i = v3f.begin() + 1;
        CHECK(v3f_i->size() == 3);
        CHECK((*v3f_i)[0] == Approx(1.0f));
        CHECK(v3f_i[2][0] == Approx(3.0f));

        // Ctor and assignment
        const auto& v3f_cref = v3f;
        auto v3f_i_ccopy = v3f_cref.begin();
        CHECK((*v3f_i_ccopy)[0] == Approx(0.0f));
        v3f_i_ccopy = v3f_i; // Copy assign
        CHECK((*v3f_i_ccopy)[0] == Approx(1.0f));
        decltype(v3f)::const_iterator v3f_ic(v3f_i); // Copy ctor
        CHECK((*v3f_ic)[0] == Approx(1.0f));

        // Check comparison
        auto v3f_ie = v3f.end();
        auto v3f_ice = v3f_cref.end();
        CHECK(v3f_ie == v3f_ice);
        CHECK(v3f_i != v3f_ie);
        CHECK(v3f_i != v3f_ice);
        CHECK(!(v3f_i > v3f_ie));
        CHECK(  v3f_i < v3f_ie );
        CHECK(!(v3f_i >= v3f_ie));
        CHECK(  v3f_i <= v3f_ie );
        CHECK(v3f_ie >= v3f_ice);
        CHECK(v3f_ie <= v3f_ice);

        // Check arithmetics
        CHECK(v3f_i - v3f_ice == -4);
        CHECK(v3f_i + 4 == v3f_ice);
        CHECK(4 + v3f_i == v3f_ice);
        CHECK(v3f_ice - 4 == v3f_i);

        // Check inc/dec
        ++v3f_i;
        CHECK((*v3f_i)[0] == Approx(2.0f));
        --v3f_i;
        CHECK((*v3f_i)[0] == Approx(1.0f));
        v3f_i += 2;
        CHECK((*v3f_i)[0] == Approx(3.0f));
        v3f_i -= 2;
        CHECK((*v3f_i)[0] == Approx(1.0f));
        auto v3f_i_copy = v3f_i++;
        CHECK((*v3f_i_copy)[0] == Approx(1.0f));
        CHECK((*v3f_i     )[0] == Approx(2.0f));
        v3f_i_copy = v3f_i--;
        CHECK((*v3f_i_copy)[0] == Approx(2.0f));
        CHECK((*v3f_i     )[0] == Approx(1.0f));
    }

    SECTION("Modifiers") {
        // Check popping
        v3f.pop_back();
        CHECK(v3f.size() == 4);
        CHECK(v3f.size_raw() == 12);

        // Check resize
        v3f.resize(5);
        REQUIRE(v3f.size() == 5);
        REQUIRE(v3f.size_raw() == 15);

    }

    SECTION("Arithmetics") {
        // Dot product
        CHECK(dot(v3f, v3f) == Approx(90.0f));

        // Increment and decrement
        VecArray< 3, float > v3f1 {{ 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f }};
        VecArray< 3, float > v3f2 {{ 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f }};
        v3f1 += v3f2;
        CHECK(v3f1[0][1] == Approx(2.2f));
        v3f1 -= v3f2;
        CHECK(v3f1[1][1] == Approx(5.0f));
    }
}
