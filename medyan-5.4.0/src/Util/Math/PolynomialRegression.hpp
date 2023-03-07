#ifndef MEDYAN_Util_Math_PolynomialRegression_hpp
#define MEDYAN_Util_Math_PolynomialRegression_hpp

#include <array>
#include <cmath>
#include <iterator>
#include <stdexcept>
#include <type_traits>

namespace medyan {

// An implementation of polynomial regression.
// The implementation is modified based upon https://gist.github.com/chrisengelsma/108f7ab0a746323beaaf7d6634cf4add as of Aug 30, 2021. The original license is attached below.
/**
 * MIT License
 * 
 * Copyright (c) 2020 Chris Engelsma
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
template< int order, typename ContainerX, typename ContainerY >
inline auto polynomialRegression(
    const ContainerX& x,
    const ContainerY& y
) {
    using namespace std;
    using ValueType = std::common_type_t< std::decay_t< decltype(x[0]) >, std::decay_t< decltype(y[0]) > >;
    using ResType = std::array< ValueType, order+1 >;

    if(std::size(x) != std::size(y)) {
        throw std::runtime_error("Size of x and y should be equal.");
    }
    if(std::size(x) <= order) {
        throw std::runtime_error("Size of input data should be greater than the order of polynomial.");
    }
    const int n = std::size(x);

    // Final result for coefficients.
    ResType res {};

    // xx is an array that stores values of sum of (xi^2n).
    array< ValueType, 2 * order + 1 > xx;
    for (int i = 0; i < 2 * order + 1; ++i) {
        xx[i] = 0;
        for (int j = 0; j < n; ++j) {
            xx[i] += pow(x[j], i);
        }
    }

    // B = normal augmented matrix that stores the equations.
    array< ValueType, (order + 1) * (order + 2) > B_ {};
    const auto B = [&B_](int i, int j) -> ValueType& {
        return B_[i * (order + 2) + j];
    };

    for (int i = 0; i <= order; ++i) 
        for (int j = 0; j <= order; ++j) 
            B(i,j) = xx[i + j];

    // yy is an array to store values of sum of (xi^n * yi).
    array< ValueType, order + 1 > yy;
    for (int i = 0; i <= order; ++i) {
        yy[i] = 0;
        for (int j = 0; j < n; ++j) {
            yy[i] += pow(x[j], i)*y[j];
        }
    }

    // Load values of yy as last column of B
    for (int i = 0; i <= order; ++i) 
        B(i, order+1) = yy[i];

    // Pivotisation of the B matrix.
    for (int i = 0; i <= order; ++i) 
        for (int k = i+1; k <= order; ++k) 
            if (B(i,i) < B(k,i)) 
                for (int j = 0; j <= order + 1; ++j)
                    std::swap(B(i,j), B(k,j));

    // Performs the Gaussian elimination.
    // (1) Make all elements below the pivot equal to zero or eliminate the variable.
    for (int i=0; i < order; ++i)
        for (int k = i+1; k <= order; ++k) {
            auto t = B(k,i) / B(i,i);
            for (int j=0; j <= order + 1; ++j)
                B(k,j) -= t * B(i,j);        // (1)
        }

    // Back substitution.
    // (1) Set the variable as the rhs of last equation
    // (2) Subtract all lhs values except the target coefficient.
    // (3) Divide rhs by coefficient of variable being calculated.
    for (int i = order; i >= 0; --i) {
        res[i] = B(i, order+1);              // (1)
        for (int j = 0; j <= order; ++j)
            if (j != i)
                res[i] -= B(i,j) * res[j];   // (2)
        res[i] /= B(i,i);                    // (3)
    }

    return res;
}

} // namespace medyan

#endif
