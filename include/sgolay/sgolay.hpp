/**
 * @file       sgolay.hpp
 * 
 * @author     Alex Skrynnyk (github.com/skrynnyk)
 *
 * @brief      This is a header-only implementation of the Slavitzky-Golay
 *             Filter, backed by Gram polynomials for convolution weight
 *             generation based on a paper by Peter A. Gorry "General
 *             Least-Squares Smoothing and Differentiation by the Convolution
 *             (Savitzky-Golay) Method".
 *
 * @copyright  Copyright 2023 Alex Skrynnyk. All rights reserved.
 * 
 * @license    This project is released under the MIT License.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in 
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef SGOLAY_HPP
#define SGOLAY_HPP

/* Includes ------------------------------------------------------------------*/
#include <array>
#include <cstdint>


/* Public typedefs -----------------------------------------------------------*/
namespace sgf {

/* "Strong-types" for various numerical coefficients/parameters */
namespace detail {
template<typename T, typename Ttag>
class StrongValue {
public:
    T v;

public: 
    constexpr explicit StrongValue(T v)
        : v(v)
    {}
};
} /* namespace detail */

using I = detail::StrongValue<float, struct Iparam>;
using K = detail::StrongValue<float, struct Kparam>;
using M = detail::StrongValue<float, struct Mparam>;
using N = detail::StrongValue<float, struct Nparam>;
using S = detail::StrongValue<float, struct Sparam>;
using T = detail::StrongValue<float, struct Tparam>;

/* More readable aliases for coefficient/parameter single-letter variables */
using HalfWidth = M;
using PolyOrder = N;
using DerivOrder = S;

/* Public functions ----------------------------------------------------------*/
constexpr std::size_t getWindowSize(M m) {
    return m.v * 2 + 1;
}

namespace detail {
/* Recursive Weight generation functions using Gram polynomials */
constexpr float gramPoly(I i, M m, K k, S s) {
    if (k.v > 0) {
        return (4 * k.v - 2)
            / (k.v * (2 * m.v - k.v + 1))
            * (i.v * gramPoly(i, m, K(k.v - 1), s) + s.v 
                * gramPoly(i, m, K(k.v - 1), S(s.v - 1)))
            - ((k.v - 1) * (2 * m.v + k.v))
            / (k.v * (2 * m.v - k.v + 1)) 
            * gramPoly(i, m, K(k.v - 2), s);
    }
    return (k.v == 0 && s.v == 0) ? 1 : 0;
}

constexpr float genFact(float a, float b) {
    float gf = 1.0f;
    for (std::int32_t j = a - b + 1; j <= a; j++) {
        gf = gf * j;
    }
    return gf;
}

constexpr float weight(I i, T t, M m, N n, S s) {
    float sum = 0;
    for (std::int32_t k = 0; k <= n.v; k++) {
        sum += (2 * k + 1) 
            * (genFact(2 * m.v, k) / genFact(2 * m.v + k + 1, k + 1))
            * gramPoly(i, m, K(k), S(0))
            * gramPoly(I(t.v), m, K(k), s);
    }
    return sum;
}
} /* namespace detail */

} /* namespace sgf */

#endif /* SGOLAY_HPP */
