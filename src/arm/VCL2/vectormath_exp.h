#pragma once

#include <cmath>

static inline Vec4f pow(Vec4f base, float exponent)
{
    alignas(16) float b[4];
    _mm_store_ps(b, base);
    return _mm_set_ps(std::pow(b[3], exponent),
                     std::pow(b[2], exponent),
                     std::pow(b[1], exponent),
                     std::pow(b[0], exponent));
}
