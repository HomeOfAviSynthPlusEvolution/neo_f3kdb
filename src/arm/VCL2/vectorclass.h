#pragma once

#include "../../sse2neon.h"
#include "instrset.h"

struct Vec4fb
{
    __m128 v;
    Vec4fb(__m128 x) : v(x) {}
    operator __m128() const
    {
         return v;
    }
};

struct Vec4f
{
    __m128 v;
    Vec4f() : v(_mm_setzero_ps()) {}
    Vec4f(float f) : v(_mm_set1_ps(f)) {}
    Vec4f(__m128 x) : v(x) {}

    operator __m128() const
    {
        return v;
    }

    Vec4f& load(void const* p)
    {
        v = _mm_loadu_ps((const float*)p);
        return *this;
    }

    Vec4fb operator <= (Vec4f const& b) const
    {
        return _mm_cmple_ps(v, b.v);
    }

    Vec4f operator * (Vec4f const& b) const
    {
        return _mm_mul_ps(v, b.v);
    }
};

static inline Vec4f select(Vec4fb const& mask, Vec4f const& a, Vec4f const& b)
{
    return _mm_blendv_ps(b.v, a.v, mask.v);
}

typedef __m128i Vec4i;
typedef __m128i Vec8s;
