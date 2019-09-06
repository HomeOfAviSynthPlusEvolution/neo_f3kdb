#pragma once

// See Intel Optimization Guide: Ch. 5.6.6.2 Clipping to an Arbitrary Unsigned Range [High, Low]
// high_add = 0xffff - high
// high_sub = 0xffff - high + low
static __m256i __forceinline high_bit_depth_pixels_clamp(__m256i pixels, __m256i high_add, __m256i high_sub, const __m256i& low)
{
    pixels = _mm256_adds_epu16(pixels, high_add);
    pixels = _mm256_subs_epu16(pixels, high_sub);
    pixels = _mm256_add_epi16(pixels, low);

    return pixels;
}


// like high_bit_depth_pixels_clamp, but all values are 8bit
static __m256i __forceinline low_bit_depth_pixels_clamp(__m256i pixels, __m256i high_add, __m256i high_sub, const __m256i& low)
{
    pixels = _mm256_adds_epu8(pixels, high_add);
    pixels = _mm256_subs_epu8(pixels, high_sub);
    pixels = _mm256_add_epi8(pixels, low);

    return pixels;
}
