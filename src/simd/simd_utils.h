#pragma once

// See Intel Optimization Guide: Ch. 5.6.6.2 Clipping to an Arbitrary Unsigned Range [High, Low]
// high_add = 0xffff - high
// high_sub = 0xffff - high + low
static SIMD::data_type __forceinline high_bit_depth_pixels_clamp(SIMD::data_type pixels, SIMD::data_type high_add, SIMD::data_type high_sub, const SIMD::data_type& low)
{
    pixels = SIMD::_adds_epu16(pixels, high_add);
    pixels = SIMD::_subs_epu16(pixels, high_sub);
    pixels = SIMD::_add_epi16(pixels, low);

    return pixels;
}


// like high_bit_depth_pixels_clamp, but all values are 8bit
static SIMD::data_type __forceinline low_bit_depth_pixels_clamp(SIMD::data_type pixels, SIMD::data_type high_add, SIMD::data_type high_sub, const SIMD::data_type& low)
{
    pixels = SIMD::_adds_epu8(pixels, high_add);
    pixels = SIMD::_subs_epu8(pixels, high_sub);
    pixels = SIMD::_add_epi8(pixels, low);

    return pixels;
}
