#pragma once

#include "core/constants.hpp"
#include "core/math_utils.hpp"

namespace neo_f3kdb::core::pixel_proc_common {

static inline int upsample(void*, unsigned char pixel)
{
    return pixel << (INTERNAL_BIT_DEPTH - 8);
}

static inline int downsample(int pixel, int pixel_min, int pixel_max, int output_depth)
{
    return clamp_pixel(pixel, pixel_min, pixel_max) >> (INTERNAL_BIT_DEPTH - output_depth);
}

static inline int avg_2(void*, int pixel1, int pixel2)
{
    return (pixel1 + pixel2 + 1) >> 1;
}

static inline int avg_4(void*, int pixel1, int pixel2, int pixel3, int pixel4)
{
    // consistent with SSE code
    int avg1 = (pixel1 + pixel2 + 1) >> 1;
    int avg2 = (pixel3 + pixel4 + 1) >> 1;
    if (avg1 > 0)
    {
        avg1 -= 1;
    }
    return (avg1 + avg2 + 1) >> 1;
}

} // namespace neo_f3kdb::core::pixel_proc_common
