#pragma once

#include <cmath>
#include <algorithm>
#include "core.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace neo_f3kdb {

static inline float saturate(float val) {
    return std::clamp(val, 0.0f, 1.0f);
}

static inline float calculate_ratio_term(float diff, float thresh) {
    if (thresh < 1e-5f)
        return (std::abs(diff) < 1e-5f) ? 1.0f : -1e6f;

    return 1.0f - std::abs(diff) / thresh;
}

template <class PixelIn>
static inline float calculate_gradient_angle(const process_plane_params& params, const PixelIn* src_base, int stride, int current_x, int current_y, int input_depth, int read_distance = 20) {
    auto get_pixel_value_at = [&](int x, int y) -> float {
        int clamped_y = std::clamp(y, 0, params.plane_height_in_pixels - 1);
        int clamped_x = std::clamp(x, 0, params.plane_width_in_pixels - 1);
        int val = static_cast<int>(src_base[clamped_y * stride + clamped_x]);
        return static_cast<float>(val << (16 - input_depth));
    };

    const float p00 = get_pixel_value_at(current_x - read_distance, current_y - read_distance);
    const float p10 = get_pixel_value_at(current_x, current_y - read_distance);
    const float p20 = get_pixel_value_at(current_x + read_distance, current_y - read_distance);
    const float p01 = get_pixel_value_at(current_x - read_distance, current_y);
    const float p21 = get_pixel_value_at(current_x + read_distance, current_y);
    const float p02 = get_pixel_value_at(current_x - read_distance, current_y + read_distance);
    const float p12 = get_pixel_value_at(current_x, current_y + read_distance);
    const float p22 = get_pixel_value_at(current_x + read_distance, current_y + read_distance);

    // Sobel-like gradient calculation
    const float gx = (p20 + 2.0f * p21 + p22) - (p00 + 2.0f * p01 + p02);
    const float gy = (p00 + 2.0f * p10 + p20) - (p02 + 2.0f * p12 + p22);

    const float scaled_epsilon_for_gx = 0.01f * (static_cast<float>(1 << (16 - input_depth)) * 3.0f);

    if (std::abs(gx) < scaled_epsilon_for_gx) {
        if (std::abs(gy) < scaled_epsilon_for_gx) {
            return 1.0f;
        }
        return 1.0f;
    }

    return std::atan(gy / gx) / static_cast<float>(M_PI) + 0.5f;
}

} // namespace neo_f3kdb
