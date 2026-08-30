#pragma once

#include "core/constants.hpp"
#include "core/plane.hpp"
#include "core/sample_math.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace neo_f3kdb::core::sample_modes {

template <class PixelIn>
inline int upsample(PixelIn pixel, int input_depth) {
  return static_cast<int>(pixel) << (16 - input_depth);
}

struct ReferenceSamples {
  int src;
  int ref_y = 0;
  int ref_x = 0;
  int ref1 = 0;
  int ref2 = 0;
  int ref3 = 0;
  int ref4 = 0;
};

template <class PixelIn>
struct TypedSampleOps {
  const process_plane_params& params;
  span2d::ReadOnlyRestrictPlane<PixelIn> src_plane;

  [[nodiscard]] int read(int row, int col) const {
    return upsample(src_plane(row, col), params.config.input_depth);
  }

  [[nodiscard]] int avg2(int pixel1, int pixel2) const {
    return sample_modes::avg2(pixel1, pixel2);
  }

  [[nodiscard]] int avg4(int pixel1, int pixel2, int pixel3, int pixel4) const {
    return sample_modes::avg4(pixel1, pixel2, pixel3, pixel4);
  }
};

template <class SampleOps>
inline int read_upsampled_sample(
  const process_plane_params& params,
  const SampleOps& ops,
  int row,
  int col
) {
  row = std::clamp(row, 0, params.plane_height() - 1);
  col = std::clamp(col, 0, params.plane_width() - 1);
  return ops.read(row, col);
}

template <class PixelIn>
inline int read_upsampled_sample(
  const process_plane_params& params,
  span2d::ReadOnlyRestrictPlane<PixelIn> src_plane,
  int row,
  int col
) {
  const TypedSampleOps<PixelIn> ops{params, src_plane};
  return read_upsampled_sample(params, ops, row, col);
}

template <class SampleOps>
inline float read_upsampled_pixel(
  const process_plane_params& params,
  const SampleOps& ops,
  int row,
  int col
) {
  return static_cast<float>(read_upsampled_sample(params, ops, row, col));
}

template <class PixelIn>
inline float read_upsampled_pixel(
  const process_plane_params& params,
  span2d::ReadOnlyRestrictPlane<PixelIn> src_plane,
  int row,
  int col
) {
  const TypedSampleOps<PixelIn> ops{params, src_plane};
  return read_upsampled_pixel(params, ops, row, col);
}

template <int kSampleMode, class SampleOps>
inline void load_reference_samples_only_ref(
  const process_plane_params& params,
  const SampleOps& ops,
  int row,
  int col,
  int& ref1,
  int& ref2,
  int& ref3,
  int& ref4
) {
  static_assert(kSampleMode >= 1 && kSampleMode <= 7);

  const auto info = params.dither_info_plane().row(row)[static_cast<std::size_t>(col)];
  const int width_subsamp = params.config.width_subsampling;
  const int height_subsamp = params.config.height_subsampling;

  if constexpr (kSampleMode == 1) {
    const int ref_y = info.ref1 >> height_subsamp;
    ref1 = read_upsampled_sample(params, ops, row + ref_y, col);
    ref2 = read_upsampled_sample(params, ops, row - ref_y, col);
  } else if constexpr (kSampleMode == 2) {
    const int ref_y1 = info.ref1 >> height_subsamp;
    const int ref_y2 = info.ref2 >> height_subsamp;
    const int ref_x1 = info.ref1 >> width_subsamp;
    const int ref_x2 = info.ref2 >> width_subsamp;
    ref1 = read_upsampled_sample(params, ops, row + ref_y2, col + ref_x1);
    ref2 = read_upsampled_sample(params, ops, row - ref_y1, col + ref_x2);
    ref3 = read_upsampled_sample(params, ops, row - ref_y2, col - ref_x1);
    ref4 = read_upsampled_sample(params, ops, row + ref_y1, col - ref_x2);
  } else if constexpr (kSampleMode == 3) {
    const int ref_x = info.ref1 >> width_subsamp;
    ref1 = read_upsampled_sample(params, ops, row, col + ref_x);
    ref2 = read_upsampled_sample(params, ops, row, col - ref_x);
  } else if constexpr (kSampleMode == 4 || kSampleMode == 5) {
    const int ref_y = info.ref1 >> height_subsamp;
    const int ref_x = info.ref1 >> width_subsamp;
    ref1 = read_upsampled_sample(params, ops, row + ref_y, col);
    ref2 = read_upsampled_sample(params, ops, row - ref_y, col);
    ref3 = read_upsampled_sample(params, ops, row, col + ref_x);
    ref4 = read_upsampled_sample(params, ops, row, col - ref_x);
  } else if constexpr (kSampleMode == 6 || kSampleMode == 7) {
    const int ref_y = info.ref1 >> height_subsamp;
    const int ref_x = info.ref1 >> width_subsamp;
    ref1 = read_upsampled_sample(params, ops, row + ref_y, col);
    ref2 = read_upsampled_sample(params, ops, row - ref_y, col);
    ref3 = read_upsampled_sample(params, ops, row, col + ref_x);
    ref4 = read_upsampled_sample(params, ops, row, col - ref_x);
  }
}

template <int kSampleMode, class SampleOps>
ReferenceSamples load_reference_samples(
  const process_plane_params& params,
  const SampleOps& ops,
  int row,
  int col
) {
  ReferenceSamples samples{.src = ops.read(row, col)};
  load_reference_samples_only_ref<kSampleMode>(params, ops, row, col, samples.ref1, samples.ref2, samples.ref3, samples.ref4);
  return samples;
}

template <int kSampleMode, class PixelIn>
ReferenceSamples load_reference_samples(
  const process_plane_params& params,
  span2d::ReadOnlyRestrictPlane<PixelIn> src_plane,
  int row,
  int col
) {
  const TypedSampleOps<PixelIn> ops{params, src_plane};
  return load_reference_samples<kSampleMode, TypedSampleOps<PixelIn>>(params, ops, row, col);
}

template <class SampleOps>
float calculate_gradient_angle(
  const process_plane_params& params,
  const SampleOps& ops,
  int current_x,
  int current_y,
  int read_distance = 20
) {
  const float p00 = read_upsampled_pixel(
    params,
    ops,
    current_y - read_distance,
    current_x - read_distance
  );
  const float p10 = read_upsampled_pixel(params, ops, current_y - read_distance, current_x);
  const float p20 = read_upsampled_pixel(
    params,
    ops,
    current_y - read_distance,
    current_x + read_distance
  );
  const float p01 = read_upsampled_pixel(params, ops, current_y, current_x - read_distance);
  const float p21 = read_upsampled_pixel(params, ops, current_y, current_x + read_distance);
  const float p02 = read_upsampled_pixel(
    params,
    ops,
    current_y + read_distance,
    current_x - read_distance
  );
  const float p12 = read_upsampled_pixel(params, ops, current_y + read_distance, current_x);
  const float p22 = read_upsampled_pixel(
    params,
    ops,
    current_y + read_distance,
    current_x + read_distance
  );

  const float gx = (p20 + 2.0f * p21 + p22) - (p00 + 2.0f * p01 + p02);
  const float gy = (p00 + 2.0f * p10 + p20) - (p02 + 2.0f * p12 + p22);
  const float scaled_epsilon_for_gx =
    0.01f * (static_cast<float>(1 << (INTERNAL_BIT_DEPTH - params.config.input_depth)) * 3.0f);

  if (std::abs(gx) < scaled_epsilon_for_gx) {
    return 1.0f;
  }

  constexpr float kPi = 3.14159265358979323846f;
  return std::atan(gy / gx) / kPi + 0.5f;
}

template <class PixelIn>
float calculate_gradient_angle(
  const process_plane_params& params,
  span2d::ReadOnlyRestrictPlane<PixelIn> src_plane,
  int current_x,
  int current_y,
  int read_distance = 20
) {
  const TypedSampleOps<PixelIn> ops{params, src_plane};
  return calculate_gradient_angle(params, ops, current_x, current_y, read_distance);
}

} // namespace neo_f3kdb::core::sample_modes
