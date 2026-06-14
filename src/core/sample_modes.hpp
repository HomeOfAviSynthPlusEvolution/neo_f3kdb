#pragma once

#include "core/constants.hpp"
#include "core/plane.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

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

inline int avg2(int pixel1, int pixel2) {
  return (pixel1 + pixel2 + 1) >> 1;
}

inline int avg4(int pixel1, int pixel2, int pixel3, int pixel4) {
  int avg1 = avg2(pixel1, pixel2);
  const int avg2_value = avg2(pixel3, pixel4);
  if (avg1 > 0) {
    avg1 -= 1;
  }
  return avg2(avg1, avg2_value);
}

inline bool over_threshold(int threshold, int diff) {
  return std::abs(diff) >= threshold;
}

inline float saturate(float value) {
  return std::clamp(value, 0.0f, 1.0f);
}

inline float calculate_ratio_term(float diff, float threshold) {
  if (threshold < 1e-5f) {
    return std::abs(diff) < 1e-5f ? 1.0f : -1e6f;
  }
  return 1.0f - std::abs(diff) / threshold;
}

template <class PixelIn>
inline int read_upsampled_sample(
  const process_plane_params& params,
  StridedPlaneView<const PixelIn> src_plane,
  int row,
  int col
) {
  row = std::clamp(row, 0, params.plane_height_in_pixels - 1);
  col = std::clamp(col, 0, params.plane_width_in_pixels - 1);
  return upsample(src_plane.row(row)[static_cast<std::size_t>(col)], params.input_depth);
}

template <class PixelIn>
inline float read_upsampled_pixel(
  const process_plane_params& params,
  StridedPlaneView<const PixelIn> src_plane,
  int row,
  int col
) {
  return static_cast<float>(read_upsampled_sample(params, src_plane, row, col));
}

template <int kSampleMode, class PixelIn>
ReferenceSamples load_reference_samples(
  const process_plane_params& params,
  StridedPlaneView<const PixelIn> src_plane,
  int row,
  int col
) {
  static_assert(kSampleMode >= 1 && kSampleMode <= 7);

  const auto info = params.dither_info_plane().row(row)[static_cast<std::size_t>(col)];
  const int width_subsamp = params.width_subsampling;
  const int height_subsamp = params.height_subsampling;

  ReferenceSamples samples{
    .src = upsample(src_plane.row(row)[static_cast<std::size_t>(col)], params.input_depth)
  };

  if constexpr (kSampleMode == 1) {
    samples.ref_y = info.ref1 >> height_subsamp;
    samples.ref1 = read_upsampled_sample(params, src_plane, row - samples.ref_y, col);
    samples.ref2 = read_upsampled_sample(params, src_plane, row + samples.ref_y, col);
  } else if constexpr (kSampleMode == 2) {
    const int ref_y1 = info.ref2 >> height_subsamp;
    const int ref_y2 = info.ref1 >> height_subsamp;
    const int ref_x1 = info.ref1 >> width_subsamp;
    const int ref_x2 = info.ref2 >> width_subsamp;
    samples.ref1 = read_upsampled_sample(params, src_plane, row - ref_y1, col - ref_x1);
    samples.ref2 = read_upsampled_sample(params, src_plane, row - ref_y2, col + ref_x2);
    samples.ref3 = read_upsampled_sample(params, src_plane, row + ref_y1, col + ref_x1);
    samples.ref4 = read_upsampled_sample(params, src_plane, row + ref_y2, col - ref_x2);
  } else if constexpr (kSampleMode == 3) {
    samples.ref_x = info.ref1 >> width_subsamp;
    samples.ref1 = read_upsampled_sample(params, src_plane, row, col - samples.ref_x);
    samples.ref2 = read_upsampled_sample(params, src_plane, row, col + samples.ref_x);
  } else if constexpr (kSampleMode == 4 || kSampleMode == 5) {
    samples.ref_y = info.ref1 >> height_subsamp;
    samples.ref_x = info.ref1 >> width_subsamp;
    samples.ref1 = read_upsampled_sample(params, src_plane, row - samples.ref_y, col);
    samples.ref2 = read_upsampled_sample(params, src_plane, row + samples.ref_y, col);
    samples.ref3 = read_upsampled_sample(params, src_plane, row, col - samples.ref_x);
    samples.ref4 = read_upsampled_sample(params, src_plane, row, col + samples.ref_x);
  } else if constexpr (kSampleMode == 6 || kSampleMode == 7) {
    samples.ref_y = info.ref1 >> height_subsamp;
    samples.ref_x = info.ref1 >> width_subsamp;
    samples.ref1 = read_upsampled_sample(params, src_plane, row + samples.ref_y, col);
    samples.ref2 = read_upsampled_sample(params, src_plane, row - samples.ref_y, col);
    samples.ref3 = read_upsampled_sample(params, src_plane, row, col + samples.ref_x);
    samples.ref4 = read_upsampled_sample(params, src_plane, row, col - samples.ref_x);
  }

  return samples;
}

template <class PixelIn>
float calculate_gradient_angle(
  const process_plane_params& params,
  StridedPlaneView<const PixelIn> src_plane,
  int current_x,
  int current_y,
  int read_distance = 20
) {
  const float p00 = read_upsampled_pixel(
    params,
    src_plane,
    current_y - read_distance,
    current_x - read_distance
  );
  const float p10 = read_upsampled_pixel(params, src_plane, current_y - read_distance, current_x);
  const float p20 = read_upsampled_pixel(
    params,
    src_plane,
    current_y - read_distance,
    current_x + read_distance
  );
  const float p01 = read_upsampled_pixel(params, src_plane, current_y, current_x - read_distance);
  const float p21 = read_upsampled_pixel(params, src_plane, current_y, current_x + read_distance);
  const float p02 = read_upsampled_pixel(
    params,
    src_plane,
    current_y + read_distance,
    current_x - read_distance
  );
  const float p12 = read_upsampled_pixel(params, src_plane, current_y + read_distance, current_x);
  const float p22 = read_upsampled_pixel(
    params,
    src_plane,
    current_y + read_distance,
    current_x + read_distance
  );

  const float gx = (p20 + 2.0f * p21 + p22) - (p00 + 2.0f * p01 + p02);
  const float gy = (p00 + 2.0f * p10 + p20) - (p02 + 2.0f * p12 + p22);
  const float scaled_epsilon_for_gx =
    0.01f * (static_cast<float>(1 << (INTERNAL_BIT_DEPTH - params.input_depth)) * 3.0f);

  if (std::abs(gx) < scaled_epsilon_for_gx) {
    return 1.0f;
  }

  constexpr float kPi = 3.14159265358979323846f;
  return std::atan(gy / gx) / kPi + 0.5f;
}

template <bool kBlurFirst>
inline int process_mode1(int src_px, int ref1, int ref2, int threshold) {
  const int avg = avg2(ref1, ref2);
  const bool use_src = kBlurFirst
    ? over_threshold(threshold, avg - src_px)
    : over_threshold(threshold, src_px - ref1) ||
        over_threshold(threshold, src_px - ref2);
  return use_src ? src_px : avg;
}

template <bool kBlurFirst>
inline int process_mode2(int src_px, int ref1, int ref2, int ref3, int ref4, int threshold) {
  const int avg = avg4(ref1, ref2, ref3, ref4);
  const bool use_src = kBlurFirst
    ? over_threshold(threshold, avg - src_px)
    : over_threshold(threshold, ref1 - src_px) ||
        over_threshold(threshold, ref2 - src_px) ||
        over_threshold(threshold, ref3 - src_px) ||
        over_threshold(threshold, ref4 - src_px);
  return use_src ? src_px : avg;
}

template <bool kBlurFirst>
inline int process_mode3(int src_px, int ref1, int ref2, int threshold) {
  return process_mode1<kBlurFirst>(src_px, ref1, ref2, threshold);
}

template <bool kBlurFirst>
inline int process_mode4(
  int src_px,
  int ref_v1,
  int ref_v2,
  int ref_h1,
  int ref_h2,
  int threshold
) {
  const int avg_v = avg2(ref_v1, ref_v2);
  const bool use_src_v = kBlurFirst
    ? over_threshold(threshold, avg_v - src_px)
    : over_threshold(threshold, src_px - ref_v1) ||
        over_threshold(threshold, src_px - ref_v2);
  const int new_v = use_src_v ? src_px : avg_v;

  const int avg_h = avg2(ref_h1, ref_h2);
  const bool use_src_h = kBlurFirst
    ? over_threshold(threshold, avg_h - src_px)
    : over_threshold(threshold, src_px - ref_h1) ||
        over_threshold(threshold, src_px - ref_h2);
  const int new_h = use_src_h ? src_px : avg_h;

  return avg2(new_v, new_h);
}

inline int process_mode5(
  int src_px,
  int ref_h1,
  int ref_h2,
  int ref_w1,
  int ref_w2,
  int threshold,
  int threshold1,
  int threshold2
) {
  const int avg = avg4(ref_h1, ref_h2, ref_w1, ref_w2);
  const int avg_dif = std::abs(avg - src_px);
  const int max_dif = std::max({
    std::abs(ref_h1 - src_px),
    std::abs(ref_h2 - src_px),
    std::abs(ref_w1 - src_px),
    std::abs(ref_w2 - src_px)
  });
  const int mid_dif1 = std::abs(ref_h1 + ref_h2 - 2 * src_px);
  const int mid_dif2 = std::abs(ref_w1 + ref_w2 - 2 * src_px);

  const bool use_src = over_threshold(threshold, avg_dif) ||
    over_threshold(threshold1, max_dif) ||
    over_threshold(threshold2, mid_dif1) ||
    over_threshold(threshold2, mid_dif2);
  return use_src ? src_px : avg;
}

inline int process_mode6(
  int src_px,
  float ref_h1,
  float ref_h2,
  float ref_w1,
  float ref_w2,
  float threshold,
  float threshold1,
  float threshold2
) {
  const float org_pix = static_cast<float>(src_px);
  const float avg_refs = (ref_h1 + ref_h2 + ref_w1 + ref_w2) * 0.25f;
  const float avg_dif = std::abs(avg_refs - org_pix);
  const float max_dif = std::max({
    std::abs(ref_h1 - org_pix),
    std::abs(ref_h2 - org_pix),
    std::abs(ref_w1 - org_pix),
    std::abs(ref_w2 - org_pix)
  });
  const float mid_dif_v = std::abs(ref_h1 + ref_h2 - 2.0f * org_pix);
  const float mid_dif_h = std::abs(ref_w1 + ref_w2 - 2.0f * org_pix);

  const float factor = std::pow(
    saturate(3.0f * calculate_ratio_term(avg_dif, threshold)) *
      saturate(3.0f * calculate_ratio_term(max_dif, threshold1)) *
      saturate(3.0f * calculate_ratio_term(mid_dif_v, threshold2)) *
      saturate(3.0f * calculate_ratio_term(mid_dif_h, threshold2)),
    0.1f
  );

  return static_cast<int>((org_pix + (avg_refs - org_pix) * factor) + 0.5f);
}

inline int process_mode7(
  const process_plane_params& params,
  int src_px,
  float ref_h1,
  float ref_h2,
  float ref_w1,
  float ref_w2,
  float angle_org,
  float angle_ref_h1,
  float angle_ref_h2,
  float angle_ref_w1,
  float angle_ref_w2
) {
  float max_angle_diff = 0.0f;
  max_angle_diff = std::max(max_angle_diff, std::abs(angle_ref_h1 - angle_org));
  max_angle_diff = std::max(max_angle_diff, std::abs(angle_ref_h2 - angle_org));
  max_angle_diff = std::max(max_angle_diff, std::abs(angle_ref_w1 - angle_org));
  max_angle_diff = std::max(max_angle_diff, std::abs(angle_ref_w2 - angle_org));

  float threshold = static_cast<float>(params.threshold);
  float threshold1 = static_cast<float>(params.threshold1);
  float threshold2 = static_cast<float>(params.threshold2);

  if (max_angle_diff <= params.max_angle) {
    threshold *= params.angle_boost;
    threshold1 *= params.angle_boost;
    threshold2 *= params.angle_boost;
  }

  return process_mode6(src_px, ref_h1, ref_h2, ref_w1, ref_w2, threshold, threshold1, threshold2);
}

template <int kSampleMode, bool kBlurFirst, class PixelIn>
int process_pixel(
  const process_plane_params& params,
  StridedPlaneView<const PixelIn> src_plane,
  int row,
  int col
) {
  const ReferenceSamples samples =
    load_reference_samples<kSampleMode>(params, src_plane, row, col);

  if constexpr (kSampleMode == 1) {
    return process_mode1<kBlurFirst>(samples.src, samples.ref1, samples.ref2, params.threshold);
  } else if constexpr (kSampleMode == 2) {
    return process_mode2<kBlurFirst>(
      samples.src,
      samples.ref1,
      samples.ref2,
      samples.ref3,
      samples.ref4,
      params.threshold
    );
  } else if constexpr (kSampleMode == 3) {
    return process_mode3<kBlurFirst>(samples.src, samples.ref1, samples.ref2, params.threshold);
  } else if constexpr (kSampleMode == 4) {
    return process_mode4<kBlurFirst>(
      samples.src,
      samples.ref1,
      samples.ref2,
      samples.ref3,
      samples.ref4,
      params.threshold
    );
  } else if constexpr (kSampleMode == 5) {
    return process_mode5(
      samples.src,
      samples.ref1,
      samples.ref2,
      samples.ref3,
      samples.ref4,
      params.threshold,
      params.threshold1,
      params.threshold2
    );
  } else if constexpr (kSampleMode == 6 || kSampleMode == 7) {
    if constexpr (kSampleMode == 6) {
      return process_mode6(
        samples.src,
        static_cast<float>(samples.ref1),
        static_cast<float>(samples.ref2),
        static_cast<float>(samples.ref3),
        static_cast<float>(samples.ref4),
        static_cast<float>(params.threshold),
        static_cast<float>(params.threshold1),
        static_cast<float>(params.threshold2)
      );
    } else {
      const float angle_org = calculate_gradient_angle(params, src_plane, col, row);
      const float angle_ref_h1 =
        calculate_gradient_angle(params, src_plane, col, row + samples.ref_y);
      const float angle_ref_h2 =
        calculate_gradient_angle(params, src_plane, col, row - samples.ref_y);
      const float angle_ref_w1 =
        calculate_gradient_angle(params, src_plane, col + samples.ref_x, row);
      const float angle_ref_w2 =
        calculate_gradient_angle(params, src_plane, col - samples.ref_x, row);
      return process_mode7(
        params,
        samples.src,
        static_cast<float>(samples.ref1),
        static_cast<float>(samples.ref2),
        static_cast<float>(samples.ref3),
        static_cast<float>(samples.ref4),
        angle_org,
        angle_ref_h1,
        angle_ref_h2,
        angle_ref_w1,
        angle_ref_w2
      );
    }
  } else {
    static_assert(kSampleMode >= 1 && kSampleMode <= 7);
  }
}

} // namespace neo_f3kdb::core::sample_modes
