#pragma once

#include "core/plane.hpp"

#include <algorithm>
#include <cmath>

namespace neo_f3kdb::core::sample_modes {

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

template <bool kBlurFirst, class SampleOps>
inline int process_mode1(const SampleOps& ops, int src_px, int ref1, int ref2, int threshold) {
  const int avg = ops.avg2(ref1, ref2);
  const bool use_src = kBlurFirst
    ? over_threshold(threshold, avg - src_px)
    : over_threshold(threshold, src_px - ref1) ||
        over_threshold(threshold, src_px - ref2);
  return use_src ? src_px : avg;
}

template <bool kBlurFirst, class SampleOps>
inline int process_mode2(
  const SampleOps& ops,
  int src_px,
  int ref1,
  int ref2,
  int ref3,
  int ref4,
  int threshold
) {
  const int avg = ops.avg4(ref1, ref2, ref3, ref4);
  const bool use_src = kBlurFirst
    ? over_threshold(threshold, avg - src_px)
    : over_threshold(threshold, ref1 - src_px) ||
        over_threshold(threshold, ref2 - src_px) ||
        over_threshold(threshold, ref3 - src_px) ||
        over_threshold(threshold, ref4 - src_px);
  return use_src ? src_px : avg;
}

template <bool kBlurFirst, class SampleOps>
inline int process_mode3(const SampleOps& ops, int src_px, int ref1, int ref2, int threshold) {
  return process_mode1<kBlurFirst>(ops, src_px, ref1, ref2, threshold);
}

template <bool kBlurFirst, class SampleOps>
inline int process_mode4(
  const SampleOps& ops,
  int src_px,
  int ref_v1,
  int ref_v2,
  int ref_h1,
  int ref_h2,
  int threshold
) {
  const int avg_v = ops.avg2(ref_v1, ref_v2);
  const bool use_src_v = kBlurFirst
    ? over_threshold(threshold, avg_v - src_px)
    : over_threshold(threshold, src_px - ref_v1) ||
        over_threshold(threshold, src_px - ref_v2);
  const int new_v = use_src_v ? src_px : avg_v;

  const int avg_h = ops.avg2(ref_h1, ref_h2);
  const bool use_src_h = kBlurFirst
    ? over_threshold(threshold, avg_h - src_px)
    : over_threshold(threshold, src_px - ref_h1) ||
        over_threshold(threshold, src_px - ref_h2);
  const int new_h = use_src_h ? src_px : avg_h;

  return ops.avg2(new_v, new_h);
}

template <class SampleOps>
inline int process_mode5(
  const SampleOps& ops,
  int src_px,
  int ref_h1,
  int ref_h2,
  int ref_w1,
  int ref_w2,
  int threshold,
  int threshold1,
  int threshold2
) {
  const int avg = ops.avg4(ref_h1, ref_h2, ref_w1, ref_w2);
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

  float threshold = static_cast<float>(params.config.threshold);
  float threshold1 = static_cast<float>(params.config.threshold1);
  float threshold2 = static_cast<float>(params.config.threshold2);

  if (max_angle_diff <= params.config.max_angle) {
    threshold *= params.config.angle_boost;
    threshold1 *= params.config.angle_boost;
    threshold2 *= params.config.angle_boost;
  }

  return process_mode6(src_px, ref_h1, ref_h2, ref_w1, ref_w2, threshold, threshold1, threshold2);
}

} // namespace neo_f3kdb::core::sample_modes
