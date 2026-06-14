#pragma once

#include "core/plane.hpp"
#include "core/sample_math.hpp"
#include "core/sample_reader.hpp"

namespace neo_f3kdb::core::sample_modes {

template <int kSampleMode, bool kBlurFirst, class SampleOps>
int process_pixel_with_ops(
  const process_plane_params& params,
  const SampleOps& ops,
  int row,
  int col
) {
  const ReferenceSamples samples =
    load_reference_samples<kSampleMode, SampleOps>(params, ops, row, col);

  if constexpr (kSampleMode == 1) {
    return process_mode1<kBlurFirst>(ops, samples.src, samples.ref1, samples.ref2, params.config.threshold);
  } else if constexpr (kSampleMode == 2) {
    return process_mode2<kBlurFirst>(
      ops,
      samples.src,
      samples.ref1,
      samples.ref2,
      samples.ref3,
      samples.ref4,
      params.config.threshold
    );
  } else if constexpr (kSampleMode == 3) {
    return process_mode3<kBlurFirst>(ops, samples.src, samples.ref1, samples.ref2, params.config.threshold);
  } else if constexpr (kSampleMode == 4) {
    return process_mode4<kBlurFirst>(
      ops,
      samples.src,
      samples.ref1,
      samples.ref2,
      samples.ref3,
      samples.ref4,
      params.config.threshold
    );
  } else if constexpr (kSampleMode == 5) {
    return process_mode5(
      ops,
      samples.src,
      samples.ref1,
      samples.ref2,
      samples.ref3,
      samples.ref4,
      params.config.threshold,
      params.config.threshold1,
      params.config.threshold2
    );
  } else if constexpr (kSampleMode == 6 || kSampleMode == 7) {
    if constexpr (kSampleMode == 6) {
      return process_mode6(
        samples.src,
        static_cast<float>(samples.ref1),
        static_cast<float>(samples.ref2),
        static_cast<float>(samples.ref3),
        static_cast<float>(samples.ref4),
        static_cast<float>(params.config.threshold),
        static_cast<float>(params.config.threshold1),
        static_cast<float>(params.config.threshold2)
      );
    } else {
      const float angle_org = calculate_gradient_angle(params, ops, col, row);
      const float angle_ref_h1 =
        calculate_gradient_angle(params, ops, col, row + samples.ref_y);
      const float angle_ref_h2 =
        calculate_gradient_angle(params, ops, col, row - samples.ref_y);
      const float angle_ref_w1 =
        calculate_gradient_angle(params, ops, col + samples.ref_x, row);
      const float angle_ref_w2 =
        calculate_gradient_angle(params, ops, col - samples.ref_x, row);
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

template <int kSampleMode, bool kBlurFirst, class PixelIn>
int process_pixel(
  const process_plane_params& params,
  StridedPlaneView<const PixelIn> src_plane,
  int row,
  int col
) {
  const TypedSampleOps<PixelIn> ops{params, src_plane};
  return process_pixel_with_ops<kSampleMode, kBlurFirst>(params, ops, row, col);
}

} // namespace neo_f3kdb::core::sample_modes
