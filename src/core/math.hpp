#pragma once

#include "core/plane.hpp"

#include <algorithm>
#include <cstdlib>

namespace neo_f3kdb::core::deband_scalar {

template <class PixelIn>
inline int upsample(PixelIn pixel, int input_depth) {
  return static_cast<int>(pixel) << (16 - input_depth);
}

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

template <int kSampleMode, bool kBlurFirst, class PixelIn>
int process_pixel(
  const process_plane_params& params,
  const PixelIn* src_base,
  const PixelIn* src_row,
  int src_stride,
  int row,
  int col
) {
  const int height = params.plane_height_in_pixels;
  const int width = params.plane_width_in_pixels;
  const int input_depth = params.input_depth;
  const int width_subsamp = params.width_subsampling;
  const int height_subsamp = params.height_subsampling;
  const auto info = params.info_ptr_base[static_cast<intptr_t>(params.info_stride) * row + col];

  const int src_px = upsample(src_row[col], input_depth);

  if constexpr (kSampleMode == 1) {
    const int ref_y = info.ref1 >> height_subsamp;
    const int ref1 = upsample(
      src_base[std::clamp(row - ref_y, 0, height - 1) * src_stride + col],
      input_depth
    );
    const int ref2 = upsample(
      src_base[std::clamp(row + ref_y, 0, height - 1) * src_stride + col],
      input_depth
    );
    return process_mode1<kBlurFirst>(src_px, ref1, ref2, params.threshold);
  } else if constexpr (kSampleMode == 2) {
    const int ref_y1 = info.ref2 >> height_subsamp;
    const int ref_y2 = info.ref1 >> height_subsamp;
    const int ref_x1 = info.ref1 >> width_subsamp;
    const int ref_x2 = info.ref2 >> width_subsamp;
    const int ref1 = upsample(
      src_base[
        std::clamp(row - ref_y1, 0, height - 1) * src_stride +
        std::clamp(col - ref_x1, 0, width - 1)
      ],
      input_depth
    );
    const int ref2 = upsample(
      src_base[
        std::clamp(row - ref_y2, 0, height - 1) * src_stride +
        std::clamp(col + ref_x2, 0, width - 1)
      ],
      input_depth
    );
    const int ref3 = upsample(
      src_base[
        std::clamp(row + ref_y1, 0, height - 1) * src_stride +
        std::clamp(col + ref_x1, 0, width - 1)
      ],
      input_depth
    );
    const int ref4 = upsample(
      src_base[
        std::clamp(row + ref_y2, 0, height - 1) * src_stride +
        std::clamp(col - ref_x2, 0, width - 1)
      ],
      input_depth
    );
    return process_mode2<kBlurFirst>(src_px, ref1, ref2, ref3, ref4, params.threshold);
  } else if constexpr (kSampleMode == 3) {
    const int ref_x = info.ref1 >> width_subsamp;
    const int ref1 = upsample(src_row[std::clamp(col - ref_x, 0, width - 1)], input_depth);
    const int ref2 = upsample(src_row[std::clamp(col + ref_x, 0, width - 1)], input_depth);
    return process_mode3<kBlurFirst>(src_px, ref1, ref2, params.threshold);
  } else if constexpr (kSampleMode == 4) {
    const int ref_y = info.ref1 >> height_subsamp;
    const int ref_x = info.ref1 >> width_subsamp;
    const int ref_v1 = upsample(
      src_base[std::clamp(row - ref_y, 0, height - 1) * src_stride + col],
      input_depth
    );
    const int ref_v2 = upsample(
      src_base[std::clamp(row + ref_y, 0, height - 1) * src_stride + col],
      input_depth
    );
    const int ref_h1 = upsample(src_row[std::clamp(col - ref_x, 0, width - 1)], input_depth);
    const int ref_h2 = upsample(src_row[std::clamp(col + ref_x, 0, width - 1)], input_depth);
    return process_mode4<kBlurFirst>(src_px, ref_v1, ref_v2, ref_h1, ref_h2, params.threshold);
  } else if constexpr (kSampleMode == 5) {
    const int ref_y = info.ref1 >> height_subsamp;
    const int ref_x = info.ref1 >> width_subsamp;
    const int ref_h1 = upsample(
      src_base[std::clamp(row - ref_y, 0, height - 1) * src_stride + col],
      input_depth
    );
    const int ref_h2 = upsample(
      src_base[std::clamp(row + ref_y, 0, height - 1) * src_stride + col],
      input_depth
    );
    const int ref_w1 = upsample(src_row[std::clamp(col - ref_x, 0, width - 1)], input_depth);
    const int ref_w2 = upsample(src_row[std::clamp(col + ref_x, 0, width - 1)], input_depth);
    return process_mode5(
      src_px,
      ref_h1,
      ref_h2,
      ref_w1,
      ref_w2,
      params.threshold,
      params.threshold1,
      params.threshold2
    );
  } else {
    static_assert(kSampleMode >= 1 && kSampleMode <= 5);
  }
}

} // namespace neo_f3kdb::core::deband_scalar
