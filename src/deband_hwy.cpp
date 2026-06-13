#include "deband_hwy.hpp"
#include "constants.h"
#include "impl_dispatch.h"
#include "deband_math.hpp"

// Standard dither headers included globally at top with robust header guards
#include "pixel_proc_c_high_ordered_dithering.h"
#include "pixel_proc_c_high_f_s_dithering.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "deband_hwy.cpp"

#include "hwy/foreach_target.h"
#include "hwy/highway.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <span>

HWY_BEFORE_NAMESPACE();
namespace neo_f3kdb::HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

template <int kSampleMode, bool kBlurFirst, int kDitherAlgo, class PixelIn, class PixelOut>
void process_plane_templated(const process_plane_params& params, void* dither_context) {
    int height = params.plane_height_in_pixels;
    int width = params.plane_width_in_pixels;
    int input_depth = params.input_depth;
    int output_depth = params.output_depth;
    int pixel_min = params.pixel_min;
    int pixel_max = params.pixel_max;
    std::uint16_t threshold = params.threshold;
    std::uint16_t threshold1 = params.threshold1;
    std::uint16_t threshold2 = params.threshold2;
    float angle_boost = params.angle_boost;
    float max_angle = params.max_angle;

    int src_stride = params.src_pitch / sizeof(PixelIn);
    int dst_stride = params.dst_pitch / sizeof(PixelOut);

    const PixelIn* src_base = reinterpret_cast<const PixelIn*>(params.src_plane_ptr);
    PixelOut* dst_base = reinterpret_cast<PixelOut*>(params.dst_plane_ptr);

    int width_subsamp = params.width_subsampling;
    int height_subsamp = params.height_subsampling;

    char fsdither_context[CONTEXT_BUFFER_SIZE];
    void* active_dither_ctx = dither_context;
    if constexpr (kDitherAlgo == 3) {
        pixel_proc_high_f_s_dithering::init_context(fsdither_context, width, output_depth);
        active_dither_ctx = fsdither_context;
    }

    for (int i = 0; i < height; i++) {
        const PixelIn* src_row = src_base + i * src_stride;
        PixelOut* dst_row = dst_base + i * dst_stride;

        const std::int16_t* grain_row = params.grain_buffer + params.grain_buffer_stride * i;
        const pixel_dither_info* info_row = params.info_ptr_base + params.info_stride * i;

        for (int j = 0; j < width; j++) {
            pixel_dither_info info = info_row[j];

            int src_px_up = static_cast<int>(src_row[j]);
            src_px_up <<= (16 - input_depth);

            int avg = 0;
            bool use_org_px_as_base = true;
            int new_pixel = src_px_up;

            if constexpr (kSampleMode == 1) {
                int ref_y = info.ref1 >> height_subsamp;
                int y1 = std::clamp(i - ref_y, 0, height - 1);
                int y2 = std::clamp(i + ref_y, 0, height - 1);

                int ref_1_up = static_cast<int>(src_base[y1 * src_stride + j]) << (16 - input_depth);
                int ref_2_up = static_cast<int>(src_base[y2 * src_stride + j]) << (16 - input_depth);

                avg = (ref_1_up + ref_2_up + 1) >> 1;

                if (kBlurFirst) {
                    int diff = avg - src_px_up;
                    use_org_px_as_base = std::abs(diff) >= threshold;
                } else {
                    int diff1 = src_px_up - ref_1_up;
                    int diff2 = src_px_up - ref_2_up;
                    use_org_px_as_base = (std::abs(diff1) >= threshold) || (std::abs(diff2) >= threshold);
                }
                new_pixel = use_org_px_as_base ? src_px_up : avg;
            }
            else if constexpr (kSampleMode == 2) {
                int ref_y1 = info.ref2 >> height_subsamp;
                int ref_y2 = info.ref1 >> height_subsamp;
                int ref_x1 = info.ref1 >> width_subsamp;
                int ref_x2 = info.ref2 >> width_subsamp;

                int y1 = std::clamp(i - ref_y1, 0, height - 1);
                int y2 = std::clamp(i - ref_y2, 0, height - 1);
                int y3 = std::clamp(i + ref_y1, 0, height - 1);
                int y4 = std::clamp(i + ref_y2, 0, height - 1);

                int x1 = std::clamp(j - ref_x1, 0, width - 1);
                int x2 = std::clamp(j + ref_x2, 0, width - 1);
                int x3 = std::clamp(j + ref_x1, 0, width - 1);
                int x4 = std::clamp(j - ref_x2, 0, width - 1);

                int ref_1_up = static_cast<int>(src_base[y1 * src_stride + x1]) << (16 - input_depth);
                int ref_2_up = static_cast<int>(src_base[y2 * src_stride + x2]) << (16 - input_depth);
                int ref_3_up = static_cast<int>(src_base[y3 * src_stride + x3]) << (16 - input_depth);
                int ref_4_up = static_cast<int>(src_base[y4 * src_stride + x4]) << (16 - input_depth);

                int avg1 = (ref_1_up + ref_2_up + 1) >> 1;
                int avg2 = (ref_3_up + ref_4_up + 1) >> 1;
                if (avg1 > 0) {
                    avg1 -= 1;
                }
                avg = (avg1 + avg2 + 1) >> 1;

                if (kBlurFirst) {
                    int diff = avg - src_px_up;
                    use_org_px_as_base = std::abs(diff) >= threshold;
                } else {
                    int diff1 = ref_1_up - src_px_up;
                    int diff2 = ref_2_up - src_px_up;
                    int diff3 = ref_3_up - src_px_up;
                    int diff4 = ref_4_up - src_px_up;
                    use_org_px_as_base = (std::abs(diff1) >= threshold) ||
                                         (std::abs(diff2) >= threshold) ||
                                         (std::abs(diff3) >= threshold) ||
                                         (std::abs(diff4) >= threshold);
                }
                new_pixel = use_org_px_as_base ? src_px_up : avg;
            }
            else if constexpr (kSampleMode == 3) {
                int ref_x = info.ref1 >> width_subsamp;
                int x1 = std::clamp(j - ref_x, 0, width - 1);
                int x2 = std::clamp(j + ref_x, 0, width - 1);

                int ref_1_up = static_cast<int>(src_row[x1]) << (16 - input_depth);
                int ref_2_up = static_cast<int>(src_row[x2]) << (16 - input_depth);

                avg = (ref_1_up + ref_2_up + 1) >> 1;

                if (kBlurFirst) {
                    int diff = avg - src_px_up;
                    use_org_px_as_base = std::abs(diff) >= threshold;
                } else {
                    int diff1 = src_px_up - ref_1_up;
                    int diff2 = src_px_up - ref_2_up;
                    use_org_px_as_base = (std::abs(diff1) >= threshold) || (std::abs(diff2) >= threshold);
                }
                new_pixel = use_org_px_as_base ? src_px_up : avg;
            }
            else if constexpr (kSampleMode == 4) {
                int ref_y = info.ref1 >> height_subsamp;
                int y1 = std::clamp(i - ref_y, 0, height - 1);
                int y2 = std::clamp(i + ref_y, 0, height - 1);
                int ref_v1 = static_cast<int>(src_base[y1 * src_stride + j]) << (16 - input_depth);
                int ref_v2 = static_cast<int>(src_base[y2 * src_stride + j]) << (16 - input_depth);
                int avg_v = (ref_v1 + ref_v2 + 1) >> 1;

                int ref_x = info.ref1 >> width_subsamp;
                int x1 = std::clamp(j - ref_x, 0, width - 1);
                int x2 = std::clamp(j + ref_x, 0, width - 1);
                int ref_h1 = static_cast<int>(src_row[x1]) << (16 - input_depth);
                int ref_h2 = static_cast<int>(src_row[x2]) << (16 - input_depth);
                int avg_h = (ref_h1 + ref_h2 + 1) >> 1;

                bool use_org_v = false;
                if (kBlurFirst) {
                    use_org_v = std::abs(avg_v - src_px_up) >= threshold;
                } else {
                    use_org_v = (std::abs(src_px_up - ref_v1) >= threshold) || (std::abs(src_px_up - ref_v2) >= threshold);
                }
                int new_v = use_org_v ? src_px_up : avg_v;

                bool use_org_h = false;
                if (kBlurFirst) {
                    use_org_h = std::abs(avg_h - src_px_up) >= threshold;
                } else {
                    use_org_h = (std::abs(src_px_up - ref_h1) >= threshold) || (std::abs(src_px_up - ref_h2) >= threshold);
                }
                int new_h = use_org_h ? src_px_up : avg_h;

                new_pixel = (new_v + new_h + 1) >> 1;
            }
            else if constexpr (kSampleMode == 5) {
                int ref_y = info.ref1 >> height_subsamp;
                int y1 = std::clamp(i - ref_y, 0, height - 1);
                int y2 = std::clamp(i + ref_y, 0, height - 1);
                int ref_1_h = static_cast<int>(src_base[y1 * src_stride + j]) << (16 - input_depth);
                int ref_2_h = static_cast<int>(src_base[y2 * src_stride + j]) << (16 - input_depth);

                int ref_x = info.ref1 >> width_subsamp;
                int x1 = std::clamp(j - ref_x, 0, width - 1);
                int x2 = std::clamp(j + ref_x, 0, width - 1);
                int ref_1_w = static_cast<int>(src_row[x1]) << (16 - input_depth);
                int ref_2_w = static_cast<int>(src_row[x2]) << (16 - input_depth);

                int avg1 = (ref_1_h + ref_2_h + 1) >> 1;
                int avg2 = (ref_1_w + ref_2_w + 1) >> 1;
                if (avg1 > 0) {
                    avg1 -= 1;
                }
                avg = (avg1 + avg2 + 1) >> 1;

                int avgDif = std::abs(avg - src_px_up);
                int maxDif = std::max({
                    std::abs(ref_1_h - src_px_up),
                    std::abs(ref_2_h - src_px_up),
                    std::abs(ref_1_w - src_px_up),
                    std::abs(ref_2_w - src_px_up)
                });
                int midDif1 = std::abs(ref_1_h + ref_2_h - 2 * src_px_up);
                int midDif2 = std::abs(ref_1_w + ref_2_w - 2 * src_px_up);

                use_org_px_as_base = (avgDif >= threshold) ||
                                     (maxDif >= threshold1) ||
                                     (midDif1 >= threshold2) ||
                                     (midDif2 >= threshold2);

                new_pixel = use_org_px_as_base ? src_px_up : avg;
            }
            else if constexpr (kSampleMode == 6 || kSampleMode == 7) {
                float org_pix_f = static_cast<float>(src_px_up);

                int ref_y = info.ref1 >> height_subsamp;
                int y1 = std::clamp(i - ref_y, 0, height - 1);
                int y2 = std::clamp(i + ref_y, 0, height - 1);
                float ref_1_h_f = static_cast<float>(static_cast<int>(src_base[y1 * src_stride + j]) << (16 - input_depth));
                float ref_2_h_f = static_cast<float>(static_cast<int>(src_base[y2 * src_stride + j]) << (16 - input_depth));

                int ref_x = info.ref1 >> width_subsamp;
                int x1 = std::clamp(j - ref_x, 0, width - 1);
                int x2 = std::clamp(j + ref_x, 0, width - 1);
                float ref_1_w_f = static_cast<float>(static_cast<int>(src_row[x1]) << (16 - input_depth));
                float ref_2_w_f = static_cast<float>(static_cast<int>(src_row[x2]) << (16 - input_depth));

                float current_thresh_avg_dif = static_cast<float>(threshold);
                float current_thresh_max_dif = static_cast<float>(threshold1);
                float current_thresh_mid_dif = static_cast<float>(threshold2);

                if constexpr (kSampleMode == 7) {
                    float angle_org = calculate_gradient_angle(params, src_base, src_stride, j, i, input_depth);
                    float angle_ref1_h = calculate_gradient_angle(params, src_base, src_stride, j, i + ref_y, input_depth);
                    float angle_ref2_h = calculate_gradient_angle(params, src_base, src_stride, j, i - ref_y, input_depth);
                    float angle_ref1_w = calculate_gradient_angle(params, src_base, src_stride, j + ref_x, i, input_depth);
                    float angle_ref2_w = calculate_gradient_angle(params, src_base, src_stride, j - ref_x, i, input_depth);

                    float max_angle_diff = 0.0f;
                    max_angle_diff = std::max(max_angle_diff, std::abs(angle_ref1_h - angle_org));
                    max_angle_diff = std::max(max_angle_diff, std::abs(angle_ref2_h - angle_org));
                    max_angle_diff = std::max(max_angle_diff, std::abs(angle_ref1_w - angle_org));
                    max_angle_diff = std::max(max_angle_diff, std::abs(angle_ref2_w - angle_org));

                    if (max_angle_diff <= max_angle) {
                        current_thresh_avg_dif *= angle_boost;
                        current_thresh_max_dif *= angle_boost;
                        current_thresh_mid_dif *= angle_boost;
                    }
                }

                float avg_refs_f = (ref_1_h_f + ref_2_h_f + ref_1_w_f + ref_2_w_f) * 0.25f;
                float avg_dif_f = std::abs(avg_refs_f - org_pix_f);
                float max_dif_f = std::max({
                    std::abs(ref_1_h_f - org_pix_f),
                    std::abs(ref_2_h_f - org_pix_f),
                    std::abs(ref_1_w_f - org_pix_f),
                    std::abs(ref_2_w_f - org_pix_f)
                });
                float mid_dif_v_f = std::abs(ref_1_h_f + ref_2_h_f - 2.0f * org_pix_f);
                float mid_dif_h_f = std::abs(ref_1_w_f + ref_2_w_f - 2.0f * org_pix_f);

                float factor = std::pow(
                    saturate(3.0f * calculate_ratio_term(avg_dif_f, current_thresh_avg_dif)) *
                    saturate(3.0f * calculate_ratio_term(max_dif_f, current_thresh_max_dif)) *
                    saturate(3.0f * calculate_ratio_term(mid_dif_v_f, current_thresh_mid_dif)) *
                    saturate(3.0f * calculate_ratio_term(mid_dif_h_f, current_thresh_mid_dif)),
                    0.1f
                );

                new_pixel = static_cast<int>((org_pix_f + (avg_refs_f - org_pix_f) * factor) + 0.5f);
            }

            new_pixel += static_cast<int>(grain_row[j]);

            if constexpr (kDitherAlgo == 1) {
                new_pixel = std::clamp(new_pixel, pixel_min, pixel_max) >> (16 - output_depth);
            }
            else if constexpr (kDitherAlgo == 2) {
                new_pixel += (pixel_proc_high_ordered_dithering::THRESHOLD_MAP[i & 15][j & 15] >> (output_depth - 8));
                new_pixel = std::clamp(new_pixel, pixel_min, pixel_max) >> (16 - output_depth);
            }
            else if constexpr (kDitherAlgo == 3) {
                new_pixel = pixel_proc_high_f_s_dithering::dither(active_dither_ctx, new_pixel, i, j);
                new_pixel = std::clamp(new_pixel, pixel_min, pixel_max) >> (16 - output_depth);
                pixel_proc_high_f_s_dithering::next_pixel(active_dither_ctx);
            }
            else if constexpr (kDitherAlgo == 4) {
                new_pixel = std::clamp(new_pixel, pixel_min, pixel_max);
            }

            dst_row[j] = static_cast<PixelOut>(new_pixel);
        }

        if constexpr (kDitherAlgo == 3) {
            pixel_proc_high_f_s_dithering::next_row(active_dither_ctx);
        }
    }

    if constexpr (kDitherAlgo == 3) {
        pixel_proc_high_f_s_dithering::destroy_context(active_dither_ctx);
    }
}

template <int kSampleMode, bool kBlurFirst, class PixelIn, class PixelOut>
void dispatch_dither(const process_plane_params& params, void* dither_context) {
    switch (params.dither_algo) {
        case 1: process_plane_templated<kSampleMode, kBlurFirst, 1, PixelIn, PixelOut>(params, dither_context); break;
        case 2: process_plane_templated<kSampleMode, kBlurFirst, 2, PixelIn, PixelOut>(params, dither_context); break;
        case 3: process_plane_templated<kSampleMode, kBlurFirst, 3, PixelIn, PixelOut>(params, dither_context); break;
        case 4: process_plane_templated<kSampleMode, kBlurFirst, 4, PixelIn, PixelOut>(params, dither_context); break;
        default: abort();
    }
}

template <int kSampleMode, class PixelIn, class PixelOut>
void dispatch_blur(const process_plane_params& params, void* dither_context) {
    if (params.blur_first) {
        dispatch_dither<kSampleMode, true, PixelIn, PixelOut>(params, dither_context);
    } else {
        dispatch_dither<kSampleMode, false, PixelIn, PixelOut>(params, dither_context);
    }
}

template <class PixelIn, class PixelOut>
void dispatch_sample_mode(const process_plane_params& params, void* dither_context) {
    switch (params.sample_mode) {
        case 1: dispatch_blur<1, PixelIn, PixelOut>(params, dither_context); break;
        case 2: dispatch_blur<2, PixelIn, PixelOut>(params, dither_context); break;
        case 3: dispatch_blur<3, PixelIn, PixelOut>(params, dither_context); break;
        case 4: dispatch_blur<4, PixelIn, PixelOut>(params, dither_context); break;
        case 5: dispatch_blur<5, PixelIn, PixelOut>(params, dither_context); break;
        case 6: dispatch_blur<6, PixelIn, PixelOut>(params, dither_context); break;
        case 7: dispatch_blur<7, PixelIn, PixelOut>(params, dither_context); break;
        default: abort();
    }
}

template <class PixelIn>
void dispatch_pixel_out(const process_plane_params& params, void* dither_context) {
    if (params.output_depth <= 8) {
        dispatch_sample_mode<PixelIn, std::uint8_t>(params, dither_context);
    } else {
        dispatch_sample_mode<PixelIn, std::uint16_t>(params, dither_context);
    }
}

void ProcessPlaneHWY(const process_plane_params& params, process_plane_context* context, [[maybe_unused]] process_plane_impl_t old_impl) {
    if (params.input_depth == 8) {
        dispatch_pixel_out<std::uint8_t>(params, context);
    } else {
        dispatch_pixel_out<std::uint16_t>(params, context);
    }
}

} // namespace neo_f3kdb::HWY_NAMESPACE
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
#include "hwy/per_target.h"

namespace neo_f3kdb {

HWY_EXPORT(ProcessPlaneHWY);

void process_plane_hwy(const process_plane_params& params, process_plane_context* context, process_plane_impl_t old_impl) {
  HWY_DYNAMIC_DISPATCH(ProcessPlaneHWY)(params, context, old_impl);
}

} // namespace neo_f3kdb
#endif
