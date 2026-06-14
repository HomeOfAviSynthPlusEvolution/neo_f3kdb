#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include "core/plane.hpp"
#include "core/pixel_proc.hpp"
#include <limits.h>

static inline bool _is_above_threshold(int threshold, int diff) {
    return std::abs(diff) >= threshold;
}

static inline bool is_above_threshold(int threshold, int diff1) {
    return _is_above_threshold(threshold, diff1);
}

static inline bool is_above_threshold(int threshold, int diff1, int diff2) {
    return _is_above_threshold(threshold, diff1) ||
           _is_above_threshold(threshold, diff2);
}

static inline bool is_above_threshold(int threshold, int diff1, int diff2, int diff3, int diff4) {
    return _is_above_threshold(threshold, diff1) ||
           _is_above_threshold(threshold, diff2) ||
           _is_above_threshold(threshold, diff3) ||
           _is_above_threshold(threshold, diff4);
}

inline int input_pixel_step(const process_plane_params& params)
{
    return params.input_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1;
}

template <int mode>
static __inline int read_pixel_value(const process_plane_params& params, void* context, const unsigned char* ptr)
{
    if (params.input_mode == LOW_BIT_DEPTH)
    {
        return neo_f3kdb::core::pixel_proc::upsample<mode>(context, *ptr);
    }

    int ret;

    ret = *(unsigned short*)ptr;
    ret <<= (INTERNAL_BIT_DEPTH - params.input_depth);
    return ret;
}

template <int mode>
static __inline int read_pixel_at(
    const process_plane_params& params,
    void* context,
    neo_f3kdb::core::StridedPlaneView<const unsigned char> src_plane,
    int row,
    int col
) {
    const auto src_row = src_plane.row(row);
    const auto* ptr = src_row.data() + static_cast<intptr_t>(col) * input_pixel_step(params);
    return read_pixel_value<mode>(params, context, ptr);
}

template <int mode>
static __inline int read_pixel_clamped(
    const process_plane_params& params,
    void* context,
    neo_f3kdb::core::StridedPlaneView<const unsigned char> src_plane,
    int row,
    int col
) {
    return read_pixel_at<mode>(
        params,
        context,
        src_plane,
        std::clamp(row, 0, params.plane_height() - 1),
        std::clamp(col, 0, params.plane_width() - 1)
    );
}

static __inline float saturate(float val)
{
    return std::clamp(val, 0.0f, 1.0f);
};

static __inline float calculate_ratio_term(float diff, float thresh)
{
    if (thresh < 1e-5f)
        return (std::abs(diff) < 1e-5f) ? 1.0f : -1e6f;

    return 1.0f - std::abs(diff) / thresh;
};

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

template <int pixel_proc_mode>
static float calculate_gradient_angle(
    const process_plane_params& params,
    void* context_pixel_proc,
    neo_f3kdb::core::StridedPlaneView<const unsigned char> src_plane,
    int current_x, int current_y, int read_distance = 20)
{
    auto get_pixel_value_at = [&](int x, int y) -> float {
        return static_cast<float>(read_pixel_clamped<pixel_proc_mode>(
            params,
            context_pixel_proc,
            src_plane,
            y,
            x
        ));
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

    const float scaled_epsilon_for_gx = 0.01f * (static_cast<float>(1 << (INTERNAL_BIT_DEPTH - params.input_depth)) * 3.0f);

    if (std::abs(gx) < scaled_epsilon_for_gx)
    {
        // gx is close to zero, gradient is predominantly vertical or area is flat.
        if (std::abs(gy) < scaled_epsilon_for_gx)
        {
            // Also flat in y direction
            return 1.0f;
        }

        // gx is negligible, gy is not. This is a vertical gradient.
        return 1.0f;
    }

    // gx is not close to zero, atan(gy/gx) is safe
    return std::atan(gy / gx) / static_cast<float>(M_PI) + 0.5f;
}

template <int sample_mode, bool blur_first, int mode, int output_mode>
static __forceinline void __cdecl process_plane_plainc_mode12_high(const process_plane_params& params, process_plane_context*)
{
    alignas(std::max_align_t) char context[CONTEXT_BUFFER_SIZE];

    unsigned short threshold = params.threshold;

    int pixel_min = params.pixel_min;
    int pixel_max = params.pixel_max;

    int width_subsamp = params.width_subsampling;

    neo_f3kdb::core::pixel_proc::init_context<mode>(context, params.plane_width(), params.output_depth);

    int dst_pixel_step = output_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1;

    int process_width = params.plane_width();

    const auto src_plane = params.src_bytes();
    auto dst_plane = params.dst_bytes();
    const auto grain_plane = params.grain_plane();
    const auto info_plane = params.dither_info_plane();

    for (int i = 0; i < params.plane_height(); i++)
    {
        auto dst_row = dst_plane.row(i);
        const auto grain_row = grain_plane.row(i);
        const auto info_row = info_plane.row(i);

        for (int j = 0; j < process_width; j++)
        {
            const auto column = static_cast<std::size_t>(j);
            unsigned char* dst_px = dst_row.data() + static_cast<intptr_t>(j) * dst_pixel_step;
            pixel_dither_info info = info_row[column];
            int src_px_up = read_pixel_at<mode>(params, context, src_plane, i, j);

            if constexpr (sample_mode == 1 || sample_mode == 2 || (sample_mode >= 4 && sample_mode <= 7))
            {
                assert(info.ref1 >= 0);
                assert((info.ref1 >> params.height_subsampling) <= i &&
                    (info.ref1 >> params.height_subsampling) + i < params.plane_height());
            }

            if constexpr (sample_mode >= 2 && sample_mode <= 7)
            {
                assert(info.ref2 >= 0);
                assert((info.ref2 >> params.height_subsampling) <= i &&
                       (info.ref2 >> params.height_subsampling) + i < params.plane_height());
            }
            int avg;
            bool use_org_px_as_base;
            int new_pixel = src_px_up, new_pixel_mode1, new_pixel_mode3;
            if constexpr (sample_mode == 1 || sample_mode == 4)
            {
                const int ref_y = info.ref1 >> params.height_subsampling;

                int ref_1_up = read_pixel_at<mode>(params, context, src_plane, i + ref_y, j);
                int ref_2_up = read_pixel_at<mode>(params, context, src_plane, i - ref_y, j);

                avg = neo_f3kdb::core::pixel_proc::avg_2<mode>(context, ref_1_up, ref_2_up);

                if (blur_first)
                {
                    int diff = avg - src_px_up;
                    use_org_px_as_base = is_above_threshold(threshold, diff);
                } else {
                    int diff = src_px_up - ref_1_up;
                    int diff_n = src_px_up - ref_2_up;
                    use_org_px_as_base = is_above_threshold(threshold, diff, diff_n);
                }
                new_pixel = new_pixel_mode1 = use_org_px_as_base ? src_px_up : avg;
            }
            if constexpr (sample_mode == 3 || sample_mode == 4)
            {
                const int ref_x = info.ref1 >> params.width_subsampling;

                int ref_1_up = read_pixel_at<mode>(params, context, src_plane, i, j + ref_x);
                int ref_2_up = read_pixel_at<mode>(params, context, src_plane, i, j - ref_x);

                avg = neo_f3kdb::core::pixel_proc::avg_2<mode>(context, ref_1_up, ref_2_up);

                if (blur_first)
                {
                    int diff = avg - src_px_up;
                    use_org_px_as_base = is_above_threshold(threshold, diff);
                } else {
                    int diff = src_px_up - ref_1_up;
                    int diff_n = src_px_up - ref_2_up;
                    use_org_px_as_base = is_above_threshold(threshold, diff, diff_n);
                }
                new_pixel = new_pixel_mode3 = use_org_px_as_base ? src_px_up : avg;
            }
            if constexpr (sample_mode == 4)
            {
                new_pixel = neo_f3kdb::core::pixel_proc::avg_2<mode>(context, new_pixel_mode1, new_pixel_mode3);
            }
            if constexpr (sample_mode == 2)
            {
                int x_multiplier = 1;
                const int ref_x1 = (info.ref1 * x_multiplier) >> width_subsamp;
                const int ref_x2 = (info.ref2 * x_multiplier) >> width_subsamp;
                const int ref_y1 = info.ref1 >> params.height_subsampling;
                const int ref_y2 = info.ref2 >> params.height_subsampling;

                assert(((info.ref1 >> width_subsamp) * x_multiplier) <= j &&
                       ((info.ref1 >> width_subsamp) * x_multiplier) + j < process_width);
                assert(((info.ref2 >> width_subsamp) * x_multiplier) <= j &&
                       ((info.ref2 >> width_subsamp) * x_multiplier) + j < process_width);

                int ref_1_up = read_pixel_at<mode>(params, context, src_plane, i + ref_y2, j + ref_x1);
                int ref_2_up = read_pixel_at<mode>(params, context, src_plane, i - ref_y1, j + ref_x2);
                int ref_3_up = read_pixel_at<mode>(params, context, src_plane, i - ref_y2, j - ref_x1);
                int ref_4_up = read_pixel_at<mode>(params, context, src_plane, i + ref_y1, j - ref_x2);

                avg = neo_f3kdb::core::pixel_proc::avg_4<mode>(context, ref_1_up, ref_2_up, ref_3_up, ref_4_up);

                if (blur_first)
                {
                    int diff = avg - src_px_up;
                    use_org_px_as_base = is_above_threshold(threshold, diff);
                } else {
                    int diff1 = ref_1_up - src_px_up;
                    int diff2 = ref_2_up - src_px_up;
                    int diff3 = ref_3_up - src_px_up;
                    int diff4 = ref_4_up - src_px_up;
                    use_org_px_as_base = is_above_threshold(threshold, diff1, diff2, diff3, diff4);
                }
                new_pixel = use_org_px_as_base ? src_px_up : avg;
            }
            if constexpr (sample_mode == 5)
            {
                const int ref_y = info.ref1 >> params.height_subsampling;
                const int ref_x = info.ref1 >> params.width_subsampling;

                int ref_1_h = read_pixel_at<mode>(params, context, src_plane, i + ref_y, j);
                int ref_2_h = read_pixel_at<mode>(params, context, src_plane, i - ref_y, j);

                int ref_1_w = read_pixel_at<mode>(params, context, src_plane, i, j + ref_x);
                int ref_2_w = read_pixel_at<mode>(params, context, src_plane, i, j - ref_x);

                const int avg = neo_f3kdb::core::pixel_proc::avg_4<mode>(context, ref_1_h, ref_2_h, ref_1_w, ref_2_w);
                const int avgDif = std::abs(avg - src_px_up);
                const int maxDif = std::max(std::abs(ref_1_h - src_px_up), std::max(std::abs(ref_2_h - src_px_up),
                    std::max(std::abs(ref_1_w - src_px_up), std::abs(ref_2_w - src_px_up))));
                const int midDif1 = std::abs(ref_1_h + ref_2_h - 2 * src_px_up);
                const int midDif2 = std::abs(ref_1_w + ref_2_w - 2 * src_px_up);
                use_org_px_as_base = is_above_threshold(threshold, avgDif) ||
                    is_above_threshold(params.threshold1, maxDif) ||
                    is_above_threshold(params.threshold2, midDif1) ||
                    is_above_threshold(params.threshold2, midDif2);

                new_pixel = use_org_px_as_base ? src_px_up : avg;
            }
            if constexpr (sample_mode == 6)
            {
                const float org_pix_f = static_cast<float>(src_px_up);
                const float thresh_avg_dif_param_f = static_cast<float>(params.threshold);
                const float thresh_max_dif_param_f = static_cast<float>(params.threshold1);
                const float thresh_mid_dif_param_f = static_cast<float>(params.threshold2);

                const int ref_y = info.ref1 >> params.height_subsampling;
                const int ref_x = info.ref1 >> params.width_subsampling;
                const float ref_1_h_f = static_cast<float>(read_pixel_at<mode>(params, context, src_plane, i + ref_y, j));
                const float ref_2_h_f = static_cast<float>(read_pixel_at<mode>(params, context, src_plane, i - ref_y, j));

                const float ref_1_w_f = static_cast<float>(read_pixel_at<mode>(params, context, src_plane, i, j + ref_x));
                const float ref_2_w_f = static_cast<float>(read_pixel_at<mode>(params, context, src_plane, i, j - ref_x));

                const float avg_refs_f = (ref_1_h_f + ref_2_h_f + ref_1_w_f + ref_2_w_f) * 0.25f;

                const float avg_dif_f = std::abs(avg_refs_f - org_pix_f);
                const float max_dif_f = std::max({ std::abs(ref_1_h_f - org_pix_f),
                                                  std::abs(ref_2_h_f - org_pix_f),
                                                  std::abs(ref_1_w_f - org_pix_f),
                                                  std::abs(ref_2_w_f - org_pix_f) });
                const float mid_dif_v_f = std::abs(ref_1_h_f + ref_2_h_f - 2.0f * org_pix_f);
                const float mid_dif_h_f = std::abs(ref_1_w_f + ref_2_w_f - 2.0f * org_pix_f);

                // Calculate the blending factor
                const float factor = std::pow(saturate(3.0f * calculate_ratio_term(avg_dif_f, thresh_avg_dif_param_f)) *
                    saturate(3.0f * calculate_ratio_term(max_dif_f, thresh_max_dif_param_f)) *
                    saturate(3.0f * calculate_ratio_term(mid_dif_v_f, thresh_mid_dif_param_f)) *
                    saturate(3.0f * calculate_ratio_term(mid_dif_h_f, thresh_mid_dif_param_f)), 0.1f);

                new_pixel = static_cast<int>((org_pix_f + (avg_refs_f - org_pix_f) * factor) + 0.5f);
            }
            if constexpr (sample_mode == 7)
            {
                // This mode is sample_mode=6 + gradient angle check

                const float org_pix_f = static_cast<float>(src_px_up);

                const int ref_y = info.ref1 >> params.height_subsampling;
                const int ref_x = info.ref1 >> params.width_subsampling;
                const float ref_1_h_f = static_cast<float>(read_pixel_at<mode>(params, context, src_plane, i + ref_y, j));
                const float ref_2_h_f = static_cast<float>(read_pixel_at<mode>(params, context, src_plane, i - ref_y, j));

                const float ref_1_w_f = static_cast<float>(read_pixel_at<mode>(params, context, src_plane, i, j + ref_x));
                const float ref_2_w_f = static_cast<float>(read_pixel_at<mode>(params, context, src_plane, i, j - ref_x));

                const float angle_org = calculate_gradient_angle<mode>(params, context, src_plane, j, i);

                const int ref1h_y_offset = (info.ref1 >> params.height_subsampling);
                const int ref1w_x_offset = (info.ref1 >> params.width_subsampling);

                const float angle_ref1_h = calculate_gradient_angle<mode>(params, context, src_plane, j, i + ref1h_y_offset);
                const float angle_ref2_h = calculate_gradient_angle<mode>(params, context, src_plane, j, i - ref1h_y_offset);
                const float angle_ref1_w = calculate_gradient_angle<mode>(params, context, src_plane, j + ref1w_x_offset, i);
                const float angle_ref2_w = calculate_gradient_angle<mode>(params, context, src_plane, j - ref1w_x_offset, i);

                float max_angle_diff = 0.0f;
                max_angle_diff = std::max(max_angle_diff, std::abs(angle_ref1_h - angle_org));
                max_angle_diff = std::max(max_angle_diff, std::abs(angle_ref2_h - angle_org));
                max_angle_diff = std::max(max_angle_diff, std::abs(angle_ref1_w - angle_org));
                max_angle_diff = std::max(max_angle_diff, std::abs(angle_ref2_w - angle_org));

                float current_thresh_avg_dif = static_cast<float>(params.threshold);
                float current_thresh_max_dif = static_cast<float>(params.threshold1);
                float current_thresh_mid_dif = static_cast<float>(params.threshold2);

                const float angle_boost_factor = params.angle_boost;
                const float max_angle_threshold = params.max_angle;

                if (max_angle_diff <= max_angle_threshold) {
                    current_thresh_avg_dif *= angle_boost_factor;
                    current_thresh_max_dif *= angle_boost_factor;
                    current_thresh_mid_dif *= angle_boost_factor;
                }

                const float avg_refs_f = (ref_1_h_f + ref_2_h_f + ref_1_w_f + ref_2_w_f) * 0.25f;
                const float avg_dif_f = std::abs(avg_refs_f - org_pix_f);
                const float max_dif_f = std::max({ std::abs(ref_1_h_f - org_pix_f),
                                                  std::abs(ref_2_h_f - org_pix_f),
                                                  std::abs(ref_1_w_f - org_pix_f),
                                                  std::abs(ref_2_w_f - org_pix_f) });
                const float mid_dif_v_f = std::abs(ref_1_h_f + ref_2_h_f - 2.0f * org_pix_f);
                const float mid_dif_h_f = std::abs(ref_1_w_f + ref_2_w_f - 2.0f * org_pix_f);

                const float factor = std::pow(
                    saturate(3.0f * calculate_ratio_term(avg_dif_f, current_thresh_avg_dif)) *
                    saturate(3.0f * calculate_ratio_term(max_dif_f, current_thresh_max_dif)) *
                    saturate(3.0f * calculate_ratio_term(mid_dif_v_f, current_thresh_mid_dif)) *
                    saturate(3.0f * calculate_ratio_term(mid_dif_h_f, current_thresh_mid_dif)),
                    0.1f);

                new_pixel = static_cast<int>((org_pix_f + (avg_refs_f - org_pix_f) * factor) + 0.5f);
            }

            new_pixel = neo_f3kdb::core::pixel_proc::downsample<mode>(
                context,
                new_pixel + grain_row[column],
                i,
                j,
                pixel_min,
                pixel_max,
                params.output_depth
            );

            switch (output_mode)
            {
            case LOW_BIT_DEPTH:
                *dst_px = (unsigned char)new_pixel;
                break;
            case HIGH_BIT_DEPTH_INTERLEAVED:
                *((unsigned short*)dst_px) = (unsigned short)(new_pixel & 0xFFFF);
                break;
            default:
                abort();
            }

            neo_f3kdb::core::pixel_proc::next_pixel<mode>(context);
        }
        neo_f3kdb::core::pixel_proc::next_row<mode>(context);
    }

    neo_f3kdb::core::pixel_proc::destroy_context<mode>(context);
}

template <int sample_mode, bool blur_first, int mode>
void __cdecl process_plane_plainc(const process_plane_params& params, process_plane_context* context)
{
    static_assert(sample_mode != 0, "No longer support sample_mode = 0");
    switch (params.output_mode)
    {
    case LOW_BIT_DEPTH:
        process_plane_plainc_mode12_high<sample_mode, blur_first, mode, LOW_BIT_DEPTH>(params, context);
        break;

    case HIGH_BIT_DEPTH_INTERLEAVED:
        process_plane_plainc_mode12_high<sample_mode, blur_first, mode, HIGH_BIT_DEPTH_INTERLEAVED>(params, context);
        break;

    default:
        abort();
    }
}
