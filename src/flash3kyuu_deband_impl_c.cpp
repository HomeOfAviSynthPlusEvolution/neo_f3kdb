#include <algorithm>
#include <cmath>
#include "core.h"
#include "pixel_proc_c.h"
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

template <int mode>
static __inline int read_pixel(const process_plane_params& params, void* context, const unsigned char* base, int offset = 0)
{
    const unsigned char* ptr = base + offset;
    if (params.input_mode == LOW_BIT_DEPTH)
    {
        return pixel_proc_upsample<mode>(context, *ptr);
    }

    int ret;

    ret = *(unsigned short*)ptr;
    ret <<= (INTERNAL_BIT_DEPTH - params.input_depth);
    return ret;
}

template <int sample_mode, bool blur_first, int mode, int output_mode>
static __forceinline void __cdecl process_plane_plainc_mode12_high(const process_plane_params& params, process_plane_context*)
{
    pixel_dither_info* info_ptr;
    char context[CONTEXT_BUFFER_SIZE];

    unsigned short threshold = params.threshold;

    int pixel_min = params.pixel_min;
    int pixel_max = params.pixel_max;

    int width_subsamp = params.width_subsampling;

    pixel_proc_init_context<mode>(context, params.plane_width_in_pixels, params.output_depth);

    int pixel_step = params.input_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1;

    int process_width = params.plane_width_in_pixels;

    for (int i = 0; i < params.plane_height_in_pixels; i++)
    {
        const unsigned char* src_px = params.src_plane_ptr + params.src_pitch * i;
        unsigned char* dst_px = params.dst_plane_ptr + params.dst_pitch * i;

        const short* grain_buffer_ptr = params.grain_buffer + params.grain_buffer_stride * i;

        info_ptr = params.info_ptr_base + params.info_stride * i;


        for (int j = 0; j < process_width; j++)
        {
            pixel_dither_info info = *info_ptr;
            int src_px_up = read_pixel<mode>(params, context, src_px);
            
            if constexpr (sample_mode == 1 || sample_mode == 2 || sample_mode == 4 || sample_mode == 5 || sample_mode == 6)
            {
                assert(info.ref1 >= 0);
                assert((info.ref1 >> params.height_subsampling) <= i && 
                    (info.ref1 >> params.height_subsampling) + i < params.plane_height_in_pixels);
            }

            if constexpr (sample_mode == 3 || sample_mode == 2 || sample_mode == 4 || sample_mode == 5 || sample_mode == 6)
            {
                assert(info.ref2 >= 0);
                assert((info.ref2 >> params.height_subsampling) <= i && 
                       (info.ref2 >> params.height_subsampling) + i < params.plane_height_in_pixels);
            }
            int avg;
            bool use_org_px_as_base;
            int ref_pos, ref_pos_2;
            int new_pixel = src_px_up, new_pixel_mode1, new_pixel_mode3;
            if constexpr (sample_mode == 1 || sample_mode == 4)
            {
                ref_pos = (info.ref1 >> params.height_subsampling) * params.src_pitch;

                int ref_1_up = read_pixel<mode>(params, context, src_px, ref_pos);
                int ref_2_up = read_pixel<mode>(params, context, src_px, -ref_pos);

                avg = pixel_proc_avg_2<mode>(context, ref_1_up, ref_2_up);

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
                ref_pos = (info.ref1 >> params.width_subsampling) * pixel_step;

                int ref_1_up = read_pixel<mode>(params, context, src_px, ref_pos);
                int ref_2_up = read_pixel<mode>(params, context, src_px, -ref_pos);

                avg = pixel_proc_avg_2<mode>(context, ref_1_up, ref_2_up);

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
                new_pixel = pixel_proc_avg_2<mode>(context, new_pixel_mode1, new_pixel_mode3);
            }
            if constexpr (sample_mode == 2)
            {
                int x_multiplier = 1;
                
                assert(((info.ref1 >> width_subsamp) * x_multiplier) <= j && 
                       ((info.ref1 >> width_subsamp) * x_multiplier) + j < process_width);
                assert(((info.ref2 >> width_subsamp) * x_multiplier) <= j && 
                       ((info.ref2 >> width_subsamp) * x_multiplier) + j < process_width);

                ref_pos = params.src_pitch * (info.ref2 >> params.height_subsampling) + 
                          ((info.ref1 * x_multiplier) >> width_subsamp) * pixel_step;

                ref_pos_2 = ((info.ref2 * x_multiplier) >> width_subsamp) * pixel_step - 
                            params.src_pitch * (info.ref1 >> params.height_subsampling);

                int ref_1_up = read_pixel<mode>(params, context, src_px, ref_pos);
                int ref_2_up = read_pixel<mode>(params, context, src_px, ref_pos_2);
                int ref_3_up = read_pixel<mode>(params, context, src_px, -ref_pos);
                int ref_4_up = read_pixel<mode>(params, context, src_px, -ref_pos_2);

                avg = pixel_proc_avg_4<mode>(context, ref_1_up, ref_2_up, ref_3_up, ref_4_up);

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
                ref_pos = (info.ref1 >> params.height_subsampling) * params.src_pitch;

                int ref_1_h = read_pixel<mode>(params, context, src_px, ref_pos);
                int ref_2_h = read_pixel<mode>(params, context, src_px, -ref_pos);

                ref_pos_2 = (info.ref1 >> params.width_subsampling) * pixel_step;

                int ref_1_w = read_pixel<mode>(params, context, src_px, ref_pos_2);
                int ref_2_w = read_pixel<mode>(params, context, src_px, -ref_pos_2);

                int avg = pixel_proc_avg_4<mode>(context, ref_1_h, ref_2_h, ref_1_w, ref_2_w);
                int avgDif = std::abs(avg - src_px_up);
                int maxDif = std::max(std::abs(ref_1_h - src_px_up), std::max(std::abs(ref_2_h - src_px_up), std::max(std::abs(ref_1_w - src_px_up), std::abs(ref_2_w - src_px_up))));
                int midDif1 = std::abs(ref_1_h + ref_2_h - 2 * src_px_up);
                int midDif2 = std::abs(ref_1_w + ref_2_w - 2 * src_px_up);
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

                const int ref_v_offset_bytes = (info.ref1 >> params.height_subsampling) * params.src_pitch;
                const float ref_1_h_f = static_cast<float>(read_pixel<mode>(params, context, src_px, ref_v_offset_bytes));
                const float ref_2_h_f = static_cast<float>(read_pixel<mode>(params, context, src_px, -ref_v_offset_bytes));

                const int ref_h_offset_bytes = (info.ref1 >> params.width_subsampling) * pixel_step;
                const float ref_1_w_f = static_cast<float>(read_pixel<mode>(params, context, src_px, ref_h_offset_bytes));
                const float ref_2_w_f = static_cast<float>(read_pixel<mode>(params, context, src_px, -ref_h_offset_bytes));

                const float avg_refs_f = (ref_1_h_f + ref_2_h_f + ref_1_w_f + ref_2_w_f) * 0.25f;

                const float avg_dif_f = std::abs(avg_refs_f - org_pix_f);
                const float max_dif_f = std::max({ std::abs(ref_1_h_f - org_pix_f),
                                                  std::abs(ref_2_h_f - org_pix_f),
                                                  std::abs(ref_1_w_f - org_pix_f),
                                                  std::abs(ref_2_w_f - org_pix_f) });
                const float mid_dif_v_f = std::abs(ref_1_h_f + ref_2_h_f - 2.0f * org_pix_f);
                const float mid_dif_h_f = std::abs(ref_1_w_f + ref_2_w_f - 2.0f * org_pix_f);

                auto saturate = [](float val) {
                    return std::clamp(val, 0.0f, 1.0f);
                    };

                auto calculate_ratio_term = [](float diff, float thresh) {
                    if (thresh < 1e-5f)
                        return (diff < 1e-5f) ? 1.0f : -1e6f;

                    return 1.0f - diff / thresh;
                    };

                // Calculate the blending factor
                float factor = std::pow(saturate(3.0f * calculate_ratio_term(avg_dif_f, thresh_avg_dif_param_f)) *
                    saturate(3.0f * calculate_ratio_term(max_dif_f, thresh_max_dif_param_f)) *
                    saturate(3.0f * calculate_ratio_term(mid_dif_v_f, thresh_mid_dif_param_f)) *
                    saturate(3.0f * calculate_ratio_term(mid_dif_h_f, thresh_mid_dif_param_f)), 0.1f);

                new_pixel = static_cast<int>((org_pix_f + (avg_refs_f - org_pix_f) * factor) + 0.5f);
            }

            new_pixel = pixel_proc_downsample<mode>(context, new_pixel + *grain_buffer_ptr, i, j, pixel_min, pixel_max, params.output_depth);

            switch (output_mode)
            {
            case LOW_BIT_DEPTH:
                *dst_px = (unsigned char)new_pixel;
                break;
            case HIGH_BIT_DEPTH_INTERLEAVED:
                *((unsigned short*)dst_px) = (unsigned short)(new_pixel & 0xFFFF);
                dst_px++;
                break;
            default:
                abort();
            }

            src_px += pixel_step;
            dst_px++;
            info_ptr++;
            grain_buffer_ptr++;
            pixel_proc_next_pixel<mode>(context);
        }
        pixel_proc_next_row<mode>(context);
    }

    pixel_proc_destroy_context<mode>(context);
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

#define DECLARE_IMPL_C
#include "impl_dispatch_decl.h"
