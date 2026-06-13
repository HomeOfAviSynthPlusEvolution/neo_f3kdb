#include <stdlib.h>
#include <mutex>

#include "impl_dispatch.h"
#include "sse_utils.h"
#include "dither_high.h"
#include "VCL2/vectorclass.h"
#include "VCL2/vectormath_exp.h"

/****************************************************************************
 * NOTE: DON'T remove static from any function in this file, it is required *
 *       for generating code in multiple SSE versions.                      *
 ****************************************************************************/

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define process_34 (sample_mode == 2 || sample_mode == 4)

typedef struct _info_cache
{
    int pitch;
    char* data_stream;
} info_cache;

static void destroy_cache(void* data)
{
    assert(data);

    info_cache* cache = (info_cache*) data;
    _aligned_free(cache->data_stream);
    free(data);
}

template <int sample_mode, int ref_part_index>
static __forceinline void process_plane_info_block(
    pixel_dither_info *&info_ptr,
    const unsigned char* src_addr_start,
    const __m128i &src_pitch_vector,
    const __m128i &minus_one,
    const __m128i &width_subsample_vector,
    const __m128i &height_subsample_vector,
    const __m128i &pixel_step_shift_bits,
    char*& info_data_stream)
{
    assert(ref_part_index <= 2);

    __m128i info_block = _mm_load_si128((__m128i*)info_ptr);

    // ref1: bit 0-7
    // left-shift & right-shift 24bits to remove other elements and preserve sign
    __m128i ref1 = info_block;
    ref1 = _mm_slli_epi32(ref1, 24); // << 24
    ref1 = _mm_srai_epi32(ref1, 24); // >> 24

    __m128i ref_offset1;
    __m128i ref_offset2 = _mm_setzero_si128();

    __m128i temp_ref1;
    switch (sample_mode)
    {
    case 0:
        // ref1 = (abs(ref1) >> height_subsampling) * (sign(ref1))
        temp_ref1 = _mm_abs_epi32(ref1);
        temp_ref1 = _mm_sra_epi32(temp_ref1, height_subsample_vector);
        temp_ref1 = _mm_mullo_epi32(temp_ref1, _mm_sign_epi32(ref1, _mm_set1_epi32(1)));
        ref_offset1 = _mm_mullo_epi32(src_pitch_vector, temp_ref1); // packed DWORD multiplication
        break;
    case 1:
        // ref1 is guarenteed to be postive
        temp_ref1 = _mm_sra_epi32(ref1, height_subsample_vector);
        ref_offset1 = _mm_mullo_epi32(src_pitch_vector, temp_ref1); // packed DWORD multiplication
        break;
    case 3:
        temp_ref1 = _mm_sra_epi32(ref1, width_subsample_vector);
        ref_offset1 = _mm_sll_epi32(temp_ref1, pixel_step_shift_bits);
        break;
    case 4:
    case 5:
    case 6:
    case 7:
        temp_ref1 = _mm_sra_epi32(ref1, height_subsample_vector);
        ref_offset1 = _mm_mullo_epi32(src_pitch_vector, temp_ref1);
        temp_ref1 = _mm_sra_epi32(ref1, width_subsample_vector);
        ref_offset2 = _mm_sll_epi32(temp_ref1, pixel_step_shift_bits);
        break;
    case 2:
        // ref2: bit 8-15
        // similar to above
        __m128i ref2;
        ref2 = info_block;
        ref2 = _mm_slli_epi32(ref2, 16); // << 16
        ref2 = _mm_srai_epi32(ref2, 24); // >> 24

        __m128i ref1_fix, ref2_fix;
        // ref_px = src_pitch * info.ref2 + info.ref1;
        ref1_fix = _mm_sra_epi32(ref1, width_subsample_vector);
        ref2_fix = _mm_sra_epi32(ref2, height_subsample_vector);
        ref_offset1 = _mm_mullo_epi32(src_pitch_vector, ref2_fix); // packed DWORD multiplication
        ref_offset1 = _mm_add_epi32(ref_offset1, _mm_sll_epi32(ref1_fix, pixel_step_shift_bits));

        // ref_px_2 = info.ref2 - src_pitch * info.ref1;
        ref1_fix = _mm_sra_epi32(ref1, height_subsample_vector);
        ref2_fix = _mm_sra_epi32(ref2, width_subsample_vector);
        ref_offset2 = _mm_mullo_epi32(src_pitch_vector, ref1_fix); // packed DWORD multiplication
        ref_offset2 = _mm_sub_epi32(_mm_sll_epi32(ref2_fix, pixel_step_shift_bits), ref_offset2);
        break;
    default:
        abort();
    }

    if (info_data_stream){
        _mm_store_si128((__m128i*)info_data_stream, ref_offset1);
        info_data_stream += 16;

        if (sample_mode == 2 || (sample_mode >= 4 && sample_mode <= 7)) {
            _mm_store_si128((__m128i*)info_data_stream, ref_offset2);
            info_data_stream += 16;
        }
    }

    info_ptr += 4;
}

static __forceinline __m128i generate_blend_mask_high(__m128i a, __m128i b, __m128i threshold)
{
    __m128i diff1 = _mm_subs_epu16(a, b);
    __m128i diff2 = _mm_subs_epu16(b, a);

    __m128i abs_diff = _mm_or_si128(diff1, diff2);

    __m128i sign_convert_vector = _mm_set1_epi16( (short)0x8000 );

    __m128i converted_diff = _mm_sub_epi16(abs_diff, sign_convert_vector);

    __m128i converted_threshold = _mm_sub_epi16(threshold, sign_convert_vector);

    // mask: if threshold >= diff, set to 0xff, otherwise 0x00
    // note that this is the opposite of low bitdepth implementation
    return _mm_cmpgt_epi16(converted_threshold, converted_diff);
}

static __forceinline __m128i generate_blend_mask_high(__m128i a, __m128i threshold)
{
    __m128i sign_convert_vector = _mm_set1_epi16((short)0x8000);

    __m128i converted_diff = _mm_sub_epi16(a, sign_convert_vector);

    __m128i converted_threshold = _mm_sub_epi16(threshold, sign_convert_vector);

    // mask: if threshold >= diff, set to 0xff, otherwise 0x00
    // note that this is the opposite of low bitdepth implementation
    return _mm_cmpgt_epi16(converted_threshold, converted_diff);
}

static __forceinline void convert_u16_to_float_x2(const __m128i val_i, __m128& val_f_lo, __m128& val_f_hi)
{
    __m128i zero = _mm_setzero_si128();

    val_f_lo = _mm_cvtepi32_ps(_mm_unpacklo_epi16(val_i, zero));
    val_f_hi = _mm_cvtepi32_ps(_mm_unpackhi_epi16(val_i, zero));
}

static __forceinline __m128i convert_float_x2_to_u16(const __m128 val_f_lo, const __m128 val_f_hi)
{
    const __m128 half = _mm_set1_ps(0.5f);
    __m128i val_i_lo_32 = _mm_cvttps_epi32(_mm_add_ps(val_f_lo, half));
    __m128i val_i_hi_32 = _mm_cvttps_epi32(_mm_add_ps(val_f_hi, half));

    return _mm_packus_epi32(val_i_lo_32, val_i_hi_32);
}

static __forceinline __m128 _mm_pow_ps_scalar_approx(__m128 base, float exponent)
{
    return pow(Vec4f(base), exponent);
}

static __forceinline __m128 abs_ps(__m128 x)
{
    return _mm_and_ps(x, _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff)));
}

static __forceinline __m128 saturate_ps(__m128 val_ps)
{
    return _mm_min_ps(_mm_max_ps(val_ps, _mm_setzero_ps()), _mm_set1_ps(1.0f));
};

static __forceinline __m128 calculate_ratio_term_ps(__m128 diff_ps, __m128 thresh_ps)
{
    __m128 ratio = _mm_div_ps(diff_ps, _mm_max_ps(thresh_ps, _mm_set1_ps(1e-5f)));

    return _mm_sub_ps(_mm_set1_ps(1.0f), ratio);
};

static __forceinline float sse_scalar_get_pixel_value_f(
    const unsigned char* src_plane_base_ptr,
    int x, int y,
    const process_plane_params& params,
    __m128i upsample_shift_simd)
{
    x = std::clamp(x, 0, params.plane_width_in_pixels - 1);
    y = std::clamp(y, 0, params.plane_height_in_pixels - 1);

    const int pixel_step_bytes = (params.input_mode == HIGH_BIT_DEPTH_INTERLEAVED) ? 2 : 1;
    const unsigned char* pixel_address = src_plane_base_ptr +
        static_cast<intptr_t>(y) * params.src_pitch +
        static_cast<intptr_t>(x) * pixel_step_bytes;

    unsigned short raw_pixel_val = (params.input_mode == LOW_BIT_DEPTH) ? *pixel_address
        : *reinterpret_cast<const unsigned short*>(pixel_address);

    return static_cast<float>(raw_pixel_val << (_mm_cvtsi128_si32(upsample_shift_simd)));
}

static __forceinline float sse_scalar_calc_gradient_angle(
    const unsigned char* src_plane_base_ptr,
    int current_x, int current_y,
    int read_distance,
    const process_plane_params& params,
    __m128i upsample_shift_simd)
{
    auto get_pixel_value_at_sse = [&](int x_coord, int y_coord) -> float {
        return sse_scalar_get_pixel_value_f(src_plane_base_ptr, x_coord, y_coord, params, upsample_shift_simd);
        };

    const float p00 = get_pixel_value_at_sse(current_x - read_distance, current_y - read_distance);
    const float p10 = get_pixel_value_at_sse(current_x, current_y - read_distance);
    const float p20 = get_pixel_value_at_sse(current_x + read_distance, current_y - read_distance);
    const float p01 = get_pixel_value_at_sse(current_x - read_distance, current_y);
    const float p21 = get_pixel_value_at_sse(current_x + read_distance, current_y);
    const float p02 = get_pixel_value_at_sse(current_x - read_distance, current_y + read_distance);
    const float p12 = get_pixel_value_at_sse(current_x, current_y + read_distance);
    const float p22 = get_pixel_value_at_sse(current_x + read_distance, current_y + read_distance);

    const float gx = (p20 + 2.0f * p21 + p22) - (p00 + 2.0f * p01 + p02);
    const float gy = (p00 + 2.0f * p10 + p20) - (p02 + 2.0f * p12 + p22);

    const float scaled_epsilon_for_gx = 0.01f * (static_cast<float>(1 << (16 - params.input_depth)) * 3.0f);

    if (std::abs(gx) < scaled_epsilon_for_gx) {
        if (std::abs(gy) < scaled_epsilon_for_gx) {
            return 1.0f; // Flat area
        }
        return 1.0f; // Predominantly vertical gradient
    }

    float angle_rad = std::atan(gy / gx);
    return angle_rad / static_cast<float>(M_PI) + 0.5f;
}

template<int sample_mode, bool blur_first>
static __m128i __forceinline process_pixels_mode12_high_part(__m128i src_pixels, __m128i threshold_vector, __m128i threshold1_vector, __m128i threshold2_vector,
    __m128i change, const __m128i& ref_pixels_1, const __m128i& ref_pixels_2, const __m128i& ref_pixels_3, const __m128i& ref_pixels_4,
    const pixel_dither_info* pdi_ptr, const process_plane_params& params, __m128i upsample_to_16_shift_bits, int row, int column)
{
    __m128i use_orig_pixel_blend_mask_12, use_orig_pixel_blend_mask_34;
    __m128i avg_12, avg_34;
    __m128i dst_pixels;

    if (sample_mode < 5)
    {
        if (true)
        {
            if (!blur_first)
            {
                // note: use AND instead of OR, because two operands are reversed
                // (different from low bit-depth mode!)
                use_orig_pixel_blend_mask_12 = _mm_and_si128(
                    generate_blend_mask_high(src_pixels, ref_pixels_1, threshold_vector),
                    generate_blend_mask_high(src_pixels, ref_pixels_2, threshold_vector));
            }

            avg_12 = _mm_avg_epu16(ref_pixels_1, ref_pixels_2);
        }

        if (process_34)
        {
            if (!blur_first)
            {
                use_orig_pixel_blend_mask_34 = _mm_and_si128(
                    generate_blend_mask_high(src_pixels, ref_pixels_3, threshold_vector),
                    generate_blend_mask_high(src_pixels, ref_pixels_4, threshold_vector));
            }

            avg_34 = _mm_avg_epu16(ref_pixels_3, ref_pixels_4);
        }

        if (sample_mode == 2)
        {
            // merge 34 into 12 for mode 2
            if (!blur_first)
                use_orig_pixel_blend_mask_12 = _mm_and_si128(use_orig_pixel_blend_mask_12, use_orig_pixel_blend_mask_34);
            avg_12 = _mm_subs_epu16(avg_12, _mm_set1_epi16(1));
            avg_12 = _mm_avg_epu16(avg_12, avg_34);
        }

        if (blur_first)
            use_orig_pixel_blend_mask_12 = generate_blend_mask_high(src_pixels, avg_12, threshold_vector);

        if (process_34 && blur_first)
            use_orig_pixel_blend_mask_34 = generate_blend_mask_high(src_pixels, avg_34, threshold_vector);

        if (sample_mode == 4)
        {
            // average dst_12 and dst_34 for mode 4
            __m128i dst_pixels_12 = _mm_blendv_epi8(src_pixels, avg_12, use_orig_pixel_blend_mask_12);
            __m128i dst_pixels_34 = _mm_blendv_epi8(src_pixels, avg_34, use_orig_pixel_blend_mask_34);
            dst_pixels = _mm_avg_epu16(dst_pixels_12, dst_pixels_34);
        }
        else
        {
            // if mask is 0xff (NOT over threshold), select second operand, otherwise select first
            // note this is different from low bitdepth code
            dst_pixels = _mm_blendv_epi8(src_pixels, avg_12, use_orig_pixel_blend_mask_12);
        }
    }
    else if (sample_mode == 5)
    {
        __m128i zero = _mm_setzero_si128();

        __m128i r1_lo = _mm_unpacklo_epi16(ref_pixels_1, zero);
        __m128i r1_hi = _mm_unpackhi_epi16(ref_pixels_1, zero);
        __m128i r2_lo = _mm_unpacklo_epi16(ref_pixels_2, zero);
        __m128i r2_hi = _mm_unpackhi_epi16(ref_pixels_2, zero);
        __m128i r3_lo = _mm_unpacklo_epi16(ref_pixels_3, zero);
        __m128i r3_hi = _mm_unpackhi_epi16(ref_pixels_3, zero);
        __m128i r4_lo = _mm_unpacklo_epi16(ref_pixels_4, zero);
        __m128i r4_hi = _mm_unpackhi_epi16(ref_pixels_4, zero);

        __m128i src_lo_32 = _mm_unpacklo_epi16(src_pixels, zero);
        __m128i src_hi_32 = _mm_unpackhi_epi16(src_pixels, zero);

        __m128i sum_avg_lo_32 = _mm_add_epi32(r1_lo, r2_lo);
        sum_avg_lo_32 = _mm_add_epi32(sum_avg_lo_32, r3_lo);
        sum_avg_lo_32 = _mm_add_epi32(sum_avg_lo_32, r4_lo);

        __m128i sum_avg_hi_32 = _mm_add_epi32(r1_hi, r2_hi);
        sum_avg_hi_32 = _mm_add_epi32(sum_avg_hi_32, r3_hi);
        sum_avg_hi_32 = _mm_add_epi32(sum_avg_hi_32, r4_hi);

        __m128i avg_lo_32 = _mm_srli_epi32(sum_avg_lo_32, 2);
        __m128i avg_hi_32 = _mm_srli_epi32(sum_avg_hi_32, 2);

        __m128i avg = _mm_packus_epi32(avg_lo_32, avg_hi_32);

        __m128i avgDif = _mm_or_si128(_mm_subs_epu16(avg, src_pixels), _mm_subs_epu16(src_pixels, avg));

        __m128i maxDif_p1 = _mm_or_si128(_mm_subs_epu16(ref_pixels_1, src_pixels), _mm_subs_epu16(src_pixels, ref_pixels_1));
        __m128i maxDif_p2 = _mm_or_si128(_mm_subs_epu16(ref_pixels_2, src_pixels), _mm_subs_epu16(src_pixels, ref_pixels_2));

        __m128i maxDif_p3 = _mm_or_si128(_mm_subs_epu16(ref_pixels_3, src_pixels), _mm_subs_epu16(src_pixels, ref_pixels_3));

        __m128i maxDif_p4 = _mm_or_si128(_mm_subs_epu16(ref_pixels_4, src_pixels), _mm_subs_epu16(src_pixels, ref_pixels_4));

        __m128i maxDif = _mm_max_epu16(_mm_max_epu16(maxDif_p1, maxDif_p2), _mm_max_epu16(maxDif_p3, maxDif_p4));

        __m128i sum_r12_lo_32 = _mm_add_epi32(r1_lo, r2_lo);
        __m128i sum_r12_hi_32 = _mm_add_epi32(r1_hi, r2_hi);
        __m128i two_src_lo_32 = _mm_slli_epi32(src_lo_32, 1);
        __m128i two_src_hi_32 = _mm_slli_epi32(src_hi_32, 1);

        __m128i midDif1_abs_lo_32 = _mm_abs_epi32(_mm_sub_epi32(sum_r12_lo_32, two_src_lo_32));
        __m128i midDif1_abs_hi_32 = _mm_abs_epi32(_mm_sub_epi32(sum_r12_hi_32, two_src_hi_32));
        __m128i midDif1_16 = _mm_packus_epi32(midDif1_abs_lo_32, midDif1_abs_hi_32);

        __m128i sum_r34_lo_32 = _mm_add_epi32(r3_lo, r4_lo);
        __m128i sum_r34_hi_32 = _mm_add_epi32(r3_hi, r4_hi);
        __m128i midDif2_abs_lo_32 = _mm_abs_epi32(_mm_sub_epi32(sum_r34_lo_32, two_src_lo_32));
        __m128i midDif2_abs_hi_32 = _mm_abs_epi32(_mm_sub_epi32(sum_r34_hi_32, two_src_hi_32));
        __m128i midDif2_16 = _mm_packus_epi32(midDif2_abs_lo_32, midDif2_abs_hi_32);

        use_orig_pixel_blend_mask_12 = generate_blend_mask_high(avgDif, threshold_vector);
        use_orig_pixel_blend_mask_12 = _mm_and_si128(use_orig_pixel_blend_mask_12, generate_blend_mask_high(maxDif, threshold1_vector));
        use_orig_pixel_blend_mask_12 = _mm_and_si128(use_orig_pixel_blend_mask_12, generate_blend_mask_high(midDif1_16, threshold2_vector));
        use_orig_pixel_blend_mask_12 = _mm_and_si128(use_orig_pixel_blend_mask_12, generate_blend_mask_high(midDif2_16, threshold2_vector));

        dst_pixels = _mm_blendv_epi8(src_pixels, avg, use_orig_pixel_blend_mask_12);
    }
    else if (sample_mode == 6 || sample_mode == 7)
    {
        __m128 final_thresh_avg_dif_f_vec_lo;
        __m128 final_thresh_avg_dif_f_vec_hi;
        __m128 final_thresh_max_dif_f_vec_lo;
        __m128 final_thresh_max_dif_f_vec_hi;
        __m128 final_thresh_mid_dif_f_vec_lo;
        __m128 final_thresh_mid_dif_f_vec_hi;

        const __m128 orig_thresh_avg_ps = _mm_set1_ps(static_cast<float>(_mm_extract_epi16(threshold_vector, 0)));
        const __m128 orig_thresh_max_ps = _mm_set1_ps(static_cast<float>(_mm_extract_epi16(threshold1_vector, 0)));
        const __m128 orig_thresh_mid_ps = _mm_set1_ps(static_cast<float>(_mm_extract_epi16(threshold2_vector, 0)));

        if (sample_mode == 7) {
            const float angle_boost_factor_val = params.angle_boost;
            const float max_angle_threshold_val = params.max_angle;
            const int grad_read_distance = 20;

            alignas(16)
                float current_pixel_max_angle_diff_buffer[4];

            for (int four_pix_group = 0; four_pix_group < 2; ++four_pix_group)
            {
                for (int k_in_group = 0; k_in_group < 4; ++k_in_group)
                {
                    const int pixel_idx_in_block = four_pix_group * 4 + k_in_group;
                    const int current_x = column + pixel_idx_in_block;

                    const float angle_org = sse_scalar_calc_gradient_angle(params.src_plane_ptr, current_x, row, grad_read_distance, params,
                        upsample_to_16_shift_bits);

                    const int ref1_val = static_cast<int>(pdi_ptr[pixel_idx_in_block].ref1);

                    const int ref1h_y_offset = ref1_val >> params.height_subsampling;
                    const int ref1w_x_offset = ref1_val >> params.width_subsampling;

                    const float angle_ref1_h = sse_scalar_calc_gradient_angle(params.src_plane_ptr, current_x, row + ref1h_y_offset,
                        grad_read_distance, params, upsample_to_16_shift_bits);
                    const float angle_ref2_h = sse_scalar_calc_gradient_angle(params.src_plane_ptr, current_x, row - ref1h_y_offset,
                        grad_read_distance, params, upsample_to_16_shift_bits);
                    const float angle_ref1_w = sse_scalar_calc_gradient_angle(params.src_plane_ptr, current_x + ref1w_x_offset, row,
                        grad_read_distance, params, upsample_to_16_shift_bits);
                    const float angle_ref2_w = sse_scalar_calc_gradient_angle(params.src_plane_ptr, current_x - ref1w_x_offset, row,
                        grad_read_distance, params, upsample_to_16_shift_bits);

                    float max_diff = 0.0f;
                    max_diff = std::max(max_diff, std::abs(angle_ref1_h - angle_org));
                    max_diff = std::max(max_diff, std::abs(angle_ref2_h - angle_org));
                    max_diff = std::max(max_diff, std::abs(angle_ref1_w - angle_org));
                    max_diff = std::max(max_diff, std::abs(angle_ref2_w - angle_org));
                    current_pixel_max_angle_diff_buffer[k_in_group] = max_diff;
                }

                Vec4f max_angle_diff_ps = Vec4f().load(current_pixel_max_angle_diff_buffer);
                Vec4fb use_boost_ps = (max_angle_diff_ps <= max_angle_threshold_val);
                Vec4f boost_factor_ps = _mm_set1_ps(angle_boost_factor_val);

                Vec4f current_thresh_avg_ps = select(use_boost_ps, static_cast<Vec4f>(orig_thresh_avg_ps) * boost_factor_ps,
                    static_cast<Vec4f>(orig_thresh_avg_ps));
                Vec4f current_thresh_max_ps = select(use_boost_ps, static_cast<Vec4f>(orig_thresh_max_ps) * boost_factor_ps,
                    static_cast<Vec4f>(orig_thresh_max_ps));
                Vec4f current_thresh_mid_ps = select(use_boost_ps, static_cast<Vec4f>(orig_thresh_mid_ps) * boost_factor_ps,
                    static_cast<Vec4f>(orig_thresh_mid_ps));

                if (four_pix_group == 0) {
                    final_thresh_avg_dif_f_vec_lo = current_thresh_avg_ps;
                    final_thresh_max_dif_f_vec_lo = current_thresh_max_ps;
                    final_thresh_mid_dif_f_vec_lo = current_thresh_mid_ps;
                }
                else {
                    final_thresh_avg_dif_f_vec_hi = current_thresh_avg_ps;
                    final_thresh_max_dif_f_vec_hi = current_thresh_max_ps;
                    final_thresh_mid_dif_f_vec_hi = current_thresh_mid_ps;
                }
            }
        }
        else {
            final_thresh_avg_dif_f_vec_lo = final_thresh_avg_dif_f_vec_hi = orig_thresh_avg_ps;
            final_thresh_max_dif_f_vec_lo = final_thresh_max_dif_f_vec_hi = orig_thresh_max_ps;
            final_thresh_mid_dif_f_vec_lo = final_thresh_mid_dif_f_vec_hi = orig_thresh_mid_ps;
        }

        const __m128 f_const_3_0 = _mm_set1_ps(3.0f);

        __m128 src_f_lo;
        __m128 src_f_hi;
        __m128 ref1_f_lo;
        __m128 ref1_f_hi;
        __m128 ref2_f_lo;
        __m128 ref2_f_hi;
        __m128 ref3_f_lo;
        __m128 ref3_f_hi;
        __m128 ref4_f_lo;
        __m128 ref4_f_hi;

        convert_u16_to_float_x2(src_pixels, src_f_lo, src_f_hi);
        convert_u16_to_float_x2(ref_pixels_1, ref1_f_lo, ref1_f_hi);
        convert_u16_to_float_x2(ref_pixels_2, ref2_f_lo, ref2_f_hi);
        convert_u16_to_float_x2(ref_pixels_3, ref3_f_lo, ref3_f_hi);
        convert_u16_to_float_x2(ref_pixels_4, ref4_f_lo, ref4_f_hi);

        __m128 blended_f_lo;
        __m128 blended_f_hi;

        for (int part = 0; part < 2; ++part)
        {
            __m128 current_thresh_avg_ps = (part == 0) ? final_thresh_avg_dif_f_vec_lo : final_thresh_avg_dif_f_vec_hi;
            __m128 current_thresh_max_ps = (part == 0) ? final_thresh_max_dif_f_vec_lo : final_thresh_max_dif_f_vec_hi;
            __m128 current_thresh_mid_ps = (part == 0) ? final_thresh_mid_dif_f_vec_lo : final_thresh_mid_dif_f_vec_hi;

            __m128 src_f = (part == 0) ? src_f_lo : src_f_hi;
            __m128 p1_f = (part == 0) ? ref1_f_lo : ref1_f_hi;
            __m128 p3_f = (part == 0) ? ref2_f_lo : ref2_f_hi;
            __m128 p2_f = (part == 0) ? ref3_f_lo : ref3_f_hi;
            __m128 p4_f = (part == 0) ? ref4_f_lo : ref4_f_hi;

            __m128 sum_refs = _mm_add_ps(_mm_add_ps(p1_f, p2_f), _mm_add_ps(p3_f, p4_f));
            __m128 avg_refs_f = _mm_mul_ps(sum_refs, _mm_set1_ps(0.25f));

            __m128 diff_avg_src = _mm_sub_ps(avg_refs_f, src_f);
            __m128 avg_dif_f = abs_ps(diff_avg_src);

            __m128 d1 = abs_ps(_mm_sub_ps(p1_f, src_f));
            __m128 d2 = abs_ps(_mm_sub_ps(p2_f, src_f));
            __m128 d3 = abs_ps(_mm_sub_ps(p3_f, src_f));
            __m128 d4 = abs_ps(_mm_sub_ps(p4_f, src_f));
            __m128 maxDif = _mm_max_ps(_mm_max_ps(d1, d3), _mm_max_ps(d2, d4));

            __m128 two_src = _mm_mul_ps(src_f, _mm_set1_ps(2.0f));
            __m128 mid_dif_v_f = abs_ps(_mm_sub_ps(_mm_add_ps(p1_f, p3_f), two_src));

            __m128 mid_dif_h_f = abs_ps(_mm_sub_ps(_mm_add_ps(p2_f, p4_f), two_src));

            __m128 comp_avg = saturate_ps(_mm_mul_ps(f_const_3_0, calculate_ratio_term_ps(avg_dif_f, current_thresh_avg_ps)));
            __m128 comp_max = saturate_ps(_mm_mul_ps(f_const_3_0, calculate_ratio_term_ps(maxDif, current_thresh_max_ps)));
            __m128 comp_mid_v = saturate_ps(_mm_mul_ps(f_const_3_0, calculate_ratio_term_ps(mid_dif_v_f, current_thresh_mid_ps)));
            __m128 comp_mid_h = saturate_ps(_mm_mul_ps(f_const_3_0, calculate_ratio_term_ps(mid_dif_h_f, current_thresh_mid_ps)));

            __m128 product_comps = _mm_mul_ps(_mm_mul_ps(comp_avg, comp_max), _mm_mul_ps(comp_mid_v, comp_mid_h));

            __m128 factor = _mm_pow_ps_scalar_approx(product_comps, 0.1f);

            __m128 blended_f = _mm_add_ps(src_f, _mm_mul_ps(diff_avg_src, factor));

            if (part == 0)
                blended_f_lo = blended_f;
            else
                blended_f_hi = blended_f;
        }

        dst_pixels = convert_float_x2_to_u16(blended_f_lo, blended_f_hi);
    }

    __m128i sign_convert_vector = _mm_set1_epi16((short)0x8000);

    // convert to signed form, since change is signed
    dst_pixels = _mm_sub_epi16(dst_pixels, sign_convert_vector);

    // saturated add
    dst_pixels = _mm_adds_epi16(dst_pixels, change);

    // convert back to unsigned
    dst_pixels = _mm_add_epi16(dst_pixels, sign_convert_vector);
    return dst_pixels;
}

template<int sample_mode, bool blur_first, int dither_algo>
static __m128i __forceinline process_pixels(
    __m128i src_pixels_0,
    __m128i threshold_vector,
    __m128i threshold1_vector,
    __m128i threshold2_vector,
    const __m128i& change_1,
    const __m128i& ref_pixels_1_0,
    const __m128i& ref_pixels_2_0,
    const __m128i& ref_pixels_3_0,
    const __m128i& ref_pixels_4_0,
    const __m128i& clamp_high_add,
    const __m128i& clamp_high_sub,
    const __m128i& clamp_low,
    bool need_clamping,
    int row,
    int column,
    void* dither_context,
    const pixel_dither_info* pdi_ptr,
    const process_plane_params& params,
    __m128i upsample_to_16_shift_bits)
{
    __m128i ret = process_pixels_mode12_high_part<sample_mode, blur_first>
        (src_pixels_0,
         threshold_vector,
         threshold1_vector,
         threshold2_vector,
         change_1,
         ref_pixels_1_0,
         ref_pixels_2_0,
         ref_pixels_3_0,
         ref_pixels_4_0,
         pdi_ptr, params,
         upsample_to_16_shift_bits,
         row,
         column);

    switch (dither_algo)
    {
    case DA_HIGH_NO_DITHERING:
    case DA_HIGH_ORDERED_DITHERING:
    case DA_HIGH_FLOYD_STEINBERG_DITHERING:
        ret = dither_high::dither<dither_algo>(dither_context, ret, row, column);
        break;
    default:
        break;
    }
    if (need_clamping)
    {
        ret = high_bit_depth_pixels_clamp(ret, clamp_high_add, clamp_high_sub, clamp_low);
    }


    return ret;
}

template <PIXEL_MODE output_mode>
static int __forceinline store_pixels(
    __m128i pixels,
    __m128i downshift_bits,
    unsigned char* dst,
    int dst_pitch,
    int height_in_pixels)
{
    switch (output_mode)
    {
    case LOW_BIT_DEPTH:
        {
            pixels = _mm_srli_epi16(pixels, 8);
            pixels = _mm_packus_epi16(pixels, pixels);
            _mm_storel_epi64((__m128i*)dst, pixels);
            return 8;
            break;
        }
    case HIGH_BIT_DEPTH_INTERLEAVED:
        pixels = _mm_srl_epi16(pixels, downshift_bits);
        _mm_store_si128((__m128i*)dst, pixels);
        return 16;
        break;
    default:
        abort();
    }
    return 0;
}


template<bool aligned>
static __m128i load_m128(const unsigned char *ptr)
{
    if (aligned)
    {
        return _mm_load_si128((const __m128i*)ptr);
    } else {
        return _mm_loadu_si128((const __m128i*)ptr);
    }
}

template<PIXEL_MODE input_mode, bool aligned>
static __m128i __forceinline read_pixels(
    const process_plane_params& params,
    const unsigned char *ptr,
    __m128i upsample_shift)
{
    __m128i ret;

    switch (input_mode)
    {
    case LOW_BIT_DEPTH:
        ret = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)ptr));
        break;
    case HIGH_BIT_DEPTH_INTERLEAVED:
        ret = load_m128<aligned>(ptr);
        break;
    default:
        abort();
    }
    ret = _mm_sll_epi16(ret, upsample_shift);
    return ret;
}

template<PIXEL_MODE input_mode>
static unsigned short __forceinline read_pixel(
    const unsigned char* base,
    int offset)
{
    const unsigned char* ptr = base + offset;

    switch (input_mode)
    {
    case LOW_BIT_DEPTH:
        return *ptr;
        break;
    case HIGH_BIT_DEPTH_INTERLEAVED:
        return *(reinterpret_cast<const unsigned short*>(ptr));
        break;
    default:
        // shouldn't happen!
        abort();
        return 0;
    }

}

template <int dither_algo>
static __m128i __forceinline load_reference_pixels(
    __m128i shift,
    const unsigned short src[8])
{
    __m128i ret = _mm_load_si128((const __m128i*)src);
    ret = _mm_sll_epi16(ret, shift);
    return ret;
}


template<int sample_mode, int dither_algo, PIXEL_MODE input_mode>
static void __forceinline read_reference_pixels(
    const process_plane_params& params,
    __m128i shift,
    const unsigned char* src_px_start,
    const char* info_data_start,
    __m128i& ref_pixels_1_0,
    __m128i& ref_pixels_2_0,
    __m128i& ref_pixels_3_0,
    __m128i& ref_pixels_4_0)
{
    alignas(16)
    unsigned short tmp_1[8];
    alignas(16)
    unsigned short tmp_2[8];
    alignas(16)
    unsigned short tmp_3[8];
    alignas(16)
    unsigned short tmp_4[8];

    // cache layout: 8 offset groups (1 or 2 offsets / group depending on sample mode) in a pack,
    //               followed by 16 bytes of change values
    // in the case of 2 offsets / group, offsets are stored like this:
    // [1 1 1 1
    //  2 2 2 2
    //  1 1 1 1
    //  2 2 2 2
    //  .. snip
    //  1 1 1 1
    //  2 2 2 2]

    int i_fix = 0;
    int i_fix_step = (input_mode != HIGH_BIT_DEPTH_INTERLEAVED ? 1 : 2);

    for (int i = 0; i < 8; i++)
    {
        const int* offset_ptr_base = reinterpret_cast<const int*>(info_data_start);

        switch (sample_mode)
        {
        case 0:
            tmp_1[i] = read_pixel<input_mode>(src_px_start, i_fix + offset_ptr_base[i]);
            break;
        case 1:
        case 3:
            tmp_1[i] = read_pixel<input_mode>(src_px_start, i_fix + offset_ptr_base[i]);
            tmp_2[i] = read_pixel<input_mode>(src_px_start, i_fix - offset_ptr_base[i]);
            break;
        case 2:
        case 4:
        case 5:
        case 6:
        case 7:
            tmp_1[i] = read_pixel<input_mode>(src_px_start, i_fix + offset_ptr_base[i + (i / 4) * 4]);
            tmp_2[i] = read_pixel<input_mode>(src_px_start, i_fix - offset_ptr_base[i + (i / 4) * 4]);
            tmp_3[i] = read_pixel<input_mode>(src_px_start, i_fix + offset_ptr_base[i + (i / 4) * 4 + 4]);
            tmp_4[i] = read_pixel<input_mode>(src_px_start, i_fix - offset_ptr_base[i + (i / 4) * 4 + 4]);
            break;
        }
        i_fix += i_fix_step;
    }

    switch (sample_mode)
    {
    case 0:
        ref_pixels_1_0 = load_reference_pixels<dither_algo>(shift, tmp_1);
        break;
    case 1:
    case 3:
        ref_pixels_1_0 = load_reference_pixels<dither_algo>(shift, tmp_1);
        ref_pixels_2_0 = load_reference_pixels<dither_algo>(shift, tmp_2);
        break;
    case 2:
    case 4:
    case 5:
    case 6:
    case 7:
        ref_pixels_1_0 = load_reference_pixels<dither_algo>(shift, tmp_1);
        ref_pixels_2_0 = load_reference_pixels<dither_algo>(shift, tmp_2);
        ref_pixels_3_0 = load_reference_pixels<dither_algo>(shift, tmp_3);
        ref_pixels_4_0 = load_reference_pixels<dither_algo>(shift, tmp_4);
        break;
    }
}

std::mutex cache_mutex;
template<int sample_mode, bool blur_first, int dither_algo, bool aligned, PIXEL_MODE output_mode>
static void __cdecl _process_plane_sse_impl(const process_plane_params& params, process_plane_context* context)
{
    assert(sample_mode > 0);

    __m128i src_pitch_vector = _mm_set1_epi32(params.src_pitch);

    __m128i threshold_vector = _mm_set1_epi16(params.threshold);
    __m128i threshold1_vector = _mm_set1_epi16(params.threshold1);
    __m128i threshold2_vector = _mm_set1_epi16(params.threshold2);

    // general-purpose constant
    __m128i minus_one = _mm_set1_epi32(-1);

    alignas(16)
    char context_buffer[DITHER_CONTEXT_BUFFER_SIZE];

    dither_high::init<dither_algo>(context_buffer, params.plane_width_in_pixels, params.output_depth);


    __m128i width_subsample_vector = _mm_set_epi32(0, 0, 0, params.width_subsampling);
    __m128i height_subsample_vector = _mm_set_epi32(0, 0, 0, params.height_subsampling);

    bool need_clamping =  INTERNAL_BIT_DEPTH < 16 ||
                          params.pixel_min > 0 ||
                          params.pixel_max < 0xffff; // (1 << INTERNAL_BIT_DEPTH) -1)
    __m128i clamp_high_add = _mm_setzero_si128();
    __m128i clamp_high_sub = _mm_setzero_si128();
    __m128i clamp_low = _mm_setzero_si128();
    if (need_clamping)
    {
        clamp_low = _mm_set1_epi16(static_cast<short>(params.pixel_min));
        clamp_high_add = _mm_sub_epi16(_mm_set1_epi16(static_cast<short>(0xFFFF)), _mm_set1_epi16(static_cast<short>(params.pixel_max)));
        clamp_high_sub = _mm_add_epi16(clamp_high_add, clamp_low);
    }

    __m128i pixel_step_shift_bits = (params.input_mode == HIGH_BIT_DEPTH_INTERLEAVED) ? _mm_set_epi32(0, 0, 0, 1) : _mm_setzero_si128();
    __m128i upsample_to_16_shift_bits = _mm_set_epi32(0, 0, 0, INTERNAL_BIT_DEPTH - params.input_depth);
    __m128i downshift_bits = _mm_set_epi32(0, 0, 0, INTERNAL_BIT_DEPTH - params.output_depth);

    bool use_cached_info = false;
    info_cache *cache = NULL;
    char* info_data_stream = NULL;

    alignas(16)
    char dummy_info_buffer[128];

    // initialize storage for pre-calculated pixel offsets
    if (context->data) {
        cache = (info_cache*) context->data;
        // we need to ensure src_pitch is the same, otherwise offsets will be completely wrong
        // also, if pitch changes, don't waste time to update the cache since it is likely to change again
        if (cache->pitch == params.src_pitch) {
            info_data_stream = cache->data_stream;
            use_cached_info = true;
            cache = NULL;
        } else {
            // info_data_stream can be NULL, in this case dummy_info_buffer will be used for temporary storage
        }
        cache = NULL;
    } else {
        // set up buffer for cache
        cache = (info_cache*)malloc(sizeof(info_cache));
        if (cache) {
            size_t cache_size = static_cast<size_t>((params.plane_width_in_pixels + 7) / 8) * params.plane_height_in_pixels *
                ((sample_mode == 2 || (sample_mode >= 4 && sample_mode <= 7)) ? 64 : 32);
            info_data_stream = (char*)_aligned_malloc(cache_size, FRAME_LUT_ALIGNMENT);
            if (info_data_stream) {
                cache->data_stream = info_data_stream;
                cache->pitch = params.src_pitch;
            }
            else {
                free(cache);
                cache = NULL;
            }
        }
    }

    const int info_cache_block_size = (sample_mode == 2 || (sample_mode >= 4 && sample_mode <= 7)) ? 64 : 32;

    int current_input_mode = params.input_mode;

    for (int row = 0; row < params.plane_height_in_pixels; row++)
    {
        const unsigned char* src_px_row_base = params.src_plane_ptr + static_cast<intptr_t>(params.src_pitch) * row;
        unsigned char* dst_px_row_base = params.dst_plane_ptr + static_cast<intptr_t>(params.dst_pitch) * row;

        pixel_dither_info* info_ptr_row_base = params.info_ptr_base + static_cast<intptr_t>(params.info_stride) * row;
        const short* grain_buffer_row_base = params.grain_buffer + static_cast<intptr_t>(params.grain_buffer_stride) * row;

        char* current_row_info_data_cache_ptr = use_cached_info ?
            (info_data_stream + static_cast<intptr_t>((params.plane_width_in_pixels + 7) / 8) * row * info_cache_block_size)
            : nullptr;
        if (use_cached_info && !info_data_stream)
            use_cached_info = false;

        char* current_row_info_data_build_ptr = (!use_cached_info && cache && info_data_stream) ?
            (info_data_stream + static_cast<intptr_t>((params.plane_width_in_pixels + 7) / 8) * row * info_cache_block_size)
            : nullptr;

        int processed_pixels = 0;
        pixel_dither_info* current_info_unit_ptr = info_ptr_row_base;

        while (processed_pixels < params.plane_width_in_pixels)
        {
            const unsigned char* current_src_px = src_px_row_base + processed_pixels *
                (current_input_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1);
            unsigned char* current_dst_px = dst_px_row_base + processed_pixels * (output_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1);
            const short* current_grain_ptr = grain_buffer_row_base + processed_pixels;
            const pixel_dither_info* pdi_for_process_pixels = info_ptr_row_base + processed_pixels;

            __m128i change_1;

            __m128i ref_pixels_1_0 = _mm_setzero_si128();
            __m128i ref_pixels_2_0 = _mm_setzero_si128();
            __m128i ref_pixels_3_0 = _mm_setzero_si128();
            __m128i ref_pixels_4_0 = _mm_setzero_si128();

            char * data_stream_for_read_refs;

            if (use_cached_info) {
                data_stream_for_read_refs = current_row_info_data_cache_ptr;
                current_row_info_data_cache_ptr += info_cache_block_size;
            }
            else {
                char* temp_info_build_ptr = current_row_info_data_build_ptr ? current_row_info_data_build_ptr : dummy_info_buffer;
                data_stream_for_read_refs = temp_info_build_ptr;

                // Process info for first 4 pixels of the 8-pixel block
                process_plane_info_block<sample_mode, 0>(current_info_unit_ptr, current_src_px, src_pitch_vector, minus_one,
                    width_subsample_vector, height_subsample_vector, pixel_step_shift_bits, temp_info_build_ptr);
                // Process info for second 4 pixels of the 8-pixel block
                process_plane_info_block<sample_mode, 1>(current_info_unit_ptr, current_src_px + 4 *
                    (current_input_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1), src_pitch_vector, minus_one, width_subsample_vector,
                    height_subsample_vector, pixel_step_shift_bits, temp_info_build_ptr);
                // current_info_unit_ptr is advanced by 8 total by these two calls.
                // temp_info_build_ptr is advanced by info_cache_block_size total.
                if (current_row_info_data_build_ptr) {
                    current_row_info_data_build_ptr += info_cache_block_size;
                }
            }

            if (current_input_mode == LOW_BIT_DEPTH) {
                read_reference_pixels<sample_mode, dither_algo, LOW_BIT_DEPTH>(
                    params, upsample_to_16_shift_bits, current_src_px, data_stream_for_read_refs,
                    ref_pixels_1_0, ref_pixels_2_0, ref_pixels_3_0, ref_pixels_4_0);
            }
            else {
                read_reference_pixels<sample_mode, dither_algo, HIGH_BIT_DEPTH_INTERLEAVED>(
                    params, upsample_to_16_shift_bits, current_src_px, data_stream_for_read_refs,
                    ref_pixels_1_0, ref_pixels_2_0, ref_pixels_3_0, ref_pixels_4_0);
            }

            __m128i src_pixels_data;
            if (current_input_mode == LOW_BIT_DEPTH) {
                src_pixels_data = read_pixels<LOW_BIT_DEPTH, aligned>(params, current_src_px, upsample_to_16_shift_bits);
            }
            else {
                src_pixels_data = read_pixels<HIGH_BIT_DEPTH_INTERLEAVED, aligned>(params, current_src_px, upsample_to_16_shift_bits);
            }

            change_1 = _mm_load_si128((__m128i*)current_grain_ptr);

            __m128i dst_pixels_data = process_pixels<sample_mode, blur_first, dither_algo>(
                                     src_pixels_data,
                                     threshold_vector,
                                     threshold1_vector,
                                     threshold2_vector,
                                     change_1,
                                     ref_pixels_1_0,
                                     ref_pixels_2_0,
                                     ref_pixels_3_0,
                                     ref_pixels_4_0,
                                     clamp_high_add,
                                     clamp_high_sub,
                                     clamp_low,
                                     need_clamping,
                                     row,
                                     processed_pixels,
                                     context_buffer,
                                     pdi_for_process_pixels,
                                     params,
                                     upsample_to_16_shift_bits);

            store_pixels<output_mode>(dst_pixels_data, downshift_bits, current_dst_px, params.dst_pitch, params.plane_height_in_pixels);

            processed_pixels += 8;
        }
        dither_high::next_row<dither_algo>(context_buffer);
    }

    dither_high::complete<dither_algo>(context_buffer);

    // for thread-safety, save context after all data is processed
    if (!use_cached_info && !context->data && cache && info_data_stream)
    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        context->destroy = destroy_cache;
        if (context->data)
            destroy_cache(cache);
        else {
            context->data = cache;
            context->destroy = destroy_cache;
        }
    }
    else if (cache && (!info_data_stream || context->data)) {
        if (info_data_stream)
            _aligned_free(info_data_stream);

        free(cache);
    }
}


template<int sample_mode, bool blur_first, int dither_algo>
static void process_plane_sse_impl(const process_plane_params& params, process_plane_context* context)
{
    switch (params.output_mode)
    {
    case LOW_BIT_DEPTH:
        _process_plane_sse_impl<sample_mode, blur_first, dither_algo, true, LOW_BIT_DEPTH>(params, context);
        break;
    case HIGH_BIT_DEPTH_INTERLEAVED:
        _process_plane_sse_impl<sample_mode, blur_first, dither_algo, true, HIGH_BIT_DEPTH_INTERLEAVED>(params, context);
        break;
    default:
        abort();
    }
}
