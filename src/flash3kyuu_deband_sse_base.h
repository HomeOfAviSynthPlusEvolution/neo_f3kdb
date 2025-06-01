#include <stdlib.h>
#include <mutex>

#include "impl_dispatch.h"
#include "sse_utils.h"
#include "dither_high.h"

/****************************************************************************
 * NOTE: DON'T remove static from any function in this file, it is required *
 *       for generating code in multiple SSE versions.                      *
 ****************************************************************************/

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
    __m128i ref_offset2;

    __m128i temp_ref1;
    switch (sample_mode)
    {
    case 0:
        // ref1 = (abs(ref1) >> height_subsampling) * (sign(ref1))
        temp_ref1 = _mm_abs_epi32(ref1);
        temp_ref1 = _mm_sra_epi32(temp_ref1, height_subsample_vector);
        temp_ref1 = _mm_mullo_epi32(temp_ref1, _mm_srai_epi32(ref1, 31));
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
    case 5:
    case 6:
    case 4:
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

        if (sample_mode == 2 || sample_mode == 4 || sample_mode == 5 || sample_mode == 6) {
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
    alignas(16) float b[4];
    alignas(16) float r[4];
    _mm_store_ps(b, base);

    for (int i = 0; i < 4; ++i)
        r[i] = std::pow(b[i], exponent);

    return _mm_load_ps(r);
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

template<int sample_mode, bool blur_first>
static __m128i __forceinline process_pixels_mode12_high_part(__m128i src_pixels, __m128i threshold_vector, __m128i threshold1_vector, __m128i threshold2_vector,
    __m128i change, const __m128i& ref_pixels_1, const __m128i& ref_pixels_2, const __m128i& ref_pixels_3, const __m128i& ref_pixels_4)
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
        __m128i avg12 = _mm_avg_epu16(ref_pixels_1, ref_pixels_2);
        __m128i avg34 = _mm_avg_epu16(ref_pixels_3, ref_pixels_4);
        __m128i avg = _mm_avg_epu16(avg12, avg34);
        __m128i avgDif = _mm_or_si128(_mm_subs_epu16(avg, src_pixels), _mm_subs_epu16(src_pixels, avg));
        
        __m128i maxDif = _mm_max_epu16(
            _mm_or_si128(_mm_subs_epu16(ref_pixels_1, src_pixels), _mm_subs_epu16(src_pixels, ref_pixels_1)),
            _mm_or_si128(_mm_subs_epu16(ref_pixels_2, src_pixels), _mm_subs_epu16(src_pixels, ref_pixels_2)));
        maxDif = _mm_max_epu16(
            _mm_or_si128(_mm_subs_epu16(ref_pixels_3, src_pixels), _mm_subs_epu16(src_pixels, ref_pixels_3)),
            maxDif);
        maxDif = _mm_max_epu16(
            _mm_or_si128(_mm_subs_epu16(ref_pixels_4, src_pixels), _mm_subs_epu16(src_pixels, ref_pixels_4)),
            maxDif);

        __m128i midDif1 = _mm_adds_epu16(ref_pixels_1, ref_pixels_2);
        midDif1 = _mm_or_si128(
            _mm_subs_epu16(midDif1, _mm_sll_epi16(src_pixels, _mm_cvtsi32_si128(1))),
            _mm_subs_epu16(_mm_sll_epi16(src_pixels, _mm_cvtsi32_si128(1)), midDif1));

        __m128i midDif2 = _mm_adds_epu16(ref_pixels_3, ref_pixels_4);
        midDif2 = _mm_or_si128(
            _mm_subs_epu16(midDif2, _mm_sll_epi16(src_pixels, _mm_cvtsi32_si128(1))),
            _mm_subs_epu16(_mm_sll_epi16(src_pixels, _mm_cvtsi32_si128(1)), midDif2));

        use_orig_pixel_blend_mask_12 = _mm_and_si128(
            generate_blend_mask_high(avgDif, threshold_vector),
            generate_blend_mask_high(maxDif, threshold1_vector));
        use_orig_pixel_blend_mask_12 = _mm_and_si128(
            use_orig_pixel_blend_mask_12,
            generate_blend_mask_high(midDif1, threshold2_vector));
        use_orig_pixel_blend_mask_12 = _mm_and_si128(
            use_orig_pixel_blend_mask_12,
            generate_blend_mask_high(midDif2, threshold2_vector));

        dst_pixels = _mm_blendv_epi8(src_pixels, avg, use_orig_pixel_blend_mask_12);       
    }
    else if (sample_mode == 6)
    {
        const __m128 f_const_3_0 = _mm_set1_ps(3.0f);

        __m128 src_f_lo, src_f_hi;
        __m128 ref1_f_lo, ref1_f_hi;
        __m128 ref2_f_lo, ref2_f_hi;
        __m128 ref3_f_lo, ref3_f_hi;
        __m128 ref4_f_lo, ref4_f_hi;

        convert_u16_to_float_x2(src_pixels, src_f_lo, src_f_hi);
        convert_u16_to_float_x2(ref_pixels_1, ref1_f_lo, ref1_f_hi);
        convert_u16_to_float_x2(ref_pixels_2, ref2_f_lo, ref2_f_hi);
        convert_u16_to_float_x2(ref_pixels_3, ref3_f_lo, ref3_f_hi);
        convert_u16_to_float_x2(ref_pixels_4, ref4_f_lo, ref4_f_hi);

        __m128 thresh_avg_dif_f_vec = _mm_set1_ps(static_cast<float>(_mm_extract_epi16(threshold_vector, 0)));
        __m128 thresh_max_dif_f_vec = _mm_set1_ps(static_cast<float>(_mm_extract_epi16(threshold1_vector, 0)));
        __m128 thresh_mid_dif_f_vec = _mm_set1_ps(static_cast<float>(_mm_extract_epi16(threshold2_vector, 0)));

        __m128 blended_f_lo;
        __m128 blended_f_hi;

        for (int part = 0; part < 2; ++part)
        {
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
            __m128 max_dif_f = _mm_max_ps(_mm_max_ps(d1, d2), _mm_max_ps(d3, d4));

            __m128 two_src = _mm_mul_ps(src_f, _mm_set1_ps(2.0f));
            __m128 mid_dif_v_f = abs_ps(_mm_sub_ps(_mm_add_ps(p1_f, p3_f), two_src));

            __m128 mid_dif_h_f = abs_ps(_mm_sub_ps(_mm_add_ps(p2_f, p4_f), two_src));

            __m128 comp_avg = saturate_ps(_mm_mul_ps(f_const_3_0, calculate_ratio_term_ps(avg_dif_f, thresh_avg_dif_f_vec)));
            __m128 comp_max = saturate_ps(_mm_mul_ps(f_const_3_0, calculate_ratio_term_ps(max_dif_f, thresh_max_dif_f_vec)));
            __m128 comp_mid_v = saturate_ps(_mm_mul_ps(f_const_3_0, calculate_ratio_term_ps(mid_dif_v_f, thresh_mid_dif_f_vec)));
            __m128 comp_mid_h = saturate_ps(_mm_mul_ps(f_const_3_0, calculate_ratio_term_ps(mid_dif_h_f, thresh_mid_dif_f_vec)));

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
    void* dither_context)
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
         ref_pixels_4_0);
    
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
        {
            __m128i zero = _mm_setzero_si128();
            ret = _mm_loadl_epi64((__m128i*)ptr);
            ret = _mm_unpacklo_epi8(zero, ret);
            return ret;
        }
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
        return *(unsigned short*)ptr;
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
        switch (sample_mode)
        {
        case 0:
            tmp_1[i] = read_pixel<input_mode>(src_px_start, i_fix + *(int*)(info_data_start + 4 * i));
            break;
        case 1:
        case 3:
            tmp_1[i] = read_pixel<input_mode>(src_px_start, i_fix + *(int*)(info_data_start + 4 * i));
            tmp_2[i] = read_pixel<input_mode>(src_px_start, i_fix + -*(int*)(info_data_start + 4 * i));
            break;
        case 2:
        case 5:
        case 6:
        case 4:
            tmp_1[i] = read_pixel<input_mode>(src_px_start, i_fix + *(int*)(info_data_start + 4 * (i + i / 4 * 4)));
            tmp_2[i] = read_pixel<input_mode>(src_px_start, i_fix + -*(int*)(info_data_start + 4 * (i + i / 4 * 4)));
            tmp_3[i] = read_pixel<input_mode>(src_px_start, i_fix + *(int*)(info_data_start + 4 * (i + i / 4 * 4 + 4)));
            tmp_4[i] = read_pixel<input_mode>(src_px_start, i_fix + -*(int*)(info_data_start + 4 * (i + i / 4 * 4 + 4)));
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
    case 5:
    case 6:
    case 4:
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

    pixel_dither_info* info_ptr = params.info_ptr_base;

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
                          params.pixel_max < 0xffff;
    __m128i clamp_high_add = _mm_setzero_si128();
    __m128i clamp_high_sub = _mm_setzero_si128();
    __m128i clamp_low = _mm_setzero_si128();
    if (need_clamping)
    {
        clamp_low = _mm_set1_epi16((short)params.pixel_min);
        clamp_high_add = _mm_sub_epi16(_mm_set1_epi16((short)0xffff), _mm_set1_epi16((short)params.pixel_max));
        clamp_high_sub = _mm_add_epi16(clamp_high_add, clamp_low);
    }
    
    __m128i pixel_step_shift_bits;
    __m128i upsample_to_16_shift_bits;

    if (params.input_mode == HIGH_BIT_DEPTH_INTERLEAVED)
    {
        pixel_step_shift_bits = _mm_set_epi32(0, 0, 0, 1);
    } else {
        pixel_step_shift_bits = _mm_setzero_si128();
    }
    upsample_to_16_shift_bits = _mm_set_epi32(0, 0, 0, 16 - params.input_depth);

    __m128i downshift_bits = _mm_set_epi32(0, 0, 0, 16 - params.output_depth);

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
        } else {
            // info_data_stream can be NULL, in this case dummy_info_buffer will be used for temporary storage
        }
        cache = NULL;
    } else {
        // set up buffer for cache
        cache = (info_cache*)malloc(sizeof(info_cache));
        // 4 offsets (2 bytes per item) + 2-byte change
        info_data_stream = (char*)_aligned_malloc(params.info_stride * (4 * 2 + 2) * params.get_src_height(), FRAME_LUT_ALIGNMENT);
        cache->data_stream = info_data_stream;
        cache->pitch = params.src_pitch;
    }

    const int info_cache_block_size = (process_34 || sample_mode == 5 || sample_mode == 6 ? 64 : 32);

    int input_mode = params.input_mode;

    for (int row = 0; row < params.plane_height_in_pixels; row++)
    {
        const unsigned char* src_px = params.src_plane_ptr + params.src_pitch * row;
        unsigned char* dst_px = params.dst_plane_ptr + params.dst_pitch * row;

        info_ptr = params.info_ptr_base + params.info_stride * row;

        const short* grain_buffer_ptr = params.grain_buffer + params.grain_buffer_stride * row;

        int processed_pixels = 0;

        while (processed_pixels < params.plane_width_in_pixels)
        {
            __m128i change_1;
            
            __m128i ref_pixels_1_0;
            __m128i ref_pixels_2_0;
            __m128i ref_pixels_3_0;
            __m128i ref_pixels_4_0;

#define READ_REFS(data_stream, inp_mode) read_reference_pixels<sample_mode, dither_algo, inp_mode>( \
                    params, \
                    upsample_to_16_shift_bits, \
                    src_px, \
                    data_stream, \
                    ref_pixels_1_0, \
                    ref_pixels_2_0, \
                    ref_pixels_3_0, \
                    ref_pixels_4_0)

            char * data_stream_block_start;

            if (use_cached_info) {
                data_stream_block_start = info_data_stream;
                info_data_stream += info_cache_block_size;
            } else {
                // we need to process the info block

                char * data_stream_ptr = info_data_stream;
                if (!data_stream_ptr)
                {
                    data_stream_ptr = dummy_info_buffer;
                }

                data_stream_block_start = data_stream_ptr;
            
    #define PROCESS_INFO_BLOCK(n) \
                process_plane_info_block<sample_mode, n>(info_ptr, src_px, src_pitch_vector, minus_one, width_subsample_vector, height_subsample_vector, pixel_step_shift_bits, data_stream_ptr);
            
                PROCESS_INFO_BLOCK(0);
                PROCESS_INFO_BLOCK(1);

    #undef PROCESS_INFO_BLOCK
                
                if (info_data_stream) {
                    info_data_stream += info_cache_block_size;
                    assert(info_data_stream == data_stream_ptr);
                }

            }

            __m128i src_pixels;
            // abuse the guard bytes on the end of frame, as long as they are present there won't be segfault
            // garbage data is not a problem
            if (input_mode == LOW_BIT_DEPTH)
            {
                READ_REFS(data_stream_block_start, LOW_BIT_DEPTH);
                src_pixels = read_pixels<LOW_BIT_DEPTH, aligned>(params, src_px, upsample_to_16_shift_bits);
            } else if (input_mode == HIGH_BIT_DEPTH_INTERLEAVED)
            {
                READ_REFS(data_stream_block_start, HIGH_BIT_DEPTH_INTERLEAVED);
                src_pixels = read_pixels<HIGH_BIT_DEPTH_INTERLEAVED, aligned>(params, src_px, upsample_to_16_shift_bits);
            } else {
                abort();
                return;
            }

            change_1 = _mm_load_si128((__m128i*)grain_buffer_ptr);

            __m128i dst_pixels = process_pixels<sample_mode, blur_first, dither_algo>(
                                     src_pixels, 
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
                                     context_buffer);

            dst_px += store_pixels<output_mode>(dst_pixels, downshift_bits, dst_px, params.dst_pitch, params.plane_height_in_pixels);
            processed_pixels += 8;
            src_px += params.input_mode != HIGH_BIT_DEPTH_INTERLEAVED ? 8 : 16;
            grain_buffer_ptr += 8;
        }
        dither_high::next_row<dither_algo>(context_buffer);
    }
    
    dither_high::complete<dither_algo>(context_buffer);

    // for thread-safety, save context after all data is processed
    if (!use_cached_info && !context->data && cache)
    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        context->destroy = destroy_cache;
        if (context->data)
            destroy_cache(cache);
        else
            context->data = cache;
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
