#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#include "../impl_dispatch.h"
#include "simd_utils.h"
#include "dither_high.h"

#include <Windows.h>
inline static void PerformanceCounter(const char* info) noexcept
{
    static FILE* fp = fopen("perf.txt", "w");
    static long long lasttime = 0;
    LARGE_INTEGER li;
    ::QueryPerformanceCounter(&li);
    auto time = li.QuadPart;
    if(info)
        fprintf(fp, "%s %" PRId64 "\n", info, time - lasttime);
    lasttime = time;
}

/****************************************************************************
 * NOTE: DON'T remove static from any function in this file, it is required *
 *       for generating code in multiple SIMD versions.                      *
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
    const SIMD::data_type &src_pitch_vector, 
    const SIMD::data_type &minus_one, 
    const SIMD::count_type &width_subsample_vector,
    const SIMD::count_type &height_subsample_vector,
    const SIMD::count_type &pixel_step_shift_bits,
    char*& info_data_stream)
{
    assert(ref_part_index <= 2);

    SIMD::data_type info_block = SIMD::_load((SIMD::data_type*)info_ptr);

    // ref1: bit 0-7
    // left-shift & right-shift 24bits to remove other elements and preserve sign
    SIMD::data_type ref1 = info_block;
    ref1 = SIMD::_slli_epi32(ref1, 24); // << 24
    ref1 = SIMD::_srai_epi32(ref1, 24); // >> 24

    SIMD::data_type ref_offset1;
    SIMD::data_type ref_offset2;

    SIMD::data_type temp_ref1;
    switch (sample_mode)
    {
    case 0:
        // ref1 = (abs(ref1) >> height_subsampling) * (sign(ref1))
        temp_ref1 = SIMD::_abs_epi32(ref1);
        temp_ref1 = SIMD::_sra_epi32(temp_ref1, height_subsample_vector);
        temp_ref1 = SIMD::_mullo_epi32(temp_ref1, SIMD::_srai_epi32(ref1, 31));
        ref_offset1 = SIMD::_mullo_epi32(src_pitch_vector, temp_ref1); // packed DWORD multiplication
        break;
    case 1:
        // ref1 is guarenteed to be postive
        temp_ref1 = SIMD::_sra_epi32(ref1, height_subsample_vector);
        ref_offset1 = SIMD::_mullo_epi32(src_pitch_vector, temp_ref1); // packed DWORD multiplication
        break;
    case 3:
        temp_ref1 = SIMD::_sra_epi32(ref1, width_subsample_vector);
        ref_offset1 = SIMD::_sll_epi32(temp_ref1, pixel_step_shift_bits);
        break;
    case 4:
        temp_ref1 = SIMD::_sra_epi32(ref1, height_subsample_vector);
        ref_offset1 = SIMD::_mullo_epi32(src_pitch_vector, temp_ref1);
        temp_ref1 = SIMD::_sra_epi32(ref1, width_subsample_vector);
        ref_offset2 = SIMD::_sll_epi32(temp_ref1, pixel_step_shift_bits);
        break;
    case 2:
        // ref2: bit 8-15
        // similar to above
        SIMD::data_type ref2;
        ref2 = info_block;
        ref2 = SIMD::_slli_epi32(ref2, 16); // << 16
        ref2 = SIMD::_srai_epi32(ref2, 24); // >> 24

        SIMD::data_type ref1_fix, ref2_fix;
        // ref_px = src_pitch * info.ref2 + info.ref1;
        ref1_fix = SIMD::_sra_epi32(ref1, width_subsample_vector);
        ref2_fix = SIMD::_sra_epi32(ref2, height_subsample_vector);
        ref_offset1 = SIMD::_mullo_epi32(src_pitch_vector, ref2_fix); // packed DWORD multiplication
        ref_offset1 = SIMD::_add_epi32(ref_offset1, SIMD::_sll_epi32(ref1_fix, pixel_step_shift_bits));

        // ref_px_2 = info.ref2 - src_pitch * info.ref1;
        ref1_fix = SIMD::_sra_epi32(ref1, height_subsample_vector);
        ref2_fix = SIMD::_sra_epi32(ref2, width_subsample_vector);
        ref_offset2 = SIMD::_mullo_epi32(src_pitch_vector, ref1_fix); // packed DWORD multiplication
        ref_offset2 = SIMD::_sub_epi32(SIMD::_sll_epi32(ref2_fix, pixel_step_shift_bits), ref_offset2);
        break;
    default:
        abort();
    }

    if (info_data_stream){
        SIMD::_store((SIMD::data_type*)info_data_stream, ref_offset1);
        info_data_stream += SIMD::width_8;

        if (sample_mode == 2 || sample_mode == 4) {
            SIMD::_store((SIMD::data_type*)info_data_stream, ref_offset2);
            info_data_stream += SIMD::width_8;
        }
    }

    info_ptr += SIMD::width_32;
}

static __forceinline SIMD::mask_type generate_blend_mask_high(SIMD::data_type a, SIMD::data_type b, SIMD::data_type threshold)
{
    SIMD::data_type diff1 = SIMD::_subs_epu16(a, b);
    SIMD::data_type diff2 = SIMD::_subs_epu16(b, a);

    SIMD::data_type abs_diff = SIMD::_or(diff1, diff2);

    SIMD::data_type sign_convert_vector = SIMD::_set1_epi16( (short)0x8000 );

    SIMD::data_type converted_diff = SIMD::_sub_epi16(abs_diff, sign_convert_vector);

    SIMD::data_type converted_threshold = SIMD::_sub_epi16(threshold, sign_convert_vector);

    // mask: if threshold >= diff, set to 0xff, otherwise 0x00
    // note that this is the opposite of low bitdepth implementation
    return SIMD::_cmpgt_epi16(converted_threshold, converted_diff);
}


template<int sample_mode, bool blur_first>
static SIMD::data_type __forceinline process_pixels_mode12_high_part(
    SIMD::data_type src_pixels,
    SIMD::data_type threshold_vector,
    SIMD::data_type change,
    const SIMD::data_type& ref_pixels_1,
    const SIMD::data_type& ref_pixels_2,
    const SIMD::data_type& ref_pixels_3,
    const SIMD::data_type& ref_pixels_4)
{
    SIMD::mask_type use_orig_pixel_blend_mask_12, use_orig_pixel_blend_mask_34;
    SIMD::data_type avg_12, avg_34;
    SIMD::data_type dst_pixels;

    if (true)
    {
        if (!blur_first)
        {
            // note: use AND instead of OR, because two operands are reversed
            // (different from low bit-depth mode!)
            use_orig_pixel_blend_mask_12 = SIMD::_and(
                generate_blend_mask_high(src_pixels, ref_pixels_1, threshold_vector),
                generate_blend_mask_high(src_pixels, ref_pixels_2, threshold_vector) );
        }

        avg_12 = SIMD::_avg_epu16(ref_pixels_1, ref_pixels_2);
    }

    if (process_34)
    {
        if (!blur_first)
        {
            use_orig_pixel_blend_mask_34 = SIMD::_and(
                generate_blend_mask_high(src_pixels, ref_pixels_3, threshold_vector),
                generate_blend_mask_high(src_pixels, ref_pixels_4, threshold_vector) );
        }

        avg_34 = SIMD::_avg_epu16(ref_pixels_3, ref_pixels_4);
    }

    if (sample_mode == 2)
    {
        // merge 34 into 12 for mode 2
        if (!blur_first)
            use_orig_pixel_blend_mask_12 = SIMD::_and(use_orig_pixel_blend_mask_12, use_orig_pixel_blend_mask_34 );
        avg_12 = SIMD::_subs_epu16(avg_12, SIMD::_set1_epi16(1));
        avg_12 = SIMD::_avg_epu16(avg_12, avg_34);
    }

    if (blur_first)
        use_orig_pixel_blend_mask_12 = generate_blend_mask_high(src_pixels, avg_12, threshold_vector);

    if (process_34 && blur_first)
        use_orig_pixel_blend_mask_34 = generate_blend_mask_high(src_pixels, avg_34, threshold_vector);
    
    if (sample_mode == 4)
    {
        // average dst_12 and dst_34 for mode 4
        SIMD::data_type dst_pixels_12 = SIMD::_blendv_epi8(src_pixels, avg_12, use_orig_pixel_blend_mask_12);
        SIMD::data_type dst_pixels_34 = SIMD::_blendv_epi8(src_pixels, avg_34, use_orig_pixel_blend_mask_34);
        dst_pixels = SIMD::_avg_epu16(dst_pixels_12, dst_pixels_34);
    }
    else
    {
        // if mask is 0xff (NOT over threshold), select second operand, otherwise select first
        // note this is different from low bitdepth code
        dst_pixels = SIMD::_blendv_epi8(src_pixels, avg_12, use_orig_pixel_blend_mask_12);
    }

    SIMD::data_type sign_convert_vector = SIMD::_set1_epi16((short)0x8000);

    // convert to signed form, since change is signed
    dst_pixels = SIMD::_sub_epi16(dst_pixels, sign_convert_vector);

    // saturated add
    dst_pixels = SIMD::_adds_epi16(dst_pixels, change);

    // convert back to unsigned
    dst_pixels = SIMD::_add_epi16(dst_pixels, sign_convert_vector);
    return dst_pixels;
}

template<int sample_mode, bool blur_first, int dither_algo>
static SIMD::data_type __forceinline process_pixels(
    SIMD::data_type src_pixels_0, 
    SIMD::data_type threshold_vector, 
    const SIMD::data_type& change_1, 
    const SIMD::data_type& ref_pixels_1_0,
    const SIMD::data_type& ref_pixels_2_0,
    const SIMD::data_type& ref_pixels_3_0,
    const SIMD::data_type& ref_pixels_4_0,
    const SIMD::data_type& clamp_high_add,
    const SIMD::data_type& clamp_high_sub,
    const SIMD::data_type& clamp_low,
    bool need_clamping,
    int row,
    int column,
    void* dither_context)
{
    SIMD::data_type ret = process_pixels_mode12_high_part<sample_mode, blur_first>
        (src_pixels_0, 
         threshold_vector, 
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
    SIMD::data_type pixels,
    SIMD::count_type downshift_bits,
    unsigned char* dst,
    int dst_pitch,
    int height_in_pixels)
{
    switch (output_mode)
    {
    case LOW_BIT_DEPTH:
        pixels = SIMD::_srli_epi16(pixels, 8);
        pixels = SIMD::_packus_epi16(pixels, pixels);
        SIMD::_storel((SIMD::half_data_type*)dst, pixels);
        return SIMD::width_16;
    case HIGH_BIT_DEPTH_INTERLEAVED:
        pixels = SIMD::_srl_epi16(pixels, downshift_bits);
        SIMD::_store((SIMD::data_type*)dst, pixels);
        return SIMD::width_16 * 2;
    default:
        abort();
    }
    return 0;
}


template<bool aligned>
static SIMD::data_type load_data(const unsigned char *ptr)
{
    if (aligned)
    {
        return SIMD::_load((const SIMD::data_type*)ptr);
    } else {
        return SIMD::_loadu((const SIMD::data_type*)ptr);
    }
}

template<PIXEL_MODE input_mode, bool aligned>
static SIMD::data_type __forceinline read_pixels(
    const process_plane_params& params,
    const unsigned char *ptr, 
    SIMD::count_type upsample_shift)
{
    SIMD::data_type ret;

    switch (input_mode)
    {
    case LOW_BIT_DEPTH:
        {
            auto zero = SIMD::_setzero();
            return SIMD::_unpacklo_epi8(zero, SIMD::_loadl((SIMD::half_data_type*)ptr));
        }
        break;
    case HIGH_BIT_DEPTH_INTERLEAVED:
        ret = load_data<aligned>(ptr);
        break;
    default:
        abort();
    }
    ret = SIMD::_sll_epi16(ret, upsample_shift);
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
static SIMD::data_type __forceinline load_reference_pixels(
    SIMD::count_type shift,
    const unsigned short src[SIMD::width_16])
{
    SIMD::data_type ret = SIMD::_load((const SIMD::data_type*)src);
    ret = SIMD::_sll_epi16(ret, shift);
    return ret;
}

template<int sample_mode, int dither_algo, PIXEL_MODE input_mode>
static void __forceinline read_reference_pixels(
    const process_plane_params& params,
    SIMD::count_type shift,
    const unsigned char* src_px_start,
    const char* info_data_start,
    SIMD::data_type& real_ret_1,
    SIMD::data_type& real_ret_2,
    SIMD::data_type& real_ret_3,
    SIMD::data_type& real_ret_4)
{
    #pragma GCC diagnostic ignored "-Wuninitialized"
    __m128i tmp_1;
    #pragma GCC diagnostic ignored "-Wuninitialized"
    __m128i tmp_2;
    #pragma GCC diagnostic ignored "-Wuninitialized"
    __m128i tmp_3;
    #pragma GCC diagnostic ignored "-Wuninitialized"
    __m128i tmp_4;
    SIMD::data_type ret_1;
    SIMD::data_type ret_2;
    SIMD::data_type ret_3;
    SIMD::data_type ret_4;

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
    
    #define UNROLLED_LOOP(i) \
    do { \
        int d1 = *(int*)(info_data_start + 4 * i); \
        int d2 = *(int*)(info_data_start + 4 * (i + i / SIMD::width_32 * SIMD::width_32)); \
        int d3 = *(int*)(info_data_start + 4 * (i + i / SIMD::width_32 * SIMD::width_32 + SIMD::width_32)); \
        switch (sample_mode) \
        { \
        case 1: \
        case 3: \
            tmp_1 = _mm_insert_epi16(tmp_1, read_pixel<input_mode>(src_px_start, i_fix + d1), i); \
            tmp_2 = _mm_insert_epi16(tmp_2, read_pixel<input_mode>(src_px_start, i_fix + -d1), i); \
            break; \
        case 2: \
        case 4: \
            tmp_1 = _mm_insert_epi16(tmp_1, read_pixel<input_mode>(src_px_start, i_fix + d2), i); \
            tmp_2 = _mm_insert_epi16(tmp_2, read_pixel<input_mode>(src_px_start, i_fix + -d2), i); \
            tmp_3 = _mm_insert_epi16(tmp_3, read_pixel<input_mode>(src_px_start, i_fix + d3), i); \
            tmp_4 = _mm_insert_epi16(tmp_4, read_pixel<input_mode>(src_px_start, i_fix + -d3), i); \
            break; \
        } \
        i_fix += i_fix_step; \
    } while (0)

    #if defined SIMD_SSE4
    UNROLLED_LOOP(0); UNROLLED_LOOP(1); UNROLLED_LOOP(2); UNROLLED_LOOP(3); UNROLLED_LOOP(4); UNROLLED_LOOP(5); UNROLLED_LOOP(6); UNROLLED_LOOP(7);
    ret_1 = tmp_1;
    ret_2 = tmp_2;
    ret_3 = tmp_3;
    ret_4 = tmp_4;

    #elif defined SIMD_AVX2
    UNROLLED_LOOP(0); UNROLLED_LOOP(1); UNROLLED_LOOP(2); UNROLLED_LOOP(3); UNROLLED_LOOP(4); UNROLLED_LOOP(5); UNROLLED_LOOP(6); UNROLLED_LOOP(7);
    ret_1 = _mm256_castsi128_si256(tmp_1);
    ret_2 = _mm256_castsi128_si256(tmp_2);
    ret_3 = _mm256_castsi128_si256(tmp_3);
    ret_4 = _mm256_castsi128_si256(tmp_4);
    UNROLLED_LOOP(0); UNROLLED_LOOP(1); UNROLLED_LOOP(2); UNROLLED_LOOP(3); UNROLLED_LOOP(4); UNROLLED_LOOP(5); UNROLLED_LOOP(6); UNROLLED_LOOP(7);
    ret_1 = _mm256_inserti128_si256(ret_1, tmp_1, 1);
    ret_2 = _mm256_inserti128_si256(ret_2, tmp_2, 1);
    ret_3 = _mm256_inserti128_si256(ret_3, tmp_3, 1);
    ret_4 = _mm256_inserti128_si256(ret_4, tmp_4, 1);

    #elif defined SIMD_AVX512
    UNROLLED_LOOP(0); UNROLLED_LOOP(1); UNROLLED_LOOP(2); UNROLLED_LOOP(3); UNROLLED_LOOP(4); UNROLLED_LOOP(5); UNROLLED_LOOP(6); UNROLLED_LOOP(7);
    ret_1 = _mm512_castsi128_si512(tmp_1);
    ret_2 = _mm512_castsi128_si512(tmp_2);
    ret_3 = _mm512_castsi128_si512(tmp_3);
    ret_4 = _mm512_castsi128_si512(tmp_4);
    UNROLLED_LOOP(0); UNROLLED_LOOP(1); UNROLLED_LOOP(2); UNROLLED_LOOP(3); UNROLLED_LOOP(4); UNROLLED_LOOP(5); UNROLLED_LOOP(6); UNROLLED_LOOP(7);
    ret_1 = _mm512_inserti32x4(ret_1, tmp_1, 1);
    ret_2 = _mm512_inserti32x4(ret_2, tmp_2, 1);
    ret_3 = _mm512_inserti32x4(ret_3, tmp_3, 1);
    ret_4 = _mm512_inserti32x4(ret_4, tmp_4, 1);
    UNROLLED_LOOP(0); UNROLLED_LOOP(1); UNROLLED_LOOP(2); UNROLLED_LOOP(3); UNROLLED_LOOP(4); UNROLLED_LOOP(5); UNROLLED_LOOP(6); UNROLLED_LOOP(7);
    ret_1 = _mm512_inserti32x4(ret_1, tmp_1, 2);
    ret_2 = _mm512_inserti32x4(ret_2, tmp_2, 2);
    ret_3 = _mm512_inserti32x4(ret_3, tmp_3, 2);
    ret_4 = _mm512_inserti32x4(ret_4, tmp_4, 2);
    UNROLLED_LOOP(0); UNROLLED_LOOP(1); UNROLLED_LOOP(2); UNROLLED_LOOP(3); UNROLLED_LOOP(4); UNROLLED_LOOP(5); UNROLLED_LOOP(6); UNROLLED_LOOP(7);
    ret_1 = _mm512_inserti32x4(ret_1, tmp_1, 3);
    ret_2 = _mm512_inserti32x4(ret_2, tmp_2, 3);
    ret_3 = _mm512_inserti32x4(ret_3, tmp_3, 3);
    ret_4 = _mm512_inserti32x4(ret_4, tmp_4, 3);

    #endif

    #undef UNROLLED_LOOP

    real_ret_1 = ret_1;
    real_ret_2 = ret_2;
    if (process_34)
    {
        real_ret_3 = ret_3;
        real_ret_4 = ret_4;
    }

    // Search vpmaxsb in asm
    // real_ret_1 = SIMD::_max_epi8(real_ret_1, real_ret_2);
}


template<int sample_mode, bool blur_first, int dither_algo, bool aligned, PIXEL_MODE output_mode>
static void __cdecl _process_plane_simd_impl(const process_plane_params& params, process_plane_context* context)
{
    assert(sample_mode > 0);

    pixel_dither_info* info_ptr = params.info_ptr_base;

    SIMD::data_type src_pitch_vector = SIMD::_set1_epi32(params.src_pitch);
           
    SIMD::data_type threshold_vector = SIMD::_set1_epi16(params.threshold);

    // general-purpose constant
    SIMD::data_type minus_one = SIMD::_set1_epi32(-1);

    alignas(16)
    char context_buffer[DITHER_CONTEXT_BUFFER_SIZE];

    dither_high::init<dither_algo>(context_buffer, params.plane_width_in_pixels, params.output_depth);

    SIMD::count_type width_subsample_vector = SIMD::count_set(params.width_subsampling);
    SIMD::count_type height_subsample_vector = SIMD::count_set(params.height_subsampling);

    bool need_clamping =  INTERNAL_BIT_DEPTH < 16 || 
                          params.pixel_min > 0 || 
                          params.pixel_max < 0xffff;
    SIMD::data_type clamp_high_add = SIMD::_setzero();
    SIMD::data_type clamp_high_sub = SIMD::_setzero();
    SIMD::data_type clamp_low = SIMD::_setzero();
    if (need_clamping)
    {
        clamp_low = SIMD::_set1_epi16((short)params.pixel_min);
        clamp_high_add = SIMD::_sub_epi16(SIMD::_set1_epi16((short)0xffff), SIMD::_set1_epi16((short)params.pixel_max));
        clamp_high_sub = SIMD::_add_epi16(clamp_high_add, clamp_low);
    }
    
    SIMD::count_type pixel_step_shift_bits;
    SIMD::count_type upsample_to_16_shift_bits;

    if (params.input_mode == HIGH_BIT_DEPTH_INTERLEAVED)
    {
        pixel_step_shift_bits = SIMD::count_set(1);
    } else {
        pixel_step_shift_bits = SIMD::count_zero();
    }
    upsample_to_16_shift_bits = SIMD::count_set(16 - params.input_depth);

    SIMD::count_type downshift_bits = SIMD::count_set(16 - params.output_depth);

    bool use_cached_info = false;
    info_cache *cache = NULL;
    char* info_data_stream = NULL;

    alignas(SIMD::align)
    char dummy_info_buffer[SIMD::width_8 * 8];

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

    const int info_cache_block_size = SIMD::width_16 * (process_34 ? 8 : 4);

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
            SIMD::data_type change_1;
            
            SIMD::data_type ref_pixels_1_0;
            SIMD::data_type ref_pixels_2_0;
            SIMD::data_type ref_pixels_3_0;
            SIMD::data_type ref_pixels_4_0;

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

            SIMD::data_type src_pixels;
            // abuse the guard bytes on the end of frame, as long as they are present there won't be segfault
            // garbage data is not a problem
            if (input_mode == LOW_BIT_DEPTH)
            {
                #if 1
                read_reference_pixels<sample_mode, dither_algo, LOW_BIT_DEPTH>( \
                    params, \
                    upsample_to_16_shift_bits, \
                    src_px, \
                    data_stream_block_start, \
                    ref_pixels_1_0, \
                    ref_pixels_2_0, \
                    ref_pixels_3_0, \
                    ref_pixels_4_0);
                #else
                const char* info_data_start = data_stream_block_start;
                alignas(SIMD::align)
                unsigned short tmp_1[SIMD::width_16];
                alignas(SIMD::align)
                unsigned short tmp_2[SIMD::width_16];
                alignas(SIMD::align)
                unsigned short tmp_3[SIMD::width_16];
                alignas(SIMD::align)
                unsigned short tmp_4[SIMD::width_16];

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
                
                for (int i = 0; i < SIMD::width_16; i++)
                {
                    int d1 = *(int*)(info_data_start + 4 * i);
                    int d2 = *(int*)(info_data_start + 4 * (i + i / SIMD::width_32 * SIMD::width_32));
                    int d3 = *(int*)(info_data_start + 4 * (i + i / SIMD::width_32 * SIMD::width_32 + SIMD::width_32));
                    switch (sample_mode)
                    {
                    case 0:
                        tmp_1[i] = read_pixel<LOW_BIT_DEPTH>(src_px, i_fix + d1);
                        break;
                    case 1:
                    case 3:
                        tmp_1[i] = read_pixel<LOW_BIT_DEPTH>(src_px, i_fix + d1);
                        tmp_2[i] = read_pixel<LOW_BIT_DEPTH>(src_px, i_fix + -d1);
                        break;
                    case 2:
                    case 4:
                        tmp_1[i] = read_pixel<LOW_BIT_DEPTH>(src_px, i_fix + d2);
                        tmp_2[i] = read_pixel<LOW_BIT_DEPTH>(src_px, i_fix + -d2);
                        tmp_3[i] = read_pixel<LOW_BIT_DEPTH>(src_px, i_fix + d3);
                        tmp_4[i] = read_pixel<LOW_BIT_DEPTH>(src_px, i_fix + -d3);
                        break;
                    }
                    i_fix += i_fix_step;
                }

                switch (sample_mode)
                {
                case 0:
                    ref_pixels_1_0 = load_reference_pixels<dither_algo>(upsample_to_16_shift_bits, tmp_1);
                    break;
                case 1:
                case 3:
                    ref_pixels_1_0 = load_reference_pixels<dither_algo>(upsample_to_16_shift_bits, tmp_1);
                    ref_pixels_2_0 = load_reference_pixels<dither_algo>(upsample_to_16_shift_bits, tmp_2);
                    break;
                case 2:
                case 4:
                    ref_pixels_1_0 = load_reference_pixels<dither_algo>(upsample_to_16_shift_bits, tmp_1);
                    ref_pixels_2_0 = load_reference_pixels<dither_algo>(upsample_to_16_shift_bits, tmp_2);
                    ref_pixels_3_0 = load_reference_pixels<dither_algo>(upsample_to_16_shift_bits, tmp_3);
                    ref_pixels_4_0 = load_reference_pixels<dither_algo>(upsample_to_16_shift_bits, tmp_4);
                    break;
                }
                #endif

                src_pixels = read_pixels<LOW_BIT_DEPTH, aligned>(params, src_px, upsample_to_16_shift_bits);
            } else if (input_mode == HIGH_BIT_DEPTH_INTERLEAVED)
            {
                READ_REFS(data_stream_block_start, HIGH_BIT_DEPTH_INTERLEAVED);
                src_pixels = read_pixels<HIGH_BIT_DEPTH_INTERLEAVED, aligned>(params, src_px, upsample_to_16_shift_bits);
            } else {
                abort();
                return;
            }

            change_1 = SIMD::_load((SIMD::data_type*)grain_buffer_ptr);

            SIMD::data_type dst_pixels = process_pixels<sample_mode, blur_first, dither_algo>(
                                     src_pixels, 
                                     threshold_vector,
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
            processed_pixels += SIMD::width_16;
            src_px += SIMD::width_16 * (params.input_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1);
            grain_buffer_ptr += SIMD::width_16;
        }
        dither_high::next_row<dither_algo>(context_buffer);
    }
    
    dither_high::complete<dither_algo>(context_buffer);

    // for thread-safety, save context after all data is processed
    if (!use_cached_info && !context->data && cache)
    {
        context->destroy = destroy_cache;
        if (_InterlockedCompareExchangePointer(&context->data, cache, NULL) != NULL)
        {
            // other thread has completed first, so we can destroy our copy
            destroy_cache(cache);
        }
    }
}


template<int sample_mode, bool blur_first, int dither_algo, bool aligned>
static void process_plane_simd_impl_stub1(const process_plane_params& params, process_plane_context* context)
{
    switch (params.output_mode)
    {
    case LOW_BIT_DEPTH:
        _process_plane_simd_impl<sample_mode, blur_first, dither_algo, aligned, LOW_BIT_DEPTH>(params, context);
        break;
    case HIGH_BIT_DEPTH_INTERLEAVED:
        _process_plane_simd_impl<sample_mode, blur_first, dither_algo, aligned, HIGH_BIT_DEPTH_INTERLEAVED>(params, context);
        break;
    default:
        abort();
    }
}

template<int sample_mode, bool blur_first, int dither_algo>
static void __cdecl process_plane_simd_impl(const process_plane_params& params, process_plane_context* context)
{
    if ( ( (intptr_t)params.src_plane_ptr & (SIMD::align - 1) ) == 0 && (params.src_pitch & (SIMD::align - 1) ) == 0 )
    {
        PerformanceCounter(NULL);
        process_plane_simd_impl_stub1<sample_mode, blur_first, dither_algo, true>(params, context);
        PerformanceCounter("aligned");
    } else {
        PerformanceCounter(NULL);
        process_plane_simd_impl_stub1<sample_mode, blur_first, dither_algo, false>(params, context);
        PerformanceCounter("unaligned");
    }
}
