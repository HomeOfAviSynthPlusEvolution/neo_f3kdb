#pragma once

#include <algorithm>
#include <mutex>
#include <stdlib.h>
#include <type_traits>

#include "impl_dispatch.h"
#include "sse_utils.h"
#include "dither_high.h"
// VCL2 headers provide a modern C++ interface for SIMD intrinsics.
#include "VCL2/vectorclass.h"
#include "VCL2/vectormath_exp.h"
#include "VCL2/vectormath_trig.h"

/****************************************************************************
 * NOTE: DON'T remove static from any function in this file, it is required *
 *       for generating code in multiple instruction sets.                  *
 ****************************************************************************/

using V_int = Vec8i;
using V_float = 
std::conditional_t<std::is_same_v<V_int, Vec4i>, Vec4f, std::conditional_t<std::is_same_v<V_int, Vec8i>, Vec8f, Vec16f>>;
using V_fbool =
std::conditional_t<std::is_same_v<V_int, Vec4i>, Vec4fb, std::conditional_t<std::is_same_v<V_int, Vec8i>, Vec8fb, Vec16fb>>;
using V_ushort =
std::conditional_t<std::is_same_v<V_int, Vec4i>, Vec8us, std::conditional_t<std::is_same_v<V_int, Vec8i>, Vec16us, Vec32us>>;
using V_short =
std::conditional_t<std::is_same_v<V_int, Vec4i>, Vec8s, std::conditional_t<std::is_same_v<V_int, Vec8i>, Vec16s, Vec32s>>;
using V_sbool =
std::conditional_t<std::is_same_v<V_int, Vec4i>, Vec8sb, std::conditional_t<std::is_same_v<V_int, Vec8i>, Vec16sb, Vec32sb>>;
using V_uchar =
std::conditional_t<std::is_same_v<V_int, Vec4i>, Vec16uc, std::conditional_t<std::is_same_v<V_int, Vec8i>, Vec32uc, Vec64uc>>;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct _info_cache_avx2
{
    int pitch;
    char* data_stream;
} info_cache_avx2;

static void destroy_cache_avx2(void* data)
{
    assert(data);
    auto* cache = reinterpret_cast<info_cache_avx2*>(data);
    _aligned_free(cache->data_stream);
    free(data);
}

template <typename V>
static auto __forceinline high_bit_depth_pixels_clamp_avx2(V pixels, V high_add, V high_sub, const V& low)
{
    pixels = add_saturated(pixels, high_add);
    pixels = sub_saturated(pixels, high_sub);
    return pixels + low;
}

template <typename V>
static __forceinline V saturate(V const& v)
{
    return max(0.0f, min(1.0f, v));
}

template <typename V, int sample_mode>
static __forceinline void process_plane_info_block_avx2_16px(
    pixel_dither_info*& info_ptr,
    const V& src_pitch_vector,
    const int width_subsample,
    const int height_subsample,
    const int pixel_step_shift_bits,
    char*& info_data_stream)
{
    // Process first block of pixels
    {
        auto ref1 = (V().load(reinterpret_cast<const int32_t*>(info_ptr)) << 24) >> 24;
        auto temp_ref1_h = ref1 >> height_subsample;
        auto ref_offset1 = src_pitch_vector * temp_ref1_h;
        auto temp_ref1_w = ref1 >> width_subsample;
        auto ref_offset2 = temp_ref1_w << pixel_step_shift_bits;

        if (info_data_stream) {
            ref_offset1.store(reinterpret_cast<int32_t*>(info_data_stream));
            info_data_stream += sizeof(V);
            ref_offset2.store(reinterpret_cast<int32_t*>(info_data_stream));
            info_data_stream += sizeof(V);
        }
    }

    // Process next block of pixels
    {
        auto ref1 = (V().load(reinterpret_cast<const int32_t*>(info_ptr + V_int().size())) << 24) >> 24;
        auto temp_ref1_h = ref1 >> height_subsample;
        auto ref_offset1 = src_pitch_vector * temp_ref1_h;
        auto temp_ref1_w = ref1 >> width_subsample;
        auto ref_offset2 = temp_ref1_w << pixel_step_shift_bits;

        if (info_data_stream) {
            ref_offset1.store(reinterpret_cast<int32_t*>(info_data_stream));
            info_data_stream += sizeof(V);
            ref_offset2.store(reinterpret_cast<int32_t*>(info_data_stream));
            info_data_stream += sizeof(V);
        }
    }

    info_ptr += V_ushort().size();
}

template <typename V, typename V_bool>
static V_bool __forceinline generate_blend_mask_high_avx2(V a, V threshold)
{
    return V_bool(a < threshold);
}

template <typename V, typename V_float>
static __forceinline V_float gather_pixel_values_avx2(const process_plane_params& params, V const y_coords, V const x_coords,
    int upsample_shift)
{
    auto clamped_y = max(V(0), min(y_coords, params.plane_height_in_pixels - 1));
    auto clamped_x = max(V(0), min(x_coords, params.plane_width_in_pixels - 1));

    V pitch(params.src_pitch);
    V pixel_bytes_v((params.input_mode == HIGH_BIT_DEPTH_INTERLEAVED) ? 2 : 1);

    auto byte_offsets = clamped_y * pitch + clamped_x * pixel_bytes_v;

    const unsigned char* base_ptr = params.src_plane_ptr;

    alignas(32)
        int32_t offsets[V_int().size()];
    byte_offsets.store(offsets);

    auto shorts_i32 = [&]() {
        if (params.input_mode == LOW_BIT_DEPTH) {
            alignas(8)
                uint8_t pixels_buf[V_int().size()];

            for (int i = 0; i < V_int().size(); ++i) {
                const unsigned char* pixel_address = base_ptr + offsets[i];
                pixels_buf[i] = *pixel_address;
            }

            if constexpr (std::is_same_v<V, Vec4i>)
                return V().load_4uc(pixels_buf);
            else if constexpr (std::is_same_v<V, Vec8i>)
                return V().load_8uc(pixels_buf);
            else
                return V().load_16uc(pixels_buf);
        }
        else {
            alignas(16)
                uint16_t pixels_buf[V_int().size()];

            for (int i = 0; i < V_int().size(); ++i) {
                const unsigned char* pixel_address = base_ptr + offsets[i];
                pixels_buf[i] = *reinterpret_cast<const uint16_t*>(pixel_address);
            }

            if constexpr (std::is_same_v<V, Vec4i>)
                return V().load_4us(pixels_buf);
            else if constexpr (std::is_same_v<V, Vec8i>)
                return V().load_8us(pixels_buf);
            else
                return V().load_16us(pixels_buf);
        }
        }();

    return to_float(shorts_i32 << upsample_shift);
}

template <typename V, typename V_float>
static __forceinline V_float calculate_gradient_angle_avx2(const process_plane_params& params, const V& y_coords, const V& x_coords,
    int read_distance, int upsample_shift)
{
    V rd(read_distance);
    auto p00 = gather_pixel_values_avx2<V, V_float> (params, y_coords - rd, x_coords - rd, upsample_shift);
    auto p10 = gather_pixel_values_avx2<V, V_float>(params, y_coords - rd, x_coords, upsample_shift);
    auto p20 = gather_pixel_values_avx2<V, V_float>(params, y_coords - rd, x_coords + rd, upsample_shift);
    auto p01 = gather_pixel_values_avx2<V, V_float>(params, y_coords, x_coords - rd, upsample_shift);
    auto p21 = gather_pixel_values_avx2<V, V_float>(params, y_coords, x_coords + rd, upsample_shift);
    auto p02 = gather_pixel_values_avx2<V, V_float>(params, y_coords + rd, x_coords - rd, upsample_shift);
    auto p12 = gather_pixel_values_avx2<V, V_float>(params, y_coords + rd, x_coords, upsample_shift);
    auto p22 = gather_pixel_values_avx2<V, V_float>(params, y_coords + rd, x_coords + rd, upsample_shift);

    auto gx = (p20 + 2.0f * p21 + p22) - (p00 + 2.0f * p01 + p02);
    auto gy = (p00 + 2.0f * p10 + p20) - (p02 + 2.0f * p12 + p22);

    const float scaled_epsilon_for_gx = 0.01f * (static_cast<float>(1 << (16 - params.input_depth)) * 3.0f);

    V_fbool gx_is_small = abs(gx) < scaled_epsilon_for_gx;

    auto angle = atan(gy / gx);
    angle = select(gx_is_small, 1.0f, angle / static_cast<float>(M_PI) + 0.5f);
    return angle;
}

template<typename V, typename V_signed, int sample_mode, bool blur_first, int dither_algo>
static auto __forceinline process_pixels_avx2(V src_pixels, V_signed change, const V& ref_pixels_1, const V& ref_pixels_2,
    const V& ref_pixels_3, const V& ref_pixels_4, const V& clamp_high_add, const V& clamp_high_sub, const V& clamp_low, bool need_clamping,
    int row, int column, void* dither_context, const pixel_dither_info* pdi_ptr, const process_plane_params& params,
    int upsample_to_16_shift_bits)
{
    const int threshold = params.threshold;
    const int threshold1 = params.threshold1;
    const int threshold2 = params.threshold2;

    V dst_pixels = zero_si256();

    if constexpr (sample_mode == 5) {
        V_int r1_lo = extend_low(ref_pixels_1);
        V_int r1_hi = extend_high(ref_pixels_1);
        V_int r2_lo = extend_low(ref_pixels_2);
        V_int r2_hi = extend_high(ref_pixels_2);
        V_int r3_lo = extend_low(ref_pixels_3);
        V_int r3_hi = extend_high(ref_pixels_3);
        V_int r4_lo = extend_low(ref_pixels_4);
        V_int r4_hi = extend_high(ref_pixels_4);

        auto sum_lo = r1_lo + r2_lo + r3_lo + r4_lo;
        auto sum_hi = r1_hi + r2_hi + r3_hi + r4_hi;

        auto avg = compress(sum_lo >> 2, sum_hi >> 2);

        auto avgDif = abs(avg - src_pixels);
        auto maxDif = max(
            max(abs(ref_pixels_1 - src_pixels), abs(ref_pixels_2 - src_pixels)),
            max(abs(ref_pixels_3 - src_pixels), abs(ref_pixels_4 - src_pixels))
        );

        auto src_lo = extend_low(src_pixels); auto src_hi = extend_high(src_pixels);
        auto two_src_lo = src_lo << 1;
        auto two_src_hi = src_hi << 1;

        auto midDif1 = compress(abs((r1_lo + r2_lo) - two_src_lo), abs((r1_hi + r2_hi) - two_src_hi));
        auto midDif2 = compress(abs((r3_lo + r4_lo) - two_src_lo), abs((r3_hi + r4_hi) - two_src_hi));

        auto use_orig_pixel_blend_mask = generate_blend_mask_high_avx2<V_short, V_sbool>(avgDif, V_short(threshold))
            & generate_blend_mask_high_avx2<V_short, V_sbool>(maxDif, V_short(threshold1))
            & generate_blend_mask_high_avx2<V_short, V_sbool>(midDif1, V_short(threshold2))
            & generate_blend_mask_high_avx2<V_short, V_sbool>(midDif2, V_short(threshold2));

        dst_pixels = select(use_orig_pixel_blend_mask, avg, src_pixels);
    }
    else { // sample_mode 6 or 7
        auto src_f_lo = to_float(extend(src_pixels.get_low()));
        auto src_f_hi = to_float(extend(src_pixels.get_high()));

        auto p1_f_lo = to_float(extend(ref_pixels_1.get_low()));
        auto p1_f_hi = to_float(extend(ref_pixels_1.get_high()));

        auto p2_f_lo = to_float(extend(ref_pixels_2.get_low()));
        auto p2_f_hi = to_float(extend(ref_pixels_2.get_high()));

        auto p3_f_lo = to_float(extend(ref_pixels_3.get_low()));
        auto p3_f_hi = to_float(extend(ref_pixels_3.get_high()));

        auto p4_f_lo = to_float(extend(ref_pixels_4.get_low()));
        auto p4_f_hi = to_float(extend(ref_pixels_4.get_high()));

        V_float current_thresh_avg_dif_lo(threshold);
        V_float current_thresh_avg_dif_hi(threshold);

        V_float current_thresh_max_dif_lo(threshold1);
        V_float current_thresh_max_dif_hi(threshold1);

        V_float current_thresh_mid_dif_lo(threshold2);
        V_float current_thresh_mid_dif_hi(threshold2);

        if constexpr (sample_mode == 7) {
            const int grad_read_distance = 20;
            const float angle_boost_factor = params.angle_boost;
            const float max_angle_threshold = params.max_angle;

            auto base_x_coords_lo = V_int(column) + V_int(0, 1, 2, 3, 4, 5, 6, 7);
            auto base_x_coords_hi = V_int(column + V_int().size()) + V_int(0, 1, 2, 3, 4, 5, 6, 7);
            V_int base_y_coords(row);

            auto angle_org_lo = calculate_gradient_angle_avx2<V_int, V_float>(params, base_y_coords, base_x_coords_lo,
                grad_read_distance, upsample_to_16_shift_bits);
            auto angle_org_hi = calculate_gradient_angle_avx2<V_int, V_float>(params, base_y_coords, base_x_coords_hi,
                grad_read_distance, upsample_to_16_shift_bits);

            alignas(32)
                int32_t ref1_buffer[V_ushort().size()];

            for (int k = 0; k < V_ushort().size(); ++k)
                ref1_buffer[k] = pdi_ptr[k].ref1;

            auto ref1_offsets_lo(V_int().load_a(ref1_buffer));
            auto ref1_offsets_hi(V_int().load_a(ref1_buffer + V_int().size()));

            auto y_offsets_h_lo = ref1_offsets_lo >> params.height_subsampling;
            auto y_offsets_h_hi = ref1_offsets_hi >> params.height_subsampling;

            auto x_offsets_w_lo = ref1_offsets_lo >> params.width_subsampling;
            auto x_offsets_w_hi = ref1_offsets_hi >> params.width_subsampling;

            auto angle_ref1_h_lo = calculate_gradient_angle_avx2<V_int, V_float>(params, base_y_coords + y_offsets_h_lo,
                base_x_coords_lo, grad_read_distance, upsample_to_16_shift_bits);
            auto angle_ref1_h_hi = calculate_gradient_angle_avx2<V_int, V_float>(params, base_y_coords + y_offsets_h_hi,
                base_x_coords_hi, grad_read_distance, upsample_to_16_shift_bits);

            auto angle_ref2_h_lo = calculate_gradient_angle_avx2<V_int, V_float>(params, base_y_coords - y_offsets_h_lo,
                base_x_coords_lo, grad_read_distance, upsample_to_16_shift_bits);
            auto angle_ref2_h_hi = calculate_gradient_angle_avx2<V_int, V_float>(params, base_y_coords - y_offsets_h_hi,
                base_x_coords_hi, grad_read_distance, upsample_to_16_shift_bits);

            auto angle_ref1_w_lo = calculate_gradient_angle_avx2<V_int, V_float>(params, base_y_coords,
                base_x_coords_lo + x_offsets_w_lo, grad_read_distance, upsample_to_16_shift_bits);
            auto angle_ref1_w_hi = calculate_gradient_angle_avx2<V_int, V_float>(params, base_y_coords,
                base_x_coords_hi + x_offsets_w_hi, grad_read_distance, upsample_to_16_shift_bits);

            auto angle_ref2_w_lo = calculate_gradient_angle_avx2<V_int, V_float>(params, base_y_coords,
                base_x_coords_lo - x_offsets_w_lo, grad_read_distance, upsample_to_16_shift_bits);
            auto angle_ref2_w_hi = calculate_gradient_angle_avx2<V_int, V_float>(params, base_y_coords,
                base_x_coords_hi - x_offsets_w_hi, grad_read_distance, upsample_to_16_shift_bits);

            auto max_angle_diff_lo = max(abs(angle_ref1_h_lo - angle_org_lo), abs(angle_ref2_h_lo - angle_org_lo));
            auto max_angle_diff_hi = max(abs(angle_ref1_h_hi - angle_org_hi), abs(angle_ref2_h_hi - angle_org_hi));

            max_angle_diff_lo = max(max_angle_diff_lo, max(abs(angle_ref1_w_lo - angle_org_lo), abs(angle_ref2_w_lo - angle_org_lo)));
            max_angle_diff_hi = max(max_angle_diff_hi, max(abs(angle_ref1_w_hi - angle_org_hi), abs(angle_ref2_w_hi - angle_org_hi)));

            decltype(max_angle_diff_lo > max_angle_diff_lo) use_boost_lo = max_angle_diff_lo <= max_angle_threshold;
            decltype(max_angle_diff_hi > max_angle_diff_hi) use_boost_hi = max_angle_diff_hi <= max_angle_threshold;

            current_thresh_avg_dif_lo = select(use_boost_lo, current_thresh_avg_dif_lo * angle_boost_factor, current_thresh_avg_dif_lo);
            current_thresh_avg_dif_hi = select(use_boost_hi, current_thresh_avg_dif_hi * angle_boost_factor, current_thresh_avg_dif_hi);

            current_thresh_max_dif_lo = select(use_boost_lo, current_thresh_max_dif_lo * angle_boost_factor, current_thresh_max_dif_lo);
            current_thresh_max_dif_hi = select(use_boost_hi, current_thresh_max_dif_hi * angle_boost_factor, current_thresh_max_dif_hi);

            current_thresh_mid_dif_lo = select(use_boost_lo, current_thresh_mid_dif_lo * angle_boost_factor, current_thresh_mid_dif_lo);
            current_thresh_mid_dif_hi = select(use_boost_hi, current_thresh_mid_dif_hi * angle_boost_factor, current_thresh_mid_dif_hi);
        }

        auto avg_refs_f_lo = (p1_f_lo + p2_f_lo + p3_f_lo + p4_f_lo) * 0.25f;
        auto avg_refs_f_hi = (p1_f_hi + p2_f_hi + p3_f_hi + p4_f_hi) * 0.25f;

        auto diff_avg_src_lo = avg_refs_f_lo - src_f_lo;
        auto diff_avg_src_hi = avg_refs_f_hi - src_f_hi;

        auto avg_dif_f_lo = abs(diff_avg_src_lo);
        auto avg_dif_f_hi = abs(diff_avg_src_hi);

        auto d1_lo = abs(p1_f_lo - src_f_lo);
        auto d1_hi = abs(p1_f_hi - src_f_hi);

        auto d2_lo = abs(p2_f_lo - src_f_lo);
        auto d2_hi = abs(p2_f_hi - src_f_hi);

        auto d3_lo = abs(p3_f_lo - src_f_lo);
        auto d3_hi = abs(p3_f_hi - src_f_hi);

        auto d4_lo = abs(p4_f_lo - src_f_lo);
        auto d4_hi = abs(p4_f_hi - src_f_hi);

        auto maxDif_lo = max(max(d1_lo, d2_lo), max(d3_lo, d4_lo));
        auto maxDif_hi = max(max(d1_hi, d2_hi), max(d3_hi, d4_hi));

        auto two_src_lo = src_f_lo * 2.0f;
        auto two_src_hi = src_f_hi * 2.0f;

        auto mid_dif_v_f_lo = abs((p1_f_lo + p2_f_lo) - two_src_lo);
        auto mid_dif_v_f_hi = abs((p1_f_hi + p2_f_hi) - two_src_hi);

        auto mid_dif_h_f_lo = abs((p3_f_lo + p4_f_lo) - two_src_lo);
        auto mid_dif_h_f_hi = abs((p3_f_hi + p4_f_hi) - two_src_hi);

        auto comp_avg_lo = saturate<V_float>(3.0f * (1.0f - avg_dif_f_lo / max(current_thresh_avg_dif_lo, 1e-5f)));
        auto comp_avg_hi = saturate<V_float>(3.0f * (1.0f - avg_dif_f_hi / max(current_thresh_avg_dif_hi, 1e-5f)));

        auto comp_max_lo = saturate<V_float>(3.0f * (1.0f - maxDif_lo / max(current_thresh_max_dif_lo, 1e-5f)));
        auto comp_max_hi = saturate<V_float>(3.0f * (1.0f - maxDif_hi / max(current_thresh_max_dif_hi, 1e-5f)));

        auto comp_mid_v_lo = saturate<V_float>(3.0f * (1.0f - mid_dif_v_f_lo / max(current_thresh_mid_dif_lo, 1e-5f)));
        auto comp_mid_v_hi = saturate<V_float>(3.0f * (1.0f - mid_dif_v_f_hi / max(current_thresh_mid_dif_hi, 1e-5f)));

        auto comp_mid_h_lo = saturate<V_float>(3.0f * (1.0f - mid_dif_h_f_lo / max(current_thresh_mid_dif_lo, 1e-5f)));
        auto comp_mid_h_hi = saturate<V_float>(3.0f * (1.0f - mid_dif_h_f_hi / max(current_thresh_mid_dif_hi, 1e-5f)));

        auto product_comps_lo = comp_avg_lo * comp_max_lo * comp_mid_v_lo * comp_mid_h_lo;
        auto product_comps_hi = comp_avg_hi * comp_max_hi * comp_mid_v_hi * comp_mid_h_hi;

        auto factor_lo = pow(product_comps_lo, 0.1f);
        auto factor_hi = pow(product_comps_hi, 0.1f);

        V_float blended_f_lo = src_f_lo + diff_avg_src_lo * factor_lo;
        V_float blended_f_hi = src_f_hi + diff_avg_src_hi * factor_hi;

        auto blended_i32_lo = truncatei(blended_f_lo + 0.5f);
        auto blended_i32_hi = truncatei(blended_f_hi + 0.5f);
        dst_pixels = compress(blended_i32_lo, blended_i32_hi);
    }

    auto sign_convert_vector = V_signed(static_cast<short>(0x8000));
    auto dst_signed = V_signed(dst_pixels) - sign_convert_vector;
    dst_signed = add_saturated(dst_signed, change);
    dst_pixels = V(dst_signed + sign_convert_vector);

    switch (dither_algo)
    {
    case DA_HIGH_NO_DITHERING:
    case DA_HIGH_ORDERED_DITHERING:
    case DA_HIGH_FLOYD_STEINBERG_DITHERING:
    {
        auto dst_lo = dither_high::dither<dither_algo>(dither_context, dst_pixels.get_low(), row, column);
        auto dst_hi = dither_high::dither<dither_algo>(dither_context, dst_pixels.get_high(), row, column + 8);
        dst_pixels = V(Vec8us(dst_lo), Vec8us(dst_hi));
    }
    break;
    default:
        break;
    }

    if (need_clamping)
        dst_pixels = high_bit_depth_pixels_clamp_avx2<V>(dst_pixels, clamp_high_add, clamp_high_sub, clamp_low);

    return dst_pixels;
}

template<PIXEL_MODE input_mode>
static unsigned short __forceinline read_pixel_avx2(const unsigned char* base, int offset)
{
    const unsigned char* ptr = base + offset;

    if constexpr (input_mode == LOW_BIT_DEPTH)
        return *ptr;
    else
        return *reinterpret_cast<const unsigned short*>(ptr);
}

template<typename V, int sample_mode, int dither_algo, PIXEL_MODE input_mode>
static void __forceinline read_reference_pixels_avx2(
    const process_plane_params& params, int shift, const unsigned char* src_px_start, const char* info_data_start,
    V& ref_pixels_1, V& ref_pixels_2, V& ref_pixels_3, V& ref_pixels_4)
{
    alignas(32)
        unsigned short tmp_1[V_ushort().size()];
    alignas(32)
        unsigned short tmp_2[V_ushort().size()];
    alignas(32)
        unsigned short tmp_3[V_ushort().size()];
    alignas(32)
        unsigned short tmp_4[V_ushort().size()];

    const int i_fix_step = (input_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1);

    const int* offsets_v1 = reinterpret_cast<const int*>(info_data_start);
    const int* offsets_h1 = offsets_v1 + V_int().size();
    const int* offsets_v2 = offsets_h1 + V_int().size();
    const int* offsets_h2 = offsets_v2 + V_int().size();

    int i_fix = 0;
    for (int i = 0; i < V_int().size(); ++i) {
        tmp_1[i] = read_pixel_avx2<input_mode>(src_px_start, i_fix + offsets_v1[i]);
        tmp_2[i] = read_pixel_avx2<input_mode>(src_px_start, i_fix - offsets_v1[i]);
        tmp_3[i] = read_pixel_avx2<input_mode>(src_px_start, i_fix + offsets_h1[i]);
        tmp_4[i] = read_pixel_avx2<input_mode>(src_px_start, i_fix - offsets_h1[i]);
        i_fix += i_fix_step;
    }

    for (int i = 0; i < V_int().size(); ++i) {
        tmp_1[i + V_int().size()] = read_pixel_avx2<input_mode>(src_px_start, i_fix + offsets_v2[i]);
        tmp_2[i + V_int().size()] = read_pixel_avx2<input_mode>(src_px_start, i_fix - offsets_v2[i]);
        tmp_3[i + V_int().size()] = read_pixel_avx2<input_mode>(src_px_start, i_fix + offsets_h2[i]);
        tmp_4[i + V_int().size()] = read_pixel_avx2<input_mode>(src_px_start, i_fix - offsets_h2[i]);
        i_fix += i_fix_step;
    }

    ref_pixels_1 = V().load(tmp_1) << shift;
    ref_pixels_2 = V().load(tmp_2) << shift;
    ref_pixels_3 = V().load(tmp_3) << shift;
    ref_pixels_4 = V().load(tmp_4) << shift;
}

std::mutex cache_mutex_avx2;
template<int sample_mode, bool blur_first, int dither_algo, bool aligned, PIXEL_MODE output_mode>
static void __cdecl _process_plane_avx2_impl(const process_plane_params& params, process_plane_context* context)
{
    auto src_pitch_vector = V_int(params.src_pitch);

    alignas(32)
        char context_buffer[DITHER_CONTEXT_BUFFER_SIZE];

    dither_high::init<dither_algo>(context_buffer, params.plane_width_in_pixels, params.output_depth);

    bool need_clamping = INTERNAL_BIT_DEPTH < 16 || params.pixel_min > 0 || params.pixel_max < 0xffff;
    auto clamp_high_add = V_ushort(0);
    auto clamp_high_sub = V_ushort(0);
    auto clamp_low = V_ushort(0);

    if (need_clamping) {
        clamp_low = V_ushort(static_cast<uint16_t>(params.pixel_min));
        clamp_high_add = V_ushort(static_cast<uint16_t>(0xFFFF)) - V_ushort(static_cast<uint16_t>(params.pixel_max));
        clamp_high_sub = clamp_high_add + clamp_low;
    }

    const int upsample_to_16_shift_bits = INTERNAL_BIT_DEPTH - params.input_depth;
    const int downshift_bits = INTERNAL_BIT_DEPTH - params.output_depth;
    const int pixel_step_shift_bits = (params.input_mode == HIGH_BIT_DEPTH_INTERLEAVED) ? 1 : 0;

    info_cache_avx2* cache = nullptr;
    char* info_data_stream = nullptr;
    bool use_cached_info = false;

    if (context->data) {
        cache = static_cast<info_cache_avx2*>(context->data);
        if (cache->pitch == params.src_pitch) {
            info_data_stream = cache->data_stream;
            use_cached_info = true;
        }
        cache = nullptr;
    }
    else {
        cache = static_cast<info_cache_avx2*>(malloc(sizeof(info_cache_avx2)));
        if (cache) {
            size_t cache_size = static_cast<size_t>((params.plane_width_in_pixels + 15) / 16) * params.plane_height_in_pixels * 128;
            info_data_stream = static_cast<char*>(_aligned_malloc(cache_size, 32));
            if (info_data_stream) {
                cache->data_stream = info_data_stream;
                cache->pitch = params.src_pitch;
            }
            else {
                free(cache); cache = nullptr;
            }
        }
    }

    const int info_cache_block_size = 128;
    const int current_input_mode = params.input_mode;

    for (int row = 0; row < params.plane_height_in_pixels; ++row) {
        const unsigned char* src_px_row_base = params.src_plane_ptr + static_cast<intptr_t>(params.src_pitch) * row;
        unsigned char* dst_px_row_base = params.dst_plane_ptr + static_cast<intptr_t>(params.dst_pitch) * row;
        pixel_dither_info* info_ptr_row_base = params.info_ptr_base + static_cast<intptr_t>(params.info_stride) * row;
        const short* grain_buffer_row_base = params.grain_buffer + static_cast<intptr_t>(params.grain_buffer_stride) * row;

        char* current_row_info_data_cache_ptr = use_cached_info ?
            (info_data_stream + static_cast<intptr_t>((params.plane_width_in_pixels + 15) / 16) * row * info_cache_block_size) : nullptr;
        char* current_row_info_data_build_ptr = (!use_cached_info && cache && info_data_stream) ?
            (info_data_stream + static_cast<intptr_t>((params.plane_width_in_pixels + 15) / 16) * row * info_cache_block_size) : nullptr;

        for (int col = 0; col < params.plane_width_in_pixels; col += V_ushort().size()) {
            const unsigned char* current_src_px = src_px_row_base + col * (current_input_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1);
            unsigned char* current_dst_px = dst_px_row_base + col * (output_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1);
            const short* current_grain_ptr = grain_buffer_row_base + col;
            pixel_dither_info* current_info_unit_ptr = info_ptr_row_base + col;

            char* data_stream_for_read_refs;
            alignas(32)
                char dummy_info_buffer[128];

            if (use_cached_info) {
                data_stream_for_read_refs = current_row_info_data_cache_ptr;
                current_row_info_data_cache_ptr += info_cache_block_size;
            }
            else {
                char* temp_info_build_ptr = current_row_info_data_build_ptr ? current_row_info_data_build_ptr : dummy_info_buffer;
                data_stream_for_read_refs = temp_info_build_ptr;
                process_plane_info_block_avx2_16px<V_int, sample_mode>(current_info_unit_ptr, src_pitch_vector, params.width_subsampling,
                    params.height_subsampling, pixel_step_shift_bits, temp_info_build_ptr);
                if (current_row_info_data_build_ptr) current_row_info_data_build_ptr += info_cache_block_size;
            }

            V_ushort ref_pixels_1 = zero_si256();
            V_ushort ref_pixels_2 = zero_si256();
            V_ushort ref_pixels_3 = zero_si256();
            V_ushort ref_pixels_4 = zero_si256();

            if (current_input_mode == LOW_BIT_DEPTH)
                read_reference_pixels_avx2<V_ushort, sample_mode, dither_algo, LOW_BIT_DEPTH>(params, upsample_to_16_shift_bits,
                    current_src_px, data_stream_for_read_refs, ref_pixels_1, ref_pixels_2, ref_pixels_3, ref_pixels_4);
            else
                read_reference_pixels_avx2<V_ushort, sample_mode, dither_algo, HIGH_BIT_DEPTH_INTERLEAVED>(params,
                    upsample_to_16_shift_bits, current_src_px, data_stream_for_read_refs, ref_pixels_1, ref_pixels_2, ref_pixels_3,
                    ref_pixels_4);

            auto src_pixels_data = (current_input_mode == LOW_BIT_DEPTH) ?
                (extend_low(V_uchar().load(current_src_px)) << upsample_to_16_shift_bits) :
                (V_ushort().load(current_src_px) << upsample_to_16_shift_bits);

            auto change = V_short().load(current_grain_ptr);

            auto dst_pixels_data = process_pixels_avx2<V_ushort, V_short, sample_mode, blur_first, dither_algo>(
                src_pixels_data, change, ref_pixels_1, ref_pixels_2, ref_pixels_3, ref_pixels_4, clamp_high_add, clamp_high_sub, clamp_low,
                need_clamping, row, col, context_buffer, info_ptr_row_base + col, params, upsample_to_16_shift_bits);

            if (output_mode == LOW_BIT_DEPTH) {
                auto p = dst_pixels_data >> downshift_bits;
                auto p_8bit = compress_saturated(p.get_low(), p.get_high());
                p_8bit.store(current_dst_px);
            }
            else {
                auto p = dst_pixels_data >> downshift_bits;
                p.store(current_dst_px);
            }
        }

        dither_high::next_row<dither_algo>(context_buffer);
    }

    dither_high::complete<dither_algo>(context_buffer);

    if (!use_cached_info && !context->data && cache && info_data_stream) {
        std::lock_guard<std::mutex> lock(cache_mutex_avx2);
        if (context->data) {
            destroy_cache_avx2(cache);
        }
        else {
            context->data = cache;
            context->destroy = destroy_cache_avx2;
        }
    }
    else if (cache && (!info_data_stream || context->data)) {
        if (info_data_stream) _aligned_free(info_data_stream);
        free(cache);
    }
}

template<int sample_mode, bool blur_first, int dither_algo>
static void process_plane_avx2_impl(const process_plane_params& params, process_plane_context* context)
{
    switch (params.output_mode)
    {
    case LOW_BIT_DEPTH:
        _process_plane_avx2_impl<sample_mode, blur_first, dither_algo, true, LOW_BIT_DEPTH>(params, context);
        break;
    case HIGH_BIT_DEPTH_INTERLEAVED:
        _process_plane_avx2_impl<sample_mode, blur_first, dither_algo, true, HIGH_BIT_DEPTH_INTERLEAVED>(params, context);
        break;
    default:
        abort();
    }
}
