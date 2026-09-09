#pragma once

#include "dither_high.h"
#include "impl_dispatch.h"
#include "sse_utils.h"

#include "VCL2/vectorclass.h"
#include "VCL2/vectormath_exp.h"
#include "VCL2/vectormath_trig.h"

#if INSTRSET >= 10 // AVX512VL
#define DEBAND_NAMESPACE ns_avx512
#elif INSTRSET >= 8 // AVX2
#define DEBAND_NAMESPACE ns_avx2
#else
#error "INSTRSET must be 8 (AVX2) or 10 (AVX512)"
#endif

namespace DEBAND_NAMESPACE {

#if INSTRSET >= 10 // AVX512VL
    using V_int = Vec16i;
    inline __m512i get_zero_int() noexcept
    {
        return zero_si512();
    }

    static constexpr int simd_align = 64;
#elif INSTRSET >= 8 // AVX2
    using V_int = Vec8i;
    inline __m256i get_zero_int() noexcept
    {
        return zero_si256();
    }

    static constexpr int simd_align = 32;
#endif

    using V_float = std::conditional_t<std::is_same_v<V_int, Vec8i>, Vec8f, Vec16f>;
    using V_fbool = std::conditional_t<std::is_same_v<V_int, Vec8i>, Vec8fb, Vec16fb>;
    using V_ushort = std::conditional_t<std::is_same_v<V_int, Vec8i>, Vec16us, Vec32us>;
    using V_short = std::conditional_t<std::is_same_v<V_int, Vec8i>, Vec16s, Vec32s>;
    using V_sbool = std::conditional_t<std::is_same_v<V_int, Vec8i>, Vec16sb, Vec32sb>;
    using V_uchar = std::conditional_t<std::is_same_v<V_int, Vec8i>, Vec32uc, Vec64uc>;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

    struct MXCSR_guard
    {
        int old_cw;
        MXCSR_guard() noexcept
            : old_cw(get_control_word())
        {
            no_subnormals();
        }
        ~MXCSR_guard() noexcept
        {
            set_control_word(old_cw);
        }
    };

    typedef struct _info_cache_avx2_avx512
    {
        int pitch;
        char* data_stream;
    } info_cache_avx2_avx512;

    static void destroy_cache_avx2_avx512(void* data)
    {
        assert(data);
        auto* cache = reinterpret_cast<info_cache_avx2_avx512*>(data);
        _aligned_free(cache->data_stream);
        free(data);
    }

    template <typename V>
    static auto __forceinline high_bit_depth_pixels_clamp_avx2_avx512(V pixels, V high_add, V high_sub, const V& low)
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

    template <typename V_float, typename V_int>
    static __forceinline V_float fast_pow01(V_float x)
    {
        auto is_zero = (x == 0.0f);
        x = select(is_zero, 1.0f, x);

        V_int xi = reinterpret_i(x);
        V_float xi_f = to_float(xi);
        V_float out_f = xi_f * 0.1f + 958817894.0f;
        V_int out_i = truncatei(out_f);
        V_float result = reinterpret_f(out_i);

        return select(is_zero, 0.0f, result);
    }

    template <typename V, int sample_mode>
    static __forceinline void process_plane_info_block_avx2_avx512_px(pixel_dither_info*& info_ptr, const V& src_pitch_vector,
        const int width_subsample, const int height_subsample, const int pixel_step_shift_bits, char*& info_data_stream)
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
    static V_bool __forceinline generate_blend_mask_high_avx2_avx512(V a, V threshold)
    {
        return V_bool(a < threshold);
    }

    template <typename V, PIXEL_MODE input_mode>
    static __forceinline V gather_pixel_values_avx2_avx512(const process_plane_params& params, V const y_coords, V const x_coords,
        int upsample_shift)
    {
        auto clamped_y = max(V(0), min(y_coords, params.plane_height_in_pixels - 1));
        auto clamped_x = max(V(0), min(x_coords, params.plane_width_in_pixels - 1));

        V pitch(params.src_pitch);
        V pixel_bytes_v((input_mode == HIGH_BIT_DEPTH_INTERLEAVED) ? 2 : 1);
        auto byte_offsets = clamped_y * pitch + clamped_x * pixel_bytes_v;

        const unsigned char* base_ptr = params.src_plane_ptr;

        if constexpr (std::is_same_v<V, Vec8i>) {
            V offsets = byte_offsets;
            V gathered = _mm256_i32gather_epi32(reinterpret_cast<const int*>(base_ptr), offsets, 1);

            if constexpr (input_mode == LOW_BIT_DEPTH) {
                V pixels_v(gathered & V(0x000000FF));
                pixels_v <<= upsample_shift;
                return pixels_v;
            }
            else {
                V pixels_v(gathered & V(0x0000FFFF));
                pixels_v <<= upsample_shift;
                return pixels_v;
            }
        }
        else {
            V offsets = byte_offsets;
            V gathered = _mm512_i32gather_epi32(offsets, base_ptr, 1);

            if constexpr (input_mode == LOW_BIT_DEPTH) {
                V pixels_v(gathered & V(0x000000FF));
                pixels_v <<= upsample_shift;
                return pixels_v;
            }
            else {
                V pixels_v(gathered & V(0x0000FFFF));
                pixels_v <<= upsample_shift;
                return pixels_v;
            }
        }
    }

    template <typename V, PIXEL_MODE input_mode>
    static __forceinline void calculate_gradient_vector_avx2_avx512(const process_plane_params& params, const V& y_coords, const V& x_coords,
        int read_distance, int upsample_shift, V_float& out_gx, V_float& out_gy)
    {
        V rd(read_distance);
        auto p00 = gather_pixel_values_avx2_avx512<V, input_mode>(params, y_coords - rd, x_coords - rd, upsample_shift);
        auto p10 = gather_pixel_values_avx2_avx512<V, input_mode>(params, y_coords - rd, x_coords, upsample_shift);
        auto p20 = gather_pixel_values_avx2_avx512<V, input_mode>(params, y_coords - rd, x_coords + rd, upsample_shift);
        auto p01 = gather_pixel_values_avx2_avx512<V, input_mode>(params, y_coords, x_coords - rd, upsample_shift);
        auto p21 = gather_pixel_values_avx2_avx512<V, input_mode>(params, y_coords, x_coords + rd, upsample_shift);
        auto p02 = gather_pixel_values_avx2_avx512<V, input_mode>(params, y_coords + rd, x_coords - rd, upsample_shift);
        auto p12 = gather_pixel_values_avx2_avx512<V, input_mode>(params, y_coords + rd, x_coords, upsample_shift);
        auto p22 = gather_pixel_values_avx2_avx512<V, input_mode>(params, y_coords + rd, x_coords + rd, upsample_shift);

        auto gx = (p20 + (p21 << 1) + p22) - (p00 + (p01 << 1) + p02);
        auto gy = (p00 + (p10 << 1) + p20) - (p02 + (p12 << 1) + p22);

        out_gx = to_float(gx);
        out_gy = to_float(gy);
    }

    struct v_thresh_params
    {
        const V_float v_tan_thresh;

        const V_float v_inv_thresh_base_avg;
        const V_float v_inv_thresh_base_max;
        const V_float v_inv_thresh_base_mid;

        const V_float v_inv_thresh_boosted_avg;
        const V_float v_inv_thresh_boosted_max;
        const V_float v_inv_thresh_boosted_mid;
    };

    template<typename V, typename V_signed, int sample_mode, bool blur_first, int dither_algo, PIXEL_MODE input_mode>
    static auto __forceinline process_pixels_avx2_avx512(V src_pixels, V_signed change, const V& ref_pixels_1, const V& ref_pixels_2,
        const V& ref_pixels_3, const V& ref_pixels_4, const V& clamp_high_add, const V& clamp_high_sub, const V& clamp_low,
        bool need_clamping, int row, int column, void* dither_context, const pixel_dither_info* pdi_ptr, const process_plane_params& params,
        int upsample_to_16_shift_bits, const v_thresh_params& thresh_params)
    {
        const int threshold = params.threshold;
        const int threshold1 = params.threshold1;
        const int threshold2 = params.threshold2;

        V dst_pixels = get_zero_int();

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

            auto avg = compress_saturated_s2u(sum_lo >> 2, sum_hi >> 2);

            auto avgDif = max(sub_saturated(avg, src_pixels), sub_saturated(src_pixels, avg));
            auto maxDif = max(
                max(max(sub_saturated(ref_pixels_1, src_pixels), sub_saturated(src_pixels, ref_pixels_1)),
                    max(sub_saturated(ref_pixels_2, src_pixels), sub_saturated(src_pixels, ref_pixels_2))),
                max(max(sub_saturated(ref_pixels_3, src_pixels), sub_saturated(src_pixels, ref_pixels_3)),
                    max(sub_saturated(ref_pixels_4, src_pixels), sub_saturated(src_pixels, ref_pixels_4)))
            );

            auto src_lo = extend_low(src_pixels);
            auto src_hi = extend_high(src_pixels);
            auto two_src_lo = src_lo << 1;
            auto two_src_hi = src_hi << 1;

            auto midDif1 = compress_saturated_s2u(abs((r1_lo + r2_lo) - two_src_lo), abs((r1_hi + r2_hi) - two_src_hi));
            auto midDif2 = compress_saturated_s2u(abs((r3_lo + r4_lo) - two_src_lo), abs((r3_hi + r4_hi) - two_src_hi));

            auto use_orig_pixel_blend_mask = generate_blend_mask_high_avx2_avx512<V_ushort, V_sbool>(avgDif, V_ushort(threshold))
                & generate_blend_mask_high_avx2_avx512<V_ushort, V_sbool>(maxDif, V_ushort(threshold1))
                & generate_blend_mask_high_avx2_avx512<V_ushort, V_sbool>(midDif1, V_ushort(threshold2))
                & generate_blend_mask_high_avx2_avx512<V_ushort, V_sbool>(midDif2, V_ushort(threshold2));

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

            V_float inv_thresh_avg_lo = thresh_params.v_inv_thresh_base_avg;
            V_float inv_thresh_avg_hi = thresh_params.v_inv_thresh_base_avg;

            V_float inv_thresh_max_lo = thresh_params.v_inv_thresh_base_max;
            V_float inv_thresh_max_hi = thresh_params.v_inv_thresh_base_max;

            V_float inv_thresh_mid_lo = thresh_params.v_inv_thresh_base_mid;
            V_float inv_thresh_mid_hi = thresh_params.v_inv_thresh_base_mid;

            if constexpr (sample_mode == 7) {
                constexpr int grad_read_distance = 20;
                const V_float v_tan_thresh = thresh_params.v_tan_thresh;

#if INSTRSET >= 10 // AVX512VL
                auto base_x_coords_lo = V_int(column) + V_int(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
                auto base_x_coords_hi = V_int(column + V_int().size()) + V_int(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
#elif INSTRSET >= 8 // AVX2
                auto base_x_coords_lo = V_int(column) + V_int(0, 1, 2, 3, 4, 5, 6, 7);
                auto base_x_coords_hi = V_int(column + V_int().size()) + V_int(0, 1, 2, 3, 4, 5, 6, 7);
#endif
                V_int base_y_coords(row);

                V_float gx_org_lo;
                V_float gy_org_lo;
                calculate_gradient_vector_avx2_avx512<V_int, input_mode>(params, base_y_coords, base_x_coords_lo, grad_read_distance,
                    upsample_to_16_shift_bits, gx_org_lo, gy_org_lo);

                V_float gx_org_hi;
                V_float gy_org_hi;
                calculate_gradient_vector_avx2_avx512<V_int, input_mode>(params, base_y_coords, base_x_coords_hi, grad_read_distance,
                    upsample_to_16_shift_bits, gx_org_hi, gy_org_hi);

                alignas(simd_align)
                    int32_t ref1_buffer[V_ushort().size()];
                for (int k = 0; k < V_ushort().size(); ++k)
                    ref1_buffer[k] = pdi_ptr[k].ref1;

                auto ref1_offsets_lo(V_int().load_a(ref1_buffer));
                auto ref1_offsets_hi(V_int().load_a(ref1_buffer + V_int().size()));

                auto y_offsets_h_lo = ref1_offsets_lo >> params.height_subsampling;
                auto y_offsets_h_hi = ref1_offsets_hi >> params.height_subsampling;
                auto x_offsets_w_lo = ref1_offsets_lo >> params.width_subsampling;
                auto x_offsets_w_hi = ref1_offsets_hi >> params.width_subsampling;

                auto check_aligned = [&](const V_float& gx1, const V_float& gy1, const V_float mag_sq1, const V_float& gx2,
                    const V_float& gy2) {
                    const auto cross = abs(gx1 * gy2 - gy1 * gx2);
                    const auto dot = abs(gx1 * gx2 + gy1 * gy2);

                    const auto mag_sq2 = gx2 * gx2 + gy2 * gy2;

                    constexpr float flat_epsilon_sq = 1.0f;

                    const auto both_flat = (mag_sq1 < flat_epsilon_sq) & (mag_sq2 < flat_epsilon_sq);
                    const auto both_active = (mag_sq1 >= flat_epsilon_sq) & (mag_sq2 >= flat_epsilon_sq);

                    return both_flat | (both_active & (cross <= v_tan_thresh * dot));
                    };

                V_float gx_ref;
                V_float gy_ref;

                auto mag_sq1 = gx_org_lo * gx_org_lo + gy_org_lo * gy_org_lo;
                calculate_gradient_vector_avx2_avx512<V_int, input_mode>(params, base_y_coords + y_offsets_h_lo, base_x_coords_lo,
                    grad_read_distance, upsample_to_16_shift_bits, gx_ref, gy_ref);
                auto use_boost_lo = check_aligned(gx_org_lo, gy_org_lo, mag_sq1, gx_ref, gy_ref);

                calculate_gradient_vector_avx2_avx512<V_int, input_mode>(params, base_y_coords - y_offsets_h_lo, base_x_coords_lo,
                    grad_read_distance, upsample_to_16_shift_bits, gx_ref, gy_ref);
                use_boost_lo &= check_aligned(gx_org_lo, gy_org_lo, mag_sq1, gx_ref, gy_ref);

                calculate_gradient_vector_avx2_avx512<V_int, input_mode>(params, base_y_coords, base_x_coords_lo + x_offsets_w_lo,
                    grad_read_distance, upsample_to_16_shift_bits, gx_ref, gy_ref);
                use_boost_lo &= check_aligned(gx_org_lo, gy_org_lo, mag_sq1,  gx_ref, gy_ref);

                calculate_gradient_vector_avx2_avx512<V_int, input_mode>(params, base_y_coords, base_x_coords_lo - x_offsets_w_lo,
                    grad_read_distance, upsample_to_16_shift_bits, gx_ref, gy_ref);
                use_boost_lo &= check_aligned(gx_org_lo, gy_org_lo, mag_sq1, gx_ref, gy_ref);

                mag_sq1 = gx_org_hi * gx_org_hi + gy_org_hi * gy_org_hi;
                calculate_gradient_vector_avx2_avx512<V_int, input_mode>(params, base_y_coords + y_offsets_h_hi, base_x_coords_hi,
                    grad_read_distance, upsample_to_16_shift_bits, gx_ref, gy_ref);
                auto use_boost_hi = check_aligned(gx_org_hi, gy_org_hi, mag_sq1, gx_ref, gy_ref);

                calculate_gradient_vector_avx2_avx512<V_int, input_mode>(params, base_y_coords - y_offsets_h_hi, base_x_coords_hi,
                    grad_read_distance, upsample_to_16_shift_bits, gx_ref, gy_ref);
                use_boost_hi &= check_aligned(gx_org_hi, gy_org_hi, mag_sq1, gx_ref, gy_ref);

                calculate_gradient_vector_avx2_avx512<V_int, input_mode>(params, base_y_coords, base_x_coords_hi + x_offsets_w_hi,
                    grad_read_distance, upsample_to_16_shift_bits, gx_ref, gy_ref);
                use_boost_hi &= check_aligned(gx_org_hi, gy_org_hi, mag_sq1, gx_ref, gy_ref);

                calculate_gradient_vector_avx2_avx512<V_int, input_mode>(params, base_y_coords, base_x_coords_hi - x_offsets_w_hi,
                    grad_read_distance, upsample_to_16_shift_bits, gx_ref, gy_ref);
                use_boost_hi &= check_aligned(gx_org_hi, gy_org_hi, mag_sq1, gx_ref, gy_ref);

                inv_thresh_avg_lo = select(use_boost_lo, thresh_params.v_inv_thresh_boosted_avg, thresh_params.v_inv_thresh_base_avg);
                inv_thresh_avg_hi = select(use_boost_hi, thresh_params.v_inv_thresh_boosted_avg, thresh_params.v_inv_thresh_base_avg);

                inv_thresh_max_lo = select(use_boost_lo, thresh_params.v_inv_thresh_boosted_max, thresh_params.v_inv_thresh_base_max);
                inv_thresh_max_hi = select(use_boost_hi, thresh_params.v_inv_thresh_boosted_max, thresh_params.v_inv_thresh_base_max);

                inv_thresh_mid_lo = select(use_boost_lo, thresh_params.v_inv_thresh_boosted_mid, thresh_params.v_inv_thresh_base_mid);
                inv_thresh_mid_hi = select(use_boost_hi, thresh_params.v_inv_thresh_boosted_mid, thresh_params.v_inv_thresh_base_mid);
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

            auto comp_avg_lo = saturate<V_float>(3.0f * (1.0f - avg_dif_f_lo * inv_thresh_avg_lo));
            auto comp_avg_hi = saturate<V_float>(3.0f * (1.0f - avg_dif_f_hi * inv_thresh_avg_hi));

            auto comp_max_lo = saturate<V_float>(3.0f * (1.0f - maxDif_lo * inv_thresh_max_lo));
            auto comp_max_hi = saturate<V_float>(3.0f * (1.0f - maxDif_hi * inv_thresh_max_hi));

            auto comp_mid_v_lo = saturate<V_float>(3.0f * (1.0f - mid_dif_v_f_lo * inv_thresh_mid_lo));
            auto comp_mid_v_hi = saturate<V_float>(3.0f * (1.0f - mid_dif_v_f_hi * inv_thresh_mid_hi));

            auto comp_mid_h_lo = saturate<V_float>(3.0f * (1.0f - mid_dif_h_f_lo * inv_thresh_mid_lo));
            auto comp_mid_h_hi = saturate<V_float>(3.0f * (1.0f - mid_dif_h_f_hi * inv_thresh_mid_hi));

            auto product_comps_lo = comp_avg_lo * comp_max_lo * comp_mid_v_lo * comp_mid_h_lo;
            auto product_comps_hi = comp_avg_hi * comp_max_hi * comp_mid_v_hi * comp_mid_h_hi;

            auto factor_lo = fast_pow01<V_float, V_int>(product_comps_lo);
            auto factor_hi = fast_pow01<V_float, V_int>(product_comps_hi);

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
#if INSTRSET >= 10 // AVX512VL
            auto dst_1 = dither_high::dither<dither_algo>(dither_context, dst_pixels.get_low().get_low(), row, column);
            auto dst_2 = dither_high::dither<dither_algo>(dither_context, dst_pixels.get_low().get_high(), row, column + 8);
            auto dst_3 = dither_high::dither<dither_algo>(dither_context, dst_pixels.get_high().get_low(), row, column + 16);
            auto dst_4 = dither_high::dither<dither_algo>(dither_context, dst_pixels.get_high().get_high(), row, column + 24);
            dst_pixels = V(Vec16us(Vec8us(dst_1), Vec8us(dst_2)), Vec16us(Vec8us(dst_3), Vec8us(dst_4)));
#elif INSTRSET >= 8 // AVX2
            auto dst_lo = dither_high::dither<dither_algo>(dither_context, dst_pixels.get_low(), row, column);
            auto dst_hi = dither_high::dither<dither_algo>(dither_context, dst_pixels.get_high(), row, column + 8);
            dst_pixels = V(Vec8us(dst_lo), Vec8us(dst_hi));
#endif
        }
        break;
        default:
            break;
        }

        if (need_clamping)
            dst_pixels = high_bit_depth_pixels_clamp_avx2_avx512<V>(dst_pixels, clamp_high_add, clamp_high_sub, clamp_low);

        return dst_pixels;
    }

    template<PIXEL_MODE input_mode, typename V = V_int>
    inline auto gather_us_vec(const unsigned char* base, V offsets1, V offsets2) {
        if constexpr (std::is_same_v<V, Vec8i>) {
            V g1 = _mm256_i32gather_epi32(reinterpret_cast<const int*>(base), offsets1, 1);
            V g2 = _mm256_i32gather_epi32(reinterpret_cast<const int*>(base), offsets2, 1);
            if constexpr (input_mode == LOW_BIT_DEPTH) {
                g1 = g1 & V(0x000000FF);
                g2 = g2 & V(0x000000FF);
            }
            else {
                g1 = g1 & V(0x0000FFFF);
                g2 = g2 & V(0x0000FFFF);
            }

            return compress_saturated_s2u(g1, g2);
        }
        else {
            V g1 = _mm512_i32gather_epi32(offsets1, base, 1);
            V g2 = _mm512_i32gather_epi32(offsets2, base, 1);

            if constexpr (input_mode == LOW_BIT_DEPTH)
            {
                g1 = g1 & V(0x000000FF);
                g2 = g2 & V(0x000000FF);
            }
            else
            {
                g1 = g1 & V(0x0000FFFF);
                g2 = g2 & V(0x0000FFFF);
            }

            auto a0 = g1.get_low();
            auto a1 = g1.get_high();

            auto b0 = g2.get_low();
            auto b1 = g2.get_high();

            auto pack_a = compress_saturated_s2u(a0, a1);
            auto pack_b = compress_saturated_s2u(b0, b1);

            return Vec32us(pack_a, pack_b);
        }
    }

    template<typename V, int sample_mode, int dither_algo, PIXEL_MODE input_mode>
    static void __forceinline read_reference_pixels_avx2_avx512(
        const process_plane_params& params, int shift, const unsigned char* src_px_start, const char* info_data_start,
        V& ref_pixels_1, V& ref_pixels_2, V& ref_pixels_3, V& ref_pixels_4)
    {
        const int i_fix_step = (input_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1);

        const int* offsets_v1 = reinterpret_cast<const int*>(info_data_start);
        const int* offsets_h1 = offsets_v1 + V_int().size();
        const int* offsets_v2 = offsets_h1 + V_int().size();
        const int* offsets_h2 = offsets_v2 + V_int().size();

#if INSTRSET >= 10 // AVX512VL
        V_int i_fix_vec1 = V_int(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15) * i_fix_step;
        V_int i_fix_vec2 = i_fix_vec1 + V_int(16 * i_fix_step);
#elif INSTRSET >= 8 // AVX2
        V_int i_fix_vec1 = V_int(0, 1, 2, 3, 4, 5, 6, 7) * i_fix_step;
        V_int i_fix_vec2 = i_fix_vec1 + V_int(8 * i_fix_step);
#endif

        auto v_plus_1 = i_fix_vec1 + V_int().load(offsets_v1);
        auto v_minus_1 = i_fix_vec1 - V_int().load(offsets_v1);
        auto h_plus_1 = i_fix_vec1 + V_int().load(offsets_h1);
        auto h_minus_1 = i_fix_vec1 - V_int().load(offsets_h1);

        auto v_plus_2 = i_fix_vec2 + V_int().load(offsets_v2);
        auto v_minus_2 = i_fix_vec2 - V_int().load(offsets_v2);
        auto h_plus_2 = i_fix_vec2 + V_int().load(offsets_h2);
        auto h_minus_2 = i_fix_vec2 - V_int().load(offsets_h2);

        ref_pixels_1 = gather_us_vec<input_mode>(src_px_start, v_plus_1, v_plus_2) << shift;
        ref_pixels_2 = gather_us_vec<input_mode>(src_px_start, v_minus_1, v_minus_2) << shift;
        ref_pixels_3 = gather_us_vec<input_mode>(src_px_start, h_plus_1, h_plus_2) << shift;
        ref_pixels_4 = gather_us_vec<input_mode>(src_px_start, h_minus_1, h_minus_2) << shift;
    }

    std::mutex cache_mutex_avx2_avx512;
    template<int sample_mode, bool blur_first, int dither_algo, bool aligned, PIXEL_MODE output_mode>
    static void __cdecl _process_plane_avx2_avx512_impl(const process_plane_params& params, process_plane_context* context)
    {
        MXCSR_guard guard;

        auto src_pitch_vector = V_int(params.src_pitch);

        alignas(simd_align)
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

        info_cache_avx2_avx512* cache = nullptr;
        char* info_data_stream = nullptr;
        bool use_cached_info = false;
        static constexpr int info_cache_block_size = sizeof(V_int) * 4;
        const size_t blocks_per_row = (params.plane_width_in_pixels + V_ushort().size() - 1) / V_ushort().size();

        if (context->data) {
            cache = static_cast<info_cache_avx2_avx512*>(context->data);
            if (cache->pitch == params.src_pitch) {
                info_data_stream = cache->data_stream;
                use_cached_info = true;
            }
            cache = nullptr;
        }
        else {
            cache = static_cast<info_cache_avx2_avx512*>(malloc(sizeof(info_cache_avx2_avx512)));
            if (cache) {
                size_t cache_size = blocks_per_row * params.plane_height_in_pixels * info_cache_block_size;
                info_data_stream = static_cast<char*>(_aligned_malloc(cache_size, simd_align));
                if (info_data_stream) {
                    cache->data_stream = info_data_stream;
                    cache->pitch = params.src_pitch;
                }
                else {
                    free(cache); cache = nullptr;
                }
            }
        }

        const int current_input_mode = params.input_mode;

        float tan_thresh = 0.0f;
        float angle_boost_factor = 1.0f;
        if constexpr (sample_mode == 7) {
            const float max_angle_rad = params.max_angle * static_cast<float>(M_PI);
            constexpr float half_pi = 1.57079632679f;
            tan_thresh = (max_angle_rad >= half_pi) ? 1e6f : std::tan(max_angle_rad);

            angle_boost_factor = params.angle_boost;
        }

        float inv_thresh_base_avg = 0.0f;
        float inv_thresh_base_max = 0.0f;
        float inv_thresh_base_mid = 0.0f;

        float inv_thresh_boosted_avg = 0.0f;
        float inv_thresh_boosted_max = 0.0f;
        float inv_thresh_boosted_mid = 0.0f;

        if constexpr (sample_mode >= 6) {
            const float base_thresh_avg = static_cast<float>(params.threshold);
            const float base_thresh_max = static_cast<float>(params.threshold1);
            const float base_thresh_mid = static_cast<float>(params.threshold2);

            const float boosted_thresh_avg = base_thresh_avg * angle_boost_factor;
            const float boosted_thresh_max = base_thresh_max * angle_boost_factor;
            const float boosted_thresh_mid = base_thresh_mid * angle_boost_factor;

            inv_thresh_base_avg = 1.0f / std::max(base_thresh_avg, 1e-5f);
            inv_thresh_base_max = 1.0f / std::max(base_thresh_max, 1e-5f);
            inv_thresh_base_mid = 1.0f / std::max(base_thresh_mid, 1e-5f);

            inv_thresh_boosted_avg = 1.0f / std::max(boosted_thresh_avg, 1e-5f);
            inv_thresh_boosted_max = 1.0f / std::max(boosted_thresh_max, 1e-5f);
            inv_thresh_boosted_mid = 1.0f / std::max(boosted_thresh_mid, 1e-5f);
        }

        v_thresh_params thresh_params = {
            V_float(tan_thresh),
            V_float(inv_thresh_base_avg),
            V_float(inv_thresh_base_max),
            V_float(inv_thresh_base_mid),
            V_float(inv_thresh_boosted_avg),
            V_float(inv_thresh_boosted_max),
            V_float(inv_thresh_boosted_mid),
        };

        for (int row = 0; row < params.plane_height_in_pixels; ++row) {
            const unsigned char* src_px_row_base = params.src_plane_ptr + static_cast<intptr_t>(params.src_pitch) * row;
            unsigned char* dst_px_row_base = params.dst_plane_ptr + static_cast<intptr_t>(params.dst_pitch) * row;
            pixel_dither_info* info_ptr_row_base = params.info_ptr_base + static_cast<intptr_t>(params.info_stride) * row;
            const short* grain_buffer_row_base = params.grain_buffer + static_cast<intptr_t>(params.grain_buffer_stride) * row;

            char* current_row_info_data_cache_ptr = use_cached_info ?
                (info_data_stream + blocks_per_row * row * info_cache_block_size) : nullptr;
            char* current_row_info_data_build_ptr = (!use_cached_info && cache && info_data_stream) ?
                (info_data_stream + blocks_per_row * row * info_cache_block_size) : nullptr;

            for (int col = 0; col < params.plane_width_in_pixels; col += V_ushort().size()) {
                const unsigned char* current_src_px = src_px_row_base + col * (current_input_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1);
                unsigned char* current_dst_px = dst_px_row_base + col * (output_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1);
                const short* current_grain_ptr = grain_buffer_row_base + col;
                pixel_dither_info* current_info_unit_ptr = info_ptr_row_base + col;

                char* data_stream_for_read_refs;
                alignas(simd_align)
                    char dummy_info_buffer[sizeof(V_int) * 4];

                if (use_cached_info) {
                    data_stream_for_read_refs = current_row_info_data_cache_ptr;
                    current_row_info_data_cache_ptr += info_cache_block_size;
                }
                else {
                    char* temp_info_build_ptr = current_row_info_data_build_ptr ? current_row_info_data_build_ptr : dummy_info_buffer;
                    data_stream_for_read_refs = temp_info_build_ptr;
                    process_plane_info_block_avx2_avx512_px<V_int, sample_mode>(current_info_unit_ptr, src_pitch_vector,
                        params.width_subsampling, params.height_subsampling, pixel_step_shift_bits, temp_info_build_ptr);
                    if (current_row_info_data_build_ptr)
                        current_row_info_data_build_ptr += info_cache_block_size;
                }

                V_ushort ref_pixels_1 = get_zero_int();
                V_ushort ref_pixels_2 = get_zero_int();
                V_ushort ref_pixels_3 = get_zero_int();
                V_ushort ref_pixels_4 = get_zero_int();

                if (current_input_mode == LOW_BIT_DEPTH)
                    read_reference_pixels_avx2_avx512<V_ushort, sample_mode, dither_algo, LOW_BIT_DEPTH>(params, upsample_to_16_shift_bits,
                        current_src_px, data_stream_for_read_refs, ref_pixels_1, ref_pixels_2, ref_pixels_3, ref_pixels_4);
                else
                    read_reference_pixels_avx2_avx512<V_ushort, sample_mode, dither_algo, HIGH_BIT_DEPTH_INTERLEAVED>(params,
                        upsample_to_16_shift_bits, current_src_px, data_stream_for_read_refs, ref_pixels_1, ref_pixels_2, ref_pixels_3,
                        ref_pixels_4);

                auto src_pixels_data = (current_input_mode == LOW_BIT_DEPTH) ?
                    (extend_low(V_uchar().load(current_src_px)) << upsample_to_16_shift_bits) :
                    (V_ushort().load(current_src_px) << upsample_to_16_shift_bits);

                auto change = V_short().load(current_grain_ptr);

                V_ushort dst_pixels_data;
                if (current_input_mode == LOW_BIT_DEPTH) {
                    dst_pixels_data = process_pixels_avx2_avx512<V_ushort, V_short, sample_mode, blur_first, dither_algo, LOW_BIT_DEPTH>(
                        src_pixels_data, change, ref_pixels_1, ref_pixels_2, ref_pixels_3, ref_pixels_4, clamp_high_add, clamp_high_sub, clamp_low,
                        need_clamping, row, col, context_buffer, info_ptr_row_base + col, params, upsample_to_16_shift_bits, thresh_params);
                }
                else {
                    dst_pixels_data = process_pixels_avx2_avx512<V_ushort, V_short, sample_mode, blur_first, dither_algo, HIGH_BIT_DEPTH_INTERLEAVED>(
                        src_pixels_data, change, ref_pixels_1, ref_pixels_2, ref_pixels_3, ref_pixels_4, clamp_high_add, clamp_high_sub, clamp_low,
                        need_clamping, row, col, context_buffer, info_ptr_row_base + col, params, upsample_to_16_shift_bits, thresh_params);
                }

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
            std::lock_guard<std::mutex> lock(cache_mutex_avx2_avx512);
            if (context->data) {
                destroy_cache_avx2_avx512(cache);
            }
            else {
                context->data = cache;
                context->destroy = destroy_cache_avx2_avx512;
            }
        }
        else if (cache && (!info_data_stream || context->data)) {
            if (info_data_stream)
                _aligned_free(info_data_stream);

            free(cache);
        }
    }

    template<int sample_mode, bool blur_first, int dither_algo>
    void process_plane_impl(const process_plane_params& params, process_plane_context* context)
    {
        switch (params.output_mode)
        {
        case LOW_BIT_DEPTH:
            _process_plane_avx2_avx512_impl<sample_mode, blur_first, dither_algo, true, LOW_BIT_DEPTH>(params, context);
            break;
        case HIGH_BIT_DEPTH_INTERLEAVED:
            _process_plane_avx2_avx512_impl<sample_mode, blur_first, dither_algo, true, HIGH_BIT_DEPTH_INTERLEAVED>(params, context);
            break;
        default:
            abort();
        }
    }
} // namespace DEBAND_NAMESPACE
