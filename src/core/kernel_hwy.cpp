#include "core/constants.hpp"
#include "core/dither_ordered.hpp"
#include "core/floyd_steinberg_dither.hpp"
#include "core/kernel.hpp"
#include "core/pixel_proc_common.hpp"
#include "core/sample_modes.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "core/kernel_hwy.cpp"

#include "hwy/foreach_target.h"
#include "hwy/highway.h"
#include "hwy/contrib/math/math-inl.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <type_traits>

HWY_BEFORE_NAMESPACE();
namespace neo_f3kdb::HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

#include "core/kernel_hwy-inl.hpp"

template <int kSampleMode, bool kBlurFirst, int kDitherAlgo, bool kUseCachedOffsets, bool kPureDeband, class PixelIn, class PixelOut>
void process_plane_rows(
    const process_plane_params& params,
    span2d::ReadOnlyRestrictPlane<PixelIn> src_plane,
    span2d::RestrictPlane<PixelOut> dst_plane,
    span2d::ReadOnlyRestrictPlane<std::int16_t> grain_plane,
    const neo_f3kdb::core::PlaneOffsetCache* cache
) {
    static_assert(!kUseCachedOffsets || (kSampleMode >= 1 && kSampleMode <= 6));

    const auto& config = params.config;
    const int height = params.plane_height();
    const int width = params.plane_width();

    // Rebind in process_block preserves the lane count, so these caps also
    // bound the uint16_t deband core and its cached reference gather.
    using BlockTag = std::conditional_t<
        kSampleMode == 2 || kSampleMode == 4 || kSampleMode == 5,
        hn::CappedTag<std::uint32_t, 8>,
        std::conditional_t<
            kSampleMode == 1 || kSampleMode == 3,
            hn::CappedTag<std::uint32_t, 16>,
            hn::ScalableTag<std::uint32_t>
        >
    >;
    const BlockTag du32;
    const hn::Rebind<std::int32_t, BlockTag> di32;
    const std::size_t lanes = hn::Lanes(di32);
    const int vec_width = width - static_cast<int>(width % lanes);

    using FloydSteinbergDither = neo_f3kdb::core::dither::FloydSteinbergDither;
    alignas(FloydSteinbergDither) char fs_context[CONTEXT_BUFFER_SIZE];
    FloydSteinbergDither* fs_dither = nullptr;
    if constexpr (kDitherAlgo == DA_HIGH_FLOYD_STEINBERG_DITHERING) {
        fs_dither = new (fs_context) FloydSteinbergDither(
            fs_context + sizeof(FloydSteinbergDither),
            CONTEXT_BUFFER_SIZE - static_cast<int>(sizeof(FloydSteinbergDither)),
            width,
            config.output_depth
        );
    }

    auto dst_cursor = dst_plane.cursor();
    auto grain_cursor = grain_plane.cursor();

    for (int row = 0; row < height; ++row) {
        auto dst_row = *dst_cursor;
        const auto grain_row = *grain_cursor;
        const auto* src_row = src_plane.row_ptr(row);

        int col = 0;
        if constexpr (kUseCachedOffsets) {
            const std::size_t row_offset = static_cast<std::size_t>(row) * width;
            const auto* off1 = cache->off1.data() + row_offset;
            const auto* off2 = cache->off2.data() + row_offset;
            const auto* off3 = cache->off3.data() + row_offset;
            const auto* off4 = cache->off4.data() + row_offset;

            for (; col < vec_width; col += static_cast<int>(lanes)) {
                deband_hwy_detail::process_block<kSampleMode, kBlurFirst, kDitherAlgo, true, kPureDeband>(
                    params,
                    fs_dither,
                    src_plane,
                    src_row,
                    dst_row,
                    grain_row,
                    row,
                    col,
                    off1 + col,
                    off2 + col,
                    off3 + col,
                    off4 + col,
                    di32,
                    du32,
                    lanes
                );
            }
        } else {
            for (; col < vec_width; col += static_cast<int>(lanes)) {
                deband_hwy_detail::process_block<kSampleMode, kBlurFirst, kDitherAlgo, false, kPureDeband>(
                    params,
                    fs_dither,
                    src_plane,
                    src_row,
                    dst_row,
                    grain_row,
                    row,
                    col,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    di32,
                    du32,
                    lanes
                );
            }
        }

        for (; col < width; ++col) {
            int pixel = neo_f3kdb::core::sample_modes::process_pixel<kSampleMode, kBlurFirst>(
                params,
                src_plane,
                row,
                col
            );
            if constexpr (kPureDeband) {
                dst_row[static_cast<std::size_t>(col)] = static_cast<PixelOut>(pixel >> 8);
            } else if constexpr (kDitherAlgo == DA_HIGH_FLOYD_STEINBERG_DITHERING) {
                pixel = deband_hwy_detail::postprocess_floyd_steinberg_pixel(
                    config,
                    *fs_dither,
                    pixel,
                    grain_row,
                    col
                );
                fs_dither->next_pixel();
            } else {
                pixel = deband_hwy_detail::postprocess_scalar_pixel<kDitherAlgo>(
                    config,
                    pixel,
                    grain_row,
                    row,
                    col
                );
            }
            dst_row[static_cast<std::size_t>(col)] = static_cast<PixelOut>(pixel);
        }

        if constexpr (kDitherAlgo == DA_HIGH_FLOYD_STEINBERG_DITHERING) {
            fs_dither->next_row();
        }

        ++dst_cursor;
        ++grain_cursor;
    }

    if constexpr (kDitherAlgo == DA_HIGH_FLOYD_STEINBERG_DITHERING) {
        fs_dither->~FloydSteinbergDither();
    }
}

template <int kSampleMode, bool kBlurFirst, bool kUseCachedOffsets, class PixelIn, class PixelOut>
void process_plane_rows_fs_wavefront(
    const process_plane_params& params,
    span2d::ReadOnlyRestrictPlane<PixelIn> src_plane,
    span2d::RestrictPlane<PixelOut> dst_plane,
    span2d::ReadOnlyRestrictPlane<std::int16_t> grain_plane,
    const neo_f3kdb::core::PlaneOffsetCache* cache
) {
    const auto& config = params.config;
    const int height = params.plane_height();
    const int width = params.plane_width();

    const hn::ScalableTag<std::uint16_t> d;
    const std::size_t N = hn::Lanes(d);

    using BlockTag = std::conditional_t<
        kSampleMode == 2 || kSampleMode == 4 || kSampleMode == 5,
        hn::CappedTag<std::uint32_t, 8>,
        std::conditional_t<
            kSampleMode == 1 || kSampleMode == 3,
            hn::CappedTag<std::uint32_t, 16>,
            hn::ScalableTag<std::uint32_t>
        >
    >;
    const hn::Rebind<std::int32_t, BlockTag> di32;
    const hn::Rebind<std::uint16_t, BlockTag> du16;
    const hn::Rebind<std::int16_t, BlockTag> d16;
    const std::size_t deband_lanes = hn::Lanes(di32);
    const int vec_width = width - static_cast<int>(width % deband_lanes);

    std::vector<std::uint16_t> stripe_intermediate(N * static_cast<std::size_t>(width), 0);
    std::vector<std::uint16_t> cur_bridge(width + 2, 0);
    std::vector<std::uint16_t> next_bridge(width + 2, 0);

    std::vector<const std::uint16_t*> stripe_src_rows(N, nullptr);
    std::vector<PixelOut*> stripe_dst_rows(N, nullptr);

    for (int y0 = 0; y0 < height; y0 += static_cast<int>(N)) {
        const int stripe_height = std::min(static_cast<int>(N), height - y0);
        std::fill(next_bridge.begin(), next_bridge.end(), static_cast<std::uint16_t>(0));

        // 1. Deband + Grain filtering for rows in this stripe into stripe_intermediate
        for (int k = 0; k < stripe_height; ++k) {
            const int row = y0 + k;
            const auto* src_row = src_plane.row_ptr(row);
            const auto grain_row = grain_plane.row(row);
            auto* inter_row = stripe_intermediate.data() + static_cast<std::size_t>(k) * width;
            stripe_src_rows[k] = inter_row;
            stripe_dst_rows[k] = dst_plane.row_ptr(row);

            int col = 0;
            if constexpr ((kSampleMode >= 1 && kSampleMode <= 5) && std::is_same_v<PixelIn, unsigned char>) {
                const hn::Rebind<unsigned char, BlockTag> d8;
                if constexpr (kUseCachedOffsets) {
                    const std::size_t row_offset = static_cast<std::size_t>(row) * width;
                    const auto* off1 = cache->off1.data() + row_offset;
                    const auto* off2 = cache->off2.data() + row_offset;
                    const auto* off3 = cache->off3.data() + row_offset;
                    const auto* off4 = cache->off4.data() + row_offset;

                    for (; col < vec_width; col += static_cast<int>(deband_lanes)) {
                        auto bytes = hn::LoadU(d8, src_row + col);
                        auto v_src_u16 = hn::ShiftLeft<8>(hn::PromoteTo(du16, bytes));
                        auto r1_v = hn::Zero(du16), r2_v = hn::Zero(du16), r3_v = hn::Zero(du16), r4_v = hn::Zero(du16);
                        deband_hwy_detail::gather_cached_reference_vectors(
                            du16, src_row + col, off1 + col, off2 + col, off3 + col, off4 + col,
                            16 - config.input_depth, deband_lanes, r1_v, r2_v, r3_v, r4_v
                        );
                        auto processed_u16 = deband_hwy_detail::process_reference_samples_u16<kSampleMode, kBlurFirst>(
                            du16, v_src_u16, r1_v, r2_v, r3_v, r4_v,
                            config.threshold, config.threshold1, config.threshold2
                        );
                        auto grain_v = hn::BitCast(du16, hn::LoadU(d16, grain_row.data() + col));
                        auto intermediate_v = hn::Add(processed_u16, grain_v);
                        hn::StoreU(intermediate_v, du16, inter_row + col);
                    }
                } else {
                    for (; col < vec_width; col += static_cast<int>(deband_lanes)) {
                        auto bytes = hn::LoadU(d8, src_row + col);
                        auto v_src_u16 = hn::ShiftLeft<8>(hn::PromoteTo(du16, bytes));
                        auto r1_v = hn::Zero(du16), r2_v = hn::Zero(du16), r3_v = hn::Zero(du16), r4_v = hn::Zero(du16);
                        deband_hwy_detail::gather_uncached_reference_vectors<kSampleMode>(
                            du16, params, src_plane, row, col, deband_lanes, r1_v, r2_v, r3_v, r4_v
                        );
                        auto processed_u16 = deband_hwy_detail::process_reference_samples_u16<kSampleMode, kBlurFirst>(
                            du16, v_src_u16, r1_v, r2_v, r3_v, r4_v,
                            config.threshold, config.threshold1, config.threshold2
                        );
                        auto grain_v = hn::BitCast(du16, hn::LoadU(d16, grain_row.data() + col));
                        auto intermediate_v = hn::Add(processed_u16, grain_v);
                        hn::StoreU(intermediate_v, du16, inter_row + col);
                    }
                }
            } else {
                for (; col < vec_width; col += static_cast<int>(deband_lanes)) {
                    decltype(hn::Zero(di32)) v_src_up;
                    if constexpr (std::is_same_v<PixelIn, unsigned char>) {
                        if (config.input_mode == LOW_BIT_DEPTH && config.input_depth == 8) {
                            const hn::Rebind<unsigned char, BlockTag> d8;
                            auto bytes = hn::LoadU(d8, src_row + col);
                            auto u16 = hn::PromoteTo(du16, bytes);
                            auto u32 = hn::PromoteTo(di32, u16);
                            v_src_up = hn::ShiftLeft<8>(u32);
                        } else {
                            alignas(64) std::int32_t src_up_buf[hn::MaxLanes(di32)];
                            for (std::size_t lane = 0; lane < deband_lanes; ++lane) {
                                src_up_buf[lane] = neo_f3kdb::core::sample_modes::upsample(
                                    src_row[col + static_cast<int>(lane)],
                                    config.input_depth
                                );
                            }
                            v_src_up = hn::LoadU(di32, src_up_buf);
                        }
                    } else {
                        alignas(64) std::int32_t src_up_buf[hn::MaxLanes(di32)];
                        for (std::size_t lane = 0; lane < deband_lanes; ++lane) {
                            src_up_buf[lane] = neo_f3kdb::core::sample_modes::upsample(
                                src_row[col + static_cast<int>(lane)],
                                config.input_depth
                            );
                        }
                        v_src_up = hn::LoadU(di32, src_up_buf);
                    }

                    auto r1_v = hn::Zero(di32), r2_v = hn::Zero(di32), r3_v = hn::Zero(di32), r4_v = hn::Zero(di32);
                    if constexpr (kUseCachedOffsets) {
                        const std::size_t row_offset = static_cast<std::size_t>(row) * width;
                        const auto* off1 = cache->off1.data() + row_offset;
                        const auto* off2 = cache->off2.data() + row_offset;
                        const auto* off3 = cache->off3.data() + row_offset;
                        const auto* off4 = cache->off4.data() + row_offset;
                        deband_hwy_detail::gather_cached_reference_vectors(
                            di32, src_row + col, off1 + col, off2 + col, off3 + col, off4 + col,
                            16 - config.input_depth, deband_lanes, r1_v, r2_v, r3_v, r4_v
                        );
                    } else {
                        deband_hwy_detail::gather_uncached_reference_vectors<kSampleMode>(
                            di32, params, src_plane, row, col, deband_lanes, r1_v, r2_v, r3_v, r4_v
                        );
                    }

                    auto processed_v = deband_hwy_detail::process_reference_samples<kSampleMode, kBlurFirst>(
                        di32, v_src_up, r1_v, r2_v, r3_v, r4_v, config.threshold, config.threshold1, config.threshold2
                    );
                    auto grain_v = hn::PromoteTo(di32, hn::LoadU(d16, grain_row.data() + col));
                    auto sum_v = hn::Add(processed_v, grain_v);
                    auto clamped_u16 = hn::DemoteTo(du16, hn::Min(hn::Max(sum_v, hn::Zero(di32)), hn::Set(di32, 65535)));
                    hn::StoreU(clamped_u16, du16, inter_row + col);
                }
            }

            // Tail columns
            for (; col < width; ++col) {
                int pixel = neo_f3kdb::core::sample_modes::process_pixel<kSampleMode, kBlurFirst>(
                    params, src_plane, row, col
                );
                pixel += static_cast<int>(grain_row[static_cast<std::size_t>(col)]);
                inter_row[col] = static_cast<std::uint16_t>(std::clamp(pixel, 0, 65535));
            }
        }

        // 2. Wavefront F-S Dither across the stripe
        deband_hwy_detail::process_stripe_fs_wavefront_hwy(
            config,
            stripe_src_rows.data(),
            stripe_dst_rows.data(),
            width,
            height,
            y0,
            stripe_height,
            cur_bridge,
            next_bridge
        );

        std::swap(cur_bridge, next_bridge);
    }
}

template <int kSampleMode, bool kBlurFirst, int kDitherAlgo, bool kPureDeband, class PixelIn, class PixelOut>
void process_plane_templated_impl(const process_plane_params& params, process_plane_context* context) {
    static_assert(kSampleMode >= 1 && kSampleMode <= 7);
    static_assert(
        kDitherAlgo == DA_HIGH_NO_DITHERING ||
        kDitherAlgo == DA_HIGH_ORDERED_DITHERING ||
        kDitherAlgo == DA_HIGH_FLOYD_STEINBERG_DITHERING ||
        kDitherAlgo == DA_16BIT_INTERLEAVED
    );

    const auto src_plane = params.src_plane<PixelIn>();
    auto dst_plane = params.dst_plane<PixelOut>();
    const auto grain_plane = params.grain_plane();

    if constexpr (kDitherAlgo == DA_HIGH_FLOYD_STEINBERG_DITHERING) {
        if constexpr (kSampleMode >= 1 && kSampleMode <= 6) {
            const auto* cache = deband_hwy_detail::prepare_offset_cache<kSampleMode>(params, src_plane, context);
            if (cache) {
                process_plane_rows_fs_wavefront<kSampleMode, kBlurFirst, true, PixelIn, PixelOut>(
                    params,
                    src_plane,
                    dst_plane,
                    grain_plane,
                    cache
                );
                return;
            }
        }

        process_plane_rows_fs_wavefront<kSampleMode, kBlurFirst, false, PixelIn, PixelOut>(
            params,
            src_plane,
            dst_plane,
            grain_plane,
            nullptr
        );
        return;
    }

    if constexpr (kSampleMode >= 1 && kSampleMode <= 6) {
        const auto* cache = deband_hwy_detail::prepare_offset_cache<kSampleMode>(params, src_plane, context);
        if (cache) {
            process_plane_rows<kSampleMode, kBlurFirst, kDitherAlgo, true, kPureDeband>(
                params,
                src_plane,
                dst_plane,
                grain_plane,
                cache
            );
            return;
        }
    }

    process_plane_rows<kSampleMode, kBlurFirst, kDitherAlgo, false, kPureDeband>(
        params,
        src_plane,
        dst_plane,
        grain_plane,
        nullptr
    );
}

template <int kSampleMode, bool kBlurFirst, int kDitherAlgo, class PixelIn, class PixelOut>
void process_plane_templated(const process_plane_params& params, process_plane_context* context) {
    const auto& config = params.config;

    if constexpr (
        kSampleMode >= 1 && kSampleMode <= 5 &&
        kDitherAlgo == DA_HIGH_NO_DITHERING &&
        std::is_same_v<PixelIn, unsigned char> &&
        std::is_same_v<PixelOut, unsigned char>
    ) {
        // Select this once per plane so the block path can omit zero-grain work.
        if (!config.has_grain && config.output_depth == 8 &&
            config.pixel_min == FULL_RANGE_Y_MIN && config.pixel_max == FULL_RANGE_Y_MAX) {
            process_plane_templated_impl<kSampleMode, kBlurFirst, kDitherAlgo, true, PixelIn, PixelOut>(params, context);
            return;
        }
    }

    process_plane_templated_impl<kSampleMode, kBlurFirst, kDitherAlgo, false, PixelIn, PixelOut>(params, context);
}

template <int kSampleMode, bool kBlurFirst, class PixelIn, class PixelOut>
void dispatch_dither(const process_plane_params& params, process_plane_context* context) {
    switch (params.config.dither_algo) {
        case DA_HIGH_NO_DITHERING:
            process_plane_templated<kSampleMode, kBlurFirst, DA_HIGH_NO_DITHERING, PixelIn, PixelOut>(params, context);
            break;
        case DA_HIGH_ORDERED_DITHERING:
            process_plane_templated<kSampleMode, kBlurFirst, DA_HIGH_ORDERED_DITHERING, PixelIn, PixelOut>(params, context);
            break;
        case DA_HIGH_FLOYD_STEINBERG_DITHERING:
            process_plane_templated<kSampleMode, kBlurFirst, DA_HIGH_FLOYD_STEINBERG_DITHERING, PixelIn, PixelOut>(
                params,
                context
            );
            break;
        case DA_16BIT_INTERLEAVED:
            process_plane_templated<kSampleMode, kBlurFirst, DA_16BIT_INTERLEAVED, PixelIn, PixelOut>(params, context);
            break;
        default:
            abort();
    }
}

template <int kSampleMode, class PixelIn, class PixelOut>
void dispatch_blur(const process_plane_params& params, process_plane_context* context) {
    if (params.config.blur_first) {
        dispatch_dither<kSampleMode, true, PixelIn, PixelOut>(params, context);
    } else {
        dispatch_dither<kSampleMode, false, PixelIn, PixelOut>(params, context);
    }
}

template <class PixelIn, class PixelOut>
void dispatch_sample_mode(const process_plane_params& params, process_plane_context* context) {
    switch (params.config.sample_mode) {
        case 1:
            dispatch_blur<1, PixelIn, PixelOut>(params, context);
            break;
        case 2:
            dispatch_blur<2, PixelIn, PixelOut>(params, context);
            break;
        case 3:
            dispatch_blur<3, PixelIn, PixelOut>(params, context);
            break;
        case 4:
            dispatch_blur<4, PixelIn, PixelOut>(params, context);
            break;
        case 5:
            dispatch_blur<5, PixelIn, PixelOut>(params, context);
            break;
        case 6:
            dispatch_blur<6, PixelIn, PixelOut>(params, context);
            break;
        case 7:
            dispatch_blur<7, PixelIn, PixelOut>(params, context);
            break;
        default:
            abort();
    }
}

template <class PixelIn>
void dispatch_pixel_out(const process_plane_params& params, process_plane_context* context) {
    if (params.config.output_depth <= 8) {
        dispatch_sample_mode<PixelIn, std::uint8_t>(params, context);
    } else {
        dispatch_sample_mode<PixelIn, std::uint16_t>(params, context);
    }
}

void ProcessPlaneHWY(const process_plane_params& params, process_plane_context* context) {
    if (params.config.input_depth == 8) {
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

namespace core {

void process_plane_highway(KernelExecution execution) {
    HWY_DYNAMIC_DISPATCH(ProcessPlaneHWY)(execution.input, execution.context);
}

} // namespace core
} // namespace neo_f3kdb
#endif
