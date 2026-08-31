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
        kSampleMode == 2 || kSampleMode == 4,
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
        kSampleMode >= 1 && kSampleMode <= 4 &&
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
