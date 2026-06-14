#include "core/constants.hpp"
#include "core/dither_floyd_steinberg.hpp"
#include "core/dither_ordered.hpp"
#include "core/kernel.hpp"
#include "core/pixel_proc_common.hpp"
#include "core/sample_modes.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "core/kernel_hwy.cpp"

#include "hwy/foreach_target.h"
#include "hwy/highway.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <type_traits>

HWY_BEFORE_NAMESPACE();
namespace neo_f3kdb::HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

#include "core/kernel_hwy-inl.hpp"

template <int kSampleMode, bool kBlurFirst, int kDitherAlgo, class PixelIn, class PixelOut>
void process_plane_templated(const process_plane_params& params) {
    static_assert(kSampleMode >= 1 && kSampleMode <= 7);
    static_assert(
        kDitherAlgo == DA_HIGH_NO_DITHERING ||
        kDitherAlgo == DA_HIGH_ORDERED_DITHERING ||
        kDitherAlgo == DA_HIGH_FLOYD_STEINBERG_DITHERING ||
        kDitherAlgo == DA_16BIT_INTERLEAVED
    );

    const int height = params.plane_height_in_pixels;
    const int width = params.plane_width_in_pixels;
    const auto src_plane = params.src_plane<PixelIn>();
    auto dst_plane = params.dst_plane<PixelOut>();
    const auto grain_plane = params.grain_plane();

    const hn::ScalableTag<std::uint32_t> du32;
    const hn::Rebind<std::int32_t, decltype(du32)> di32;
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
            params.output_depth
        );
    }

    for (int row = 0; row < height; ++row) {
        auto dst_row = dst_plane.row(row);
        const auto grain_row = grain_plane.row(row);

        int col = 0;
        for (; col < vec_width; col += static_cast<int>(lanes)) {
            deband_hwy_detail::process_block<kSampleMode, kBlurFirst, kDitherAlgo>(
                params,
                fs_dither,
                src_plane,
                dst_row,
                grain_row,
                row,
                col,
                di32,
                du32,
                lanes
            );
        }

        for (; col < width; ++col) {
            int pixel = neo_f3kdb::core::sample_modes::process_pixel<kSampleMode, kBlurFirst>(
                params,
                src_plane,
                row,
                col
            );
            if constexpr (kDitherAlgo == DA_HIGH_FLOYD_STEINBERG_DITHERING) {
                pixel = deband_hwy_detail::postprocess_floyd_steinberg_pixel(
                    params,
                    *fs_dither,
                    pixel,
                    grain_row,
                    col
                );
                fs_dither->next_pixel();
            } else {
                pixel = deband_hwy_detail::postprocess_scalar_pixel<kDitherAlgo>(
                    params,
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
    }

    if constexpr (kDitherAlgo == DA_HIGH_FLOYD_STEINBERG_DITHERING) {
        fs_dither->~FloydSteinbergDither();
    }
}

template <int kSampleMode, bool kBlurFirst, class PixelIn, class PixelOut>
void dispatch_dither(const process_plane_params& params) {
    switch (params.dither_algo) {
        case DA_HIGH_NO_DITHERING:
            process_plane_templated<kSampleMode, kBlurFirst, DA_HIGH_NO_DITHERING, PixelIn, PixelOut>(params);
            break;
        case DA_HIGH_ORDERED_DITHERING:
            process_plane_templated<kSampleMode, kBlurFirst, DA_HIGH_ORDERED_DITHERING, PixelIn, PixelOut>(params);
            break;
        case DA_HIGH_FLOYD_STEINBERG_DITHERING:
            process_plane_templated<kSampleMode, kBlurFirst, DA_HIGH_FLOYD_STEINBERG_DITHERING, PixelIn, PixelOut>(
                params
            );
            break;
        case DA_16BIT_INTERLEAVED:
            process_plane_templated<kSampleMode, kBlurFirst, DA_16BIT_INTERLEAVED, PixelIn, PixelOut>(params);
            break;
        default:
            abort();
    }
}

template <int kSampleMode, class PixelIn, class PixelOut>
void dispatch_blur(const process_plane_params& params) {
    if (params.blur_first) {
        dispatch_dither<kSampleMode, true, PixelIn, PixelOut>(params);
    } else {
        dispatch_dither<kSampleMode, false, PixelIn, PixelOut>(params);
    }
}

template <class PixelIn, class PixelOut>
void dispatch_sample_mode(const process_plane_params& params) {
    switch (params.sample_mode) {
        case 1:
            dispatch_blur<1, PixelIn, PixelOut>(params);
            break;
        case 2:
            dispatch_blur<2, PixelIn, PixelOut>(params);
            break;
        case 3:
            dispatch_blur<3, PixelIn, PixelOut>(params);
            break;
        case 4:
            dispatch_blur<4, PixelIn, PixelOut>(params);
            break;
        case 5:
            dispatch_blur<5, PixelIn, PixelOut>(params);
            break;
        case 6:
            dispatch_blur<6, PixelIn, PixelOut>(params);
            break;
        case 7:
            dispatch_blur<7, PixelIn, PixelOut>(params);
            break;
        default:
            abort();
    }
}

template <class PixelIn>
void dispatch_pixel_out(const process_plane_params& params) {
    if (params.output_depth <= 8) {
        dispatch_sample_mode<PixelIn, std::uint8_t>(params);
    } else {
        dispatch_sample_mode<PixelIn, std::uint16_t>(params);
    }
}

void ProcessPlaneHWY(const process_plane_params& params, [[maybe_unused]] process_plane_context* context) {
    if (params.input_depth == 8) {
        dispatch_pixel_out<std::uint8_t>(params);
    } else {
        dispatch_pixel_out<std::uint16_t>(params);
    }
}

} // namespace neo_f3kdb::HWY_NAMESPACE
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
#include "hwy/per_target.h"

namespace neo_f3kdb {

HWY_EXPORT(ProcessPlaneHWY);

namespace core {

void process_plane_highway(const PlaneJob& job) {
    HWY_DYNAMIC_DISPATCH(ProcessPlaneHWY)(job.params, job.context);
}

} // namespace core
} // namespace neo_f3kdb
#endif
