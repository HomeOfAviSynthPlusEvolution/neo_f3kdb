#include "constants.h"
#include "core/kernel.hpp"
#include "core/math.hpp"
#include "core/dither_ordered.hpp"

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
    static_assert(kSampleMode >= 1 && kSampleMode <= 5);
    static_assert(
        kDitherAlgo == DA_HIGH_NO_DITHERING ||
        kDitherAlgo == DA_HIGH_ORDERED_DITHERING ||
        kDitherAlgo == DA_16BIT_INTERLEAVED
    );

    const int height = params.plane_height_in_pixels;
    const int width = params.plane_width_in_pixels;
    const int src_stride = params.src_pitch / sizeof(PixelIn);
    const int dst_stride = params.dst_pitch / sizeof(PixelOut);

    const PixelIn* src_base = reinterpret_cast<const PixelIn*>(params.src_plane_ptr);
    PixelOut* dst_base = reinterpret_cast<PixelOut*>(params.dst_plane_ptr);

    const hn::ScalableTag<std::uint32_t> du32;
    const hn::Rebind<std::int32_t, decltype(du32)> di32;
    const std::size_t lanes = hn::Lanes(di32);
    const int vec_width = width - static_cast<int>(width % lanes);

    for (int row = 0; row < height; ++row) {
        const PixelIn* src_row = src_base + row * src_stride;
        PixelOut* dst_row = dst_base + row * dst_stride;
        const std::int16_t* grain_row = params.grain_buffer + params.grain_buffer_stride * row;
        const pixel_dither_info* info_row = params.info_ptr_base + params.info_stride * row;

        int col = 0;
        for (; col < vec_width; col += static_cast<int>(lanes)) {
            deband_hwy_detail::process_block<kSampleMode, kBlurFirst, kDitherAlgo>(
                params,
                src_base,
                src_row,
                dst_row,
                grain_row,
                info_row,
                src_stride,
                row,
                col,
                di32,
                du32,
                lanes
            );
        }

        for (; col < width; ++col) {
            int pixel = neo_f3kdb::core::deband_scalar::process_pixel<kSampleMode, kBlurFirst>(
                params,
                src_base,
                src_row,
                src_stride,
                row,
                col
            );
            pixel = deband_hwy_detail::postprocess_scalar_pixel<kDitherAlgo>(params, pixel, grain_row, row, col);
            dst_row[col] = static_cast<PixelOut>(pixel);
        }
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
