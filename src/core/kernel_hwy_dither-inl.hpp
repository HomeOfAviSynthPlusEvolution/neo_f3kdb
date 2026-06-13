// Intentionally no include guard: included by kernel_hwy-inl.hpp once per Highway target.

template <int kDitherAlgo>
int postprocess_scalar_pixel(const process_plane_params& params, int pixel, const std::int16_t* grain_row, int row, int col) {
    pixel += static_cast<int>(grain_row[col]);

    if constexpr (kDitherAlgo == DA_HIGH_NO_DITHERING) {
        return std::clamp(pixel, params.pixel_min, params.pixel_max) >> (16 - params.output_depth);
    } else if constexpr (kDitherAlgo == DA_HIGH_ORDERED_DITHERING) {
        pixel += pixel_proc_high_ordered_dithering::THRESHOLD_MAP[row & 15][col & 15] >> (params.output_depth - 8);
        return std::clamp(pixel, params.pixel_min, params.pixel_max) >> (16 - params.output_depth);
    } else {
        static_assert(kDitherAlgo == DA_16BIT_INTERLEAVED);
        return std::clamp(pixel, params.pixel_min, params.pixel_max);
    }
}

template <int kDitherAlgo, class D32, class V32>
V32 apply_dither_and_grain(
    const process_plane_params& params,
    D32 d32,
    V32 pixel,
    const std::int16_t* grain_row,
    int row,
    int col,
    std::size_t lanes
) {
    const hn::Rebind<std::int16_t, D32> d16;
    pixel = hn::Add(pixel, hn::PromoteTo(d32, hn::LoadU(d16, grain_row + col)));

    if constexpr (kDitherAlgo == DA_HIGH_ORDERED_DITHERING) {
        alignas(64) std::int32_t dither[hn::MaxLanes(d32)];
        for (std::size_t lane = 0; lane < lanes; ++lane) {
            dither[lane] = pixel_proc_high_ordered_dithering::THRESHOLD_MAP[row & 15][(col + static_cast<int>(lane)) & 15] >> (params.output_depth - 8);
        }
        pixel = hn::Add(pixel, hn::LoadU(d32, dither));
    }

    pixel = hn::Min(hn::Max(pixel, hn::Set(d32, params.pixel_min)), hn::Set(d32, params.pixel_max));
    if constexpr (kDitherAlgo != DA_16BIT_INTERLEAVED) {
        pixel = hn::ShiftRightSame(pixel, 16 - params.output_depth);
    }
    return pixel;
}
