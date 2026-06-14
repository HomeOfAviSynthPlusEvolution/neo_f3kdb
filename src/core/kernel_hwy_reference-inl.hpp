// Intentionally no include guard: included by kernel_hwy-inl.hpp once per Highway target.

template <int kSampleMode, class PixelIn>
void gather_reference_pixels(
    const process_plane_params& params,
    const PixelIn* src_base,
    const PixelIn* src_row,
    const pixel_dither_info* info_row,
    int src_stride,
    int row,
    int col,
    std::size_t lanes,
    std::int32_t* src_up,
    std::int32_t* ref1,
    std::int32_t* ref2,
    std::int32_t* ref3,
    std::int32_t* ref4
) {
    static_assert(kSampleMode >= 1 && kSampleMode <= 5);
    (void)info_row;

    for (std::size_t lane = 0; lane < lanes; ++lane) {
        const int pixel_col = col + static_cast<int>(lane);
        const auto samples = neo_f3kdb::core::sample_modes::load_reference_samples<kSampleMode>(
            params,
            src_base,
            src_row,
            src_stride,
            row,
            pixel_col
        );
        src_up[lane] = samples.src;
        ref1[lane] = samples.ref1;
        ref2[lane] = samples.ref2;
        ref3[lane] = samples.ref3;
        ref4[lane] = samples.ref4;
    }
}
