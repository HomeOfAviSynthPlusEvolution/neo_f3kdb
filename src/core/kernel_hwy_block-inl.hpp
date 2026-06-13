// Intentionally no include guard: included by kernel_hwy-inl.hpp once per Highway target.

template <int kSampleMode, bool kBlurFirst, int kDitherAlgo, class PixelIn, class PixelOut, class D32, class DU32>
void process_block(
    const process_plane_params& params,
    const PixelIn* src_base,
    const PixelIn* src_row,
    PixelOut* dst_row,
    const std::int16_t* grain_row,
    const pixel_dither_info* info_row,
    int src_stride,
    int row,
    int col,
    D32 d32,
    DU32 du32,
    std::size_t lanes
) {
    alignas(64) std::int32_t src_up[hn::MaxLanes(d32)];
    alignas(64) std::int32_t ref1[hn::MaxLanes(d32)];
    alignas(64) std::int32_t ref2[hn::MaxLanes(d32)];
    alignas(64) std::int32_t ref3[hn::MaxLanes(d32)];
    alignas(64) std::int32_t ref4[hn::MaxLanes(d32)];

    if constexpr (kSampleMode == 6 || kSampleMode == 7) {
        for (std::size_t lane = 0; lane < lanes; ++lane) {
            src_up[lane] = neo_f3kdb::core::sample_modes::process_pixel<kSampleMode, kBlurFirst>(
                params,
                src_base,
                src_row,
                src_stride,
                row,
                col + static_cast<int>(lane)
            );
        }
    } else {
        gather_reference_pixels<kSampleMode>(
            params,
            src_base,
            src_row,
            info_row,
            src_stride,
            row,
            col,
            lanes,
            src_up,
            ref1,
            ref2,
            ref3,
            ref4
        );

        hn::StoreU(
            process_refs<kSampleMode, kBlurFirst>(
                d32,
                hn::LoadU(d32, src_up),
                ref1,
                ref2,
                ref3,
                ref4,
                params.threshold,
                params.threshold1,
                params.threshold2
            ),
            d32,
            src_up
        );
    }

    auto pixel = hn::LoadU(d32, src_up);
    pixel = apply_dither_and_grain<kDitherAlgo>(params, d32, pixel, grain_row, row, col, lanes);
    store_pixels_from_u32(du32, dst_row + col, hn::BitCast(du32, pixel));
}
