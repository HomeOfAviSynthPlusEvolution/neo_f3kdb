// Intentionally no include guard: included by kernel_hwy-inl.hpp once per Highway target.

template <int kSampleMode, bool kBlurFirst, int kDitherAlgo, bool kUseCachedOffsets, bool kPureDeband, class PixelIn, class PixelOut, class D32, class DU32>
HWY_INLINE void process_block(
    const process_plane_params& params,
    neo_f3kdb::core::dither::FloydSteinbergDither* fs_dither,
    span2d::ReadOnlyRestrictPlane<PixelIn> src_plane,
    const PixelIn* src_row,
    span2d::RestrictRow<PixelOut> dst_row,
    span2d::ReadOnlyRestrictRow<std::int16_t> grain_row,
    int row,
    int col,
    const std::int32_t* off1,
    const std::int32_t* off2,
    const std::int32_t* off3,
    const std::int32_t* off4,
    D32 d32,
    DU32 du32,
    std::size_t lanes
) {
    static_assert(!kUseCachedOffsets || (kSampleMode >= 1 && kSampleMode <= 6));
    static_assert(!kPureDeband || (
        kSampleMode >= 1 && kSampleMode <= 4 &&
        kDitherAlgo == DA_HIGH_NO_DITHERING &&
        std::is_same_v<PixelIn, unsigned char> &&
        std::is_same_v<PixelOut, unsigned char>
    ));

    const auto& config = params.config;
    decltype(hn::Zero(d32)) v_src_up;
    decltype(hn::Zero(d32)) processed_v;

    if constexpr ((kSampleMode >= 1 && kSampleMode <= 4) && std::is_same_v<PixelIn, unsigned char> && kDitherAlgo != DA_HIGH_FLOYD_STEINBERG_DITHERING) {
        const hn::Rebind<std::uint16_t, D32> d16;
        const hn::Rebind<unsigned char, D32> d8;

        auto bytes = hn::LoadU(d8, src_row + col);
        auto v_src_u16 = hn::ShiftLeft<8>(hn::PromoteTo(d16, bytes));

        auto r1_v = hn::Zero(d16);
        auto r2_v = hn::Zero(d16);
        auto r3_v = hn::Zero(d16);
        auto r4_v = hn::Zero(d16);

        if constexpr (kUseCachedOffsets) {
            gather_cached_reference_vectors(
                d16,
                src_row + col,
                off1,
                off2,
                off3,
                off4,
                16 - config.input_depth,
                lanes,
                r1_v,
                r2_v,
                r3_v,
                r4_v
            );
        } else {
            gather_uncached_reference_vectors<kSampleMode>(
                d16,
                params,
                src_plane,
                row,
                col,
                lanes,
                r1_v,
                r2_v,
                r3_v,
                r4_v
            );
        }

        auto processed_u16 = process_reference_samples_u16<kSampleMode, kBlurFirst>(
            d16,
            v_src_u16,
            r1_v,
            r2_v,
            r3_v,
            r4_v,
            config.threshold
        );

        if constexpr (kPureDeband) {
            const auto output = hn::DemoteTo(d8, hn::ShiftRight<8>(processed_u16));
            hn::StoreU(output, d8, dst_row.data() + col);
        } else {
            auto processed_i32 = hn::PromoteTo(d32, processed_u16);
            auto pixel = apply_dither_and_grain<kDitherAlgo>(config, d32, processed_i32, grain_row, row, col, lanes);
            store_pixels_from_u32(du32, dst_row.data() + col, hn::BitCast(du32, pixel));
        }
    } else {
        if constexpr (std::is_same_v<PixelIn, unsigned char>) {
            if (config.input_mode == LOW_BIT_DEPTH && config.input_depth == 8) {
                const hn::Rebind<unsigned char, D32> d8;
                const hn::Rebind<std::uint16_t, D32> du16;
                auto bytes = hn::LoadU(d8, src_row + col);
                auto u16 = hn::PromoteTo(du16, bytes);
                auto u32 = hn::PromoteTo(d32, u16);
                v_src_up = hn::ShiftLeft<8>(u32);
            } else {
                alignas(64) std::int32_t src_up_buf[hn::MaxLanes(d32)];
                for (std::size_t lane = 0; lane < lanes; ++lane) {
                    src_up_buf[lane] = neo_f3kdb::core::sample_modes::upsample(
                        src_row[col + static_cast<int>(lane)],
                        config.input_depth
                    );
                }
                v_src_up = hn::LoadU(d32, src_up_buf);
            }
        } else {
            alignas(64) std::int32_t src_up_buf[hn::MaxLanes(d32)];
            for (std::size_t lane = 0; lane < lanes; ++lane) {
                src_up_buf[lane] = neo_f3kdb::core::sample_modes::upsample(
                    src_row[col + static_cast<int>(lane)],
                    config.input_depth
                );
            }
            v_src_up = hn::LoadU(d32, src_up_buf);
        }

        if constexpr (kSampleMode == 7) {
            alignas(64) std::int32_t src_up_buf[hn::MaxLanes(d32)];
            for (std::size_t lane = 0; lane < lanes; ++lane) {
                src_up_buf[lane] = neo_f3kdb::core::sample_modes::process_pixel<kSampleMode, kBlurFirst>(
                    params,
                    src_plane,
                    row,
                    col + static_cast<int>(lane)
                );
            }
            processed_v = hn::LoadU(d32, src_up_buf);
        } else {
            auto r1_v = hn::Zero(d32);
            auto r2_v = hn::Zero(d32);
            auto r3_v = hn::Zero(d32);
            auto r4_v = hn::Zero(d32);

            if constexpr (kUseCachedOffsets) {
                gather_cached_reference_vectors(
                    d32,
                    src_row + col,
                    off1,
                    off2,
                    off3,
                    off4,
                    16 - config.input_depth,
                    lanes,
                    r1_v,
                    r2_v,
                    r3_v,
                    r4_v
                );
            } else {
                gather_uncached_reference_vectors<kSampleMode>(
                    d32,
                    params,
                    src_plane,
                    row,
                    col,
                    lanes,
                    r1_v,
                    r2_v,
                    r3_v,
                    r4_v
                );
            }

            processed_v = process_reference_samples<kSampleMode, kBlurFirst>(
                d32,
                v_src_up,
                r1_v,
                r2_v,
                r3_v,
                r4_v,
                config.threshold,
                config.threshold1,
                config.threshold2
            );
        }

        if constexpr (kDitherAlgo == DA_HIGH_FLOYD_STEINBERG_DITHERING) {
            alignas(64) std::int32_t fs_buf[hn::MaxLanes(d32)];
            hn::StoreU(processed_v, d32, fs_buf);
            for (std::size_t lane = 0; lane < lanes; ++lane) {
                const int column = col + static_cast<int>(lane);
                const int pixel = postprocess_floyd_steinberg_pixel(
                    config,
                    *fs_dither,
                    fs_buf[lane],
                    grain_row,
                    column
                );
                dst_row[static_cast<std::size_t>(column)] = static_cast<PixelOut>(pixel);
                fs_dither->next_pixel();
            }
        } else {
            auto pixel = apply_dither_and_grain<kDitherAlgo>(config, d32, processed_v, grain_row, row, col, lanes);
            store_pixels_from_u32(du32, dst_row.data() + col, hn::BitCast(du32, pixel));
        }
    }
}
