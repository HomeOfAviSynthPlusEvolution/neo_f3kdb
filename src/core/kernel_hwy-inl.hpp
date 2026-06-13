// Intentionally no include guard: this file is included once per Highway target.

namespace deband_hwy_detail {

template <class D32, class Pixel, class V32>
void store_pixels_from_u32(D32 d32, Pixel* dest_ptr, V32 values) {
    if constexpr (std::is_same_v<Pixel, std::uint8_t>) {
        const hn::Rebind<std::uint16_t, D32> d16;
        const hn::Rebind<std::uint8_t, D32> d8;
        const auto out16 = hn::DemoteTo(d16, values);
        const auto out8 = hn::DemoteTo(d8, out16);
        hn::StoreU(out8, d8, dest_ptr);
    } else {
        const hn::Rebind<std::uint16_t, D32> d16;
        const auto out16 = hn::DemoteTo(d16, values);
        hn::StoreU(out16, d16, dest_ptr);
    }
}

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
    const int height = params.plane_height_in_pixels;
    const int width = params.plane_width_in_pixels;
    const int input_depth = params.input_depth;
    const int width_subsamp = params.width_subsampling;
    const int height_subsamp = params.height_subsampling;

    for (std::size_t lane = 0; lane < lanes; ++lane) {
        const int pixel_col = col + static_cast<int>(lane);
        const pixel_dither_info& info = info_row[pixel_col];
        src_up[lane] = static_cast<std::int32_t>(src_row[pixel_col]) << (16 - input_depth);

        if constexpr (kSampleMode == 1 || kSampleMode == 3 || kSampleMode == 4 || kSampleMode == 5) {
            const int ref_y = info.ref1 >> height_subsamp;
            const int ref_x = info.ref1 >> width_subsamp;

            if constexpr (kSampleMode == 1) {
                ref1[lane] = static_cast<std::int32_t>(
                    src_base[std::clamp(row - ref_y, 0, height - 1) * src_stride + pixel_col]
                ) << (16 - input_depth);
                ref2[lane] = static_cast<std::int32_t>(
                    src_base[std::clamp(row + ref_y, 0, height - 1) * src_stride + pixel_col]
                ) << (16 - input_depth);
            } else if constexpr (kSampleMode == 3) {
                ref1[lane] = static_cast<std::int32_t>(
                    src_row[std::clamp(pixel_col - ref_x, 0, width - 1)]
                ) << (16 - input_depth);
                ref2[lane] = static_cast<std::int32_t>(
                    src_row[std::clamp(pixel_col + ref_x, 0, width - 1)]
                ) << (16 - input_depth);
            } else {
                ref1[lane] = static_cast<std::int32_t>(
                    src_base[std::clamp(row - ref_y, 0, height - 1) * src_stride + pixel_col]
                ) << (16 - input_depth);
                ref2[lane] = static_cast<std::int32_t>(
                    src_base[std::clamp(row + ref_y, 0, height - 1) * src_stride + pixel_col]
                ) << (16 - input_depth);
                ref3[lane] = static_cast<std::int32_t>(
                    src_row[std::clamp(pixel_col - ref_x, 0, width - 1)]
                ) << (16 - input_depth);
                ref4[lane] = static_cast<std::int32_t>(
                    src_row[std::clamp(pixel_col + ref_x, 0, width - 1)]
                ) << (16 - input_depth);
            }
        } else if constexpr (kSampleMode == 2) {
            const int ref_y1 = info.ref2 >> height_subsamp;
            const int ref_y2 = info.ref1 >> height_subsamp;
            const int ref_x1 = info.ref1 >> width_subsamp;
            const int ref_x2 = info.ref2 >> width_subsamp;

            ref1[lane] = static_cast<std::int32_t>(
                src_base[std::clamp(row - ref_y1, 0, height - 1) * src_stride + std::clamp(pixel_col - ref_x1, 0, width - 1)]
            ) << (16 - input_depth);
            ref2[lane] = static_cast<std::int32_t>(
                src_base[std::clamp(row - ref_y2, 0, height - 1) * src_stride + std::clamp(pixel_col + ref_x2, 0, width - 1)]
            ) << (16 - input_depth);
            ref3[lane] = static_cast<std::int32_t>(
                src_base[std::clamp(row + ref_y1, 0, height - 1) * src_stride + std::clamp(pixel_col + ref_x1, 0, width - 1)]
            ) << (16 - input_depth);
            ref4[lane] = static_cast<std::int32_t>(
                src_base[std::clamp(row + ref_y2, 0, height - 1) * src_stride + std::clamp(pixel_col - ref_x2, 0, width - 1)]
            ) << (16 - input_depth);
        }
    }
}

template <int kSampleMode, bool kBlurFirst, class D32, class V32>
V32 process_refs(D32 d32, V32 src, const std::int32_t* ref1, const std::int32_t* ref2, const std::int32_t* ref3, const std::int32_t* ref4, int threshold, int threshold1, int threshold2) {
    const auto one = hn::Set(d32, 1);
    const auto threshold_v = hn::Set(d32, static_cast<std::int32_t>(threshold));

    if constexpr (kSampleMode == 1 || kSampleMode == 3) {
        const auto r1 = hn::LoadU(d32, ref1);
        const auto r2 = hn::LoadU(d32, ref2);
        const auto avg = hn::ShiftRight<1>(hn::Add(hn::Add(r1, r2), one));
        const auto use_src = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg, src), threshold_v)
            : hn::Or(hn::Ge(hn::AbsDiff(src, r1), threshold_v), hn::Ge(hn::AbsDiff(src, r2), threshold_v));
        return hn::IfThenElse(use_src, src, avg);
    } else if constexpr (kSampleMode == 2) {
        const auto r1 = hn::LoadU(d32, ref1);
        const auto r2 = hn::LoadU(d32, ref2);
        const auto r3 = hn::LoadU(d32, ref3);
        const auto r4 = hn::LoadU(d32, ref4);
        auto avg1 = hn::ShiftRight<1>(hn::Add(hn::Add(r1, r2), one));
        const auto avg2 = hn::ShiftRight<1>(hn::Add(hn::Add(r3, r4), one));
        avg1 = hn::IfThenElse(hn::Gt(avg1, hn::Zero(d32)), hn::Sub(avg1, one), avg1);
        const auto avg = hn::ShiftRight<1>(hn::Add(hn::Add(avg1, avg2), one));
        const auto use_src = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg, src), threshold_v)
            : hn::Or(
                hn::Or(hn::Ge(hn::AbsDiff(r1, src), threshold_v), hn::Ge(hn::AbsDiff(r2, src), threshold_v)),
                hn::Or(hn::Ge(hn::AbsDiff(r3, src), threshold_v), hn::Ge(hn::AbsDiff(r4, src), threshold_v))
            );
        return hn::IfThenElse(use_src, src, avg);
    } else if constexpr (kSampleMode == 4) {
        const auto ref_v1 = hn::LoadU(d32, ref1);
        const auto ref_v2 = hn::LoadU(d32, ref2);
        const auto ref_h1 = hn::LoadU(d32, ref3);
        const auto ref_h2 = hn::LoadU(d32, ref4);
        const auto avg_v = hn::ShiftRight<1>(hn::Add(hn::Add(ref_v1, ref_v2), one));
        const auto avg_h = hn::ShiftRight<1>(hn::Add(hn::Add(ref_h1, ref_h2), one));
        const auto use_src_v = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg_v, src), threshold_v)
            : hn::Or(hn::Ge(hn::AbsDiff(src, ref_v1), threshold_v), hn::Ge(hn::AbsDiff(src, ref_v2), threshold_v));
        const auto use_src_h = kBlurFirst
            ? hn::Ge(hn::AbsDiff(avg_h, src), threshold_v)
            : hn::Or(hn::Ge(hn::AbsDiff(src, ref_h1), threshold_v), hn::Ge(hn::AbsDiff(src, ref_h2), threshold_v));
        const auto new_v = hn::IfThenElse(use_src_v, src, avg_v);
        const auto new_h = hn::IfThenElse(use_src_h, src, avg_h);
        return hn::ShiftRight<1>(hn::Add(hn::Add(new_v, new_h), one));
    } else {
        static_assert(kSampleMode == 5);
        const auto r1 = hn::LoadU(d32, ref1);
        const auto r2 = hn::LoadU(d32, ref2);
        const auto r3 = hn::LoadU(d32, ref3);
        const auto r4 = hn::LoadU(d32, ref4);
        auto avg1 = hn::ShiftRight<1>(hn::Add(hn::Add(r1, r2), one));
        const auto avg2 = hn::ShiftRight<1>(hn::Add(hn::Add(r3, r4), one));
        avg1 = hn::IfThenElse(hn::Gt(avg1, hn::Zero(d32)), hn::Sub(avg1, one), avg1);
        const auto avg = hn::ShiftRight<1>(hn::Add(hn::Add(avg1, avg2), one));
        const auto threshold1_v = hn::Set(d32, static_cast<std::int32_t>(threshold1));
        const auto threshold2_v = hn::Set(d32, static_cast<std::int32_t>(threshold2));
        const auto max_diff = hn::Max(
            hn::AbsDiff(r1, src),
            hn::Max(hn::AbsDiff(r2, src), hn::Max(hn::AbsDiff(r3, src), hn::AbsDiff(r4, src)))
        );
        const auto use_src = hn::Or(
            hn::Or(hn::Ge(hn::AbsDiff(avg, src), threshold_v), hn::Ge(max_diff, threshold1_v)),
            hn::Or(
                hn::Ge(hn::AbsDiff(hn::Add(r1, r2), hn::ShiftLeft<1>(src)), threshold2_v),
                hn::Ge(hn::AbsDiff(hn::Add(r3, r4), hn::ShiftLeft<1>(src)), threshold2_v)
            )
        );
        return hn::IfThenElse(use_src, src, avg);
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

    auto pixel = process_refs<kSampleMode, kBlurFirst>(
        d32,
        hn::LoadU(d32, src_up),
        ref1,
        ref2,
        ref3,
        ref4,
        params.threshold,
        params.threshold1,
        params.threshold2
    );
    pixel = apply_dither_and_grain<kDitherAlgo>(params, d32, pixel, grain_row, row, col, lanes);
    store_pixels_from_u32(du32, dst_row + col, hn::BitCast(du32, pixel));
}

} // namespace deband_hwy_detail
