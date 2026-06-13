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
