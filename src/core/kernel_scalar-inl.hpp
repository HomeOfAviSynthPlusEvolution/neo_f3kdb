#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include "core/plane.hpp"
#include "core/pixel_proc.hpp"
#include "core/sample_modes.hpp"

inline int input_pixel_step(const process_plane_params& params)
{
    return params.config.input_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1;
}

template <int mode>
static __inline int read_pixel_value(const process_plane_params& params, void* context, const unsigned char* ptr)
{
    if (params.config.input_mode == LOW_BIT_DEPTH)
    {
        return neo_f3kdb::core::pixel_proc::upsample<mode>(context, *ptr);
    }

    int ret;

    ret = *(unsigned short*)ptr;
    ret <<= (INTERNAL_BIT_DEPTH - params.config.input_depth);
    return ret;
}

template <int mode>
static __inline int read_pixel_at(
    const process_plane_params& params,
    void* context,
    span2d::ReadOnlyRestrictPlane<unsigned char> src_plane,
    int row,
    int col
) {
    const auto* ptr = src_plane.row_ptr(row) + static_cast<intptr_t>(col) * input_pixel_step(params);
    return read_pixel_value<mode>(params, context, ptr);
}

template <int mode>
struct ScalarSampleOps {
    const process_plane_params& params;
    void* context;
    span2d::ReadOnlyRestrictPlane<unsigned char> src_plane;

    [[nodiscard]] int read(int row, int col) const {
        return read_pixel_at<mode>(params, context, src_plane, row, col);
    }

    [[nodiscard]] int avg2(int pixel1, int pixel2) const {
        return neo_f3kdb::core::pixel_proc::avg_2<mode>(context, pixel1, pixel2);
    }

    [[nodiscard]] int avg4(int pixel1, int pixel2, int pixel3, int pixel4) const {
        return neo_f3kdb::core::pixel_proc::avg_4<mode>(context, pixel1, pixel2, pixel3, pixel4);
    }
};

template <int sample_mode, bool blur_first, int mode, int output_mode>
static __forceinline void __cdecl process_plane_plainc_mode12_high(const process_plane_params& params, process_plane_context*)
{
    alignas(std::max_align_t) char context[CONTEXT_BUFFER_SIZE];

    neo_f3kdb::core::pixel_proc::init_context<mode>(context, params.plane_width(), params.config.output_depth);

    int dst_pixel_step = output_mode == HIGH_BIT_DEPTH_INTERLEAVED ? 2 : 1;
    int process_width = params.plane_width();

    const auto src_plane = params.src_bytes();
    auto dst_plane = params.dst_bytes();
    const auto grain_plane = params.grain_plane();
    const ScalarSampleOps<mode> sample_ops{params, context, src_plane};

    auto dst_cursor = dst_plane.cursor();
    auto grain_cursor = grain_plane.cursor();

    for (int i = 0; i < params.plane_height(); i++)
    {
        auto dst_row = *dst_cursor;
        const auto grain_row = *grain_cursor;

        for (int j = 0; j < process_width; j++)
        {
            const auto column = static_cast<std::size_t>(j);
            unsigned char* dst_px = dst_row.data() + static_cast<intptr_t>(j) * dst_pixel_step;
            int new_pixel = neo_f3kdb::core::sample_modes::process_pixel_with_ops<sample_mode, blur_first>(
                params,
                sample_ops,
                i,
                j
            );

            new_pixel = neo_f3kdb::core::pixel_proc::downsample<mode>(
                context,
                new_pixel + grain_row[column],
                i,
                j,
                params.config.pixel_min,
                params.config.pixel_max,
                params.config.output_depth
            );

            switch (output_mode)
            {
            case LOW_BIT_DEPTH:
                *dst_px = (unsigned char)new_pixel;
                break;
            case HIGH_BIT_DEPTH_INTERLEAVED:
                *((unsigned short*)dst_px) = (unsigned short)(new_pixel & 0xFFFF);
                break;
            default:
                abort();
            }

            neo_f3kdb::core::pixel_proc::next_pixel<mode>(context);
        }
        neo_f3kdb::core::pixel_proc::next_row<mode>(context);
        ++dst_cursor;
        ++grain_cursor;
    }

    neo_f3kdb::core::pixel_proc::destroy_context<mode>(context);
}

template <int sample_mode, bool blur_first, int mode>
void __cdecl process_plane_plainc(const process_plane_params& params, process_plane_context* context)
{
    static_assert(sample_mode != 0, "No longer support sample_mode = 0");
    switch (params.config.output_mode)
    {
    case LOW_BIT_DEPTH:
        process_plane_plainc_mode12_high<sample_mode, blur_first, mode, LOW_BIT_DEPTH>(params, context);
        break;

    case HIGH_BIT_DEPTH_INTERLEAVED:
        process_plane_plainc_mode12_high<sample_mode, blur_first, mode, HIGH_BIT_DEPTH_INTERLEAVED>(params, context);
        break;

    default:
        abort();
    }
}
