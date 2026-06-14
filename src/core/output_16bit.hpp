#include <assert.h>

#include "core/constants.hpp"
#include "core/pixel_proc_common.hpp"

namespace neo_f3kdb::core::pixel_proc_detail::output_16bit {
    
    static inline void init_context(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth)
    {
        // sanity check only
        assert(output_depth == 16);
    }

    static inline void destroy_context(void* context)
    {
        // nothing to do
    }

    static inline void next_pixel(void* context)
    {
        // nothing to do
    }

    static inline void next_row(void* context)
    {
        // nothing to do
    }

    static inline int dither(void* context, int pixel, int row, int column)
    {
        return pixel;
    }

    static inline int upsample(void* context, unsigned char pixel)
    {
        return neo_f3kdb::core::pixel_proc_common::upsample(context, pixel);
    }

    static inline int downsample(void* context, int pixel, int row, int column, int pixel_min, int pixel_max, int output_depth)
    {
        assert(output_depth == 16);
        // I know the method name is totally wrong...
        return clamp_pixel(pixel, pixel_min, pixel_max) << (output_depth - INTERNAL_BIT_DEPTH);
    }

    static inline int avg_2(void* context, int pixel1, int pixel2)
    {
        return neo_f3kdb::core::pixel_proc_common::avg_2(context, pixel1, pixel2);
    }

    static inline int avg_4(void* context, int pixel1, int pixel2, int pixel3, int pixel4)
    {
        return neo_f3kdb::core::pixel_proc_common::avg_4(context, pixel1, pixel2, pixel3, pixel4);
    }

} // namespace neo_f3kdb::core::pixel_proc_detail::output_16bit
