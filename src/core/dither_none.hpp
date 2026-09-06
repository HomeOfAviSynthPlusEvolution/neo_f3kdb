

#include "core/kernel.hpp"
#include "core/pixel_proc_common.hpp"

namespace neo_f3kdb::core::pixel_proc_detail::none {
	
	inline void init_context(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth)
	{
		// nothing to do
	}

	inline void destroy_context(void* context)
	{
		// nothing to do
	}

	inline void next_pixel(void* context)
	{
		// nothing to do
	}

	inline void next_row(void* context)
	{
		// nothing to do
	}

	inline int dither(void* context, int pixel, int row, int column)
	{
		return pixel;
	}

    inline int upsample(void* context, unsigned char pixel)
    {
        return neo_f3kdb::core::pixel_proc_common::upsample(context, pixel);
    }

    inline int downsample(void* context, int pixel, int row, int column, int pixel_min, int pixel_max, int output_depth)
    {
        pixel = dither(context, pixel, row, column);
        return neo_f3kdb::core::pixel_proc_common::downsample(pixel, pixel_min, pixel_max, output_depth);
    }

    inline int avg_2(void* context, int pixel1, int pixel2)
    {
        return neo_f3kdb::core::pixel_proc_common::avg_2(context, pixel1, pixel2);
    }

    inline int avg_4(void* context, int pixel1, int pixel2, int pixel3, int pixel4)
    {
        return neo_f3kdb::core::pixel_proc_common::avg_4(context, pixel1, pixel2, pixel3, pixel4);
    }

} // namespace neo_f3kdb::core::pixel_proc_detail::none
