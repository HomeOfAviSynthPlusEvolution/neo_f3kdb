#pragma once

#include "core/kernel.hpp"

#include "core/dither_none.hpp"
#include "core/dither_ordered.hpp"
#include "core/dither_floyd_steinberg.hpp"

#include "core/output_16bit.hpp"

namespace neo_f3kdb::core::pixel_proc {

template <int mode>
struct Impl;

template <>
struct Impl<DA_HIGH_NO_DITHERING> {
	static inline void init_context(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth) {
		neo_f3kdb::core::pixel_proc_detail::none::init_context(context_buffer, frame_width, output_depth);
	}

	static inline void destroy_context(void* context) {
		neo_f3kdb::core::pixel_proc_detail::none::destroy_context(context);
	}

	static inline void next_pixel(void* context) {
		neo_f3kdb::core::pixel_proc_detail::none::next_pixel(context);
	}

	static inline void next_row(void* context) {
		neo_f3kdb::core::pixel_proc_detail::none::next_row(context);
	}

	static inline int upsample(void* context, unsigned char pixel) {
		return neo_f3kdb::core::pixel_proc_detail::none::upsample(context, pixel);
	}

	static inline int downsample(
		void* context,
		int pixel,
		int row,
		int column,
		int pixel_min,
		int pixel_max,
		int output_depth
	) {
		return neo_f3kdb::core::pixel_proc_detail::none::downsample(
			context,
			pixel,
			row,
			column,
			pixel_min,
			pixel_max,
			output_depth
		);
	}

	static inline int avg_2(void* context, int pixel1, int pixel2) {
		return neo_f3kdb::core::pixel_proc_detail::none::avg_2(context, pixel1, pixel2);
	}

	static inline int avg_4(void* context, int pixel1, int pixel2, int pixel3, int pixel4) {
		return neo_f3kdb::core::pixel_proc_detail::none::avg_4(context, pixel1, pixel2, pixel3, pixel4);
	}
};

template <>
struct Impl<DA_HIGH_ORDERED_DITHERING> {
	static inline void init_context(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth) {
		neo_f3kdb::core::pixel_proc_detail::ordered::init_context(context_buffer, frame_width, output_depth);
	}

	static inline void destroy_context(void* context) {
		neo_f3kdb::core::pixel_proc_detail::ordered::destroy_context(context);
	}

	static inline void next_pixel(void* context) {
		neo_f3kdb::core::pixel_proc_detail::ordered::next_pixel(context);
	}

	static inline void next_row(void* context) {
		neo_f3kdb::core::pixel_proc_detail::ordered::next_row(context);
	}

	static inline int upsample(void* context, unsigned char pixel) {
		return neo_f3kdb::core::pixel_proc_detail::ordered::upsample(context, pixel);
	}

	static inline int downsample(
		void* context,
		int pixel,
		int row,
		int column,
		int pixel_min,
		int pixel_max,
		int output_depth
	) {
		return neo_f3kdb::core::pixel_proc_detail::ordered::downsample(
			context,
			pixel,
			row,
			column,
			pixel_min,
			pixel_max,
			output_depth
		);
	}

	static inline int avg_2(void* context, int pixel1, int pixel2) {
		return neo_f3kdb::core::pixel_proc_detail::ordered::avg_2(context, pixel1, pixel2);
	}

	static inline int avg_4(void* context, int pixel1, int pixel2, int pixel3, int pixel4) {
		return neo_f3kdb::core::pixel_proc_detail::ordered::avg_4(context, pixel1, pixel2, pixel3, pixel4);
	}
};

template <>
struct Impl<DA_HIGH_FLOYD_STEINBERG_DITHERING> {
	static inline void init_context(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth) {
		neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::init_context(context_buffer, frame_width, output_depth);
	}

	static inline void destroy_context(void* context) {
		neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::destroy_context(context);
	}

	static inline void next_pixel(void* context) {
		neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::next_pixel(context);
	}

	static inline void next_row(void* context) {
		neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::next_row(context);
	}

	static inline int upsample(void* context, unsigned char pixel) {
		return neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::upsample(context, pixel);
	}

	static inline int downsample(
		void* context,
		int pixel,
		int row,
		int column,
		int pixel_min,
		int pixel_max,
		int output_depth
	) {
		return neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::downsample(
			context,
			pixel,
			row,
			column,
			pixel_min,
			pixel_max,
			output_depth
		);
	}

	static inline int avg_2(void* context, int pixel1, int pixel2) {
		return neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::avg_2(context, pixel1, pixel2);
	}

	static inline int avg_4(void* context, int pixel1, int pixel2, int pixel3, int pixel4) {
		return neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::avg_4(context, pixel1, pixel2, pixel3, pixel4);
	}
};

template <>
struct Impl<DA_16BIT_INTERLEAVED> {
	static inline void init_context(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth) {
		neo_f3kdb::core::pixel_proc_detail::output_16bit::init_context(context_buffer, frame_width, output_depth);
	}

	static inline void destroy_context(void* context) {
		neo_f3kdb::core::pixel_proc_detail::output_16bit::destroy_context(context);
	}

	static inline void next_pixel(void* context) {
		neo_f3kdb::core::pixel_proc_detail::output_16bit::next_pixel(context);
	}

	static inline void next_row(void* context) {
		neo_f3kdb::core::pixel_proc_detail::output_16bit::next_row(context);
	}

	static inline int upsample(void* context, unsigned char pixel) {
		return neo_f3kdb::core::pixel_proc_detail::output_16bit::upsample(context, pixel);
	}

	static inline int downsample(
		void* context,
		int pixel,
		int row,
		int column,
		int pixel_min,
		int pixel_max,
		int output_depth
	) {
		return neo_f3kdb::core::pixel_proc_detail::output_16bit::downsample(
			context,
			pixel,
			row,
			column,
			pixel_min,
			pixel_max,
			output_depth
		);
	}

	static inline int avg_2(void* context, int pixel1, int pixel2) {
		return neo_f3kdb::core::pixel_proc_detail::output_16bit::avg_2(context, pixel1, pixel2);
	}

	static inline int avg_4(void* context, int pixel1, int pixel2, int pixel3, int pixel4) {
		return neo_f3kdb::core::pixel_proc_detail::output_16bit::avg_4(context, pixel1, pixel2, pixel3, pixel4);
	}
};

template <int mode>
static inline void init_context(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth)
{
	Impl<mode>::init_context(context_buffer, frame_width, output_depth);
}

template <int mode>
static inline void destroy_context(void* context)
{
	Impl<mode>::destroy_context(context);
}

template <int mode>
static inline void next_pixel(void* context)
{
	Impl<mode>::next_pixel(context);
}

template <int mode>
static inline void next_row(void* context)
{
	Impl<mode>::next_row(context);
}

template <int mode>
static inline int upsample(void* context, unsigned char pixel)
{
	return Impl<mode>::upsample(context, pixel);
}

template <int mode>
static inline int downsample(void* context, int pixel, int row, int column, int pixel_min, int pixel_max, int output_depth)
{
	return Impl<mode>::downsample(context, pixel, row, column, pixel_min, pixel_max, output_depth);
}

template <int mode>
static inline int avg_2(void* context, int pixel1, int pixel2)
{
	return Impl<mode>::avg_2(context, pixel1, pixel2);
}

template <int mode>
static inline int avg_4(void* context, int pixel1, int pixel2, int pixel3, int pixel4)
{
	return Impl<mode>::avg_4(context, pixel1, pixel2, pixel3, pixel4);
}

} // namespace neo_f3kdb::core::pixel_proc
