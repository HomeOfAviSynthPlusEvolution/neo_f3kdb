#pragma once

#include "core/kernel.hpp"

#include "core/dither_none.hpp"
#include "core/dither_ordered.hpp"
#include "core/dither_floyd_steinberg.hpp"

#include "core/output_16bit.hpp"

template <int mode>
struct pixel_proc_impl;

template <>
struct pixel_proc_impl<DA_HIGH_NO_DITHERING> {
	static inline void init_context(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth) {
		pixel_proc_high_no_dithering::init_context(context_buffer, frame_width, output_depth);
	}

	static inline void destroy_context(void* context) {
		pixel_proc_high_no_dithering::destroy_context(context);
	}

	static inline void next_pixel(void* context) {
		pixel_proc_high_no_dithering::next_pixel(context);
	}

	static inline void next_row(void* context) {
		pixel_proc_high_no_dithering::next_row(context);
	}

	static inline int upsample(void* context, unsigned char pixel) {
		return pixel_proc_high_no_dithering::upsample(context, pixel);
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
		return pixel_proc_high_no_dithering::downsample(
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
		return pixel_proc_high_no_dithering::avg_2(context, pixel1, pixel2);
	}

	static inline int avg_4(void* context, int pixel1, int pixel2, int pixel3, int pixel4) {
		return pixel_proc_high_no_dithering::avg_4(context, pixel1, pixel2, pixel3, pixel4);
	}
};

template <>
struct pixel_proc_impl<DA_HIGH_ORDERED_DITHERING> {
	static inline void init_context(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth) {
		pixel_proc_high_ordered_dithering::init_context(context_buffer, frame_width, output_depth);
	}

	static inline void destroy_context(void* context) {
		pixel_proc_high_ordered_dithering::destroy_context(context);
	}

	static inline void next_pixel(void* context) {
		pixel_proc_high_ordered_dithering::next_pixel(context);
	}

	static inline void next_row(void* context) {
		pixel_proc_high_ordered_dithering::next_row(context);
	}

	static inline int upsample(void* context, unsigned char pixel) {
		return pixel_proc_high_ordered_dithering::upsample(context, pixel);
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
		return pixel_proc_high_ordered_dithering::downsample(
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
		return pixel_proc_high_ordered_dithering::avg_2(context, pixel1, pixel2);
	}

	static inline int avg_4(void* context, int pixel1, int pixel2, int pixel3, int pixel4) {
		return pixel_proc_high_ordered_dithering::avg_4(context, pixel1, pixel2, pixel3, pixel4);
	}
};

template <>
struct pixel_proc_impl<DA_HIGH_FLOYD_STEINBERG_DITHERING> {
	static inline void init_context(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth) {
		pixel_proc_high_f_s_dithering::init_context(context_buffer, frame_width, output_depth);
	}

	static inline void destroy_context(void* context) {
		pixel_proc_high_f_s_dithering::destroy_context(context);
	}

	static inline void next_pixel(void* context) {
		pixel_proc_high_f_s_dithering::next_pixel(context);
	}

	static inline void next_row(void* context) {
		pixel_proc_high_f_s_dithering::next_row(context);
	}

	static inline int upsample(void* context, unsigned char pixel) {
		return pixel_proc_high_f_s_dithering::upsample(context, pixel);
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
		return pixel_proc_high_f_s_dithering::downsample(
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
		return pixel_proc_high_f_s_dithering::avg_2(context, pixel1, pixel2);
	}

	static inline int avg_4(void* context, int pixel1, int pixel2, int pixel3, int pixel4) {
		return pixel_proc_high_f_s_dithering::avg_4(context, pixel1, pixel2, pixel3, pixel4);
	}
};

template <>
struct pixel_proc_impl<DA_16BIT_INTERLEAVED> {
	static inline void init_context(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth) {
		pixel_proc_16bit::init_context(context_buffer, frame_width, output_depth);
	}

	static inline void destroy_context(void* context) {
		pixel_proc_16bit::destroy_context(context);
	}

	static inline void next_pixel(void* context) {
		pixel_proc_16bit::next_pixel(context);
	}

	static inline void next_row(void* context) {
		pixel_proc_16bit::next_row(context);
	}

	static inline int upsample(void* context, unsigned char pixel) {
		return pixel_proc_16bit::upsample(context, pixel);
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
		return pixel_proc_16bit::downsample(
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
		return pixel_proc_16bit::avg_2(context, pixel1, pixel2);
	}

	static inline int avg_4(void* context, int pixel1, int pixel2, int pixel3, int pixel4) {
		return pixel_proc_16bit::avg_4(context, pixel1, pixel2, pixel3, pixel4);
	}
};

template <int mode>
static inline void pixel_proc_init_context(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth)
{
	pixel_proc_impl<mode>::init_context(context_buffer, frame_width, output_depth);
}

template <int mode>
static inline void pixel_proc_destroy_context(void* context)
{
	pixel_proc_impl<mode>::destroy_context(context);
}

template <int mode>
static inline void pixel_proc_next_pixel(void* context)
{
	pixel_proc_impl<mode>::next_pixel(context);
}

template <int mode>
static inline void pixel_proc_next_row(void* context)
{
	pixel_proc_impl<mode>::next_row(context);
}

template <int mode>
static inline int pixel_proc_upsample(void* context, unsigned char pixel)
{
	return pixel_proc_impl<mode>::upsample(context, pixel);
}

template <int mode>
static inline int pixel_proc_downsample(void* context, int pixel, int row, int column, int pixel_min, int pixel_max, int output_depth)
{
	return pixel_proc_impl<mode>::downsample(context, pixel, row, column, pixel_min, pixel_max, output_depth);
}

template <int mode>
static inline int pixel_proc_avg_2(void* context, int pixel1, int pixel2)
{
	return pixel_proc_impl<mode>::avg_2(context, pixel1, pixel2);
}

template <int mode>
static inline int pixel_proc_avg_4(void* context, int pixel1, int pixel2, int pixel3, int pixel4)
{
	return pixel_proc_impl<mode>::avg_4(context, pixel1, pixel2, pixel3, pixel4);
}
