#pragma once

#include "core/kernel.hpp"

#include "core/dither_none.hpp"
#include "core/dither_ordered.hpp"
#include "core/dither_floyd_steinberg.hpp"

#include "core/output_16bit.hpp"

namespace neo_f3kdb::core::pixel_proc {

namespace adapter {

template <
	void (*InitContext)(char*, int, int),
	void (*DestroyContext)(void*),
	void (*NextPixel)(void*),
	void (*NextRow)(void*),
	int (*Upsample)(void*, unsigned char),
	int (*Downsample)(void*, int, int, int, int, int, int),
	int (*Avg2)(void*, int, int),
	int (*Avg4)(void*, int, int, int, int)>
struct FunctionTable {
	static inline void init_context(char context_buffer[CONTEXT_BUFFER_SIZE], int frame_width, int output_depth) {
		InitContext(context_buffer, frame_width, output_depth);
	}

	static inline void destroy_context(void* context) {
		DestroyContext(context);
	}

	static inline void next_pixel(void* context) {
		NextPixel(context);
	}

	static inline void next_row(void* context) {
		NextRow(context);
	}

	static inline int upsample(void* context, unsigned char pixel) {
		return Upsample(context, pixel);
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
		return Downsample(context, pixel, row, column, pixel_min, pixel_max, output_depth);
	}

	static inline int avg_2(void* context, int pixel1, int pixel2) {
		return Avg2(context, pixel1, pixel2);
	}

	static inline int avg_4(void* context, int pixel1, int pixel2, int pixel3, int pixel4) {
		return Avg4(context, pixel1, pixel2, pixel3, pixel4);
	}
};

} // namespace adapter

template <int mode>
struct Impl;

template <>
struct Impl<DA_HIGH_NO_DITHERING> : adapter::FunctionTable<
	neo_f3kdb::core::pixel_proc_detail::none::init_context,
	neo_f3kdb::core::pixel_proc_detail::none::destroy_context,
	neo_f3kdb::core::pixel_proc_detail::none::next_pixel,
	neo_f3kdb::core::pixel_proc_detail::none::next_row,
	neo_f3kdb::core::pixel_proc_detail::none::upsample,
	neo_f3kdb::core::pixel_proc_detail::none::downsample,
	neo_f3kdb::core::pixel_proc_detail::none::avg_2,
	neo_f3kdb::core::pixel_proc_detail::none::avg_4> {};

template <>
struct Impl<DA_HIGH_ORDERED_DITHERING> : adapter::FunctionTable<
	neo_f3kdb::core::pixel_proc_detail::ordered::init_context,
	neo_f3kdb::core::pixel_proc_detail::ordered::destroy_context,
	neo_f3kdb::core::pixel_proc_detail::ordered::next_pixel,
	neo_f3kdb::core::pixel_proc_detail::ordered::next_row,
	neo_f3kdb::core::pixel_proc_detail::ordered::upsample,
	neo_f3kdb::core::pixel_proc_detail::ordered::downsample,
	neo_f3kdb::core::pixel_proc_detail::ordered::avg_2,
	neo_f3kdb::core::pixel_proc_detail::ordered::avg_4> {};

template <>
struct Impl<DA_HIGH_FLOYD_STEINBERG_DITHERING> : adapter::FunctionTable<
	neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::init_context,
	neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::destroy_context,
	neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::next_pixel,
	neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::next_row,
	neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::upsample,
	neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::downsample,
	neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::avg_2,
	neo_f3kdb::core::pixel_proc_detail::floyd_steinberg::avg_4> {};

template <>
struct Impl<DA_16BIT_INTERLEAVED> : adapter::FunctionTable<
	neo_f3kdb::core::pixel_proc_detail::output_16bit::init_context,
	neo_f3kdb::core::pixel_proc_detail::output_16bit::destroy_context,
	neo_f3kdb::core::pixel_proc_detail::output_16bit::next_pixel,
	neo_f3kdb::core::pixel_proc_detail::output_16bit::next_row,
	neo_f3kdb::core::pixel_proc_detail::output_16bit::upsample,
	neo_f3kdb::core::pixel_proc_detail::output_16bit::downsample,
	neo_f3kdb::core::pixel_proc_detail::output_16bit::avg_2,
	neo_f3kdb::core::pixel_proc_detail::output_16bit::avg_4> {};

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
