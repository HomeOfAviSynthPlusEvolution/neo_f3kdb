#include "core.h"

#define IMPL_DISPATCH_IMPORT_DECLARATION

#include "impl_dispatch_decl.h"

const process_plane_impl_t* process_plane_impl_high_precision_no_dithering[] = {
    process_plane_impl_c_high_no_dithering,
    process_plane_impl_c_high_no_dithering,
    process_plane_impl_c_high_no_dithering,
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
    process_plane_impl_sse4_high_no_dithering,
    process_plane_impl_avx2_high_no_dithering,
    process_plane_impl_avx512_high_no_dithering,
#elif defined(__arm__) || defined(__aarch64__) || defined(_M_ARM)
    process_plane_impl_sse4_high_no_dithering,
    nullptr,
    nullptr,
#else
    nullptr,
    nullptr,
    nullptr,
#endif
};

const process_plane_impl_t* process_plane_impl_high_precision_ordered_dithering[] = {
    process_plane_impl_c_high_ordered_dithering,
    process_plane_impl_c_high_ordered_dithering,
    process_plane_impl_c_high_ordered_dithering,
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
    process_plane_impl_sse4_high_ordered_dithering,
    process_plane_impl_avx2_high_ordered_dithering,
    process_plane_impl_avx512_high_ordered_dithering,
#elif defined(__arm__) || defined(__aarch64__) || defined(_M_ARM)
    process_plane_impl_sse4_high_ordered_dithering,
    nullptr,
    nullptr,
#else
    nullptr,
    nullptr,
    nullptr,
#endif
};

const process_plane_impl_t* process_plane_impl_high_precision_floyd_steinberg_dithering[] = {
    process_plane_impl_c_high_floyd_steinberg_dithering,
    process_plane_impl_c_high_floyd_steinberg_dithering,
    process_plane_impl_c_high_floyd_steinberg_dithering,
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
    process_plane_impl_sse4_high_floyd_steinberg_dithering,
    process_plane_impl_avx2_high_floyd_steinberg_dithering,
    process_plane_impl_avx512_high_floyd_steinberg_dithering,
#elif defined(__arm__) || defined(__aarch64__) || defined(_M_ARM)
    process_plane_impl_sse4_high_floyd_steinberg_dithering,
    nullptr,
    nullptr,
#else
    nullptr,
    nullptr,
    nullptr,
#endif
};

const process_plane_impl_t* process_plane_impl_16bit_interleaved[] = {
    process_plane_impl_c_16bit_interleaved,
    process_plane_impl_c_16bit_interleaved,
    process_plane_impl_c_16bit_interleaved,
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
    process_plane_impl_sse4_16bit_interleaved,
    process_plane_impl_avx2_16bit_interleaved,
    process_plane_impl_avx512_16bit_interleaved,
#elif defined(__arm__) || defined(__aarch64__) || defined(_M_ARM)
    process_plane_impl_sse4_16bit_interleaved,
    nullptr,
    nullptr,
#else
    nullptr,
    nullptr,
    nullptr,
#endif
};


const process_plane_impl_t** process_plane_impls[] = {
	nullptr, // process_plane_impl_low_precision has been removed,
	process_plane_impl_high_precision_no_dithering,
	process_plane_impl_high_precision_ordered_dithering,
	process_plane_impl_high_precision_floyd_steinberg_dithering,
    process_plane_impl_16bit_interleaved
};
