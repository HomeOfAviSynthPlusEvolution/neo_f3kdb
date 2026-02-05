#include "core.h"

#define IMPL_DISPATCH_IMPORT_DECLARATION

#include "impl_dispatch_decl.h"

#ifdef HAS_SSE4
  #define SSE4_IMPL(name) process_plane_impl_sse4_##name
#else
  #define SSE4_IMPL(name) process_plane_impl_c_##name
#endif

#ifdef HAS_AVX2
  #define AVX2_IMPL(name) process_plane_impl_avx2_##name
#else
  #define AVX2_IMPL(name) process_plane_impl_c_##name
#endif

#ifdef HAS_AVX512
  #define AVX512_IMPL(name) process_plane_impl_avx512_##name
#else
  #define AVX512_IMPL(name) process_plane_impl_c_##name
#endif

const process_plane_impl_t* process_plane_impl_high_precision_no_dithering[] = {
    process_plane_impl_c_high_no_dithering,
    process_plane_impl_c_high_no_dithering,
    process_plane_impl_c_high_no_dithering,
    SSE4_IMPL(high_no_dithering),
    AVX2_IMPL(high_no_dithering),
    AVX512_IMPL(high_no_dithering)
};

const process_plane_impl_t* process_plane_impl_high_precision_ordered_dithering[] = {
    process_plane_impl_c_high_ordered_dithering,
    process_plane_impl_c_high_ordered_dithering,
    process_plane_impl_c_high_ordered_dithering,
    SSE4_IMPL(high_ordered_dithering),
    AVX2_IMPL(high_ordered_dithering),
    AVX512_IMPL(high_ordered_dithering)
};

const process_plane_impl_t* process_plane_impl_high_precision_floyd_steinberg_dithering[] = {
    process_plane_impl_c_high_floyd_steinberg_dithering,
    process_plane_impl_c_high_floyd_steinberg_dithering,
    process_plane_impl_c_high_floyd_steinberg_dithering,
    SSE4_IMPL(high_floyd_steinberg_dithering),
    AVX2_IMPL(high_floyd_steinberg_dithering),
    AVX512_IMPL(high_floyd_steinberg_dithering)
};

const process_plane_impl_t* process_plane_impl_16bit_interleaved[] = {
    process_plane_impl_c_16bit_interleaved,
    process_plane_impl_c_16bit_interleaved,
    process_plane_impl_c_16bit_interleaved,
    SSE4_IMPL(16bit_interleaved),
    AVX2_IMPL(16bit_interleaved),
    AVX512_IMPL(16bit_interleaved)
};


const process_plane_impl_t** process_plane_impls[] = {
	nullptr, // process_plane_impl_low_precision has been removed,
	process_plane_impl_high_precision_no_dithering,
	process_plane_impl_high_precision_ordered_dithering,
	process_plane_impl_high_precision_floyd_steinberg_dithering,
    process_plane_impl_16bit_interleaved
};
