#pragma once

#include "core.h"

#ifdef IMPL_DISPATCH_IMPORT_DECLARATION

#define DEFINE_IMPL(n, \
					nullptr, \
					impl_func_mode1_blur, \
					impl_func_mode1_noblur, \
					impl_func_mode2_blur, \
					impl_func_mode2_noblur, \
					impl_func_mode3_blur, \
					impl_func_mode3_noblur, \
					impl_func_mode4_blur, \
					impl_func_mode4_noblur, \
					impl_func_mode5_blur, \
					impl_func_mode5_noblur, \
					impl_func_mode6_blur, \
					impl_func_mode6_noblur, \
					impl_func_mode7_blur, \
					impl_func_mode7_noblur) \
	extern const process_plane_impl_t process_plane_impl_##n [];

#else

#define DEFINE_IMPL(n, \
					nullptr, \
					impl_func_mode1_blur, \
					impl_func_mode1_noblur, \
					impl_func_mode2_blur, \
					impl_func_mode2_noblur, \
					impl_func_mode3_blur, \
					impl_func_mode3_noblur, \
					impl_func_mode4_blur, \
					impl_func_mode4_noblur, \
					impl_func_mode5_blur, \
					impl_func_mode5_noblur, \
					impl_func_mode6_blur, \
				    impl_func_mode6_noblur, \
					impl_func_mode7_blur, \
				    impl_func_mode7_noblur) \
	extern const process_plane_impl_t process_plane_impl_##n [] = { \
					nullptr, \
					impl_func_mode1_blur, \
					impl_func_mode1_noblur, \
					impl_func_mode2_blur, \
					impl_func_mode2_noblur, \
					impl_func_mode3_blur, \
					impl_func_mode3_noblur, \
					impl_func_mode4_blur, \
					impl_func_mode4_noblur, \
					impl_func_mode5_blur, \
					impl_func_mode5_noblur, \
					impl_func_mode6_blur, \
					impl_func_mode6_noblur, \
					impl_func_mode7_blur, \
					impl_func_mode7_noblur};

#endif


#define DEFINE_TEMPLATE_IMPL(name, impl_func, ...) \
	DEFINE_IMPL(name, \
				(nullptr), \
				(&impl_func<1, true, __VA_ARGS__>), \
				(&impl_func<1, false, __VA_ARGS__>), \
				(&impl_func<2, true, __VA_ARGS__>), \
				(&impl_func<2, false, __VA_ARGS__>), \
				(&impl_func<3, true, __VA_ARGS__>), \
				(&impl_func<3, false, __VA_ARGS__>), \
				(&impl_func<4, true, __VA_ARGS__>), \
				(&impl_func<4, false, __VA_ARGS__>), \
				(&impl_func<5, true, __VA_ARGS__>), \
				(&impl_func<5, false, __VA_ARGS__>), \
				(&impl_func<6, true, __VA_ARGS__>), \
				(&impl_func<6, false, __VA_ARGS__>), \
				(&impl_func<7, true, __VA_ARGS__>), \
				(&impl_func<7, false, __VA_ARGS__>) );

#define DEFINE_SSE_IMPL(name, ...) \
	DEFINE_TEMPLATE_IMPL(name, process_plane_sse_impl, __VA_ARGS__);


#if defined(IMPL_DISPATCH_IMPORT_DECLARATION) || defined(DECLARE_IMPL_C)
	DEFINE_TEMPLATE_IMPL(c_high_no_dithering, process_plane_plainc, DA_HIGH_NO_DITHERING);
	DEFINE_TEMPLATE_IMPL(c_high_ordered_dithering, process_plane_plainc, DA_HIGH_ORDERED_DITHERING);
	DEFINE_TEMPLATE_IMPL(c_high_floyd_steinberg_dithering, process_plane_plainc, DA_HIGH_FLOYD_STEINBERG_DITHERING);
	DEFINE_TEMPLATE_IMPL(c_16bit_interleaved, process_plane_plainc, DA_16BIT_INTERLEAVED);
#endif


#if defined(IMPL_DISPATCH_IMPORT_DECLARATION) || defined(DECLARE_IMPL_SSE4)
	DEFINE_SSE_IMPL(sse4_high_no_dithering, DA_HIGH_NO_DITHERING);
	DEFINE_SSE_IMPL(sse4_high_ordered_dithering, DA_HIGH_ORDERED_DITHERING);
	DEFINE_SSE_IMPL(sse4_high_floyd_steinberg_dithering, DA_HIGH_FLOYD_STEINBERG_DITHERING);
	DEFINE_SSE_IMPL(sse4_16bit_interleaved, DA_16BIT_INTERLEAVED);
#endif

#if defined(IMPL_DISPATCH_IMPORT_DECLARATION) || defined(DECLARE_IMPL_AVX2)
#define DEFINE_AVX2_IMPL(name, ...) \
		DEFINE_TEMPLATE_IMPL(name, process_plane_avx2_impl, __VA_ARGS__);
	DEFINE_AVX2_IMPL(avx2_high_no_dithering, DA_HIGH_NO_DITHERING);
	DEFINE_AVX2_IMPL(avx2_high_ordered_dithering, DA_HIGH_ORDERED_DITHERING);
	DEFINE_AVX2_IMPL(avx2_high_floyd_steinberg_dithering, DA_HIGH_FLOYD_STEINBERG_DITHERING);
	DEFINE_AVX2_IMPL(avx2_16bit_interleaved, DA_16BIT_INTERLEAVED);
#endif

#if defined(IMPL_DISPATCH_IMPORT_DECLARATION) || defined(DECLARE_IMPL_AVX512)
#define DEFINE_AVX512_IMPL(name, ...) \
		DEFINE_TEMPLATE_IMPL(name, process_plane_avx512_impl, __VA_ARGS__);
	DEFINE_AVX512_IMPL(avx512_high_no_dithering, DA_HIGH_NO_DITHERING);
	DEFINE_AVX512_IMPL(avx512_high_ordered_dithering, DA_HIGH_ORDERED_DITHERING);
	DEFINE_AVX512_IMPL(avx512_high_floyd_steinberg_dithering, DA_HIGH_FLOYD_STEINBERG_DITHERING);
	DEFINE_AVX512_IMPL(avx512_16bit_interleaved, DA_16BIT_INTERLEAVED);
#endif
