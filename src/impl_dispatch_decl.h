#pragma once

#include "core.h"

#ifdef IMPL_DISPATCH_IMPORT_DECLARATION
#define DEFINE_TEMPLATE_IMPL_C(name, impl_func, ...) \
        extern const process_plane_impl_t process_plane_impl_##name [];
#else

#define DEFINE_TEMPLATE_IMPL_C(name, impl_func, ...) \
        extern const process_plane_impl_t process_plane_impl_##name [] = { \
            nullptr, \
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
            (&impl_func<7, false, __VA_ARGS__>)};

#endif

#if defined(IMPL_DISPATCH_IMPORT_DECLARATION) || defined(DECLARE_IMPL_C)
DEFINE_TEMPLATE_IMPL_C(c_high_no_dithering, process_plane_plainc, DA_HIGH_NO_DITHERING);
DEFINE_TEMPLATE_IMPL_C(c_high_ordered_dithering, process_plane_plainc, DA_HIGH_ORDERED_DITHERING);
DEFINE_TEMPLATE_IMPL_C(c_high_floyd_steinberg_dithering, process_plane_plainc, DA_HIGH_FLOYD_STEINBERG_DITHERING);
DEFINE_TEMPLATE_IMPL_C(c_16bit_interleaved, process_plane_plainc, DA_16BIT_INTERLEAVED);
#endif

#if defined(IMPL_DISPATCH_IMPORT_DECLARATION) || defined(DECLARE_IMPL_SSE4)
#if defined(HAS_SSE4) || !defined(IMPL_DISPATCH_IMPORT_DECLARATION)
#define DEFINE_SSE_IMPL(name, ...) \
            DEFINE_TEMPLATE_IMPL_C(sse4_##name, process_plane_sse_impl, __VA_ARGS__);
DEFINE_SSE_IMPL(high_no_dithering, DA_HIGH_NO_DITHERING);
DEFINE_SSE_IMPL(high_ordered_dithering, DA_HIGH_ORDERED_DITHERING);
DEFINE_SSE_IMPL(high_floyd_steinberg_dithering, DA_HIGH_FLOYD_STEINBERG_DITHERING);
DEFINE_SSE_IMPL(16bit_interleaved, DA_16BIT_INTERLEAVED);
#endif
#endif
