#pragma once

#include "core/kernel.hpp"

#ifdef NEO_F3KDB_KERNEL_SCALAR_IMPORT_DECLARATION

#define DEFINE_IMPL(n, \
                    null_entry, \
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
  extern const deband_kernel_proc_t process_plane_impl_##n [];

#else

#define DEFINE_IMPL(n, \
                    null_entry, \
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
  extern const deband_kernel_proc_t process_plane_impl_##n [] = { \
    null_entry, \
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
    impl_func_mode7_noblur \
  };

#endif

#define DEFINE_TEMPLATE_IMPL(name, impl_func, ...) \
  DEFINE_IMPL(name, \
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
              (&impl_func<7, false, __VA_ARGS__>))

#if defined(NEO_F3KDB_KERNEL_SCALAR_IMPORT_DECLARATION) || defined(NEO_F3KDB_KERNEL_SCALAR_DEFINE_TABLES)
DEFINE_TEMPLATE_IMPL(scalar_high_no_dithering, process_plane_plainc, DA_HIGH_NO_DITHERING);
DEFINE_TEMPLATE_IMPL(scalar_high_ordered_dithering, process_plane_plainc, DA_HIGH_ORDERED_DITHERING);
DEFINE_TEMPLATE_IMPL(
  scalar_high_floyd_steinberg_dithering,
  process_plane_plainc,
  DA_HIGH_FLOYD_STEINBERG_DITHERING
);
DEFINE_TEMPLATE_IMPL(scalar_16bit_interleaved, process_plane_plainc, DA_16BIT_INTERLEAVED);
#endif

#undef DEFINE_TEMPLATE_IMPL
#undef DEFINE_IMPL
