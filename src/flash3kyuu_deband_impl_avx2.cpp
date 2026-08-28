#if defined(HAS_AVX2) && (defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86))

#include <immintrin.h>
#include "flash3kyuu_deband_avx_base.h"

#define INSTANTIATE_IMPL(mode, blur, dither) \
    template void ns_avx2::process_plane_impl<mode, blur, dither>( \
        const process_plane_params&, process_plane_context*);

#define INSTANTIATE_DITHER(dither) \
    INSTANTIATE_IMPL(1, true, dither) INSTANTIATE_IMPL(1, false, dither) \
    INSTANTIATE_IMPL(2, true, dither) INSTANTIATE_IMPL(2, false, dither) \
    INSTANTIATE_IMPL(3, true, dither) INSTANTIATE_IMPL(3, false, dither) \
    INSTANTIATE_IMPL(4, true, dither) INSTANTIATE_IMPL(4, false, dither) \
    INSTANTIATE_IMPL(5, true, dither) INSTANTIATE_IMPL(5, false, dither) \
    INSTANTIATE_IMPL(6, true, dither) INSTANTIATE_IMPL(6, false, dither) \
    INSTANTIATE_IMPL(7, true, dither) INSTANTIATE_IMPL(7, false, dither)

INSTANTIATE_DITHER(DA_HIGH_NO_DITHERING)
INSTANTIATE_DITHER(DA_HIGH_ORDERED_DITHERING)
INSTANTIATE_DITHER(DA_HIGH_FLOYD_STEINBERG_DITHERING)
INSTANTIATE_DITHER(DA_16BIT_INTERLEAVED)

#undef INSTANTIATE_DITHER
#undef INSTANTIATE_IMPL

#endif // HAS_AVX2 && x86
