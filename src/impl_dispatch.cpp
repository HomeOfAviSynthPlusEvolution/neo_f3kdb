#include "core.h"
#include "avx_dispatch.h"

#define IMPL_DISPATCH_IMPORT_DECLARATION

#include "impl_dispatch_decl.h"

#ifdef HAS_SSE4
#define SSE4_IMPL(name) process_plane_impl_sse4_##name
#else
#define SSE4_IMPL(name) process_plane_impl_c_##name
#endif

const process_plane_impl_t* process_plane_impl_high_precision_no_dithering[] = {
    process_plane_impl_c_high_no_dithering,
    process_plane_impl_c_high_no_dithering,
    process_plane_impl_c_high_no_dithering,
    SSE4_IMPL(high_no_dithering),
};

const process_plane_impl_t* process_plane_impl_high_precision_ordered_dithering[] = {
    process_plane_impl_c_high_ordered_dithering,
    process_plane_impl_c_high_ordered_dithering,
    process_plane_impl_c_high_ordered_dithering,
    SSE4_IMPL(high_ordered_dithering),
};

const process_plane_impl_t* process_plane_impl_high_precision_floyd_steinberg_dithering[] = {
    process_plane_impl_c_high_floyd_steinberg_dithering,
    process_plane_impl_c_high_floyd_steinberg_dithering,
    process_plane_impl_c_high_floyd_steinberg_dithering,
    SSE4_IMPL(high_floyd_steinberg_dithering),
};

const process_plane_impl_t* process_plane_impl_16bit_interleaved[] = {
    process_plane_impl_c_16bit_interleaved,
    process_plane_impl_c_16bit_interleaved,
    process_plane_impl_c_16bit_interleaved,
    SSE4_IMPL(16bit_interleaved),
};

const process_plane_impl_t** process_plane_impls[] = {
    nullptr,
    process_plane_impl_high_precision_no_dithering,
    process_plane_impl_high_precision_ordered_dithering,
    process_plane_impl_high_precision_floyd_steinberg_dithering,
    process_plane_impl_16bit_interleaved
};

#define BUILD_NS_ARRAY(ns, dither) { \
    nullptr, \
    &ns::process_plane_impl<1, true,  dither>, &ns::process_plane_impl<1, false, dither>, \
    &ns::process_plane_impl<2, true,  dither>, &ns::process_plane_impl<2, false, dither>, \
    &ns::process_plane_impl<3, true,  dither>, &ns::process_plane_impl<3, false, dither>, \
    &ns::process_plane_impl<4, true,  dither>, &ns::process_plane_impl<4, false, dither>, \
    &ns::process_plane_impl<5, true,  dither>, &ns::process_plane_impl<5, false, dither>, \
    &ns::process_plane_impl<6, true,  dither>, &ns::process_plane_impl<6, false, dither>, \
    &ns::process_plane_impl<7, true,  dither>, &ns::process_plane_impl<7, false, dither>  \
}

#ifdef HAS_AVX2
static const process_plane_impl_t avx2_no_dither[] = BUILD_NS_ARRAY(ns_avx2, DA_HIGH_NO_DITHERING);
static const process_plane_impl_t avx2_ordered_dither[] = BUILD_NS_ARRAY(ns_avx2, DA_HIGH_ORDERED_DITHERING);
static const process_plane_impl_t avx2_floyd_steinberg[] = BUILD_NS_ARRAY(ns_avx2, DA_HIGH_FLOYD_STEINBERG_DITHERING);
static const process_plane_impl_t avx2_16bit[] = BUILD_NS_ARRAY(ns_avx2, DA_16BIT_INTERLEAVED);
#endif

#ifdef HAS_AVX512
static const process_plane_impl_t avx512_no_dither[] = BUILD_NS_ARRAY(ns_avx512, DA_HIGH_NO_DITHERING);
static const process_plane_impl_t avx512_ordered_dither[] = BUILD_NS_ARRAY(ns_avx512, DA_HIGH_ORDERED_DITHERING);
static const process_plane_impl_t avx512_floyd_steinberg[] = BUILD_NS_ARRAY(ns_avx512, DA_HIGH_FLOYD_STEINBERG_DITHERING);
static const process_plane_impl_t avx512_16bit[] = BUILD_NS_ARRAY(ns_avx512, DA_16BIT_INTERLEAVED);
#endif

#undef BUILD_NS_ARRAY

const process_plane_impl_t* select_avx_impl(int opt, int cpu_flags, int dither_algo) {
#if !defined(HAS_AVX2) && !defined(HAS_AVX512)
    (void)opt; (void)cpu_flags; (void)dither_algo;
    return nullptr;
#else
    const int AVX512_REQUIRED = CPUF_AVX512F | CPUF_AVX512BW | CPUF_AVX512DQ | CPUF_AVX512VL | CPUF_AVX512CD;

#ifdef HAS_AVX512
    const bool has_avx512 = ((cpu_flags & AVX512_REQUIRED) == AVX512_REQUIRED);
#else
    const bool has_avx512 = false;
#endif

#ifdef HAS_AVX2
    const bool has_avx2 = (cpu_flags & CPUF_AVX2) != 0;
#else
    const bool has_avx2 = false;
#endif

    bool use_avx512 = false;
    bool use_avx2 = false;

    if (opt == IMPL_AVX512) {
        use_avx512 = has_avx512;
    }
    else if (opt == IMPL_AVX2) {
        use_avx2 = has_avx2;
    }
    else if (opt == IMPL_AUTO_DETECT) {
        if (has_avx512)
            use_avx512 = true;
        else if (has_avx2)
            use_avx2 = true;
    }

    auto get_array = [](int dither, const process_plane_impl_t* no_dith, const process_plane_impl_t* ord_dith,
        const process_plane_impl_t* fs_dith, const process_plane_impl_t* int16_dith) -> const process_plane_impl_t*
        {
            switch (dither) {
            case DA_HIGH_NO_DITHERING: return no_dith;
            case DA_HIGH_ORDERED_DITHERING: return ord_dith;
            case DA_HIGH_FLOYD_STEINBERG_DITHERING: return fs_dith;
            case DA_16BIT_INTERLEAVED: return int16_dith;
            default: return nullptr;
            }
        };

#ifdef HAS_AVX512
    if (use_avx512)
        return get_array(dither_algo, avx512_no_dither, avx512_ordered_dither, avx512_floyd_steinberg, avx512_16bit);
#endif

#ifdef HAS_AVX2
    if (use_avx2)
        return get_array(dither_algo, avx2_no_dither, avx2_ordered_dither, avx2_floyd_steinberg, avx2_16bit);
#endif

    return nullptr;
#endif // HAS_AVX2 || HAS_AVX512
}
