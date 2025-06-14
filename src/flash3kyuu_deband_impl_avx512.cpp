#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
#include <immintrin.h>
#include "flash3kyuu_deband_avx512_base.h"

#define DECLARE_IMPL_AVX512
#include "impl_dispatch_decl.h"
#endif
