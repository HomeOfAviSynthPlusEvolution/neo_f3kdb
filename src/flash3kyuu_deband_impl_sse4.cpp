#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
#include <smmintrin.h>
#elif defined(__arm__) || defined(__aarch64__) || defined(_M_ARM)
#include "sse2neon.h"
#endif
#include "flash3kyuu_deband_sse_base.h"

#define DECLARE_IMPL_SSE4
#include "impl_dispatch_decl.h"
