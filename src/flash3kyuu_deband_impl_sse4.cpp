#include <smmintrin.h>

#include "simd/T_SSE4.hpp"
using SIMD = T_SSE4;
#include "simd/engine.hpp"

#define DECLARE_IMPL_SSE4
#include "impl_dispatch_decl.h"
