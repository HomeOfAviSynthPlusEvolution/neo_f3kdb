// #include <intrin.h>

#define SIMD_PREFIX _mm_
#define SIMD_ANY_SUFFIX _si128

#define REG_ANY(n) REG_ANY_TYPE_INTERNAL1(SIMD_PREFIX, n, SIMD_ANY_SUFFIX)

#define REG(n)  REG_INTERNAL1(SIMD_PREFIX, n)
#define REG_I08(n) REG_TYPE_INTERNAL1(SIMD_PREFIX, n, _epi8)
#define REG_I16(n) REG_TYPE_INTERNAL1(SIMD_PREFIX, n, _epi16)
#define REG_I32(n) REG_TYPE_INTERNAL1(SIMD_PREFIX, n, _epi32)
#define REG_I64(n) REG_TYPE_INTERNAL1(SIMD_PREFIX, n, _epi64)
#define REG_U08(n) REG_TYPE_INTERNAL1(SIMD_PREFIX, n, _epu8)
#define REG_U16(n) REG_TYPE_INTERNAL1(SIMD_PREFIX, n, _epu16)
#define REG_U32(n) REG_TYPE_INTERNAL1(SIMD_PREFIX, n, _epu32)
#define REG_U64(n) REG_TYPE_INTERNAL1(SIMD_PREFIX, n, _epu64)

#define REG_I(n)   REG_I08(n) REG_I16(n) REG_I32(n)
#define REG_U(n)   REG_U08(n) REG_U16(n) REG_U32(n)

#define REG_L(n)   REG_I08(n) REG_I16(n) \
                   REG_U08(n) REG_U16(n)
#define REG_A(n)   REG_I(n) REG_U(n)

#define REG_ANY_TYPE_INTERNAL1(prefix, n, suffix) REG_ANY_TYPE_INTERNAL2(prefix, n, suffix)

#define REG_TYPE_INTERNAL1(prefix, n, suffix) REG_TYPE_INTERNAL2(prefix, n, suffix)

#define REG_INTERNAL1(prefix, n) REG_INTERNAL2(prefix, n)

#define REG_ANY_TYPE_INTERNAL2(prefix, n, suffix) \
  static constexpr auto _##n = prefix##n##suffix;

#define REG_TYPE_INTERNAL2(prefix, n, suffix) \
  static constexpr auto _##n##suffix = prefix##n##suffix;

#define REG_INTERNAL2(prefix, n) \
  static constexpr auto _##n = prefix##n;

class T_SSE4
{
  public:
  using data_type = __m128i;
  using half_data_type = __m128i;
  using count_type = __m128i;
  using mask_type = __m128i;
  static const int width_8 = 16;
  static constexpr int align = width_8;
  static constexpr int width_16 = width_8 >> 1;
  static constexpr int width_32 = width_8 >> 2;
  static constexpr int width_64 = width_8 >> 3;
  #include "T_shared.inc"
  #include "T_shared_low.inc"
  static constexpr auto _storel = _mm_storel_epi64;
  static constexpr auto _loadl = _mm_loadl_epi64;
  static constexpr auto count_zero = _mm_setzero_si128;
  static data_type cast_full(half_data_type reg)
  {
    return reg;
  }
  static count_type count_set(int cnt)
  {
    return _mm_set_epi32(0, 0, 0, cnt);
  }
};
