// #include <intrin.h>

#define SIMD_PREFIX _mm256_
#define SIMD_ANY_SUFFIX _si256

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

class T_AVX2
{
  public:
  using data_type = __m256i;
  using half_data_type = __m128i;
  using count_type = __m128i;
  using mask_type = __m256i;
  static const int width_8 = 32;
  static constexpr int align = width_8;
  static constexpr int width_16 = width_8 >> 1;
  static constexpr int width_32 = width_8 >> 2;
  static constexpr int width_64 = width_8 >> 3;
  #include "T_shared.inc"
  #include "T_shared_low.inc"
  static constexpr auto count_zero = _mm_setzero_si128;
  static data_type cast_full(half_data_type reg)
  {
    return _mm256_castsi128_si256(reg);
  }
  static void _storel(half_data_type *mem_addr, data_type reg)
  {
    _mm_store_si128(mem_addr, _mm256_castsi256_si128(reg));
  }
  static half_data_type _loadl(half_data_type const *mem_addr)
  {
    return _mm_load_si128(mem_addr);
  }
  static count_type count_set(int cnt)
  {
    return _mm_set_epi32(0, 0, 0, cnt);
  }
};
