#include "core/random.hpp"

#include <cassert>
#include <cmath>
#include <cstdint>

typedef double (*rand_impl_t)(int& seed, double param);

double rand_old(int& seed, double param);
double rand_uniform(int& seed, double param);
double rand_gaussian(int& seed, double param);

static const rand_impl_t rand_algorithms[] = {
  rand_old,
  rand_uniform,
  rand_gaussian
};

inline double round(double r) {
  return (r > 0.0) ? std::floor(r + 0.5) : std::ceil(r - 0.5);
}

int random(RANDOM_ALGORITHM algo, int& seed, int range, double param) {
  assert(algo >= 0 && algo < RANDOM_ALGORITHM_COUNT);

  const double num = rand_algorithms[algo](seed, param);
  assert(num >= -1.0 && num <= 1.0);
  return static_cast<int>(round(num * range));
}

double rand_to_double(int rand_num) {
  union {
    std::uint64_t itemp;
    double result;
  };
  itemp = static_cast<std::uint64_t>(rand_num) & 0xffffffffULL;
  itemp = itemp << 20 | itemp >> 12;
  itemp |= 0x3ff0000000000000ULL;
  return (result - 1.0) * 2 - 1.0;
}

double rand_old(int& seed, double) {
  const int seed_tmp = (((seed << 13) ^ static_cast<unsigned int>(seed)) >> 17) ^
    (seed << 13) ^ seed;
  seed = 32 * seed_tmp ^ seed_tmp;
  return rand_to_double(seed);
}

double rand_uniform(int& seed, double) {
  seed = 1664525 * seed + 1013904223;
  return rand_to_double(seed);
}

double rand_gaussian(int& seed, double param) {
  double ret;
  double x;
  double y;
  double r2;

  do {
    do {
      x = rand_uniform(seed, param);
      y = rand_uniform(seed, param);
      r2 = x * x + y * y;
    } while (r2 > 1.0 || r2 == 0);

    ret = param * y * std::sqrt(-2.0 * std::log(r2) / r2);
  } while (ret <= -1.0 || ret >= 1.0);

  return ret;
}
