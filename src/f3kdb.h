#pragma once

typedef enum _PIXEL_MODE : int {
  DEFAULT_PIXEL_MODE = -1,
  LOW_BIT_DEPTH = 0,
  INVALID_OPTION1,
  HIGH_BIT_DEPTH_INTERLEAVED,
  PIXEL_MODE_COUNT
} PIXEL_MODE;

typedef enum _DITHER_ALGORITHM : int {
  // _DEPRECATED_DA_LOW = 0,
  DA_HIGH_NO_DITHERING = 1,
  DA_HIGH_ORDERED_DITHERING,
  DA_HIGH_FLOYD_STEINBERG_DITHERING,
  DA_16BIT_INTERLEAVED,

  DA_COUNT,
  DA_USER_PARAM_MAX = DA_HIGH_FLOYD_STEINBERG_DITHERING
} DITHER_ALGORITHM;

typedef enum _RANDOM_ALGORITHM : int {
  RANDOM_ALGORITHM_OLD = 0,
  RANDOM_ALGORITHM_UNIFORM,
  RANDOM_ALGORITHM_GAUSSIAN,
  RANDOM_ALGORITHM_COUNT
} RANDOM_ALGORITHM;

typedef enum _OPTIMIZATION_MODE : int {
  IMPL_AUTO_DETECT = -1,
  IMPL_C = 0,
  IMPL_SSE2,
  IMPL_SSSE3,
  IMPL_SSE4,
  IMPL_AVX2,
  IMPL_AVX512,

  IMPL_COUNT
} OPTIMIZATION_MODE;

typedef struct _f3kdb_params_t {
  int range {15};
  int Y {64};
  int Cb {64};
  int Cr {64};
  int grainY {64};
  int grainC {64};
  int sample_mode {2};
  int seed {0};
  bool blur_first {true};
  bool dynamic_grain {false};
  DITHER_ALGORITHM dither_algo {DA_HIGH_FLOYD_STEINBERG_DITHERING};
  bool keep_tv_range {false};
  int output_depth {-1};
  RANDOM_ALGORITHM random_algo_ref {RANDOM_ALGORITHM_UNIFORM};
  RANDOM_ALGORITHM random_algo_grain {RANDOM_ALGORITHM_UNIFORM};
  double random_param_ref {1.0f};
  double random_param_grain {1.0f};
  int Y_1 {-1};
  int Cb_1 {-1};
  int Cr_1 {-1};
  int Y_2 {-1};
  int Cb_2 {-1};
  int Cr_2 {-1};
  double angle_boost {1.5};
  double max_angle {0.15};
} f3kdb_params_t;
