#pragma once

#include "f3kdb.h"

#include <dualsynth/param.hpp>
#include <dualsynth/video_filter.hpp>

#include <cstdint>
#include <string>
#include <string_view>

namespace neo_f3kdb::core {

enum class Backend : std::uint8_t {
  Scalar,
  Highway,
};

struct DebandParameters {
  int range = 15;
  int y = 64;
  int cb = 64;
  int cr = 64;
  int grain_y = 64;
  int grain_c = 64;
  int sample_mode = 2;
  int seed = 0;
  bool blur_first = true;
  bool dynamic_grain = false;
  DITHER_ALGORITHM dither_algo = DA_HIGH_FLOYD_STEINBERG_DITHERING;
  bool keep_tv_range = false;
  int output_depth = -1;
  RANDOM_ALGORITHM random_algo_ref = RANDOM_ALGORITHM_UNIFORM;
  RANDOM_ALGORITHM random_algo_grain = RANDOM_ALGORITHM_UNIFORM;
  double random_param_ref = 1.0;
  double random_param_grain = 1.0;
  int y_1 = -1;
  int cb_1 = -1;
  int cr_1 = -1;
  int y_2 = -1;
  int cb_2 = -1;
  int cr_2 = -1;
  double angle_boost = 1.5;
  double max_angle = 0.15;
};

struct ParsedParameters {
  DebandParameters params{};
  bool mt = true;
  int opt = -1;
  bool scale = false;
};

[[nodiscard]] ds::Result<ParsedParameters> parse_parameters(const ds::ParamValues* values);

void apply_preset(DebandParameters& params, std::string_view preset, bool scale);
void resolve_parameter_defaults(DebandParameters& params, int input_depth);
void normalize_parameters(DebandParameters& params, int input_depth, bool scale);
void validate_parameters(
  const DebandParameters& params,
  const ds::VideoInputInfo& input,
  bool scale
);

[[nodiscard]] ds::VideoOutputInfo make_output_info(
  const ds::VideoInputInfo& input,
  const DebandParameters& params
);

[[nodiscard]] bool supports_highway(const DebandParameters& params);
[[nodiscard]] Backend select_backend(const DebandParameters& params, int opt);

} // namespace neo_f3kdb::core
