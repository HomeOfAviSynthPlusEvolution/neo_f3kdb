#include "core/parameters.hpp"

#include "core/constants.hpp"

#include <dualsynth/format.hpp>

#include <cstdio>
#include <exception>
#include <sstream>
#include <stdexcept>
#include <string>

namespace neo_f3kdb::core {

namespace {

template <class T>
T get_value(const ds::Result<T>& res, T default_val) {
  return res.has_value() ? res.value() : default_val;
}

void invalid_param_if(bool condition, const char* expression) {
  if (condition) {
    throw std::invalid_argument(std::string("Invalid parameter condition: ") + expression);
  }
}

void check_param(const char* name, int value, int lower_bound, int upper_bound) {
  if (value < lower_bound || value > upper_bound) {
    char err[256];
    std::snprintf(
      err,
      sizeof(err),
      "Invalid parameter %s, must be between %d and %d",
      name,
      lower_bound,
      upper_bound
    );
    throw std::invalid_argument(err);
  }
}

} // namespace

ds::Result<ParsedParameters> parse_parameters(const ds::ParamValues* values) {
  try {
    ParsedParameters parsed{};
    std::string preset;

    if (values != nullptr) {
      preset = get_value(values->get_string("preset", ""), std::string(""));
      parsed.scale = get_value(values->get_bool("scale", false), false);
      parsed.mt = get_value(values->get_bool("mt", true), true);
      parsed.opt = get_value(values->get_int("opt", -1), -1);
    }

    apply_preset(parsed.params, preset, parsed.scale);

    if (values != nullptr) {
      auto& p = parsed.params;
      p.range = get_value(values->get_int("range", p.range), p.range);
      p.y = get_value(values->get_int("y", p.y), p.y);
      p.cb = get_value(values->get_int("cb", p.cb), p.cb);
      p.cr = get_value(values->get_int("cr", p.cr), p.cr);
      p.grain_y = get_value(values->get_int("grainy", p.grain_y), p.grain_y);
      p.grain_c = get_value(values->get_int("grainc", p.grain_c), p.grain_c);
      p.sample_mode = get_value(values->get_int("sample_mode", p.sample_mode), p.sample_mode);
      p.seed = get_value(values->get_int("seed", p.seed), p.seed);
      p.blur_first = get_value(values->get_bool("blur_first", p.blur_first), p.blur_first);
      p.dynamic_grain =
        get_value(values->get_bool("dynamic_grain", p.dynamic_grain), p.dynamic_grain);
      p.keep_tv_range =
        get_value(values->get_bool("keep_tv_range", p.keep_tv_range), p.keep_tv_range);
      p.output_depth = get_value(values->get_int("output_depth", p.output_depth), p.output_depth);
      p.random_algo_ref = static_cast<RANDOM_ALGORITHM>(
        get_value(
          values->get_int("random_algo_ref", static_cast<int>(p.random_algo_ref)),
          static_cast<int>(p.random_algo_ref)
        )
      );
      p.random_algo_grain = static_cast<RANDOM_ALGORITHM>(
        get_value(
          values->get_int("random_algo_grain", static_cast<int>(p.random_algo_grain)),
          static_cast<int>(p.random_algo_grain)
        )
      );
      p.random_param_ref =
        get_value(values->get_double("random_param_ref", p.random_param_ref), p.random_param_ref);
      p.random_param_grain = get_value(
        values->get_double("random_param_grain", p.random_param_grain),
        p.random_param_grain
      );
      p.y_1 = get_value(values->get_int("y_1", p.y_1), p.y_1);
      p.cb_1 = get_value(values->get_int("cb_1", p.cb_1), p.cb_1);
      p.cr_1 = get_value(values->get_int("cr_1", p.cr_1), p.cr_1);
      p.y_2 = get_value(values->get_int("y_2", p.y_2), p.y_2);
      p.cb_2 = get_value(values->get_int("cb_2", p.cb_2), p.cb_2);
      p.cr_2 = get_value(values->get_int("cr_2", p.cr_2), p.cr_2);
      p.angle_boost = get_value(values->get_double("angle_boost", p.angle_boost), p.angle_boost);
      p.max_angle = get_value(values->get_double("max_angle", p.max_angle), p.max_angle);
      p.dither_algo = static_cast<DITHER_ALGORITHM>(
        get_value(
          values->get_int("dither_algo", static_cast<int>(p.dither_algo)),
          static_cast<int>(p.dither_algo)
        )
      );
    }

    return ds::Result<ParsedParameters>::success(std::move(parsed));
  } catch (const std::exception& error) {
    return ds::Result<ParsedParameters>::failure(
      ds::Error{ds::ErrorCode::InvalidArgument, error.what()}
    );
  }
}

void apply_preset(DebandParameters& p, std::string_view preset, bool scale) {
  std::istringstream stream{std::string(preset)};
  while (!stream.eof()) {
    std::string token;
    std::getline(stream, token, '/');
    if (token == "depth") {
      p.y = p.cb = p.cr = p.grain_y = p.grain_c = p.y_1 = p.cb_1 = p.cr_1 = p.y_2 =
        p.cb_2 = p.cr_2 = 0;
    } else if (token == "low") {
      p.y = p.cb = p.cr = p.grain_y = p.grain_c = p.y_1 = p.cb_1 = p.cr_1 = p.y_2 =
        p.cb_2 = p.cr_2 = scale ? 128 : 32;
    } else if (token == "medium") {
      p.y = p.cb = p.cr = p.grain_y = p.grain_c = p.y_1 = p.cb_1 = p.cr_1 = p.y_2 =
        p.cb_2 = p.cr_2 = scale ? 192 : 48;
    } else if (token == "high") {
      p.y = p.cb = p.cr = p.grain_y = p.grain_c = p.y_1 = p.cb_1 = p.cr_1 = p.y_2 =
        p.cb_2 = p.cr_2 = scale ? 256 : 64;
    } else if (token == "veryhigh") {
      p.y = p.cb = p.cr = p.grain_y = p.grain_c = p.y_1 = p.cb_1 = p.cr_1 = p.y_2 =
        p.cb_2 = p.cr_2 = scale ? 320 : 80;
    } else if (token == "nograin") {
      p.grain_y = p.grain_c = 0;
    } else if (token == "luma") {
      p.cb = p.cr = p.grain_c = 0;
    } else if (token == "chroma") {
      p.y = p.grain_y = 0;
    }
  }
}

void validate_parameters(const DebandParameters& p, const ds::VideoInputInfo& input, bool scale) {
  invalid_param_if(input.format.color_family != ds::ColorFamily::Yuv, "input is not YUV");
  invalid_param_if(input.width < 16, "input.width < 16");
  invalid_param_if(input.height < 16, "input.height < 16");
  invalid_param_if(input.num_frames <= 0, "input.num_frames <= 0");

  const int bits = ds::bits_per_sample(input.format.sample_format);
  invalid_param_if(bits < 8 || bits > INTERNAL_BIT_DEPTH, "bits out of range");
  invalid_param_if(input.format.sample_format == ds::SampleFormat::Float32, "float input");

  const int y_threshold_upper_limit = scale ? 65535 : 511;
  const int cb_threshold_upper_limit = scale ? 65535 : 511;
  const int cr_threshold_upper_limit = scale ? 65535 : 511;
  constexpr int dither_upper_limit = 4096;

  check_param("range", p.range, 0, 255);
  check_param("Y", p.y, 0, y_threshold_upper_limit);
  check_param("Cb", p.cb, 0, cb_threshold_upper_limit);
  check_param("Cr", p.cr, 0, cr_threshold_upper_limit);
  check_param("grainY", p.grain_y, 0, dither_upper_limit);
  check_param("grainC", p.grain_c, 0, dither_upper_limit);
  check_param("sample_mode", p.sample_mode, 1, 7);
  check_param("dither_algo", p.dither_algo, DA_HIGH_NO_DITHERING, DA_COUNT - 1);
  check_param("random_algo_ref", p.random_algo_ref, RANDOM_ALGORITHM_OLD, RANDOM_ALGORITHM_COUNT - 1);
  check_param(
    "random_algo_grain",
    p.random_algo_grain,
    RANDOM_ALGORITHM_OLD,
    RANDOM_ALGORITHM_COUNT - 1
  );
  check_param("Y_1", p.y_1, 0, y_threshold_upper_limit);
  check_param("Cb_1", p.cb_1, 0, cb_threshold_upper_limit);
  check_param("Cr_1", p.cr_1, 0, cr_threshold_upper_limit);
  check_param("Y_2", p.y_2, 0, y_threshold_upper_limit);
  check_param("Cb_2", p.cb_2, 0, cb_threshold_upper_limit);
  check_param("Cr_2", p.cr_2, 0, cr_threshold_upper_limit);

  if (p.angle_boost < 0.0) {
    throw std::invalid_argument("invalid parameter angle_boost, must be positive value");
  }
  if (p.max_angle < 0.0 || p.max_angle > 1.0) {
    throw std::invalid_argument("invalid parameter max_angle, must be between 0.0 and 1.0");
  }
}

void resolve_parameter_defaults(DebandParameters& p, int input_depth) {
  p.y_1 = p.y_1 == -1 ? p.y : p.y_1;
  p.cb_1 = p.cb_1 == -1 ? p.cb : p.cb_1;
  p.cr_1 = p.cr_1 == -1 ? p.cr : p.cr_1;
  p.y_2 = p.y_2 == -1 ? p.y : p.y_2;
  p.cb_2 = p.cb_2 == -1 ? p.cb : p.cb_2;
  p.cr_2 = p.cr_2 == -1 ? p.cr : p.cr_2;

  if (p.output_depth < 0) {
    p.output_depth = input_depth;
  }
  if (p.output_depth == 16) {
    p.dither_algo = DA_16BIT_INTERLEAVED;
  }
}

void normalize_parameters(DebandParameters& p, int input_depth, bool scale) {
  (void)input_depth;
  p.y = scale ? p.y : p.y << 2;
  p.cb = scale ? p.cb : p.cb << 2;
  p.cr = scale ? p.cr : p.cr << 2;
  p.y_1 = scale ? p.y_1 : p.y_1 << 2;
  p.cb_1 = scale ? p.cb_1 : p.cb_1 << 2;
  p.cr_1 = scale ? p.cr_1 : p.cr_1 << 2;
  p.y_2 = scale ? p.y_2 : p.y_2 << 2;
  p.cb_2 = scale ? p.cb_2 : p.cb_2 << 2;
  p.cr_2 = scale ? p.cr_2 : p.cr_2 << 2;
  p.grain_y <<= 2;
  p.grain_c <<= 2;
}

ds::VideoOutputInfo make_output_info(
  const ds::VideoInputInfo& input,
  const DebandParameters& params
) {
  return ds::VideoOutputInfo{
    .width = input.width,
    .height = input.height,
    .num_frames = input.num_frames,
    .format = ds::VideoFormat{
      .color_family = input.format.color_family,
      .sample_format = params.output_depth == 8 ? ds::SampleFormat::UInt8 : ds::SampleFormat::UInt16,
      .plane_count = input.format.plane_count,
      .subsampling_w = input.format.subsampling_w,
      .subsampling_h = input.format.subsampling_h,
    },
    .fps = input.fps
  };
}

bool supports_highway(const DebandParameters& params) {
  return params.dither_algo != DA_HIGH_FLOYD_STEINBERG_DITHERING &&
    params.sample_mode >= 1 &&
    params.sample_mode <= 5;
}

Backend select_backend(const DebandParameters& params, int opt) {
  if (opt == 0) {
    return Backend::Scalar;
  }
  return supports_highway(params) ? Backend::Highway : Backend::Scalar;
}

} // namespace neo_f3kdb::core
