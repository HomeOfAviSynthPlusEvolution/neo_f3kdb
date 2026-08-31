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
      p.y_raw = get_value(values->get_double("y", p.y_raw), p.y_raw);
      p.cb_raw = get_value(values->get_double("cb", p.cb_raw), p.cb_raw);
      p.cr_raw = get_value(values->get_double("cr", p.cr_raw), p.cr_raw);
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
      p.y_1_raw = get_value(values->get_double("y_1", p.y_1_raw), p.y_1_raw);
      p.cb_1_raw = get_value(values->get_double("cb_1", p.cb_1_raw), p.cb_1_raw);
      p.cr_1_raw = get_value(values->get_double("cr_1", p.cr_1_raw), p.cr_1_raw);
      p.y_2_raw = get_value(values->get_double("y_2", p.y_2_raw), p.y_2_raw);
      p.cb_2_raw = get_value(values->get_double("cb_2", p.cb_2_raw), p.cb_2_raw);
      p.cr_2_raw = get_value(values->get_double("cr_2", p.cr_2_raw), p.cr_2_raw);
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
      p.y_raw = p.cb_raw = p.cr_raw = p.grain_y = p.grain_c = p.y_1_raw = p.cb_1_raw = p.cr_1_raw = p.y_2_raw =
        p.cb_2_raw = p.cr_2_raw = 0.0;
    } else if (token == "low") {
      p.y_raw = p.cb_raw = p.cr_raw = p.y_1_raw = p.cb_1_raw = p.cr_1_raw = p.y_2_raw =
        p.cb_2_raw = p.cr_2_raw = scale ? 128.0 : 32.0;
      p.grain_y = p.grain_c = scale ? 128 : 32;
    } else if (token == "medium") {
      p.y_raw = p.cb_raw = p.cr_raw = p.y_1_raw = p.cb_1_raw = p.cr_1_raw = p.y_2_raw =
        p.cb_2_raw = p.cr_2_raw = scale ? 192.0 : 48.0;
      p.grain_y = p.grain_c = scale ? 192 : 48;
    } else if (token == "high") {
      p.y_raw = p.cb_raw = p.cr_raw = p.y_1_raw = p.cb_1_raw = p.cr_1_raw = p.y_2_raw =
        p.cb_2_raw = p.cr_2_raw = scale ? 256.0 : 64.0;
      p.grain_y = p.grain_c = scale ? 256 : 64;
    } else if (token == "veryhigh") {
      p.y_raw = p.cb_raw = p.cr_raw = p.y_1_raw = p.cb_1_raw = p.cr_1_raw = p.y_2_raw =
        p.cb_2_raw = p.cr_2_raw = scale ? 320.0 : 80.0;
      p.grain_y = p.grain_c = scale ? 320 : 80;
    } else if (token == "nograin") {
      p.grain_y = p.grain_c = 0;
    } else if (token == "luma") {
      p.cb_raw = p.cr_raw = p.cb_1_raw = p.cr_1_raw = p.cb_2_raw = p.cr_2_raw = 0.0;
      p.grain_c = 0;
    } else if (token == "chroma") {
      p.y_raw = p.y_1_raw = p.y_2_raw = 0.0;
      p.grain_y = 0;
    }
  }
}

void validate_parameters(const DebandParameters& p, const ds::VideoInputInfo& input, bool scale) {
  invalid_param_if(
    input.format.color_family != ds::ColorFamily::Yuv && input.format.color_family != ds::ColorFamily::Gray,
    "input is not YUV or Gray"
  );
  invalid_param_if(input.width < 16, "input.width < 16");
  invalid_param_if(input.height < 16, "input.height < 16");
  invalid_param_if(input.num_frames <= 0, "input.num_frames <= 0");

  const int bits = ds::bits_per_sample(input.format.sample_format);
  invalid_param_if(bits < 8 || bits > INTERNAL_BIT_DEPTH, "bits out of range");
  invalid_param_if(input.format.sample_format == ds::SampleFormat::Float32, "float input");

  const double threshold_upper_limit = scale ? 65535.0 : 511.0;
  constexpr int dither_upper_limit = 4096;

  check_param("range", p.range, 0, 255);
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

  auto validate_thresh = [&](const char* name, double val) {
    if (val < 0.0 || val > threshold_upper_limit) {
      char err[256];
      std::snprintf(err, sizeof(err), "Invalid parameter %s, must be between 0.0 and %.1f", name, threshold_upper_limit);
      throw std::invalid_argument(err);
    }
  };

  validate_thresh("Y", p.y_raw);
  validate_thresh("Cb", p.cb_raw);
  validate_thresh("Cr", p.cr_raw);
  if (p.y_1_raw >= 0.0) validate_thresh("Y_1", p.y_1_raw);
  if (p.cb_1_raw >= 0.0) validate_thresh("Cb_1", p.cb_1_raw);
  if (p.cr_1_raw >= 0.0) validate_thresh("Cr_1", p.cr_1_raw);
  if (p.y_2_raw >= 0.0) validate_thresh("Y_2", p.y_2_raw);
  if (p.cb_2_raw >= 0.0) validate_thresh("Cb_2", p.cb_2_raw);
  if (p.cr_2_raw >= 0.0) validate_thresh("Cr_2", p.cr_2_raw);

  if (p.angle_boost < 0.0) {
    throw std::invalid_argument("invalid parameter angle_boost, must be positive value");
  }
  if (p.max_angle < 0.0 || p.max_angle > 1.0) {
    throw std::invalid_argument("invalid parameter max_angle, must be between 0.0 and 1.0");
  }
}

void resolve_parameter_defaults(DebandParameters& p, int input_depth) {
  p.y_1_raw = p.y_1_raw < 0.0 ? p.y_raw : p.y_1_raw;
  p.cb_1_raw = p.cb_1_raw < 0.0 ? p.cb_raw : p.cb_1_raw;
  p.cr_1_raw = p.cr_1_raw < 0.0 ? p.cr_raw : p.cr_1_raw;
  p.y_2_raw = p.y_2_raw < 0.0 ? p.y_raw : p.y_2_raw;
  p.cb_2_raw = p.cb_2_raw < 0.0 ? p.cb_raw : p.cb_2_raw;
  p.cr_2_raw = p.cr_2_raw < 0.0 ? p.cr_raw : p.cr_2_raw;

  if (p.output_depth < 0) {
    p.output_depth = input_depth;
  }
  if (p.output_depth == 16) {
    p.dither_algo = DA_16BIT_INTERLEAVED;
  }
}

static inline int normalize_single_threshold(double val, bool scale) {
  if (val <= 0.0) {
    return 0;
  }
  if (!scale && val <= 1.0) {
    return static_cast<int>(val * 65535.0 + 0.5);
  } else if (scale) {
    return static_cast<int>(val + 0.5);
  } else {
    return static_cast<int>(val + 0.5) << 2;
  }
}

void normalize_parameters(DebandParameters& p, int input_depth, bool scale) {
  (void)input_depth;
  p.y = normalize_single_threshold(p.y_raw, scale);
  p.cb = normalize_single_threshold(p.cb_raw, scale);
  p.cr = normalize_single_threshold(p.cr_raw, scale);
  p.y_1 = normalize_single_threshold(p.y_1_raw, scale);
  p.cb_1 = normalize_single_threshold(p.cb_1_raw, scale);
  p.cr_1 = normalize_single_threshold(p.cr_1_raw, scale);
  p.y_2 = normalize_single_threshold(p.y_2_raw, scale);
  p.cb_2 = normalize_single_threshold(p.cb_2_raw, scale);
  p.cr_2 = normalize_single_threshold(p.cr_2_raw, scale);
  p.grain_y <<= 2;
  p.grain_c <<= 2;
}

ds::VideoOutputInfo make_output_info(
  const ds::VideoInputInfo& input,
  const DebandParameters& params
) {
  ds::VideoOutputInfo output{};
  output.width = input.width;
  output.height = input.height;
  output.num_frames = input.num_frames;
  output.format.color_family = input.format.color_family;
  output.format.sample_format = params.output_depth == 8 ? ds::SampleFormat::UInt8 : ds::SampleFormat::UInt16;
  output.format.plane_count = input.format.plane_count;
  output.format.subsampling_w = input.format.subsampling_w;
  output.format.subsampling_h = input.format.subsampling_h;
  output.fps = input.fps;
  return output;
}

Backend select_backend(const DebandParameters&, int opt) {
  if (opt == 0) {
    return Backend::Scalar;
  }
  return Backend::Highway;
}

} // namespace neo_f3kdb::core
