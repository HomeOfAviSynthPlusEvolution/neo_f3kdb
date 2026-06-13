#include "plugin/f3kdb_filter.hpp"

#include "constants.h"
#include "core.h"

int GetCPUFlags();

#include <climits>
#include <cstdint>
#include <exception>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace neo_f3kdb {

namespace {

template <class T>
T get_value(const ds::Result<T>& res, T default_val) {
  return res.has_value() ? res.value() : default_val;
}

DSFormat make_ds_format(ds::VideoFormat format) {
  DSFormat df;
  df.Planes = format.plane_count;
  df.IsFamilyYUV = format.color_family == ds::ColorFamily::Yuv || format.color_family == ds::ColorFamily::Gray;
  df.SSW = format.subsampling_w;
  df.SSH = format.subsampling_h;
  df.BitsPerSample = ds::bits_per_sample(format.sample_format);
  return df;
}

DSVideoInfo make_ds_video_info(ds::VideoInputInfo vi) {
  DSVideoInfo dvi;
  dvi.Format = make_ds_format(vi.format);
  dvi.Width = vi.width;
  dvi.Height = vi.height;
  dvi.Frames = vi.num_frames;
  return dvi;
}

#define INVALID_PARAM_IF(cond) \
do { if (cond) { throw std::invalid_argument("Invalid parameter condition: " #cond); } } while (0)

#define CHECK_PARAM(value, lower_bound, upper_bound) \
do { if (static_cast<int>(value) < static_cast<int>(lower_bound) || static_cast<int>(value) > static_cast<int>(upper_bound)) { \
    char err[256]; \
    snprintf(err, sizeof(err), "Invalid parameter %s, must be between %d and %d", #value, static_cast<int>(lower_bound), static_cast<int>(upper_bound)); \
    throw std::invalid_argument(err); \
} } while(0)

} // namespace

F3KDBFilterCore::State::State(std::unique_ptr<f3kdb_core_t> engine, f3kdb_params_t params, bool mt)
  : engine(std::move(engine)), params(params), mt(mt) {}

F3KDBFilterCore::State::State(State&&) noexcept = default;
F3KDBFilterCore::State& F3KDBFilterCore::State::operator=(State&&) noexcept = default;

ds::Result<ds::VideoInitStateResult<F3KDBFilterCore::State>> F3KDBFilterCore::init(ds::VideoInitContext& context) {
  if (context.inputs.size() != 1) {
    return ds::Result<ds::VideoInitStateResult<State>>::failure(
      ds::Error{ds::ErrorCode::InvalidArgument, "Neo-F3KDB: exactly one video input is required"}
    );
  }

  const ds::VideoInputInfo& input_vi = context.inputs[0];
  f3kdb_params_t ep;

  std::string preset = "";
  bool scale = false;
  bool mt = true;
  int opt_in = -1;

  if (context.params) {
    preset = get_value(context.params->get_string("preset", ""), std::string(""));
    scale = get_value(context.params->get_bool("scale", false), false);
    mt = get_value(context.params->get_bool("mt", true), true);
    opt_in = get_value(context.params->get_int("opt", -1), -1);
  }

  if (!preset.empty()) {
    std::istringstream piss(preset);
    while(!piss.eof()) {
      std::string piss1;
      std::getline(piss, piss1, '/');
      if (piss1 == "depth")
          ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = ep.Y_1 = ep.Cb_1 = ep.Cr_1 = ep.Y_2 = ep.Cb_2 = ep.Cr_2 = 0;
      else if (piss1 == "low")
          ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = ep.Y_1 = ep.Cb_1 = ep.Cr_1 = ep.Y_2 = ep.Cb_2 = ep.Cr_2 = (scale) ? 128 : 32;
      else if (piss1 == "medium")
          ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = ep.Y_1 = ep.Cb_1 = ep.Cr_1 = ep.Y_2 = ep.Cb_2 = ep.Cr_2 = (scale) ? 192 : 48;
      else if (piss1 == "high")
          ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = ep.Y_1 = ep.Cb_1 = ep.Cr_1 = ep.Y_2 = ep.Cb_2 = ep.Cr_2 = (scale) ? 256 : 64;
      else if (piss1 == "veryhigh")
          ep.Y = ep.Cb = ep.Cr = ep.grainY = ep.grainC = ep.Y_1 = ep.Cb_1 = ep.Cr_1 = ep.Y_2 = ep.Cb_2 = ep.Cr_2 = (scale) ? 320 : 80;
      else if (piss1 == "nograin")
        ep.grainY = ep.grainC = 0;
      else if (piss1 == "luma")
        ep.Cb = ep.Cr = ep.grainC = 0;
      else if (piss1 == "chroma")
        ep.Y = ep.grainY = 0;
    }
  }

  if (context.params) {
    ep.range = get_value(context.params->get_int("range", ep.range), ep.range);
    ep.Y = get_value(context.params->get_int("y", ep.Y), ep.Y);
    ep.Cb = get_value(context.params->get_int("cb", ep.Cb), ep.Cb);
    ep.Cr = get_value(context.params->get_int("cr", ep.Cr), ep.Cr);
    ep.grainY = get_value(context.params->get_int("grainy", ep.grainY), ep.grainY);
    ep.grainC = get_value(context.params->get_int("grainc", ep.grainC), ep.grainC);
    ep.sample_mode = get_value(context.params->get_int("sample_mode", ep.sample_mode), ep.sample_mode);
    ep.seed = get_value(context.params->get_int("seed", ep.seed), ep.seed);
    ep.blur_first = get_value(context.params->get_bool("blur_first", ep.blur_first), ep.blur_first);
    ep.dynamic_grain = get_value(context.params->get_bool("dynamic_grain", ep.dynamic_grain), ep.dynamic_grain);
    ep.keep_tv_range = get_value(context.params->get_bool("keep_tv_range", ep.keep_tv_range), ep.keep_tv_range);
    ep.output_depth = get_value(context.params->get_int("output_depth", ep.output_depth), ep.output_depth);
    ep.random_algo_ref = static_cast<RANDOM_ALGORITHM>(
      get_value(
        context.params->get_int("random_algo_ref", static_cast<int>(ep.random_algo_ref)),
        static_cast<int>(ep.random_algo_ref)
      )
    );
    ep.random_algo_grain = static_cast<RANDOM_ALGORITHM>(
      get_value(
        context.params->get_int("random_algo_grain", static_cast<int>(ep.random_algo_grain)),
        static_cast<int>(ep.random_algo_grain)
      )
    );
    ep.random_param_ref = get_value(context.params->get_double("random_param_ref", ep.random_param_ref), ep.random_param_ref);
    ep.random_param_grain = get_value(context.params->get_double("random_param_grain", ep.random_param_grain), ep.random_param_grain);
    ep.Y_1 = get_value(context.params->get_int("y_1", ep.Y_1), ep.Y_1);
    ep.Cb_1 = get_value(context.params->get_int("cb_1", ep.Cb_1), ep.Cb_1);
    ep.Cr_1 = get_value(context.params->get_int("cr_1", ep.Cr_1), ep.Cr_1);
    ep.Y_2 = get_value(context.params->get_int("y_2", ep.Y_2), ep.Y_2);
    ep.Cb_2 = get_value(context.params->get_int("cb_2", ep.Cb_2), ep.Cb_2);
    ep.Cr_2 = get_value(context.params->get_int("cr_2", ep.Cr_2), ep.Cr_2);
    ep.angle_boost = get_value(context.params->get_double("angle_boost", ep.angle_boost), ep.angle_boost);
    ep.max_angle = get_value(context.params->get_double("max_angle", ep.max_angle), ep.max_angle);
    ep.dither_algo = static_cast<DITHER_ALGORITHM>(
      get_value(
        context.params->get_int("dither_algo", static_cast<int>(ep.dither_algo)),
        static_cast<int>(ep.dither_algo)
      )
    );
  }

  ep.Y_1 = ep.Y_1 == -1 ? ep.Y : ep.Y_1;
  ep.Cb_1 = ep.Cb_1 == -1 ? ep.Cb : ep.Cb_1;
  ep.Cr_1 = ep.Cr_1 == -1 ? ep.Cr : ep.Cr_1;
  ep.Y_2 = ep.Y_2 == -1 ? ep.Y : ep.Y_2;
  ep.Cb_2 = ep.Cb_2 == -1 ? ep.Cb : ep.Cb_2;
  ep.Cr_2 = ep.Cr_2 == -1 ? ep.Cr : ep.Cr_2;

  // CPU dispatch
  OPTIMIZATION_MODE opt = [&]() {
      const int CPUFlags = GetCPUFlags();
      if (ep.sample_mode >= 5 && ep.sample_mode <= 7) {
          const int AVX512_REQUIRED_FLAGS = CPUF_AVX512F | CPUF_AVX512BW | CPUF_AVX512DQ | CPUF_AVX512VL | CPUF_AVX512CD;
          if (((CPUFlags & AVX512_REQUIRED_FLAGS) == AVX512_REQUIRED_FLAGS) && (opt_in == 3 || opt_in < 0))
              return IMPL_AVX512;
          if ((CPUFlags & CPUF_AVX2) && (opt_in == 2 || opt_in < 0))
              return IMPL_AVX2;
      }
      if ((CPUFlags & CPUF_SSE4_1) && (opt_in > 0 || opt_in < 0))
          return IMPL_SSE4;
      return IMPL_C;
  }();

  // Perform original parameter validation
  try {
    INVALID_PARAM_IF(input_vi.format.color_family != ds::ColorFamily::Yuv);
    INVALID_PARAM_IF(input_vi.width < 16);
    INVALID_PARAM_IF(input_vi.height < 16);
    INVALID_PARAM_IF(input_vi.num_frames <= 0);

    int bits = ds::bits_per_sample(input_vi.format.sample_format);
    INVALID_PARAM_IF(bits < 8 || bits > INTERNAL_BIT_DEPTH);
    INVALID_PARAM_IF(input_vi.format.sample_format == ds::SampleFormat::Float32);

    if (ep.output_depth < 0)
      ep.output_depth = bits;
    if (ep.output_depth == 16)
        ep.dither_algo = DA_16BIT_INTERLEAVED;

    const int y_threshold_upper_limit = scale ? 65535 : 511;
    const int cb_threshold_upper_limit = scale ? 65535 : 511;
    const int cr_threshold_upper_limit = scale ? 65535 : 511;
    constexpr int dither_upper_limit = 4096;

    CHECK_PARAM(ep.range, 0, 255);
    CHECK_PARAM(ep.Y, 0, y_threshold_upper_limit);
    CHECK_PARAM(ep.Cb, 0, cb_threshold_upper_limit);
    CHECK_PARAM(ep.Cr, 0, cr_threshold_upper_limit);
    CHECK_PARAM(ep.grainY, 0, dither_upper_limit);
    CHECK_PARAM(ep.grainC, 0, dither_upper_limit);
    CHECK_PARAM(ep.sample_mode, 1, 7);
    CHECK_PARAM(ep.dither_algo, DA_HIGH_NO_DITHERING, (DA_COUNT - 1));
    CHECK_PARAM(ep.random_algo_ref, 0, (RANDOM_ALGORITHM_COUNT - 1));
    CHECK_PARAM(ep.random_algo_grain, 0, (RANDOM_ALGORITHM_COUNT - 1));
    CHECK_PARAM(ep.Y_1, 0, y_threshold_upper_limit);
    CHECK_PARAM(ep.Cb_1, 0, cb_threshold_upper_limit);
    CHECK_PARAM(ep.Cr_1, 0, cr_threshold_upper_limit);
    CHECK_PARAM(ep.Y_2, 0, y_threshold_upper_limit);
    CHECK_PARAM(ep.Cb_2, 0, cb_threshold_upper_limit);
    CHECK_PARAM(ep.Cr_2, 0, cr_threshold_upper_limit);

    if (ep.angle_boost < 0.0f)
        throw std::invalid_argument("invalid parameter angle_boost, must be positive value");
    if (ep.max_angle < 0.0f || ep.max_angle > 1.0f)
        throw std::invalid_argument("invalid parameter max_angle, must be between 0.0 and 1.0");

    ep.Y = scale ? ep.Y : ep.Y << 2;
    ep.Cb = scale ? ep.Cb : ep.Cb << 2;
    ep.Cr = scale ? ep.Cr : ep.Cr << 2;
    ep.Y_1 = scale ? ep.Y_1 : ep.Y_1 << 2;
    ep.Cb_1 = scale ? ep.Cb_1 : ep.Cb_1 << 2;
    ep.Cr_1 = scale ? ep.Cr_1 : ep.Cr_1 << 2;
    ep.Y_2 = scale ? ep.Y_2 : ep.Y_2 << 2;
    ep.Cb_2 = scale ? ep.Cb_2 : ep.Cb_2 << 2;
    ep.Cr_2 = scale ? ep.Cr_2 : ep.Cr_2 << 2;
    ep.grainY <<= 2;
    ep.grainC <<= 2;

  } catch (const std::exception& e) {
    return ds::Result<ds::VideoInitStateResult<State>>::failure(
      ds::Error{ds::ErrorCode::InvalidArgument, std::string("Neo-F3KDB Parameter Validation Error: ") + e.what()}
    );
  }

  ds::VideoOutputInfo output_vi{
    .width = input_vi.width,
    .height = input_vi.height,
    .num_frames = input_vi.num_frames,
    .format = ds::VideoFormat{
      .color_family = input_vi.format.color_family,
      .sample_format = ep.output_depth == 8 ? ds::SampleFormat::UInt8 : ds::SampleFormat::UInt16,
      .plane_count = input_vi.format.plane_count,
      .subsampling_w = input_vi.format.subsampling_w,
      .subsampling_h = input_vi.format.subsampling_h,
    },
    .fps = input_vi.fps
  };

  DSVideoInfo old_vi = make_ds_video_info(input_vi);
  std::unique_ptr<f3kdb_core_t> engine;
  try {
    engine = std::make_unique<f3kdb_core_t>(old_vi, ep, opt);
  } catch (const std::exception& e) {
    return ds::Result<ds::VideoInitStateResult<State>>::failure(
      ds::Error{ds::ErrorCode::InternalError, std::string("Neo-F3KDB: memory allocation failed: ") + e.what()}
    );
  }

  return ds::Result<ds::VideoInitStateResult<State>>::success(
    ds::VideoInitStateResult<State>{
      .output = output_vi,
      .state = State(std::move(engine), ep, mt)
    }
  );
}

ds::Result<ds::VideoRequestResult> F3KDBFilterCore::request(ds::VideoRequestContext& context) {
  context.request_frame(0, context.output_frame);
  return ds::Result<ds::VideoRequestResult>::success(ds::VideoRequestResult{});
}

ds::Result<ds::VideoProcessResult> F3KDBFilterCore::process(ds::VideoProcessContext& context) {
  auto& filter_state = context.state<State>();

  auto frame_res = context.frames.get(0, context.output_frame);
  if (!frame_res.has_value()) {
    return ds::Result<ds::VideoProcessResult>::failure(frame_res.error());
  }
  const auto src_frame = frame_res.value().frame;

  for (int p = 0; p < src_frame.plane_count; ++p) {
    const auto& src_plane = src_frame.plane(p);
    auto& dst_plane = context.dst.plane(p);

    const int src_stride = static_cast<int>(src_plane.stride_bytes);
    const int dst_stride = static_cast<int>(dst_plane.stride_bytes);

    const unsigned char* src_ptr = reinterpret_cast<const unsigned char*>(src_plane.data);
    unsigned char* dst_ptr = reinterpret_cast<unsigned char*>(dst_plane.data);

    filter_state.engine->process_plane(
      context.output_frame,
      p,
      dst_ptr,
      dst_stride,
      src_ptr,
      src_stride
    );
  }

  return ds::Result<ds::VideoProcessResult>::success(ds::VideoProcessResult{});
}

bool F3KDBBridge::accepts_video_format(ds::VideoFormat format) {
  return format.color_family == ds::ColorFamily::Yuv &&
         (format.sample_format == ds::SampleFormat::UInt8 || format.sample_format == ds::SampleFormat::UInt16);
}

ds::FilterDescriptor F3KDBBridge::descriptor() {
  return ds::FilterDescriptor{
    .name = "Deband",
    .params = {
      ds::ParamSpec{"clip", ds::ParamType::Clip, {}, true},
      ds::ParamSpec{"range", ds::ParamType::Integer, 15},
      ds::ParamSpec{"y", ds::ParamType::Integer, 64},
      ds::ParamSpec{"cb", ds::ParamType::Integer, 64},
      ds::ParamSpec{"cr", ds::ParamType::Integer, 64},
      ds::ParamSpec{"grainy", ds::ParamType::Integer, 64},
      ds::ParamSpec{"grainc", ds::ParamType::Integer, 64},
      ds::ParamSpec{"sample_mode", ds::ParamType::Integer, 2},
      ds::ParamSpec{"seed", ds::ParamType::Integer, 0},
      ds::ParamSpec{"blur_first", ds::ParamType::Boolean, true},
      ds::ParamSpec{"dynamic_grain", ds::ParamType::Boolean, false},
      ds::ParamSpec{"opt", ds::ParamType::Integer, -1},
      ds::ParamSpec{"mt", ds::ParamType::Boolean, true},
      ds::ParamSpec{"dither_algo", ds::ParamType::Integer, 3},
      ds::ParamSpec{"keep_tv_range", ds::ParamType::Boolean, false},
      ds::ParamSpec{"output_depth", ds::ParamType::Integer, -1},
      ds::ParamSpec{"random_algo_ref", ds::ParamType::Integer, 1},
      ds::ParamSpec{"random_algo_grain", ds::ParamType::Integer, 1},
      ds::ParamSpec{"random_param_ref", ds::ParamType::Float, 1.0},
      ds::ParamSpec{"random_param_grain", ds::ParamType::Float, 1.0},
      ds::ParamSpec{"preset", ds::ParamType::String, ""},
      ds::ParamSpec{"y_1", ds::ParamType::Integer, -1},
      ds::ParamSpec{"cb_1", ds::ParamType::Integer, -1},
      ds::ParamSpec{"cr_1", ds::ParamType::Integer, -1},
      ds::ParamSpec{"y_2", ds::ParamType::Integer, -1},
      ds::ParamSpec{"cb_2", ds::ParamType::Integer, -1},
      ds::ParamSpec{"cr_2", ds::ParamType::Integer, -1},
      ds::ParamSpec{"scale", ds::ParamType::Boolean, false},
      ds::ParamSpec{"angle_boost", ds::ParamType::Float, 1.5},
      ds::ParamSpec{"max_angle", ds::ParamType::Float, 0.15}
    }
  };
}

} // namespace neo_f3kdb
