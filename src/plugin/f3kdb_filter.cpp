#include "plugin/f3kdb_filter.hpp"

#include <exception>
#include <string>
#include <utility>
#include <vector>

namespace neo_f3kdb {

F3KDBFilterCore::State::State(
  std::unique_ptr<core::DebandProcessor> input_processor,
  core::DebandParameters input_params,
  bool input_mt
) : processor(std::move(input_processor)),
    params(input_params),
    mt(input_mt) {}

F3KDBFilterCore::State::State(State&&) noexcept = default;
F3KDBFilterCore::State& F3KDBFilterCore::State::operator=(State&&) noexcept = default;

ds::Result<ds::VideoInitStateResult<F3KDBFilterCore::State>> F3KDBFilterCore::init(ds::VideoInitContext& context) {
  if (context.inputs.size() != 1) {
    return ds::Result<ds::VideoInitStateResult<State>>::failure(
      ds::Error{ds::ErrorCode::InvalidArgument, "Neo-F3KDB: exactly one video input is required"}
    );
  }

  const ds::VideoInputInfo& input_vi = context.inputs[0];
  auto parsed_result = core::parse_parameters(context.params);
  if (!parsed_result.has_value()) {
    return ds::Result<ds::VideoInitStateResult<State>>::failure(parsed_result.error());
  }

  auto parsed = parsed_result.value();
  const int bits = ds::bits_per_sample(input_vi.format.sample_format);
  try {
    core::resolve_parameter_defaults(parsed.params, bits);
    core::validate_parameters(parsed.params, input_vi, parsed.scale);
    core::normalize_parameters(parsed.params, bits, parsed.scale);
  } catch (const std::exception& e) {
    return ds::Result<ds::VideoInitStateResult<State>>::failure(
      ds::Error{ds::ErrorCode::InvalidArgument, std::string("Neo-F3KDB Parameter Validation Error: ") + e.what()}
    );
  }

  std::unique_ptr<core::DebandProcessor> processor;
  try {
    core::DebandProcessorConfig config{
      parsed.params,
      input_vi,
      core::select_backend(parsed.params, parsed.opt)
    };
    processor = std::make_unique<core::DebandProcessor>(config);
  } catch (const std::exception& e) {
    return ds::Result<ds::VideoInitStateResult<State>>::failure(
      ds::Error{ds::ErrorCode::InternalError, std::string("Neo-F3KDB: memory allocation failed: ") + e.what()}
    );
  }

  ds::VideoInitStateResult<State> init_state{};
  init_state.output = processor->output_info();
  init_state.state = State(std::move(processor), parsed.params, parsed.mt);

  return ds::Result<ds::VideoInitStateResult<State>>::success(
    std::move(init_state)
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

  filter_state.processor->process(context.dst, src_frame, context.output_frame);

  return ds::Result<ds::VideoProcessResult>::success(ds::VideoProcessResult{});
}

bool F3KDBBridge::accepts_video_format(ds::VideoFormat format) {
  return (format.color_family == ds::ColorFamily::Yuv || format.color_family == ds::ColorFamily::Gray) &&
         (format.sample_format == ds::SampleFormat::UInt8 ||
          format.sample_format == ds::SampleFormat::UInt16 ||
          format.sample_format == ds::SampleFormat::Float32);
}

ds::FilterDescriptor F3KDBBridge::descriptor() {
  ds::FilterDescriptor desc{};
  desc.name = "Deband";
  desc.params = {
    ds::ParamSpec{"clip", ds::ParamType::Clip, {}, true},
    ds::ParamSpec{"range", ds::ParamType::Integer, 15},
    ds::ParamSpec{"y", ds::ParamType::Float, 64.0},
    ds::ParamSpec{"cb", ds::ParamType::Float, 64.0},
    ds::ParamSpec{"cr", ds::ParamType::Float, 64.0},
    ds::ParamSpec{"grainy", ds::ParamType::Float, 64.0},
    ds::ParamSpec{"grainc", ds::ParamType::Float, 64.0},
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
    ds::ParamSpec{"y_1", ds::ParamType::Float, -1.0},
    ds::ParamSpec{"cb_1", ds::ParamType::Float, -1.0},
    ds::ParamSpec{"cr_1", ds::ParamType::Float, -1.0},
    ds::ParamSpec{"y_2", ds::ParamType::Float, -1.0},
    ds::ParamSpec{"cb_2", ds::ParamType::Float, -1.0},
    ds::ParamSpec{"cr_2", ds::ParamType::Float, -1.0},
    ds::ParamSpec{"scale", ds::ParamType::Boolean, false},
    ds::ParamSpec{"angle_boost", ds::ParamType::Float, 1.5},
    ds::ParamSpec{"max_angle", ds::ParamType::Float, 0.15}
  };
  return desc;
}

} // namespace neo_f3kdb
