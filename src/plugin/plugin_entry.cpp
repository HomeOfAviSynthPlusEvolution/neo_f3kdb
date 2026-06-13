#include <avisynth.h>
#include <VapourSynth4.h>

#include <dualsynth/avisynth/video_bridge.hpp>
#include <dualsynth/vapoursynth/video_bridge.hpp>
#include <dualsynth/video_bridge.hpp>

#include "plugin/f3kdb_filter.hpp"

#include <array>
#include <exception>
#include <string>
#include <utility>

#if defined(_WIN32)
#define F3KDB_AVS_PLUGIN_EXPORT extern "C" __declspec(dllexport)
#elif defined(__clang__) || defined(__GNUC__)
#define F3KDB_AVS_PLUGIN_EXPORT extern "C" __attribute__((visibility("default")))
#else
#define F3KDB_AVS_PLUGIN_EXPORT extern "C"
#endif

const AVS_Linkage* AVS_linkage = nullptr;

namespace {

bool set_avisynth_host_var(void* user, const char* name, const ds::ParamValue& value) {
  if (user == nullptr || name == nullptr) {
    return false;
  }

  auto* env = static_cast<IScriptEnvironment*>(user);
  AVSValue avs_value;
  if (!ds::avisynth::assign_avisynth_value(env, value, avs_value)) {
    return false;
  }

  env->SetVar(env->SaveString(name), avs_value);
  return true;
}

ds::HostVariableCallbacks avisynth_host_variable_callbacks(IScriptEnvironment* env) {
  if (env == nullptr) {
    return {};
  }
  return ds::HostVariableCallbacks{.user = env, .set = &set_avisynth_host_var};
}

template <class Bridge>
const char* avs_signature() {
  static const std::string signature = [] {
    const auto generated = ds::make_avisynth_signature(Bridge::descriptor());
    if (!generated.has_value()) {
      return std::string{
        "c"
        "[range]i[y]i[cb]i[cr]i[grainy]i[grainc]i[sample_mode]i[seed]i"
        "[blur_first]b[dynamic_grain]b[opt]i[mt]b[dither_algo]i[keep_tv_range]b"
        "[output_depth]i[random_algo_ref]i[random_algo_grain]i"
        "[random_param_ref]f[random_param_grain]f[preset]s"
        "[y_1]i[cb_1]i[cr_1]i[y_2]i[cb_2]i[cr_2]i[scale]b"
        "[angle_boost]f[max_angle]f"
      };
    }
    return generated.value();
  }();
  return signature.c_str();
}

const char* vs_signature() {
  return
    "clip:vnode;"
    "range:int:opt;"
    "y:int:opt;"
    "cb:int:opt;"
    "cr:int:opt;"
    "grainy:int:opt;"
    "grainc:int:opt;"
    "sample_mode:int:opt;"
    "seed:int:opt;"
    "blur_first:int:opt;"
    "dynamic_grain:int:opt;"
    "opt:int:opt;"
    "mt:int:opt;"
    "dither_algo:int:opt;"
    "keep_tv_range:int:opt;"
    "output_depth:int:opt;"
    "random_algo_ref:int:opt;"
    "random_algo_grain:int:opt;"
    "random_param_ref:float:opt;"
    "random_param_grain:float:opt;"
    "preset:data:opt;"
    "y_1:int:opt;"
    "cb_1:int:opt;"
    "cr_1:int:opt;"
    "y_2:int:opt;"
    "cb_2:int:opt;"
    "cr_2:int:opt;"
    "scale:int:opt;"
    "angle_boost:float:opt;"
    "max_angle:float:opt;";
}

template <class Bridge>
void VS_CC create_vapoursynth_filter(
  const VSMap* in,
  VSMap* out,
  void* /*user_data*/,
  VSCore* core,
  const VSAPI* vsapi
) {
  ds::vapoursynth::create_video_filter_bridge<Bridge>(in, out, core, vsapi);
}

template <class Bridge>
class AvisynthCompatVideoFilter final : public IClip {
public:
  using Filter = typename Bridge::Core;

  AvisynthCompatVideoFilter(
    const PClip& clip,
    ds::VideoInputInfo input_info,
    ds::VideoFilterState<Filter> state,
    ds::avisynth::MtMode mt_mode
  ) : clip_(clip),
      input_info_(input_info),
      state_(std::move(state)),
      mt_mode_(mt_mode),
      vi_(clip_->GetVideoInfo()) {}

  PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env) override {
    try {
      PVideoFrame src = clip_->GetFrame(n, env);
      PVideoFrame dst = src;
      bool inplace = env->MakeWritable(&dst);
      if (!inplace) {
        dst = env->NewVideoFrame(vi_);
        ds::avisynth::copy_video_frame_pixels(src, dst, input_info_.format);
      }

      const auto src_frame_view = ds::avisynth::make_video_frame_view(src, input_info_.format);
      std::array<ds::RequestedVideoFrame, 1> frames{
        ds::RequestedVideoFrame{.input_index = 0, .frame_number = n, .frame = src_frame_view}
      };
      ds::RequestedVideoFrameProvider provider(frames);

      const auto result = ds::process_video_filter<Filter>(
        n,
        provider,
        ds::avisynth::make_mutable_video_frame_view(dst, input_info_.format),
        state_
      );
      if (!result.has_value()) {
        env->ThrowError(result.error().message.c_str());
      }
      return dst;
    } catch (const AvisynthError&) {
      throw;
    } catch (const std::exception& error) {
      env->ThrowError(error.what());
    } catch (...) {
      env->ThrowError("Neo-F3KDB: unhandled exception in AviSynth video wrapper");
    }
    return {};
  }

  bool __stdcall GetParity(int n) override {
    return clip_->GetParity(n);
  }

  void __stdcall GetAudio(void* buf, int64_t start, int64_t count, IScriptEnvironment* env) override {
    clip_->GetAudio(buf, start, count, env);
  }

  int __stdcall SetCacheHints(int cachehints, int frame_range) override {
    return ds::cache_hints_video_filter<Filter>(
      cachehints,
      frame_range,
      ds::avisynth::cache_hint_response(cachehints, frame_range, mt_mode_),
      state_
    );
  }

  const VideoInfo& __stdcall GetVideoInfo() override {
    return vi_;
  }

private:
  PClip clip_;
  ds::VideoInputInfo input_info_;
  ds::VideoFilterState<Filter> state_;
  ds::avisynth::MtMode mt_mode_;
  VideoInfo vi_{};
};

template <class Bridge>
AVSValue __cdecl create_avisynth_filter(
  AVSValue args,
  void* /*user_data*/,
  IScriptEnvironment* env
) {
  using Filter = typename Bridge::Core;

  try {
    PClip clip = args[0].AsClip();
    const VideoInfo& vi = clip->GetVideoInfo();
    const auto format = ds::avisynth::make_video_format(vi);
    if (!vi.HasVideo() || !format.has_value() || !Bridge::accepts_video_format(format.value())) {
      env->ThrowError(Bridge::avs_format_error);
    }

    std::array<ds::VideoInputInfo, 1> input_infos{
      ds::VideoInputInfo{
        .width = vi.width,
        .height = vi.height,
        .num_frames = vi.num_frames,
        .format = format.value(),
        .fps = ds::FrameRate{.numerator = vi.fps_numerator, .denominator = vi.fps_denominator}
      }
    };

    auto params = ds::avisynth::read_params(args, Bridge::descriptor());
    if (!params.has_value()) {
      env->ThrowError(params.error().message.c_str());
    }

    auto init_result = ds::init_video_filter_instance<Filter>(
      input_infos,
      &params.value(),
      ds::avisynth::host_global_lock_callbacks(env),
      avisynth_host_variable_callbacks(env)
    );
    if (!init_result.has_value()) {
      env->ThrowError(init_result.error().message.c_str());
    }

    const ds::VideoOutputInfo& output = init_result.value().output;
    if (
      output.width != vi.width ||
      output.height != vi.height ||
      output.num_frames != vi.num_frames ||
      output.format != format.value()
    ) {
      env->ThrowError("Neo-F3KDB: AviSynth compatibility bridge requires source-shaped output");
    }

    return new AvisynthCompatVideoFilter<Bridge>(
      std::move(clip),
      input_infos[0],
      std::move(init_result.value().state),
      ds::avisynth::bridge_mt_mode<Bridge>()
    );
  } catch (const AvisynthError&) {
    throw;
  } catch (const std::exception& error) {
    env->ThrowError(error.what());
  } catch (...) {
    env->ThrowError("Neo-F3KDB: unhandled exception in AviSynth creation");
  }

  return {};
}

template <class Bridge>
void register_avisynth_filter(IScriptEnvironment* env, bool register_mt_mode) {
  env->AddFunction(
    Bridge::avs_name,
    avs_signature<Bridge>(),
    create_avisynth_filter<Bridge>,
    nullptr
  );
  if (register_mt_mode) {
    ds::avisynth::set_video_filter_mt_mode<Bridge>(env);
  }
}

const char* register_avisynth_filters(IScriptEnvironment* env, bool register_mt_mode) {
  register_avisynth_filter<neo_f3kdb::F3KDBBridge>(env, register_mt_mode);
  return neo_f3kdb::Plugin::Description;
}

} // namespace

VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin* plugin, const VSPLUGINAPI* vspapi) {
  vspapi->configPlugin(
    neo_f3kdb::Plugin::Identifier,
    neo_f3kdb::Plugin::Namespace,
    neo_f3kdb::Plugin::Description,
    VS_MAKE_VERSION(10, 0),
    VAPOURSYNTH_API_VERSION,
    0,
    plugin
  );

  vspapi->registerFunction(
    neo_f3kdb::F3KDBBridge::vs_name,
    vs_signature(),
    "clip:vnode;",
    create_vapoursynth_filter<neo_f3kdb::F3KDBBridge>,
    nullptr,
    plugin
  );
}

F3KDB_AVS_PLUGIN_EXPORT const char* __stdcall AvisynthPluginInit2(IScriptEnvironment* env) {
  AVS_linkage = env->GetAVSLinkage();
  return register_avisynth_filters(env, false);
}

F3KDB_AVS_PLUGIN_EXPORT const char* __stdcall AvisynthPluginInit3(
  IScriptEnvironment* env,
  const AVS_Linkage* const vectors
) {
  AVS_linkage = vectors;
  return register_avisynth_filters(env, true);
}
