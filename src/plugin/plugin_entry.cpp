#include <avisynth.h>
#include <VapourSynth4.h>

#include <dualsynth/avisynth/video_bridge.hpp>
#include <dualsynth/vapoursynth/video_bridge.hpp>
#include <dualsynth/video_bridge.hpp>

#include "plugin/f3kdb_filter.hpp"

#include <string>

#if defined(_WIN32)
#define F3KDB_AVS_PLUGIN_EXPORT extern "C" __declspec(dllexport)
#elif defined(__clang__) || defined(__GNUC__)
#define F3KDB_AVS_PLUGIN_EXPORT extern "C" __attribute__((visibility("default")))
#else
#define F3KDB_AVS_PLUGIN_EXPORT extern "C"
#endif

#if defined(_MSC_VER) && !defined(_M_ARM64) && !defined(__aarch64__)
const AVS_Linkage* AVS_linkage = nullptr;
#endif

namespace {

const char* vs_signature() {
  return
    "clip:vnode;"
    "range:int:opt;"
    "y:float:opt;"
    "cb:float:opt;"
    "cr:float:opt;"
    "grainy:float:opt;"
    "grainc:float:opt;"
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
    "y_1:float:opt;"
    "cb_1:float:opt;"
    "cr_1:float:opt;"
    "y_2:float:opt;"
    "cb_2:float:opt;"
    "cr_2:float:opt;"
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

#if defined(_MSC_VER) && !defined(_M_ARM64) && !defined(__aarch64__)
template <class Bridge>
AVSValue __cdecl create_avisynth_filter(
  AVSValue args,
  void* /*user_data*/,
  IScriptEnvironment* env
) {
  return ds::avisynth::create_video_filter_bridge<Bridge>(args, env);
}

template <class Bridge>
void register_avisynth_filter(IScriptEnvironment* env, bool register_mt_mode) {
  env->AddFunction(
    Bridge::avs_name,
    ds::avisynth::bridge_avs_signature<Bridge>(),
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
#endif

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

#if defined(_MSC_VER) && !defined(_M_ARM64) && !defined(__aarch64__)
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
#endif

F3KDB_AVS_PLUGIN_EXPORT const char* AVSC_CC avisynth_c_plugin_init2(
  AVS_ScriptEnvironment* env
) {
  ds::avisynth::c::register_video_filter<neo_f3kdb::F3KDBBridge>(env);
  return neo_f3kdb::Plugin::Description;
}
