#pragma once

#include <dualsynth/avisynth/video_bridge.hpp>
#include <dualsynth/format.hpp>
#include <dualsynth/param.hpp>
#include <dualsynth/video_bridge.hpp>
#include <dualsynth/video_filter.hpp>

#include "core.h"
#include "f3kdb.h"
#include "version.hpp"

#include <memory>

namespace neo_f3kdb {

struct F3KDBFilterCore {
  static constexpr const char* name = "Deband";
  static constexpr int input_count = 1;
  static constexpr ds::OutputOrigin output_origin = ds::OutputOrigin::fresh();

  struct State {
    State() = default;
    State(std::unique_ptr<f3kdb_core_t> engine, f3kdb_params_t params, bool mt);
    ~State() = default;

    State(State&&) noexcept;
    State& operator=(State&&) noexcept;
    State(const State&) = delete;
    State& operator=(const State&) = delete;

    std::unique_ptr<f3kdb_core_t> engine;
    f3kdb_params_t params;
    bool mt = true;
  };

  static ds::Result<ds::VideoInitStateResult<State>> init(ds::VideoInitContext& context);
  static ds::Result<ds::VideoRequestResult> request(ds::VideoRequestContext& context);
  static ds::Result<ds::VideoProcessResult> process(ds::VideoProcessContext& context);
};

struct F3KDBBridge : ds::SingleInputVideoBridgeDefaults<F3KDBFilterCore> {
  static constexpr const char* vs_name = "Deband";
  static constexpr const char* avs_name = "neo_f3kdb";
  static constexpr const char* missing_input_error = "Neo-F3KDB: missing required video clip";
  static constexpr const char* vs_format_error =
    "Neo-F3KDB: only integer YUV format is supported";
  static constexpr const char* avs_format_error =
    "Neo-F3KDB: only integer YUV format is supported";
  static constexpr ds::avisynth::MtMode avs_mt_mode = ds::avisynth::MtMode::NiceFilter;

  static bool accepts_video_format(ds::VideoFormat format);
  static ds::FilterDescriptor descriptor();
};

namespace Plugin {
inline constexpr const char* Identifier = "in.7086.neo_f3kdb";
inline constexpr const char* Namespace = "neo_f3kdb";
inline constexpr const char* Description = "Neo F3KDB Deband Filter " PLUGIN_VERSION;
} // namespace Plugin

} // namespace neo_f3kdb
