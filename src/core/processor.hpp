#pragma once

#include "core/frame_state.hpp"
#include "core/parameters.hpp"

#include <dualsynth/frame.hpp>
#include <dualsynth/video_filter.hpp>

namespace neo_f3kdb::core {

struct DebandProcessorConfig {
  DebandParameters params{};
  ds::VideoInputInfo input{};
  Backend requested_backend = Backend::Highway;
};

class DebandProcessor {
public:
  explicit DebandProcessor(const DebandProcessorConfig& config);

  [[nodiscard]] ds::VideoOutputInfo output_info() const;

  void process(
    ds::MutableVideoFrameView& dst,
    const ds::VideoFrameView& src,
    int frame_number
  );

private:
  [[nodiscard]] PlaneJob make_plane_job(
    int frame_number,
    int plane,
    unsigned char* dst_frame_ptr,
    int dst_pitch,
    const unsigned char* src_frame_ptr,
    int src_pitch
  );

  [[nodiscard]] bool should_copy_plane(const PlaneJob& job, int grain_setting) const;
  void copy_plane(const PlaneJob& job) const;

  DebandParameters params_;
  ds::VideoInputInfo input_;
  ds::VideoOutputInfo output_;
  FrameState frame_state_;
  Backend requested_backend_ = Backend::Highway;
};

} // namespace neo_f3kdb::core
