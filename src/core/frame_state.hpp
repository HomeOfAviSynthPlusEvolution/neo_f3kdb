#pragma once

#include "core/aligned_buffer.hpp"
#include "core/parameters.hpp"
#include "core/plane.hpp"

#include <dualsynth/video_filter.hpp>

#include <array>
#include <cstdint>
#include <span>
#include <vector>

namespace neo_f3kdb::core {

class FrameState {
public:
  FrameState() = default;
  FrameState(const ds::VideoInputInfo& input, const DebandParameters& params);
  ~FrameState();

  FrameState(const FrameState&) = delete;
  FrameState& operator=(const FrameState&) = delete;
  FrameState(FrameState&&) noexcept = default;
  FrameState& operator=(FrameState&&) noexcept = default;

  void reset(const ds::VideoInputInfo& input, const DebandParameters& params);

  [[nodiscard]] int info_stride(int plane) const noexcept;
  [[nodiscard]] int grain_stride(int plane) const noexcept;

  [[nodiscard]] std::span<pixel_dither_info> dither_info(int plane) noexcept;
  [[nodiscard]] std::span<std::int16_t> grain_buffer(int plane) noexcept;
  [[nodiscard]] std::int16_t* grain_row_base(int plane, int frame_number, int frame_count) noexcept;

  [[nodiscard]] process_plane_context* context(int plane) noexcept;

private:
  void destroy_contexts() noexcept;

  AlignedBuffer<pixel_dither_info, 128> y_info_;
  AlignedBuffer<pixel_dither_info, 128> cb_info_;
  AlignedBuffer<pixel_dither_info, 128> cr_info_;

  std::array<process_plane_context, 3> contexts_{};

  AlignedBuffer<std::int16_t, 128> grain_buffer_y_;
  AlignedBuffer<std::int16_t, 128> grain_buffer_c_;
  std::vector<std::int32_t> grain_buffer_offsets_;

  int y_info_stride_ = 0;
  int c_info_stride_ = 0;
  int y_grain_stride_ = 0;
  int c_grain_stride_ = 0;
};

[[nodiscard]] int frame_lut_stride(int width_in_pixels);

} // namespace neo_f3kdb::core
