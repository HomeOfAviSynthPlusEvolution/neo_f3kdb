#pragma once

#include "compiler_compat.h"
#include "core/kernel_types.hpp"

#include <cstdint>
#include <cstddef>
#include <span>

namespace neo_f3kdb::core {

template <class T>
struct StridedPlaneView {
  T* data = nullptr;
  int width = 0;
  int height = 0;
  int stride = 0;

  [[nodiscard]] std::span<T> row(int y) const noexcept {
    return std::span<T>{
      data + static_cast<std::ptrdiff_t>(y) * stride,
      static_cast<std::size_t>(width)
    };
  }
};

} // namespace neo_f3kdb::core

typedef struct _pixel_dither_info {
  alignas(4) std::int8_t ref1;
  std::int8_t ref2;
  std::int16_t change;
} pixel_dither_info;

static_assert(sizeof(pixel_dither_info) == 4, "Something wrong in pixel_dither_info");

namespace neo_f3kdb::core {

struct PlaneGeometry {
  int width_pixels = 0;
  int height_pixels = 0;
};

struct PlaneBuffers {
  StridedPlaneView<const unsigned char> src_bytes{};
  StridedPlaneView<unsigned char> dst_bytes{};
  StridedPlaneView<const pixel_dither_info> dither_info{};
  StridedPlaneView<const std::int16_t> grain{};
};

} // namespace neo_f3kdb::core

typedef void (*destroy_data_t)(void* data);

typedef struct _process_plane_context {
  void* data;
  destroy_data_t destroy;
} process_plane_context;

void destroy_context(process_plane_context* context);
void init_context(process_plane_context* context);

typedef struct _process_plane_params {
  neo_f3kdb::core::PlaneGeometry geometry{};
  neo_f3kdb::core::PlaneBuffers buffers{};

  PIXEL_MODE input_mode;
  int input_depth;
  PIXEL_MODE output_mode;
  int output_depth;
  DITHER_ALGORITHM dither_algo;
  bool blur_first;
  int sample_mode;

  std::uint16_t threshold;
  std::uint16_t threshold1;
  std::uint16_t threshold2;
  float angle_boost;
  float max_angle;

  int plane;

  unsigned char width_subsampling;
  unsigned char height_subsampling;

  int pixel_max;
  int pixel_min;

  inline void set_geometry(int width_pixels, int height_pixels) noexcept {
    geometry = {width_pixels, height_pixels};
  }

  inline int plane_width() const noexcept {
    return geometry.width_pixels;
  }

  inline int plane_height() const noexcept {
    return geometry.height_pixels;
  }

  inline void set_frame_planes(
    const unsigned char* src_plane,
    int src_pitch,
    unsigned char* dst_plane,
    int dst_pitch
  ) noexcept {
    buffers.src_bytes = {src_plane, get_src_width(), get_src_height(), src_pitch};
    buffers.dst_bytes = {dst_plane, get_dst_width(), get_dst_height(), dst_pitch};
  }

  inline void set_dither_info_plane(std::span<pixel_dither_info> info, int stride) noexcept {
    buffers.dither_info = {info.data(), plane_width(), plane_height(), stride};
  }

  inline void set_grain_plane(std::int16_t* grain, int stride) noexcept {
    buffers.grain = {grain, plane_width(), plane_height(), stride};
  }

  inline int get_dst_width() const {
    return output_mode == HIGH_BIT_DEPTH_INTERLEAVED ? plane_width() * 2 : plane_width();
  }

  inline int get_dst_height() const {
    return plane_height();
  }

  inline int get_src_width() const {
    return input_mode == HIGH_BIT_DEPTH_INTERLEAVED ? plane_width() * 2 : plane_width();
  }

  inline int get_src_height() const {
    return plane_height();
  }

  inline int copy_line_size_bytes() const noexcept {
    return get_src_width();
  }

  inline bool has_contiguous_byte_planes_for_copy() const noexcept {
    const int line_size = copy_line_size_bytes();
    return buffers.src_bytes.stride == line_size && buffers.dst_bytes.stride == line_size;
  }

  template <class Pixel>
  [[nodiscard]] neo_f3kdb::core::StridedPlaneView<const Pixel> src_plane() const noexcept {
    return {
      reinterpret_cast<const Pixel*>(buffers.src_bytes.data),
      plane_width(),
      plane_height(),
      buffers.src_bytes.stride / static_cast<int>(sizeof(Pixel))
    };
  }

  template <class Pixel>
  [[nodiscard]] neo_f3kdb::core::StridedPlaneView<Pixel> dst_plane() const noexcept {
    return {
      reinterpret_cast<Pixel*>(buffers.dst_bytes.data),
      plane_width(),
      plane_height(),
      buffers.dst_bytes.stride / static_cast<int>(sizeof(Pixel))
    };
  }

  [[nodiscard]] neo_f3kdb::core::StridedPlaneView<const unsigned char> src_bytes() const noexcept {
    return buffers.src_bytes;
  }

  [[nodiscard]] neo_f3kdb::core::StridedPlaneView<unsigned char> dst_bytes() const noexcept {
    return buffers.dst_bytes;
  }

  [[nodiscard]] neo_f3kdb::core::StridedPlaneView<const pixel_dither_info> dither_info_plane() const noexcept {
    return buffers.dither_info;
  }

  [[nodiscard]] neo_f3kdb::core::StridedPlaneView<const std::int16_t> grain_plane() const noexcept {
    return buffers.grain;
  }
} process_plane_params;

namespace neo_f3kdb::core {

struct PlaneJob {
  process_plane_params params{};
  process_plane_context* context = nullptr;
  std::span<pixel_dither_info> dither_info{};
  std::span<std::int16_t> grain{};
};

} // namespace neo_f3kdb::core
