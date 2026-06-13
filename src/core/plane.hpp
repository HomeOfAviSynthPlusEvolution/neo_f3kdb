#pragma once

#include "compiler_compat.h"
#include "core/kernel_types.hpp"

#include <cstdint>
#include <span>

typedef struct _pixel_dither_info {
  alignas(4) std::int8_t ref1;
  std::int8_t ref2;
  std::int16_t change;
} pixel_dither_info;

static_assert(sizeof(pixel_dither_info) == 4, "Something wrong in pixel_dither_info");

typedef void (*destroy_data_t)(void* data);

typedef struct _process_plane_context {
  void* data;
  destroy_data_t destroy;
} process_plane_context;

void destroy_context(process_plane_context* context);
void init_context(process_plane_context* context);

typedef struct _process_plane_params {
  const unsigned char* src_plane_ptr;
  int src_pitch;

  unsigned char* dst_plane_ptr;
  int dst_pitch;

  int plane_width_in_pixels;
  int plane_height_in_pixels;

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
  pixel_dither_info* info_ptr_base;
  int info_stride;

  std::int16_t* grain_buffer;
  int grain_buffer_stride;

  int plane;

  unsigned char width_subsampling;
  unsigned char height_subsampling;

  int pixel_max;
  int pixel_min;

  inline int get_dst_width() const {
    return output_mode == HIGH_BIT_DEPTH_INTERLEAVED ? plane_width_in_pixels * 2 : plane_width_in_pixels;
  }

  inline int get_dst_height() const {
    return plane_height_in_pixels;
  }

  inline int get_src_width() const {
    return input_mode == HIGH_BIT_DEPTH_INTERLEAVED ? plane_width_in_pixels * 2 : plane_width_in_pixels;
  }

  inline int get_src_height() const {
    return plane_height_in_pixels;
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
