#include "core/frame_state.hpp"

#include "core/random.hpp"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <initializer_list>

namespace neo_f3kdb::core {

namespace {

int min_multi(std::initializer_list<int> values) {
  auto it = values.begin();
  int ret = *it;
  for (; it != values.end(); ++it) {
    if (*it >= 0 && *it < ret) {
      ret = *it;
    }
  }
  return ret;
}

AlignedBuffer<std::int16_t, 128> generate_grain_buffer(
  std::size_t item_count,
  RANDOM_ALGORITHM algo,
  int& seed,
  double param,
  int range
) {
  AlignedBuffer<std::int16_t, 128> buffer(item_count);
  for (std::size_t i = 0; i < item_count; i++) {
    buffer[i] = random(algo, seed, range, param);
  }
  return buffer;
}

} // namespace

int frame_lut_stride(int width_in_pixels) {
  return (((width_in_pixels - 1) | (FRAME_LUT_ALIGNMENT - 1)) + 1);
}

FrameState::FrameState(const ds::VideoInputInfo& input, const DebandParameters& params) {
  reset(input, params);
}

FrameState::~FrameState() {
  destroy_contexts();
}

void FrameState::reset(const ds::VideoInputInfo& input, const DebandParameters& params) {
  destroy_contexts();
  for (auto& context : contexts_) {
    init_context(&context);
  }

  grain_buffer_offsets_.clear();

  int seed = 0x92D68CA2 - params.seed;
  seed ^= (input.width << 16) ^ input.height;
  seed ^= (input.num_frames << 16) ^ input.num_frames;

  const int height_in_pixels = input.height;
  const int width_in_pixels = input.width;
  const int width_subsampling = input.format.subsampling_w;
  const int height_subsampling = input.format.subsampling_h;

  y_info_stride_ = frame_lut_stride(width_in_pixels);
  const int y_size = y_info_stride_ * height_in_pixels;
  y_info_.reset(y_size);

  c_info_stride_ = frame_lut_stride(width_in_pixels >> width_subsampling);
  const int c_size = c_info_stride_ * (height_in_pixels >> height_subsampling);
  cb_info_.reset(c_size);
  cr_info_.reset(c_size);

  for (int y = 0; y < height_in_pixels; y++) {
    pixel_dither_info* y_info_ptr = y_info_.data() + y * y_info_stride_;
    pixel_dither_info* cb_info_ptr = cb_info_.data() + (y >> height_subsampling) * c_info_stride_;
    pixel_dither_info* cr_info_ptr = cr_info_.data() + (y >> height_subsampling) * c_info_stride_;

    for (int x = 0; x < width_in_pixels; x++) {
      pixel_dither_info info_y = {0, 0, 0};
      info_y.change = random(
        params.random_algo_grain,
        seed,
        params.grain_y,
        params.random_param_grain
      );

      const int x_range = min_multi({params.range, x, width_in_pixels - x - 1});
      const int y_range = min_multi({params.range, y, height_in_pixels - y - 1});
      const int cur_range = [&]() {
        switch (params.sample_mode) {
        case 1:
          return y_range;
        case 3:
          return x_range;
        case 2:
        case 4:
        case 5:
        case 6:
        case 7:
          return min_multi({x_range, y_range});
        default:
          return 0;
        }
      }();

      if (cur_range > 0) {
        info_y.ref1 = static_cast<signed char>(
          random(params.random_algo_ref, seed, cur_range, params.random_param_ref)
        );
        if (params.sample_mode == 2) {
          info_y.ref2 = static_cast<signed char>(
            random(params.random_algo_ref, seed, cur_range, params.random_param_ref)
          );
        }
        if (params.sample_mode > 0) {
          info_y.ref1 = static_cast<std::int8_t>(std::abs(info_y.ref1));
          info_y.ref2 = static_cast<std::int8_t>(std::abs(info_y.ref2));
        }
      }

      *y_info_ptr = info_y;

      const bool should_set_c =
        ((x & ((1 << width_subsampling) - 1)) == 0 &&
         (y & ((1 << height_subsampling) - 1)) == 0);

      if (should_set_c) {
        pixel_dither_info info_cb = info_y;
        pixel_dither_info info_cr = info_cb;

        info_cb.change = random(
          params.random_algo_grain,
          seed,
          params.grain_c,
          params.random_param_grain
        );
        info_cr.change = random(
          params.random_algo_grain,
          seed,
          params.grain_c,
          params.random_param_grain
        );

        *cb_info_ptr = info_cb;
        *cr_info_ptr = info_cr;
        cb_info_ptr++;
        cr_info_ptr++;
      }
      y_info_ptr++;
    }
  }

  const int multiplier = params.dynamic_grain ? 3 : 1;
  int item_count = width_in_pixels;
  item_count += 255;
  item_count &= 0xffffff80;
  item_count *= height_in_pixels;

  grain_buffer_y_ = generate_grain_buffer(
    static_cast<std::size_t>(item_count * multiplier),
    params.random_algo_grain,
    seed,
    params.random_param_grain,
    params.grain_y
  );
  grain_buffer_c_ = generate_grain_buffer(
    static_cast<std::size_t>(item_count * multiplier),
    params.random_algo_grain,
    seed,
    params.random_param_grain,
    params.grain_c
  );

  y_grain_stride_ = frame_lut_stride(width_in_pixels);
  c_grain_stride_ = frame_lut_stride(width_in_pixels >> width_subsampling);

  if (params.dynamic_grain) {
    grain_buffer_offsets_.resize(input.num_frames);
    for (int i = 0; i < input.num_frames; i++) {
      int offset = item_count + random(
        RANDOM_ALGORITHM_UNIFORM,
        seed,
        item_count,
        DEFAULT_RANDOM_PARAM
      );
      offset &= 0xfffffff0;
      assert(offset >= 0);
      grain_buffer_offsets_[static_cast<std::size_t>(i)] = offset;
    }
  }
}

int FrameState::info_stride(int plane) const noexcept {
  return plane == 0 ? y_info_stride_ : c_info_stride_;
}

int FrameState::grain_stride(int plane) const noexcept {
  return plane == 0 ? y_grain_stride_ : c_grain_stride_;
}

span2d::Span<pixel_dither_info> FrameState::dither_info(int plane) noexcept {
  switch (plane) {
  case 0:
    return y_info_.span();
  case 1:
    return cb_info_.span();
  case 2:
    return cr_info_.span();
  default:
    std::abort();
  }
}

span2d::Span<std::int16_t> FrameState::grain_buffer(int plane) noexcept {
  return plane == 0 ? grain_buffer_y_.span() : grain_buffer_c_.span();
}

std::int16_t* FrameState::grain_row_base(int plane, int frame_number, int frame_count) noexcept {
  std::int16_t* grain = grain_buffer(plane).data();
  if (!grain_buffer_offsets_.empty()) {
    grain += grain_buffer_offsets_[static_cast<std::size_t>(frame_number % frame_count)];
  }
  return grain;
}

process_plane_context* FrameState::context(int plane) noexcept {
  return &contexts_[static_cast<std::size_t>(plane)];
}

void FrameState::destroy_contexts() noexcept {
  for (auto& context : contexts_) {
    destroy_context(&context);
  }
}

} // namespace neo_f3kdb::core
