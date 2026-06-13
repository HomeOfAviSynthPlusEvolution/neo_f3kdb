#include "core/processor.hpp"

#include "core/constants.hpp"
#include "core/kernel.hpp"

#include <cstdlib>
#include <cstring>

namespace neo_f3kdb::core {

DebandProcessor::DebandProcessor(const DebandProcessorConfig& config)
  : params_(config.params),
    input_(config.input),
    output_(make_output_info(input_, params_)),
    frame_state_(input_, params_),
    requested_backend_(config.requested_backend) {}

ds::VideoOutputInfo DebandProcessor::output_info() const {
  return output_;
}

void DebandProcessor::process(
  ds::MutableVideoFrameView& dst,
  const ds::VideoFrameView& src,
  int frame_number
) {
  for (int p = 0; p < src.plane_count; ++p) {
    const auto& src_plane = src.plane(p);
    auto& dst_plane = dst.plane(p);

    auto* dst_ptr = reinterpret_cast<unsigned char*>(dst_plane.data);
    const auto* src_ptr = reinterpret_cast<const unsigned char*>(src_plane.data);

    auto job = make_plane_job(
      frame_number,
      p,
      dst_ptr,
      static_cast<int>(dst_plane.stride_bytes),
      src_ptr,
      static_cast<int>(src_plane.stride_bytes)
    );

    const int grain_setting = p == 0 ? params_.grain_y : params_.grain_c;
    if (should_copy_plane(job, grain_setting)) {
      copy_plane(job);
      continue;
    }

    if (requested_backend_ == Backend::Highway && supports_highway(job.params)) {
      process_plane_highway(job);
    } else {
      process_plane_scalar(job);
    }
  }
}

PlaneJob DebandProcessor::make_plane_job(
  int frame_number,
  int plane,
  unsigned char* dst_frame_ptr,
  int dst_pitch,
  const unsigned char* src_frame_ptr,
  int src_pitch
) {
  PlaneJob job{};
  auto& p = job.params;

  p.src_plane_ptr = src_frame_ptr;
  p.src_pitch = src_pitch;
  p.dst_plane_ptr = dst_frame_ptr;
  p.dst_pitch = dst_pitch;

  const int input_depth = ds::bits_per_sample(input_.format.sample_format);
  p.input_mode = input_depth == 8 ? LOW_BIT_DEPTH : HIGH_BIT_DEPTH_INTERLEAVED;
  p.input_depth = input_depth;
  p.output_mode = params_.output_depth <= 8 ? LOW_BIT_DEPTH : HIGH_BIT_DEPTH_INTERLEAVED;
  p.output_depth = params_.output_depth;
  p.dither_algo = params_.dither_algo;
  p.blur_first = params_.blur_first;
  p.sample_mode = params_.sample_mode;
  p.angle_boost = static_cast<float>(params_.angle_boost);
  p.max_angle = static_cast<float>(params_.max_angle);

  p.plane = plane;
  p.width_subsampling = static_cast<unsigned char>(plane == 0 ? 0 : input_.format.subsampling_w);
  p.height_subsampling = static_cast<unsigned char>(plane == 0 ? 0 : input_.format.subsampling_h);
  p.plane_width_in_pixels = plane == 0 ? input_.width : (input_.width >> input_.format.subsampling_w);
  p.plane_height_in_pixels = plane == 0 ? input_.height : (input_.height >> input_.format.subsampling_h);

  p.info_stride = frame_state_.info_stride(plane);
  p.grain_buffer_stride = frame_state_.grain_stride(plane);

  job.dither_info = frame_state_.dither_info(plane);
  job.grain = frame_state_.grain_buffer(plane);
  p.info_ptr_base = job.dither_info.data();
  p.grain_buffer = frame_state_.grain_row_base(plane, frame_number, input_.num_frames);
  job.context = frame_state_.context(plane);

  switch (plane) {
  case 0:
    p.threshold = static_cast<std::uint16_t>(params_.y);
    p.threshold1 = static_cast<std::uint16_t>(params_.y_1);
    p.threshold2 = static_cast<std::uint16_t>(params_.y_2);
    p.pixel_max = params_.keep_tv_range ? TV_RANGE_Y_MAX : FULL_RANGE_Y_MAX;
    p.pixel_min = params_.keep_tv_range ? TV_RANGE_Y_MIN : FULL_RANGE_Y_MIN;
    break;
  case 1:
    p.threshold = static_cast<std::uint16_t>(params_.cb);
    p.threshold1 = static_cast<std::uint16_t>(params_.cb_1);
    p.threshold2 = static_cast<std::uint16_t>(params_.cb_2);
    p.pixel_max = params_.keep_tv_range ? TV_RANGE_C_MAX : FULL_RANGE_C_MAX;
    p.pixel_min = params_.keep_tv_range ? TV_RANGE_C_MIN : FULL_RANGE_C_MIN;
    break;
  case 2:
    p.threshold = static_cast<std::uint16_t>(params_.cr);
    p.threshold1 = static_cast<std::uint16_t>(params_.cr_1);
    p.threshold2 = static_cast<std::uint16_t>(params_.cr_2);
    p.pixel_max = params_.keep_tv_range ? TV_RANGE_C_MAX : FULL_RANGE_C_MAX;
    p.pixel_min = params_.keep_tv_range ? TV_RANGE_C_MIN : FULL_RANGE_C_MIN;
    break;
  default:
    std::abort();
  }

  return job;
}

bool DebandProcessor::should_copy_plane(const PlaneJob& job, int grain_setting) const {
  return ds::bits_per_sample(input_.format.sample_format) == params_.output_depth &&
    grain_setting == 0 &&
    job.params.threshold == 0 &&
    job.params.threshold1 == 0 &&
    job.params.threshold2 == 0;
}

void DebandProcessor::copy_plane(const PlaneJob& job) const {
  const auto& p = job.params;
  const int line_size = p.get_src_width();
  auto src = p.src_plane_ptr;
  auto dst = p.dst_plane_ptr;
  if (line_size == p.src_pitch && p.src_pitch == p.dst_pitch) {
    std::memcpy(dst, src, static_cast<std::size_t>(line_size * p.get_src_height()));
    return;
  }

  for (int row = 0; row < p.get_src_height(); row++) {
    std::memcpy(dst, src, static_cast<std::size_t>(line_size));
    src += p.src_pitch;
    dst += p.dst_pitch;
  }
}

bool supports_highway(const process_plane_params& params) {
  return params.dither_algo != DA_HIGH_FLOYD_STEINBERG_DITHERING &&
    params.sample_mode >= 1 &&
    params.sample_mode <= 5;
}

} // namespace neo_f3kdb::core
