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

    if (requested_backend_ == Backend::Highway) {
      process_plane_highway(job.kernel_execution());
    } else {
      process_plane_scalar(job.kernel_execution());
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

  const int input_depth = ds::bits_per_sample(input_.format.sample_format);
  p.config.input_mode = input_depth == 8 ? LOW_BIT_DEPTH : HIGH_BIT_DEPTH_INTERLEAVED;
  p.config.input_depth = input_depth;
  p.config.output_mode = params_.output_depth <= 8 ? LOW_BIT_DEPTH : HIGH_BIT_DEPTH_INTERLEAVED;
  p.config.output_depth = params_.output_depth;
  p.config.dither_algo = params_.dither_algo;
  p.config.blur_first = params_.blur_first;
  p.config.sample_mode = params_.sample_mode;
  p.config.angle_boost = static_cast<float>(params_.angle_boost);
  p.config.max_angle = static_cast<float>(params_.max_angle);

  p.config.plane_index = plane;
  p.config.width_subsampling = static_cast<unsigned char>(plane == 0 ? 0 : input_.format.subsampling_w);
  p.config.height_subsampling = static_cast<unsigned char>(plane == 0 ? 0 : input_.format.subsampling_h);
  p.set_geometry(
    plane == 0 ? input_.width : (input_.width >> input_.format.subsampling_w),
    plane == 0 ? input_.height : (input_.height >> input_.format.subsampling_h)
  );
  p.set_frame_planes(src_frame_ptr, src_pitch, dst_frame_ptr, dst_pitch);

  job.dither_info = frame_state_.dither_info(plane);
  job.grain = frame_state_.grain_buffer(plane);
  p.set_dither_info_plane(job.dither_info, frame_state_.info_stride(plane));
  p.set_grain_plane(
    frame_state_.grain_row_base(plane, frame_number, input_.num_frames),
    frame_state_.grain_stride(plane)
  );
  job.context = frame_state_.context(plane);

  switch (plane) {
  case 0:
    p.config.threshold = static_cast<std::uint16_t>(params_.y);
    p.config.threshold1 = static_cast<std::uint16_t>(params_.y_1);
    p.config.threshold2 = static_cast<std::uint16_t>(params_.y_2);
    p.config.pixel_max = params_.keep_tv_range ? TV_RANGE_Y_MAX : FULL_RANGE_Y_MAX;
    p.config.pixel_min = params_.keep_tv_range ? TV_RANGE_Y_MIN : FULL_RANGE_Y_MIN;
    break;
  case 1:
    p.config.threshold = static_cast<std::uint16_t>(params_.cb);
    p.config.threshold1 = static_cast<std::uint16_t>(params_.cb_1);
    p.config.threshold2 = static_cast<std::uint16_t>(params_.cb_2);
    p.config.pixel_max = params_.keep_tv_range ? TV_RANGE_C_MAX : FULL_RANGE_C_MAX;
    p.config.pixel_min = params_.keep_tv_range ? TV_RANGE_C_MIN : FULL_RANGE_C_MIN;
    break;
  case 2:
    p.config.threshold = static_cast<std::uint16_t>(params_.cr);
    p.config.threshold1 = static_cast<std::uint16_t>(params_.cr_1);
    p.config.threshold2 = static_cast<std::uint16_t>(params_.cr_2);
    p.config.pixel_max = params_.keep_tv_range ? TV_RANGE_C_MAX : FULL_RANGE_C_MAX;
    p.config.pixel_min = params_.keep_tv_range ? TV_RANGE_C_MIN : FULL_RANGE_C_MIN;
    break;
  default:
    std::abort();
  }

  return job;
}

bool DebandProcessor::should_copy_plane(const PlaneJob& job, int grain_setting) const {
  return ds::bits_per_sample(input_.format.sample_format) == params_.output_depth &&
    grain_setting == 0 &&
    job.params.config.threshold == 0 &&
    job.params.config.threshold1 == 0 &&
    job.params.config.threshold2 == 0;
}

void DebandProcessor::copy_plane(const PlaneJob& job) const {
  const auto& p = job.params;
  const int line_size = p.copy_line_size_bytes();
  const auto src_plane = p.src_bytes();
  auto dst_plane = p.dst_bytes();
  if (p.has_contiguous_byte_planes_for_copy()) {
    std::memcpy(
      dst_plane.data,
      src_plane.data,
      static_cast<std::size_t>(line_size * p.get_src_height())
    );
    return;
  }

  for (int row = 0; row < p.get_src_height(); row++) {
    std::memcpy(
      dst_plane.row(row).data(),
      src_plane.row(row).data(),
      static_cast<std::size_t>(line_size)
    );
  }
}

} // namespace neo_f3kdb::core
