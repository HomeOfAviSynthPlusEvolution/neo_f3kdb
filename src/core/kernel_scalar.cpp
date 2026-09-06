#include "core/kernel.hpp"

#define NEO_F3KDB_KERNEL_SCALAR_IMPORT_DECLARATION
#include "core/kernel_scalar_dispatch.hpp"

#include <array>
#include <cassert>
#include <cstdlib>

namespace neo_f3kdb::core {

namespace {

int select_impl_index(int sample_mode, bool blur_first) {
  assert(sample_mode != 0);
  return sample_mode * 2 + (blur_first ? 0 : 1) - 1;
}

const deband_kernel_proc_t* scalar_table(DITHER_ALGORITHM dither_algo) {
  switch (dither_algo) {
  case DA_HIGH_NO_DITHERING:
    return process_plane_impl_scalar_high_no_dithering;
  case DA_HIGH_ORDERED_DITHERING:
    return process_plane_impl_scalar_high_ordered_dithering;
  case DA_HIGH_FLOYD_STEINBERG_DITHERING:
    return process_plane_impl_scalar_high_floyd_steinberg_dithering;
  case DA_16BIT_INTERLEAVED:
    return process_plane_impl_scalar_16bit_interleaved;
  default:
    std::abort();
  }
}

} // namespace

void process_plane_scalar(KernelExecution execution) {
  const auto& params = execution.input;
  const auto* table = scalar_table(params.config.dither_algo);
  const deband_kernel_proc_t kernel =
    table[select_impl_index(params.config.sample_mode, params.config.blur_first)];
  kernel(params, execution.context);
}

} // namespace neo_f3kdb::core
