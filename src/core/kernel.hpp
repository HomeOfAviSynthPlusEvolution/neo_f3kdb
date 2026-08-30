#pragma once

#include "core/parameters.hpp"
#include "core/plane.hpp"

using deband_kernel_proc_t = void (__cdecl *)(
  const process_plane_params& params,
  process_plane_context* context
);

#define DITHER_CONTEXT_BUFFER_SIZE 8192
#define CONTEXT_BUFFER_SIZE DITHER_CONTEXT_BUFFER_SIZE

namespace neo_f3kdb::core {

void process_plane_scalar(KernelExecution execution);
void process_plane_highway(KernelExecution execution);

} // namespace neo_f3kdb::core
