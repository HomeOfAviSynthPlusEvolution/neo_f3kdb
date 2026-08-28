#pragma once

#include "core.h"

namespace ns_avx2 {
    template<int sample_mode, bool blur_first, int dither_algo>
    void process_plane_impl(const process_plane_params& params, process_plane_context* context);
}

namespace ns_avx512 {
    template<int sample_mode, bool blur_first, int dither_algo>
    void process_plane_impl(const process_plane_params& params, process_plane_context* context);
}

const process_plane_impl_t* select_avx_impl(int opt, int cpu_flags, int dither_algo);
