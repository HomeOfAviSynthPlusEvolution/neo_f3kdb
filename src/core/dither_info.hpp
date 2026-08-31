#pragma once

#include "core/plane.hpp"
#include <dualsynth/span2d.hpp>

namespace neo_f3kdb::core {

struct DitherInfoPlane {
  span2d::Span<pixel_dither_info> pixels{};
  int stride = 0;
};

} // namespace neo_f3kdb::core
