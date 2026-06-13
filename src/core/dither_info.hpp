#pragma once

#include "core/plane.hpp"

#include <span>

namespace neo_f3kdb::core {

struct DitherInfoPlane {
  std::span<pixel_dither_info> pixels{};
  int stride = 0;
};

} // namespace neo_f3kdb::core
