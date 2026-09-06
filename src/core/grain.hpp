#pragma once

#include <cstdint>
#include <dualsynth/span2d.hpp>

namespace neo_f3kdb::core {

struct GrainPlane {
  span2d::Span<std::int16_t> samples{};
  int stride = 0;
};

} // namespace neo_f3kdb::core
