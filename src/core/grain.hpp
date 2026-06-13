#pragma once

#include <cstdint>
#include <span>

namespace neo_f3kdb::core {

struct GrainPlane {
  std::span<std::int16_t> samples{};
  int stride = 0;
};

} // namespace neo_f3kdb::core
