#pragma once

#include <cstdint>
#include <cassert>

inline uint16_t
set_msb(const uint16_t x, const std::size_t s) {
  assert(16 >= s);
  uint16_t all_1 = ((1 << 15) - 1) | (1 << 15);
  if (16 == s) {
    return all_1;
  }
  uint16_t lower_1 = (1 << (16 - s)) - 1;
  uint16_t mask = all_1 ^ lower_1;
  return x | mask;
}
