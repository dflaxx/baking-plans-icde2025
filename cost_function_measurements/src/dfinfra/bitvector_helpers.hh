#pragma once

#include "standard_includes.hh"

namespace df::infra {

template <typename Tuint>
inline std::string
bv_to_string(const Tuint aBv, const uint aNoBits) {
  std::string lRes(aNoBits, '0');
  for (uint i = 0; i < aNoBits; ++i) {
    if (0x1 == (0x1 & (aBv >> i))) {  // test for '1' bit at curr pos
      lRes[aNoBits - 1 - i] = '1';    // put a '1' in the string
    }
  }
  return lRes;
}

} // end namespace df::infra
