#pragma once

#include "standard_includes.hh"

namespace df::infra {

// print vector as: 0 1 2 3
template <typename T>
inline void
vec_print(const std::vector<T> v, std::ostream& os = std::cout) {
  for (const auto& elem : v) {
    os << elem << " ";
  }
  os << "\n";
}

// output operator for vector, e.g., [0, 1, 2, 3]
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
  os << "[";
  if (!v.empty()) {
    os << v.at(0);
    for (size_t i = 1; i < v.size(); ++i) {
      os << ", " << v.at(i);
    }
  }
  os << "]";
  return os;
}

// fill vector with 'step' many values: start, start + stepsize, start + 2*stepsize...
template <typename Tint>
inline void
vec_fill_step(std::vector<Tint>& v, const Tint start, const Tint steps, const Tint stepsize = 1) {
  assert(steps > 0);
  for (Tint i = 0; i < steps; ++i) {
    v.push_back(start + (i * stepsize));
  }
}

// duplicate each element of v n-times.
template <typename T>
inline void
vec_inflate(std::vector<T>& v, const size_t n) {
  if (0 == v.size()) { return; }
  std::vector<T> u;
  u.reserve(v.size() * n);
  for (const auto& elem : v) {
    for (size_t i = 0; i < n; ++i) {
      u.push_back(elem);
    }
  }
  std::swap(v, u);
  assert((u.size() * n) == v.size());
}


}
