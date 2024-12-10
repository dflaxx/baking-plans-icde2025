#pragma once

#include "standard_includes.hh"
#include <cctype>
#include <algorithm>

namespace df::infra {

// trim from start (in place)
inline void
ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
    return !std::isspace(ch);
  }));
}

// trim from end (in place)
inline void
rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
    return !std::isspace(ch);
  }).base(), s.end());
  // reverse_iterator::base() converts a reverse iterator into the corresponding forward one.
  // Note that dereferencing a reverse iterator is tricky, see CPP-Reference and Stackoverflow.
}

// trim from left and right (in place)
inline void
trim(std::string &s) {
  ltrim(s);
  rtrim(s);
}

inline std::string&
to_lower(std::string &s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [] (unsigned char c) { return std::tolower(c); });
  return s;
}

inline std::string&
to_upper(std::string &s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [] (unsigned char c) { return std::toupper(c); });
  return s;
}

} // namespace df::infra
