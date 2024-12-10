#pragma once

#include <stdexcept>
#include <source_location>

namespace df::infra {

class NotImplementedException : public std::logic_error {
  public:
    NotImplementedException(const std::source_location& s = std::source_location::current())
      : std::logic_error(
          "In file " + std::string(s.file_name() + std::string(":") + std::to_string(s.line())) + std::string(":\n")
          + std::string(_what_indent, ' ') + std::string(s.function_name()) + ":\n"
          + std::string(_what_indent, ' ') + "This function has not been implemented yet. "
        ) {}
  private:
    inline static size_t _what_indent = 2 + 4 + 2 + 1 + 2; // '  what():  '
};

}
