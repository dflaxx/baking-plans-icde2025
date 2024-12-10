#pragma once

#include <iostream>
#include <string_view>	
#include <source_location>

// Macros for printing variables, variable names, types (poor man's reflection)
// Source: Magnus

#define COUT_LOG std::cout << __FILE__ << ":" << __LINE__ << ": "

#define STRINGIFY(X) #X
#define EXPAND_AND_STRINGIFY(X) STRINGIFY(X)
#define PRINT_REFLECTION(X) (#X) << ": " << (X)
#define TRACE_AND_EXECUTE(X) COUT_LOG << EXPAND_AND_STRINGIFY(X) << std::endl; X;

// Raise breaokpoint in debugger.
// Source: Roland Leissa, https://github.com/leissa/fe/blob/acade5d220f96f68db07933bb6b327955305b45b/include/fe/assert.h#L22
#if (defined(__clang__) || defined(__GNUC__)) && (defined(__x86_64__) || defined(__i386__))
inline void breakpoint() { asm("int3"); }
#else
inline void breakpoint() {
  volatile int* p = nullptr;
  *p              = 42;
}
#endif

// Get typename as string in human-readable form
// usage:  int x;
//         std::cout << type_name<decltype(x)>() << std::endl;  // prints 'int'
// source: https://stackoverflow.com/questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c/
template <typename T>
constexpr auto type_name() noexcept {
  std::string_view name = "Error: unsupported compiler", prefix, suffix;
#ifdef __clang__
  name = __PRETTY_FUNCTION__;
  prefix = "auto type_name() [T = ";
  suffix = "]";
#elif defined(__GNUC__)
  name = __PRETTY_FUNCTION__;
  prefix = "constexpr auto type_name() [with T = ";
  suffix = "]";
#elif defined(_MSC_VER)
  name = __FUNCSIG__;
  prefix = "auto __cdecl type_name<";
  suffix = ">(void) noexcept";
#endif
  name.remove_prefix(prefix.size());
  name.remove_suffix(suffix.size());
  return name;
}

enum class src_loc_output_et {
  file_line,        // debugging_helpers.hh:43
  function_fq,      // (fully qualified) C::f(T x) [with T = long int]
  full              // debugging_helpers.hh:43 : C::f(T x) [...]
};

inline std::string
src_loc_str(src_loc_output_et aKind,
    const std::source_location loc = std::source_location::current()) {
  std::string str;
  switch (aKind) {
    case src_loc_output_et::file_line:
      str += loc.file_name();
      str += ":";
      str += loc.line();
      break;
    case src_loc_output_et::function_fq:
      str += loc.function_name();
      break;
    case src_loc_output_et::full:
      str += loc.file_name();
      str += ":";
      str += loc.line();
      str += " : ";
      str += loc.function_name();
      break;
  }
  return str;
}

inline std::string
src_file_str(const std::source_location loc = std::source_location::current()) {
  return src_loc_str(src_loc_output_et::file_line, loc);
}

inline std::string
src_func_str(const std::source_location loc = std::source_location::current()) {
  return src_loc_str(src_loc_output_et::function_fq, loc);
}

inline void
msg(const std::string_view message,
    std::ostream& os = std::clog,
    const std::source_location location = std::source_location::current()) {
  os << "file: "
    << location.file_name() << '('
    << location.line() << ':'
    << location.column() << ") `"
    << location.function_name() << "`: "
    << message << '\n';
}
