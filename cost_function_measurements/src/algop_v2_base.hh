#pragma once

#include "dfinfra/standard_includes.hh"
#include "algebra_types_v2.hh"

#include <concepts>


/* ==========
 *  Concepts
 * ========== */

#if __cpp_concepts >= 201907L

class AlgOpBase2;  // fwd decl

template <typename T>
concept alg2_operator_c =
  std::derived_from<T, AlgOpBase2> &&
  requires (T t, const T tc) {
    typename T::input_t;
    typename T::output_t;
    typename T::consumer_t;  // might be void
    //{ t.consumer() } -> std::same_as<typename T::consumer_t*>;
    { tc.consumer() } -> std::same_as<const typename T::consumer_t*>;
    { tc.consumer_poly() } -> std::same_as<const AlgOpBase2*>;  // polymorphism
  };

template <typename T>
concept alg2_consumer_c =
  alg2_operator_c<T> &&
  requires (T t) {
    { t.mem_alloc() } -> std::same_as<bool>;
    { t.mem_init() } -> std::same_as<bool>;
    { t.init() } -> std::same_as<void>;
    { t.step(static_cast<const T::input_t*>(nullptr)) } -> std::same_as<void>;
    { t.fin() } -> std::same_as<void>;
    { t.clear() } -> std::same_as<void>;
    { t.mem_free() } -> std::same_as<void>;
  };

template <typename T>
concept alg2_producer_c =
  alg2_operator_c<T> &&
  requires (T t) {
    { t.run() } -> std::same_as<void>;
  };

template <typename T>
concept alg2_operator_build_c =
  alg2_operator_c<T> &&
  requires (T t) {
    typename T::hashtable_t;
    { t.hashtable() } -> std::same_as<const typename T::hashtable_t&>;
  };

#else
#error Compiling with C++ < 20, concepts not supported.
#endif  // __cpp_concepts >= 201907L


/* =======================
 *  Algebra Base Operator
 * ======================= */

class AlgOpBase2 {
  protected:  // make base class not instantiable
    explicit AlgOpBase2(const algop_et aKind)
      : _count(0), _ok(true), _kind(aKind), _runs(0) {}
  public:
    inline void     reset() { _count = 0; _ok = true; ++_runs; }
    inline void     inc()   { ++_count; }
    inline algop_et kind()  const { return _kind; }
    inline uint64_t runs()  const { return _runs; }
    inline uint64_t count() const { return _count; }
    inline bool     ok()    const { return _ok; }
    inline bool     ok(const bool b) { _ok = b; return _ok; }
  protected:
    uint64_t     _count;     // e.g. number of step() calls for consumers
    bool         _ok;        // e.g. to indicate error states
    algop_et     _kind;      // operator type (enum) like "top operator" (kAlgOpTop)
    uint64_t     _runs;      // number of times the operator was re-run (e.g. for repeated experiments)
};
