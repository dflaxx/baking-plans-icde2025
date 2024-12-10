#pragma once

#include "algop_v2_base.hh"
#include "RelationRS.hh"

#include <functional>
#include <unordered_set>
#include <type_traits>

template <typename Tinput, bool Tmat, bool Tprint>
class AlgOpTop2 : public AlgOpBase2 {
  public:
    using input_t = Tinput;
    using output_t = void;  // top operator is the root and has no consumer
    using consumer_t = void;
    using print_fun_t = std::function<void (const input_t*, std::ostream&)>;
  public:
    using join_res_t = Tinput;
    using attr_t = typename join_res_t::attr_t;
    using key_pair_t = key_pair_tt<attr_t, join_res_t>;
    using key_pair_st = typename std::unordered_set<key_pair_t>;
    using key_pair_ost = typename std::conditional<Tmat, key_pair_st, int>::type;
  public:
    AlgOpTop2() : AlgOpTop2(std::cout) {}
    AlgOpTop2(std::ostream& aOs)
      : AlgOpBase2(algop_et::kAlgOpTop), _os(aOs) {}
    AlgOpTop2(std::ostream& aOs, print_fun_t aPrintFunction)
      : AlgOpBase2(algop_et::kAlgOpTop), _os(aOs), _print_fun(aPrintFunction), _res() {}
  public:
    inline bool mem_alloc() { return true; }
    inline bool mem_init() { return true; }
    inline void init() { reset(); assert(0 == count()); }
    inline void fin() {}
    inline void clear() {}
    inline void mem_free() {}
           void step(const input_t*);
    inline const consumer_t* consumer() const { return nullptr; }
    inline const AlgOpBase2* consumer_poly() const { return nullptr; }
  public:
    constexpr inline bool has_mat() const { return Tmat; }
    inline const key_pair_ost& result() const { return _res; }
  private:
    std::ostream& _os;
    print_fun_t   _print_fun = [] (const input_t* aInput, std::ostream& aOs) { aOs << aInput; };
    key_pair_ost  _res;
};

template <typename Tinput, bool Tmat, bool Tprint>
void
AlgOpTop2<Tinput, Tmat, Tprint>::step(const input_t* aTuple) {
  inc();
  if constexpr (Tprint) {
    // only print result the first time the operator is run
    if (1 == runs()) {
      _print_fun(aTuple, _os);
      _os << std::endl;
    }
  }
  if constexpr (Tmat) {
    // must make sure that tuples are buns of type [T k, T a]
    // with T = uint32_t || uint64_t
    _res.insert(key_pair_t(*aTuple));
  }
}
