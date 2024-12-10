#pragma once

#include "dfinfra/standard_includes.hh"

#include <concepts>
#include <exception>
#include <format>
#include <random>
#include <stdexcept>
#include <utility>

#include "dfinfra/flex_vector.hh"
#include "dfinfra/math.hh"
#include "hasht.hh"
#include "zipf_distribution.hh"

/*
 * This files contains the implementation of the row store relation used in the physical algebra.
 * There are two versions, v1 and v2.
 * Further, there are utility functions to generate key and foreign key relations with two attributes k and a
 * (for use in the runtime experiment(s)).
 * Also, this file provides structs with an eval() function to implement hash functions, join predicates,
 * equality comparisons, and concat functions (to be used by algebra operators).
 *
 * Conventions:
 * Suffixes:
 *   _c  = concept (C++ 20)
 *   _t  = type
 *   _tt = template type
 *
 * Relation names
 *   R(k, a): key relation -> join key: k
 *   S(k, a): foreign key relation -> join key a
 *
 * XXX All relations can only have unsigned integer attributes (uint32_t, uint64_t, enforced by concept).
 *
 * Adapted from Guido Moerkotte.
 */


/* ==========
 *  Concepts
 * ========== */

// at this point, the algebra can only process uint32_t and uint64_t
template <typename T>
concept attr_type_c = requires (T t) {
  requires std::same_as<T, uint32_t> || std::same_as<T, uint64_t>;
};

template <typename T>
concept bun_c = requires(T t) {
  typename T::attr_t;
  t.k;
  t.a;
};


/* =============
 *  Tuple Types
 * ============= */

// tuple (binary unit) from binary association table (BAT), cf. CWI/Hollaender & column store
// [k, a] (k: key, a: foreign key)
template <attr_type_c Tattr>
struct bun_tt {
  using attr_t = Tattr;
  attr_t k;
  attr_t a;
};


using bun32_t = bun_tt<uint32_t>;
using bun64_t = bun_tt<uint64_t>;

template <attr_type_c Tattr>
inline std::ostream& print(std::ostream& os, const bun_tt<Tattr>& t, const size_t aIndent, const size_t aWidth) {
  std::cout << std::setw(aIndent) << "";
  if (aWidth > 0) {
    std::cout << std::setw(aWidth) << t.k << " " << std::setw(aWidth) << t.a;
  } else {
    std::cout << "[" << t.k << "," << t.a << "]";
  }
  return os;
}

template <attr_type_c Tattr>
inline std::ostream& operator<<(std::ostream& os, const bun_tt<Tattr>& t) {
  std::cout << "[" << t.k << "," << t.a << "]";
  return os;
}

// join result tuple: [r*, s*] (pointers to base buns from R and S)
template<bun_c Tbun>
struct join_res_tt {
  using attr_t = typename Tbun::attr_t;
  using bun_t = Tbun;
  const bun_t* _r;
  const bun_t* _s;
};

template<bun_c Tbun>
std::ostream&
operator<<(std::ostream& os, const join_res_tt<Tbun>& x) {
  os << "jr[" << x._r->k << ',' << x._s->k << ']';
  return os;
}


/* =============================
 *  Hash Functions & Predicates
 * ============================= */

// Hash

// R: key relation -> hash key attribute k
template<bun_c Tbun>
struct hash_R_tt {
  using input_t  = Tbun;
  using output_t = typename Tbun::attr_t;
  using attr_t   = typename Tbun::attr_t;
  static inline output_t eval(const input_t* r) {
    return ht::murmur_hash<attr_t>(r->k);
  }
};

// S: foreign key relation -> hash foreign key attribute a
template<bun_c Tbun>
struct hash_S_tt {
  using input_t  = Tbun;
  using output_t = typename Tbun::attr_t;
  using attr_t   = typename Tbun::attr_t;
  static inline output_t eval(const input_t* s) {
    return ht::murmur_hash<attr_t>(s->a);
  }
};

// Predicates

// pred: left.k == right.a <==> R.k == S.a
template<bun_c Tleft, bun_c Tright>
struct join_pred_RS {
  using left_t  = Tleft;
  using right_t = Tright;
  static inline bool eval(const left_t* x, const right_t* y) {
    return (x->k == y->a);
  }
};

// pred: left.a == right.k <==> S.a == R.k
template<bun_c Tleft, bun_c Tright>
struct join_pred_SR {
  using left_t  = Tleft;
  using right_t = Tright;
  static inline bool eval(const left_t* x, const right_t* y) {
    return (x->a == y->k);
  }
};

// pred: compare two R buns for equality (i.e., self-join predicate)
template<bun_c Tbun>
struct eq_join_attr_R {
  using left_t  = Tbun;
  using right_t = Tbun;
  static inline bool eval(const Tbun* x, const Tbun* y) {
    return (x->k == y->k);
  }
};

// pred: compare two S buns for equality (i.e., self-join predicate)
template<bun_c Tbun>
struct eq_join_attr_S {
  using left_t  = Tbun;
  using right_t = Tbun;
  static inline bool eval(const Tbun* x, const Tbun* y) {
    return (x->a == y->a);
  }
};

// Concat functions

template<bun_c Tbun>
struct concat_RS {
  using bun_t = Tbun;
  using left_t  = Tbun;
  using right_t = Tbun;
  using output_t = join_res_tt<bun_t>;
  // binary eval: deprecated
  //static inline void eval(output_t& x, const left_t* r, const right_t* s) {
  //  x._r = r;
  //  x._s = s;
  //}
  static inline void eval_prb(output_t& x, const bun_t* r) {
    x._r = r;
  }
  static inline void eval_bld(output_t& x, const bun_t* s) {
    x._s = s;
  }
};

template<bun_c Tbun>
struct concat_SR {
  using bun_t = Tbun;
  using left_t  = Tbun;
  using right_t = Tbun;
  using output_t = join_res_tt<bun_t>;
  // binary eval: deprecated
  //static inline void eval(output_t& x, const left_t* s, const right_t* r) {
  //  x._r = r;
  //  x._s = s;
  //}
  static inline void eval_prb(output_t& x, const bun_t* s) {
    x._s = s;
  }
  static inline void eval_bld(output_t& x, const bun_t* r) {
    x._r = r;
  }
};


/* =======================
 *  Testing: Top Operator
 * ======================= */

template<typename Tuint, typename Tjoinres>
struct key_pair_tt {
  Tuint _key_r;
  Tuint _key_s;
  key_pair_tt() : _key_r(0), _key_s(0) {}
  key_pair_tt(const Tjoinres& x) : _key_r(x._r->k), _key_s(x._s->k) {}

  inline Tuint hash() const {
    return ht::hash_combine<Tuint>( ht::murmur_hash<Tuint>(_key_r), ht::murmur_hash<Tuint>(_key_s));
  }

  inline bool operator==(const key_pair_tt& x) const {
    return ((_key_r == x._key_r) && (_key_s == x._key_s));
  }
};

template<typename Tuint, typename Tjoinres>
struct std::hash< key_pair_tt<Tuint, Tjoinres> > {
  inline size_t operator()(const key_pair_tt<Tuint, Tjoinres>& x) const { return x.hash(); }
};

template<typename Tuint, typename Tjoinres>
std::ostream&
operator<<(std::ostream& os, const key_pair_tt<Tuint, Tjoinres>& x) {
  os << "kp[" << x._key_r << ',' << x._key_s << ']';
  return os;
}


/* ============
 *  RelationRS
 * ============ */

// row store relation v2
template <typename Ttuple>
struct RelationRS2 {
  public:
    using tuple_t  = Ttuple;
    using tuple_vt = df::infra::FlexVectorNC<tuple_t>;

  public:
    inline RelationRS2() : _tuples() {}

  public:
    inline size_t card() const { return _tuples.size(); }
    inline size_t size() const { return _tuples.size(); }
    inline size_t capacity() const { return _tuples.capacity(); }

    inline void resize(const size_t aNewSize) { _tuples.resize(aNewSize); }
    inline void mem_alloc(const size_t aCapacity) { _tuples.mem_alloc(aCapacity); }
    inline void mem_init() { _tuples.mem_init(); }
    inline void mem_free() { _tuples.mem_free(); }

    inline const tuple_vt& tuples() const { return _tuples; }
    inline       tuple_vt& tuples()       { return _tuples; }

    inline const tuple_t& operator[](const size_t i) const { return _tuples[i]; }
    inline       tuple_t& operator[](const size_t i)       { return _tuples[i]; }

    inline const tuple_t* ptr(const size_t i) const { return &_tuples[i]; }
    inline       tuple_t* ptr(const size_t i)       { return &_tuples[i]; }

  private:
    template <typename Ttuple2>
    friend std::ostream& operator<<(std::ostream& os, const RelationRS2<Ttuple2>& aRel);

  private:
    tuple_vt _tuples;
};

template <typename Ttuple>
std::ostream& operator<<(std::ostream& os, const RelationRS2<Ttuple>& aRel) {
  for (auto const& t : aRel.tuples()) {
    os << t << "\n";
  }
  return os;
}

// row store relation v1
template <typename Ttuple>
struct RelationRS1 {
  using tuple_t  = Ttuple;
  using tuple_vt = std::vector<tuple_t>;

  tuple_vt _tuples;

  inline size_t card() const { return _tuples.size(); }
};

template <typename Ttuple>
std::ostream& operator<<(std::ostream& os, const RelationRS1<Ttuple>& aRel) {
  for (const auto& t : aRel._tuples) {
    os << t << "\n";
  }
  return os;
}


/* =================
 *  Data Generation
 * ================= */

// parameters:
// - aCard: relation cardinality
// - aDomSize: foreign key domain size, generate foreign keys from domain [0, aDomSize - 1]
// - aSkew: 0 = uniform distribution, 1 = standard Zipf distribution (alpha = 1.0)
// - aRng: random number generator engine (e.g., std::mt19937)

// build key relation R
template<bun_c Tbun>
bool
build_rel_key(RelationRS2<Tbun>& R, const size_t aCard) {
  R.resize(aCard);
  for (uint64_t i = 0; i < R.capacity(); ++i) {
    R[i].k = i;
    R[i].a = i;
  }
  return true;
}

// build foreign key relation S
template<bun_c Tbun, typename Trng>
bool 
build_rel_fk(RelationRS2<Tbun>& S, const size_t aCard, const size_t aDomSize, const uint aSkew, Trng& aRng) {
  S.resize(aCard);
  if (0 == aSkew) {
    return build_rel_fk_uni<Tbun>(S, aCard, aDomSize, aRng);
  } else if (1 == aSkew) {
    return build_rel_fk_zipf<Tbun>(S, aCard, aDomSize, aRng);
  } else {
    assert(0 == 1);
    throw std::invalid_argument("Skew must be 0 or 1, but was " + std::to_string(aSkew) + ".");
  }
  return true;
}

// build foreign key relation S: FK S.a according to uniform distribution
template<bun_c Tbun, typename Trng>
bool 
build_rel_fk_uni(RelationRS2<Tbun>& S, const size_t aCard, const size_t aDomSize, Trng& aRng) {
  using attr_t = typename Tbun::attr_t;
  assert(0 < aDomSize);
  if (1 == aDomSize) {
    for (size_t i = 0; i < aCard; ++i) {
      S[i].k = i;
      S[i].a = 0;
    }
  } else {
    std::uniform_int_distribution<attr_t> lDist(0, aDomSize - 1);
    for (size_t i = 0; i < aCard; ++i) {
      S[i].k = i;
      S[i].a = lDist(aRng);
    }
  }
  return true;
}

// build foreign key relation S: FK S.a according to standard Zipf distribution
template<bun_c Tbun, typename Trng>
bool 
build_rel_fk_zipf(RelationRS2<Tbun>& S, const size_t aCard, const size_t aDomSize, Trng& aRng) {
  using attr_t = typename Tbun::attr_t;
  assert(0 < aDomSize);
  if (1 == aDomSize) {
    for (size_t i = 0; i < aCard; ++i) {
      S[i].k = i;
      S[i].a = 0;
    }
  } else {
    zipf_distribution<attr_t> lDist(aDomSize, 1);  // parameters: (n, alpha)
    for (size_t i = 0; i < aCard; ++i) {
      S[i].k = i;
      S[i].a = lDist(aRng) - 1; // NOTE: zipf_distribution generates numbers in [1,n], must subtract one
    }
  }
  return true;
}

/* =================
 *  Print Functions
 * ================= */

template<bun_c Tbun>
std::ostream&
print_RS(std::ostream& os, const RelationRS2<Tbun>& R, const RelationRS2<Tbun>& S) {
  const size_t lIndent = 2;
  const size_t lMaxCard = std::max(R.card(), S.card());
  size_t lColWidth = df::infra::number_of_digits(lMaxCard);
  lColWidth += (~(lColWidth & 0x1) & 0x1);  // make odd so that k and a can be centered across their column
  const size_t lRelColWidth = 2*lColWidth + 3;  // 2 attributes + 3 chars separator

  os << std::format("{3:{4}}| {0:^{2}} | {1:^{2}} |\n", "R", "S", lRelColWidth, "", lIndent);
  os << std::format("{3:{4}}| {0:^{2}}   {1:^{2}} | {0:^{2}}   {1:^{2}} |\n", "k", "a", lColWidth, "", lIndent);
  os << std::format("{3:{4}}| {0:{2}}   {0:{2}} | {0:{2}}   {0:{2}} |\n", std::string(lColWidth, '-'), "", lColWidth, "", lIndent);

  const uint n = std::max(R.card(), S.card());
  for(uint i = 0; i < n; ++i) {
    if(i < R.card()) {
      os << std::format("{3:{4}}| {0:>{2}}   {1:>{2}} | ", R[i].k, R[i].a, lColWidth, "", lIndent);
    } else {
      os << std::format("{3:{4}}| {0:>{2}}   {1:>{2}} | ", "", "", lColWidth, "", lIndent);
    }
    if(i < S.card()) {
      os << std::format("{0:>{2}}   {1:>{2}} |", S[i].k, S[i].a, lColWidth, "", lIndent);
    } else {
      os << std::format("{0:>{2}}   {1:>{2}} |", "", "", lColWidth, "", lIndent);
    }
    os << std::endl;
  }
  return os;
}

template<bun_c Tbun, typename Thashfun>
std::ostream&
print_RS_wh(std::ostream& os, const RelationRS2<Tbun>& R, const RelationRS2<Tbun>& S,
            const uint aDirSizeR, const uint aDirSizeS) {
  os << "     R        |      S        |" << std::endl;
  os << "  k     a   h |   k     a   h |" << std::endl;
  const uint n = std::max(R.card(), S.card());
  for(uint i = 0; i < n; ++i) {
    if(i < R.card()) {
      const uint lHashR = Thashfun()(R[i].k) % aDirSizeR;
      os << std::setw(3) << R[i] << ' ' << std::setw(3) << lHashR << " | ";
    } else {
      os << std::string(13, ' ') << " | ";
    }
    if(i < S.card()) {
      const uint lHashS = Thashfun()(S[i].a) % aDirSizeS;
      os << std::setw(3) << S[i] << ' ' << std::setw(3) << lHashS << " | ";
    } else {
      os << std::string(13, ' ') << " | ";
    }
    os << std::endl;
  }
  return os;
}
