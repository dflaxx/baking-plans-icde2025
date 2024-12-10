#ifndef HT_HASHJOIN_NODE_CONTENT_HH
#define HT_HASHJOIN_NODE_CONTENT_HH
#pragma once

#include "dfinfra/standard_includes.hh"
#include "pointer_array.hh"

/*
 * Author: Guido Moerkotte
 */

/* ==========
 *  Concepts
 * ========== */

// concept for content of hashtable that is common to both chained and nested hashtable
template <typename T>
concept hj_common_content_c = requires(T t, const T tc) {
  // typedefs
  typename T::tuple_t;
  // getters
  { tc.is_empty() } -> std::same_as<bool>;
  { tc.is_full() } -> std::same_as<bool>;
  { tc.begin() } -> std::same_as<uint>;
  { tc.end() } -> std::same_as<uint>;
  // mutators
  { t.clear() } -> std::same_as<void>;
  { t.do_insert(static_cast<const T::tuple_t*>(nullptr)) } -> std::same_as<bool>;
};

// concept for content of chained hashtable
template <typename T>
concept hj_ch_content_c =
  hj_common_content_c<T>
  && requires(T t, const T tc) {
    typename T::hashval_t;
    { tc.tuple(static_cast<uint>(7)) } -> std::same_as<const typename T::tuple_t*>;
    { tc.hashval(static_cast<uint>(7)) } -> std::same_as<typename T::hashval_t>;
    { tc.is_full() } -> std::same_as<bool>;
    { tc.is_empty() } -> std::same_as<bool>;
    { t.do_insert(static_cast<const T::tuple_t*>(nullptr), static_cast<T::hashval_t>(0)) }
      -> std::same_as<bool>;
  };

// concept for main content of nested hashtable
template <typename T>
concept hj_3d_content_main_c =
  hj_common_content_c<T>
  && requires(T t, const T tc) {
    typename T::hashval_t;
    { tc.tuple(static_cast<uint>(7)) } -> std::same_as<const typename T::tuple_t*>;
    { tc.hashval() } -> std::same_as<typename T::hashval_t>;
  };

// concept for sub content of nested hashtable
template <typename T>
concept hj_3d_content_sub_c =
  hj_common_content_c<T>
  && requires(T t, const T tc) {
    true; // dummy, cannot have empty concept
  };


/* =================
 *  Content Structs
 * ================= */

/*
 *  structs for contents for nodes in hashtables
 *  next pointer in hashtable_xxx::entry_t
 *  sub pointer in hashtable_nested::node_sub_t
 */

template<typename Ttuple, typename Thashval, uint n>
struct ht_content_ch_tt {
  using tuple_t   = Ttuple;
  using array_t   = pointer_array_tt<const Ttuple*, n>;
  using hashval_t = Thashval; 

  // we use the first element of _tuple to hold an idx to the next free element in _tuple
  // (as long as the node is not full)
  static constexpr uint cnt_idx = array_t::cnt_idx;
  static constexpr uint capacity = array_t::capacity;

  inline __attribute__((always_inline))       uint      begin() const { return _tuple.begin(); }
  inline __attribute__((always_inline))       uint      end() const { return _tuple.end(); }
  inline __attribute__((always_inline))       uint      size() const { return _tuple.size(); }
  inline __attribute__((always_inline)) const tuple_t*  tuple(const uint i) const { return _tuple[i]; }
  inline __attribute__((always_inline))       hashval_t hashval(const uint i) const { return _hashval[i]; }

  inline __attribute__((always_inline)) bool is_empty() const { return _tuple.is_empty(); }
  inline __attribute__((always_inline)) bool is_full()  const { return _tuple.is_full(); }

  // inserts into this, returns true if not full, returns false if full
  inline __attribute__((always_inline)) bool do_insert(const tuple_t* aTuple) { return _tuple.insert(aTuple); }
  inline __attribute__((always_inline)) bool do_insert(const tuple_t* aTuple, const hashval_t aHashVal) { 
    const bool lRes = do_insert(aTuple); 
    if(lRes) {
      _hashval[begin()] = aHashVal;
    }
    return lRes;
  }

  inline __attribute__((always_inline)) void clear() { _tuple.clear(); }

  array_t   _tuple;
  hashval_t _hashval[n];
};

/*
 *  struct ht_node_3d_[main|sub]_tt
 */

template<typename Thashval, typename Tsub, uint n>
struct ht_content_3d_main_tt {
  using tuple_t   = typename Tsub::tuple_t;
  using hashval_t = Thashval;
  using content_main_t = ht_content_3d_main_tt;
  using array_t   = pointer_array_tt<const tuple_t*, n>;

  hashval_t   _hashval;
  array_t     _tuple;

  // we use the first element of _tuple to hold an idx to the next free element in _tuple
  // (as long as the node is not full)
  static constexpr uint cnt_idx = array_t::cnt_idx;
  static constexpr uint capacity = array_t::capacity;

  inline __attribute__((always_inline))       uint        max_no_tuple() const { return n; }
  inline __attribute__((always_inline))       hashval_t   hashval() const { return _hashval; }
  inline __attribute__((always_inline))       void        hashval(const hashval_t aHashval) { _hashval = aHashval; }

  inline __attribute__((always_inline))       uint     begin() const { return _tuple.begin(); }
  inline __attribute__((always_inline))       uint     end() const { return _tuple.end(); }
  inline __attribute__((always_inline))       uint     size() const { return _tuple.size(); }
  inline __attribute__((always_inline)) const tuple_t* tuple(const uint i) const { return _tuple[i]; }
  inline __attribute__((always_inline)) const tuple_t* get_first_tuple() const { return _tuple.get_first_ptr(); }

  inline __attribute__((always_inline)) bool is_empty() const { return _tuple.is_empty(); }
  inline __attribute__((always_inline)) bool is_full()  const { return _tuple.is_full(); }

  // inserts into this, returns true if not full, returns false if full
  inline __attribute__((always_inline)) bool do_insert(const tuple_t* aTuple) { return _tuple.insert(aTuple); }

  inline __attribute__((always_inline)) void clear() {
    _hashval = 0; 
    _tuple.clear();
  }
};


template<typename Ttuple, uint n>
struct ht_content_3d_sub_tt {
  using tuple_t = Ttuple;
  using array_t = pointer_array_tt<const Ttuple*, n>;

  array_t _tuple;

  // we use the first element of _tuple to hold an idx to the next free element in _tuple
  // (as long as the node is not full)
  static constexpr uint cnt_idx = array_t::cnt_idx;
  static constexpr uint capacity = array_t::capacity;

  inline __attribute__((always_inline)) bool  is_empty() const { return _tuple.is_empty(); }
  inline __attribute__((always_inline)) bool  is_full() const { return _tuple.is_full(); }
  inline __attribute__((always_inline)) bool  do_insert(const tuple_t* aTuple) { return _tuple.insert(aTuple); }
  inline __attribute__((always_inline))       uint     begin() const { return _tuple.begin(); }
  inline __attribute__((always_inline))       uint     end() const { return _tuple.end(); }
  inline __attribute__((always_inline))       uint     size() const { return _tuple.size(); }
  inline __attribute__((always_inline)) const tuple_t* tuple(const uint i) const { return _tuple[i]; }
  inline __attribute__((always_inline)) const tuple_t* get_first_tuple() const { return _tuple.get_first_ptr(); }
  inline __attribute__((always_inline)) void  clear() { _tuple.clear(); }
};


/* ================
 *  Static Asserts
 * ================ */

static_assert(hj_common_content_c<ht_content_ch_tt<void, uint64_t, 1>>);

static_assert(hj_ch_content_c<ht_content_ch_tt<void, uint64_t, 1>>);
static_assert(hj_3d_content_sub_c<ht_content_3d_sub_tt<void, 1>>);
static_assert(hj_3d_content_main_c<ht_content_3d_main_tt<uint64_t, ht_content_3d_sub_tt<void, 1>, 1>>);

static_assert(hj_ch_content_c<ht_content_ch_tt<void, uint64_t, 3>>);
static_assert(hj_3d_content_sub_c<ht_content_3d_sub_tt<void, 3>>);
static_assert(hj_3d_content_main_c<ht_content_3d_main_tt<uint64_t, ht_content_3d_sub_tt<void, 3>, 3>>);

#endif
