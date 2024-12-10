#pragma once

/*
 * A templated hashtable that stores Tcontentmain and Tcontentsub objects in a two-level linked-list structure.
 *
 * On the first level of collision chain nodes, there are main nodes (entry_main_t) that each encapsulate a Tcontentmain object.
 * Each of those nodes contains pointers to tuples with the same value on the join attribute(s),
 * so one main node represents on distinct value.
 * Each main node contains
 * - one or more tuple pointers
 * - a hash value (that is the same for all tuples of this main node)
 * - a pointer to the next main node
 * - a pointer to this main node's sub chain
 *
 * If the space for tuple pointers in one main node does not suffice,
 * the other tuples pointers are stored in the sub chain of this main node, in sub nodes (entry_sub_t).
 *
 * The assumptions about the available member of Tcontentmain and Tcontentsub
 * are enforced by concepts hj_3d_content_main_c and hj_3d_content_sub_c.
 *
 */

#include "dfinfra/math.hh"
#include "dfinfra/standard_includes.hh"
#include "dfinfra/fixed_vector.hh"

#include "globals.hh"
#include "hj_htnode_content.hh"
#include "reservoirnc.hh"

#include <format>
#include <limits>


template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
class HashtableNested2 {
  public:
    struct entry_main_t;  // fwd decl
    struct entry_sub_t;
    using content_main_t = Tcontentmain;
    using content_sub_t  = Tcontentsub;
    using tuple_t = content_main_t::tuple_t;
    using hashval_t = content_main_t::hashval_t;
    using rsv_main_node_t = ReservoirNC<entry_main_t>;
    using rsv_sub_node_t = ReservoirNC<entry_sub_t>;
    using dir_t = df::infra::FixedVectorNC<entry_main_t, true>;
    static_assert(std::is_same_v<typename content_main_t::tuple_t, typename content_sub_t::tuple_t>);
  public:
    // sub chain node
    struct entry_sub_t {
      entry_sub_t*  _next;
      content_sub_t _content;

      inline const entry_sub_t* next() const { return _next; }
      inline       entry_sub_t* next()       { return _next; }
      inline       void next(entry_sub_t* p) { _next = p; }

      inline const content_sub_t& content() const { return _content; }
      inline       content_sub_t& content()       { return _content; }

      entry_sub_t* clear();

      std::ostream& print_chain(std::ostream& os = std::cout) const;
      std::ostream& print_one(std::ostream& os = std::cout) const;
      std::ostream& print_chain_occupancy(std::ostream& os = std::cout) const;
      std::ostream& print_one_occupancy(std::ostream& os = std::cout) const;
    };

    // main collision chain node
    struct entry_main_t {
      entry_main_t*  _next;
      entry_sub_t*   _sub;
      content_main_t _content;

      inline const entry_main_t* next() const { return _next; }
      inline       entry_main_t* next()       { return _next; }
      inline       void next(entry_main_t* p) { _next = p; }

      inline const entry_sub_t* sub() const { return _sub; }
      inline       entry_sub_t* sub()       { return _sub; }
      inline       void sub(entry_sub_t* p) { _sub = p; }

      inline const content_main_t& content() const { return _content; }
      inline       content_main_t& content()       { return _content; }

      entry_main_t* clear();

      std::ostream& print_chain(std::ostream& os = std::cout, const size_t lIndent = 0) const;
      std::ostream& print_one(std::ostream& os = std::cout) const;
      std::ostream& print_chain_occupancy(std::ostream& os = std::cout, const size_t lIndent = 0) const;
      std::ostream& print_one_occupancy(std::ostream& os = std::cout) const;
    };
  public:
    HashtableNested2(const size_t aDirSize, const uint aLog2ChunkSizeMain, const uint aLog2ChunkSizeSub)
      : _dir(), _rsv_main(aLog2ChunkSizeMain), _rsv_sub(aLog2ChunkSizeSub), _dir_size(aDirSize) {}
    HashtableNested2(const size_t aDirSize, const uint aLog2ChunkSize)
      : HashtableNested2(aDirSize, aLog2ChunkSize, aLog2ChunkSize) {}
    ~HashtableNested2() = default;
  public:
    inline size_t dir_size() const { return _dir_size; }
    inline bool is_allocated() const { return (dir_size() > 0) && _dir.is_allocated(); }
  public:
    bool mem_alloc() { return _dir.mem_alloc(dir_size()); }
    bool mem_init() { return _dir.mem_init(); }
    void clear() { _rsv_sub.clear(); _rsv_main.clear(); }  // TODO also clear directory?
    bool mem_free() { _rsv_sub.erase(); _rsv_main.erase(); return _dir.mem_free(); }
  public:
    const entry_main_t* get_entry(const hashval_t aHashval) const;
          entry_main_t* get_entry(const hashval_t aHashval);
    const entry_main_t* get_entry_by_idx(const size_t aIdx) const;
          entry_main_t* get_entry_by_idx(const size_t aIdx);
  public:
    std::ostream& print(std::ostream& os = std::cout) const;
    std::ostream& print_occupancy(std::ostream& os = std::cout) const;
    std::tuple<size_t, size_t, size_t> get_stat() const;  // number of tuples, non-empty main entries, non-empty sub entries
    inline __attribute__((always_inline)) entry_main_t* new_main_entry() { return _rsv_main.get_new_entry()->clear(); }  // make sure to init!
    inline __attribute__((always_inline)) entry_sub_t*  new_sub_entry()  { return _rsv_sub.get_new_entry()->clear(); }   // otherwise, mem contains garbage!

    // testing only
    template <typename Tcmp>
    void insert(const hashval_t aHashval, const tuple_t* aTuple);
    template <typename Tprobetuple, typename Tcmp>
    const entry_main_t* probe(const hashval_t aHashval, const Tprobetuple* aProbeTuple) const;
  public:
    inline size_t mem_consumption() const { return mem_consumption_directory() + mem_consumption_reservoir(); }
    inline size_t rsv_main_card() const { return _rsv_main.cardinality(); }
    inline size_t rsv_sub_card() const { return _rsv_sub.cardinality(); }
    inline double dir_fill_factor() const;  // fraction of non-empty directory entries
  private:
    inline size_t mem_consumption_directory() const { return dir_size() * sizeof(entry_main_t); }  // in Bytes
    inline size_t mem_consumption_reservoir() const { return _rsv_main.totalsize() + _rsv_sub.totalsize(); }  // in Bytes
    inline size_t hashval_to_dir_idx(const hashval_t aHashval) const { return aHashval % dir_size(); }
  private:
    static constexpr size_t DIR_ALIGNMENT = CACHELINE_SIZE;
  private:
    dir_t _dir;                 // hash directory: fixed-size vector
    rsv_main_node_t _rsv_main;  // reservoirs for main and sub chain nodes
    rsv_sub_node_t _rsv_sub;
    size_t _dir_size;
};

/*
 * Out-of-line Function Definitions
 */

/* HashtableNested2 */

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
const HashtableNested2<Tcontentmain, Tcontentsub>::entry_main_t*
HashtableNested2<Tcontentmain, Tcontentsub>::get_entry(const hashval_t aHashval) const {
  return &(_dir[hashval_to_dir_idx(aHashval)]);
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
HashtableNested2<Tcontentmain, Tcontentsub>::entry_main_t*
HashtableNested2<Tcontentmain, Tcontentsub>::get_entry(const hashval_t aHashval) {
  return &(_dir[hashval_to_dir_idx(aHashval)]);
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
const HashtableNested2<Tcontentmain, Tcontentsub>::entry_main_t*
HashtableNested2<Tcontentmain, Tcontentsub>::get_entry_by_idx(const size_t aIdx) const {
  return &(_dir[aIdx]);
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
HashtableNested2<Tcontentmain, Tcontentsub>::entry_main_t*
HashtableNested2<Tcontentmain, Tcontentsub>::get_entry_by_idx(const size_t aIdx) {
  return &(_dir[aIdx]);
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
std::ostream&
HashtableNested2<Tcontentmain, Tcontentsub>::print(std::ostream& os) const {
  const size_t lWMargin = 2;
  const size_t lWBucketIdx = df::infra::number_of_digits(dir_size());

  for (size_t i = 0; i < _dir.size(); ++i) {
    os << std::string(lWMargin, ' ') << std::setw(lWBucketIdx) << i << " ";
    const entry_main_t* e = get_entry_by_idx(i);
    if (!e->_content.is_empty()) {
      const size_t lIndent = lWMargin + lWBucketIdx + 1;
      e->print_chain(os, lIndent);
    } else {
      os << "<empty>";
    }
    os << std::endl;
  }
  auto [lNumTuples, lNumNodesMain, lNumNodesSub] = get_stat();
  std::cout << "  " << "number of main nodes (non-empty) = " << lNumNodesMain << std::endl;
  std::cout << "  " << "number of sub nodes (non-empty)  = " << lNumNodesSub << std::endl;
  std::cout << "  " << "number of tuples                 = " << lNumTuples << std::endl;
  return os;
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
std::ostream&
HashtableNested2<Tcontentmain, Tcontentsub>::print_occupancy(std::ostream& os) const {
  const size_t lWMargin = 2;
  const size_t lWBucketIdx = df::infra::number_of_digits(dir_size());

  for (size_t i = 0; i < _dir.size(); ++i) {
    os << std::string(lWMargin, ' ') << std::setw(lWBucketIdx) << i << " ";
    const entry_main_t* e = get_entry_by_idx(i);
    if (!e->_content.is_empty()) {
      const size_t lIndent = lWMargin + lWBucketIdx + 1;
      e->print_chain_occupancy(os, lIndent);
    } else {
      os << "<empty>";
    }
    os << std::endl;
  }
  return os;
}

/* 
 * number of tuples, non-empty main entries, non-empty sub entrie
 */
template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
std::tuple<size_t, size_t, size_t>
HashtableNested2<Tcontentmain, Tcontentsub>::get_stat() const {
  size_t lNumNodesMain = 0, lNumNodesSub = 0, lNumTuples = 0;
  for (size_t i = 0; i < _dir.size(); ++i) {
    const entry_main_t* e = get_entry_by_idx(i);
    while (nullptr != e) {
      if (!e->content().is_empty()) {
        lNumTuples += e->content().size();
        ++lNumNodesMain;
        const entry_sub_t* s = e->sub();
        while (nullptr != s) {
          assert(!s->content().is_empty());
          lNumTuples += s->content().size();
          ++lNumNodesSub;
          s = s->next();
        }
      } else {
        assert(nullptr == e->next() && nullptr == e->sub());
      }

      e = e->next();
    }
  }
  return {lNumTuples, lNumNodesMain, lNumNodesSub};
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
double
HashtableNested2<Tcontentmain, Tcontentsub>::dir_fill_factor() const {
  size_t lNoNotEmpty = 0;
  for (size_t i = 0; i < _dir.size(); ++i) {
    lNoNotEmpty += !(_dir.at(i).content().is_empty());
  }
  return (static_cast<double>(lNoNotEmpty) / dir_size());
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
template <typename Tcmp>
void
HashtableNested2<Tcontentmain, Tcontentsub>::insert(const hashval_t aHashval, const tuple_t* aTuple) {
  /*
   * Insertion of tuples is supposed to be performed by the build operators of the physical algebra
   * using, e.g., get_entry, new_main_entry, new_sub_entry member functions.
   * This function is for testing only.
   *
   * To insert aHashval and aTuple...
   * - go to directory entryd at aHashval % dir_size()
   * - check if d's content is empty.
   *   - empty: add aTuple to d's content and set d's hashvalue to aHashval. end.
   *   - not empty: Walk the main collision chain nodes e until you either
   *                (a) find a node with matching hashvalue and Tcmp returns true, or
   *                (b) reach the end (e == nullptr).
   *     - (a): Check if matching node e is full
   *       - e is full: Check if sub chain node s exists and if it still has space.
   *         - true: Insert into s.
   *         - false: Create a new subchain node s and insert it at the head of e's sub chain.
   *       - e is not full: Insert into e.
   *     - (b): Create a new main node e2.
   *            Insert aHashval and aTuple into e2 and insert at head of collision chain: d -> e2 -> ...
   */

  entry_main_t* d = get_entry(aHashval);
  [[maybe_unused]] bool b = false;

  if (d->content().is_empty()) {
    // directory entry is empty, can just insert here
    assert(0 == d->content().size());
    b = d->content().do_insert(aTuple);
    d->content().hashval(aHashval);
    assert(b);
    assert(1 == d->content().size());
  } else {
    // walk the main collision chain
    entry_main_t* e = d;
    while ((nullptr != e) &&
           (e->content().hashval() != aHashval) &&
           (!Tcmp::eval(e->content().get_first_tuple(), aTuple))) {
      e = e->next();
    }
    if (nullptr == e) {
      // didn't find a main node with matching hashvalue and join key,
      // insert a new one.
      // situation before: d -> e0 -> e1 -> ...
      // situation after:  d -> e_new -> e0 -> e1 -> ...
      entry_main_t* e_new = new_main_entry();
      assert(0 == e_new->content().size());
      b = e_new->content().do_insert(aTuple);
      e_new->content().hashval(aHashval);
      assert(b);
      assert(1 == e_new->content().size());
      e_new->next(d->next());
      d->next(e_new);
    } else {
      // found a node e with matching hashvalue and join key --> tuple belongs here
      if (e->content().is_full()) {
        // main node full, must insert into sub node
        entry_sub_t* s = e->sub();
        if (nullptr == s || s->content().is_full()) {
          // no sub chain yet or first subchain node is full. create a new node.
          s = new_sub_entry();
          assert(s->content().is_empty());
          b = s->content().do_insert(aTuple);
          assert(b);
          assert(1 == s->content().size());
          s->next(e->sub());
          e->sub(s);
          assert(s->content().size() == 1);
        } else {
          // s exists and is not yet full. can safely insert.
          b = s->content().do_insert(aTuple);
          assert(b);
          assert(s->content().size() >= 2);
        }
      } else {
        b = e->content().do_insert(aTuple);
        assert(b);
        assert(!e->content().is_empty());
        assert(aHashval == e->content().hashval());
      }
    }
  }
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
template <typename Tprobetuple, typename Tcmp>
const HashtableNested2<Tcontentmain, Tcontentsub>::entry_main_t*
HashtableNested2<Tcontentmain, Tcontentsub>::probe(const hashval_t aHashval, const Tprobetuple* aTuple) const {
  /*
   * Finding a matching entry given a probe tuple using a comparison operator (join predicate).
   *
   * - Calculate the directory entry associated with the aHashval.
   * - If directory entry does not contain any tuples, return no match.
   * - Walk the main collision chain and compare each node by hashvalue and join predicate (Tcmp).
   * - Terminate if matching node is found, or no more main nodes to examine (no match).
   */
  const entry_main_t* d = get_entry(aHashval);
  // must catch the case where directory entry is empty
  if (d->content().is_empty()) {
    return nullptr;
  }
  while ((nullptr != d) &&
         (d->content().hashval() != aHashval) &&
         (!Tcmp::eval(aTuple, d->content().get_first_tuple()))) {
    d = d->next();
  }
  return d;
}



/* HashtableNested2::entry_main_t */

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
HashtableNested2<Tcontentmain, Tcontentsub>::entry_main_t*
HashtableNested2<Tcontentmain, Tcontentsub>::entry_main_t::clear() {
  _next = nullptr;
  _sub = nullptr;
  _content.clear();
  return this;
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
std::ostream&
HashtableNested2<Tcontentmain, Tcontentsub>::entry_main_t::print_chain(std::ostream& os, const size_t aIndent) const {
  const std::string lIndent(aIndent, ' ');
  const entry_main_t* e = this;
  e->print_one(os);
  e = e->next();
  while (nullptr != e) {
    os << std::endl << lIndent;
    e->print_one(os);
    e = e->next();
  }
  return os;
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
std::ostream&
HashtableNested2<Tcontentmain, Tcontentsub>::entry_main_t::print_one(std::ostream& os) const {
  size_t lHexWidth = df::infra::number_of_digits(std::numeric_limits<hashval_t>::max(), static_cast<hashval_t>(16)) + 2;  // +2 for 0x
  os << "[[ " << std::format("{0:#0{1}x}", content().hashval(), lHexWidth) << ": ";
  os << "{ M";
  for (size_t i = content().begin(); i < content().end(); ++i) {
    os << " " << (*(content().tuple(i)));
  }
  os << " }";
  const entry_sub_t* s = sub();
  if (nullptr != s) {
    os << " ";
    s->print_chain(os);
  }
  os << " ]]";
  return os;
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
std::ostream&
HashtableNested2<Tcontentmain, Tcontentsub>::entry_main_t::print_chain_occupancy(std::ostream& os, const size_t aIndent) const {
  const std::string lIndent(aIndent, ' ');
  const entry_main_t* e = this;
  e->print_one_occupancy(os);
  e = e->next();
  while (nullptr != e) {
    os << std::endl << lIndent;
    e->print_one_occupancy(os);
    e = e->next();
  }
  return os;
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
std::ostream&
HashtableNested2<Tcontentmain, Tcontentsub>::entry_main_t::print_one_occupancy(std::ostream& os) const {
  /*
   * Print one entry_main_t. Format: 
   * [[ <hashval in hex>: { M: <main node occupancy> } <list of sub nodes> ]]
   */

  // number of digits needed for hashval in hex (base 16)
  size_t lHexWidth = df::infra::number_of_digits(std::numeric_limits<hashval_t>::max(), static_cast<hashval_t>(16)) + 2; // +2 for 0x

  os << "[[ " << std::format("{0:#0{1}x}", content().hashval(), lHexWidth) << ": ";
  os << "{ M: " << content().size() << "/" << content().capacity << " }";
  const entry_sub_t* s = sub();
  if (nullptr != s) {
    os << " ";
    s->print_chain_occupancy(os);
  }
  os << " ]]";
  return os;
}

/* HashtableNested2::entry_sub_t */

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
HashtableNested2<Tcontentmain, Tcontentsub>::entry_sub_t*
HashtableNested2<Tcontentmain, Tcontentsub>::entry_sub_t::clear() {
  _next = nullptr;
  _content.clear();
  return this;
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
std::ostream&
HashtableNested2<Tcontentmain, Tcontentsub>::entry_sub_t::print_chain(std::ostream& os) const {
  const entry_sub_t* s = this;
  s->print_one(os);
  s = s->next();
  while (nullptr != s) {
    os << " ";
    s->print_one(os);
    s = s->next();
  }
  return os;
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
std::ostream&
HashtableNested2<Tcontentmain, Tcontentsub>::entry_sub_t::print_one(std::ostream& os) const {
  os << "{ S";
  for (size_t i = content().begin(); i < content().end(); ++i) {
    os << " " << (*(content().tuple(i)));
  }
  os << " }";
  return os;
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
std::ostream&
HashtableNested2<Tcontentmain, Tcontentsub>::entry_sub_t::print_chain_occupancy(std::ostream& os) const {
  const entry_sub_t* s = this;
  s->print_one_occupancy(os);
  s = s->next();
  while (nullptr != s) {
    os << " ";
    s->print_one_occupancy(os);
    s = s->next();
  }
  return os;
}

template <hj_3d_content_main_c Tcontentmain, hj_3d_content_sub_c Tcontentsub>
std::ostream&
HashtableNested2<Tcontentmain, Tcontentsub>::entry_sub_t::print_one_occupancy(std::ostream& os) const {
  /*
   * Print one entry_sub_t. Format:
   * { S: <sub node occupancy> }
   */
  os << "{ S: " << content().size() << "/" << content().capacity << " }";
  return os;
}
