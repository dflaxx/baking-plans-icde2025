#pragma once

/*
 * A templated hashtable that stores Tcontent objects in buckets of linked lists.
 *
 * The hashtable is used for a hash join, so Tcontent is assumed to contain tuple pointers
 * and pre-computed hash values. This is enforced by concept hj_ch_content_c.
 *
 * The hashtable directory has a size of _dir_size (= number of buckets).
 *
 * The hashtable directory and collision chains contain entry_t objects.
 * An entry_t object stores a pointer to the next collision chain node and a Tcontent object.
 *
 * Supported (public) operations:
 * - get a new entry_t from the reservoir
 * - get pointer to directory entry by index
 * - get pointer to directory entry by hash value
 *
 * Insert and find operations need to be implemented outside the hash table,
 * e.g., inside the algebra operators for build and probe.
 *
 */

#include "dfinfra/math.hh"
#include "dfinfra/standard_includes.hh"
#include "dfinfra/fixed_vector.hh"

#include "globals.hh"
#include "hj_htnode_content.hh"
#include "reservoirnc.hh"


template <hj_ch_content_c Tcontent>
class HashtableChained2 {
  public:
    using content_t = Tcontent;
    /*
     * Collision chain node
     */
    struct entry_t {
      entry_t* _next;
      content_t _content;

      inline const entry_t* next() const { return _next; }
      inline       entry_t* next()       { return _next; }
      inline       void next(entry_t* p) { _next = p; }

      inline const content_t& content() const { return _content; }
      inline       content_t& content()       { return _content; }

      entry_t* clear();
    };
    using node_t = entry_t;
    using tuple_t = content_t::tuple_t;
    using hashval_t = content_t::hashval_t;
    using rsv_node_t = ReservoirNC<node_t>;
    using dir_t = df::infra::FixedVectorNC<node_t, true>;
  public:
    HashtableChained2(const size_t aDirSize, const uint aLog2ChunkSize)
      : _dir(), _rsv(aLog2ChunkSize), _dir_size(aDirSize) {}
    ~HashtableChained2() = default;
  public:
    inline size_t dir_size() const { return _dir_size; }
    inline bool is_allocated() const { return (dir_size() > 0) && _dir.is_allocated(); }
  public:
    bool mem_alloc() { return _dir.mem_alloc(dir_size()); }
    bool mem_init() { return _dir.mem_init(); }
    void clear() { _rsv.clear(); }  // TODO also clear directory?
    bool mem_free() { _rsv.erase(); return _dir.mem_free(); }
  public:
    const entry_t* get_entry(const hashval_t aHashval) const;
          entry_t* get_entry(const hashval_t aHashval);
    const entry_t* get_entry_by_idx(const size_t aIdx) const;
          entry_t* get_entry_by_idx(const size_t aIdx);
  public:
    std::ostream& print(std::ostream& os = std::cout) const;
    std::ostream& print_occupancy(std::ostream& os = std::cout) const;
    inline __attribute__((always_inline)) entry_t* new_entry() { return _rsv.get_new_entry()->clear(); }
      // make sure to init/clear, otherwise mem contains garbage!

    // testing only
    void insert(const hashval_t aHashval, const tuple_t* aTuple);
  public:
    inline size_t mem_consumption() const { return mem_consumption_directory() + mem_consumption_reservoir(); }
    inline size_t rsv_card() const { return _rsv.cardinality(); }
    inline double dir_fill_factor() const;  // fraction of non-empty directory entries
  private:
    inline size_t mem_consumption_directory() const { return dir_size() * sizeof(entry_t); }  // in Bytes
    inline size_t mem_consumption_reservoir() const { return _rsv.totalsize(); }  // in Bytes
  private:
    inline size_t hashval_to_dir_idx(const hashval_t aHashval) const { return aHashval % dir_size(); }
  private:
    static constexpr size_t DIR_ALIGNMENT = CACHELINE_SIZE;
  private:
    dir_t _dir;       // hash directory: fixed-size vector
    rsv_node_t _rsv;  // reservoir for collision chain nodes
    size_t _dir_size;
};

/*
 * Out-of-line Function Definitions
 */

/* HashtableChained2 */

template <hj_ch_content_c Tcontent>
const HashtableChained2<Tcontent>::entry_t*
HashtableChained2<Tcontent>::get_entry(const hashval_t aHashval) const {
  return &(_dir[hashval_to_dir_idx(aHashval)]);
}

template <hj_ch_content_c Tcontent>
HashtableChained2<Tcontent>::entry_t*
HashtableChained2<Tcontent>::get_entry(const hashval_t aHashval) {
  return &(_dir[hashval_to_dir_idx(aHashval)]);
}

template <hj_ch_content_c Tcontent>
const HashtableChained2<Tcontent>::entry_t*
HashtableChained2<Tcontent>::get_entry_by_idx(const size_t aIdx) const {
  return &(_dir[aIdx]);
}

template <hj_ch_content_c Tcontent>
HashtableChained2<Tcontent>::entry_t*
HashtableChained2<Tcontent>::get_entry_by_idx(const size_t aIdx) {
  return &(_dir[aIdx]);
}

template <hj_ch_content_c Tcontent>
std::ostream&
HashtableChained2<Tcontent>::print(std::ostream& os) const {
  size_t lNumNodes = 0, lNumTuples = 0;
  size_t lWIdx = df::infra::number_of_digits(dir_size());
  for (size_t i = 0; i < _dir.size(); ++i) {
    std::cout << "  " << std::setw(lWIdx) << i;
    const entry_t* e = &(_dir[i]);
    while (nullptr != e) {
      std::cout << " [";
      for (uint i = e->content().begin(); i < e->content().end(); ++i) {
        os << ' ' << (*(e->content().tuple(i)));
        ++lNumTuples;
      }
      std::cout << " ]";
      lNumNodes += !(e->content().is_empty());
      e = e->next();
    }
    std::cout << std::endl;
  }
  std::cout << "  " << "number of nodes (non-empty)  = " << lNumNodes << std::endl;
  std::cout << "  " << "number of tuples             = " << lNumTuples << std::endl;
  return os;
}

template <hj_ch_content_c Tcontent>
std::ostream&
HashtableChained2<Tcontent>::print_occupancy(std::ostream& os) const {
  size_t lWIdx = df::infra::number_of_digits(dir_size());
  for (size_t i = 0; i < _dir.size(); ++i) {
    os << "  " << std::setw(lWIdx) << i;
    const entry_t* e = &(_dir[i]);
    while (nullptr != e) {
      os << " [";
      os << e->content().size() << "/" << e->content().capacity;
      os << "]";
      e = e->next();
      if (nullptr != e) {
        os << " ->";
      }
    }
    os << std::endl;
  }
  return os;
}

template <hj_ch_content_c Tcontent>
double
HashtableChained2<Tcontent>::dir_fill_factor() const {
  size_t lNoNotEmpty = 0;
  for (size_t i = 0; i < _dir.size(); ++i) {
    lNoNotEmpty += !(_dir.at(i).content().is_empty());
  }
  return (static_cast<double>(lNoNotEmpty) / dir_size());
}


template <hj_ch_content_c Tcontent>
void
HashtableChained2<Tcontent>::insert(const hashval_t aHashval, const tuple_t* aTuple) {
  /*
   * Insertion is supposed to be performed by the build operators of the physical algebra.
   * This function is just for testing.
   *
   * To insert...
   * - go to directory entry d at aHashval % dir_size()
   * - check if d's content is full.
   *   - not full: insert. end.
   *   - full: check if d has a next entry e.
   *     - no next entry e: create this entry, set next pointer and insert into e. end.
   *     - e already exists: check if full.
   *       - not full: insert into e. end.
   *       - full: create new entry f s.t. d -> f -> e and insert into f. end.
   */
  entry_t* e0 = get_entry(aHashval);
  [[maybe_unused]] bool b;
  if (!(e0->content().is_full())) {
    b = e0->content().do_insert(aTuple, aHashval);
    assert(b);
  } else {
    // before insert: e0 -> e1 (e1 may be nullptr)
    entry_t* e1 = e0->next();
    if (nullptr == e1 || e1->content().is_full()) {
      // after insert:  e0 -> e2 -> e1 (e1 may be nullptr)
      entry_t* e2 = new_entry();
      assert(e2->content().is_empty());
      b = e2->content().do_insert(aTuple, aHashval);
      assert(b);
      e2->next(e1);
      e0->next(e2);
    } else {
      b = e1->content().do_insert(aTuple, aHashval);
      assert(b);
    }
  }
}



/* HashtableChained2::entry_t */

template <hj_ch_content_c Tcontent>
HashtableChained2<Tcontent>::entry_t*
HashtableChained2<Tcontent>::entry_t::clear() {
  _next = nullptr;
  _content.clear();
  return this;
}

