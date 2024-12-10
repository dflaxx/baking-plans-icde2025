#pragma once

#include "dfinfra/standard_includes.hh"
#include "dfinfra/exceptions.hh"

#include "algop_v2_base.hh"
#include "hashtable_nested_v2.hh"
#include "hj_htnode_content.hh"
#include "concepts.hh"

/*
 * 3D hash join build operator using the nested hash table v2
 * with Rolling Prefetching
 *
 * In the 3D hash table (HashtableNested2),
 * directory entries and main collision chain entries are of type entry_main_t,
 * sub chain entries are of type entry_sub_t.
 * All entry_main_t objects store only pointers to build tuples with equal join attribute value(s).
 * They have a next pointer (next main entry with different join attribute value),
 * a sub chain pointer (more tuple pointers with the same join attribute value),
 * and a payload (content).
 * The payload is a hj_3d_content_main_c which consists of a build hash value (for a cheap pre-test)
 * and an array of tuple pointers with a counter.
 * The number of tuple pointers is directed by this class's template parameter Tnoptrmain.
 * All entry_sub_t store a pointer to the next sub chain entry and a payload
 * (a pointer/counter array, similar to entry_main_t's payload).
 * The number of pointers in the sub entries is directed by this class's template parameter Tnoptrsub.
 *
 * Template parameters:
 * - Thashfun: build hash function mapping tuples to hash values (usually uint64_t) using the eval() function
 * - Tequalfun: equality comparison function that compares two tuples for join attribute equality
 *   using its eval() function
 * - Tnoptrmain, Tnoptrsub: maximum number of tuple pointers in main and sub chain entries, respectively
 */
template <alg_hashfun_c Thashfun,
          alg_binary_predicate_c Tequalfun,
          uint Tnoptrmain,
          uint Tnoptrsub>
class AlgOpJoinHash3dRpBuild2 : public AlgOpBase2 {
  public:
    using hashfun_t   = Thashfun;
    using hashval_t   = typename hashfun_t::output_t;
    using input_t     = typename hashfun_t::input_t;
    using tuple_t     = input_t;
    using output_t    = void;
    using consumer_t  = void;
    using eqfun_t     = Tequalfun;
    using ht_content_sub_t = ht_content_3d_sub_tt<input_t, Tnoptrsub>;
    using ht_content_main_t = ht_content_3d_main_tt<hashval_t, ht_content_sub_t, Tnoptrmain>;
    using hashtable_t = HashtableNested2<ht_content_main_t, ht_content_sub_t>;
    using entry_main_t = typename hashtable_t::entry_main_t;
    using entry_sub_t = typename hashtable_t::entry_sub_t;
  public:
    AlgOpJoinHash3dRpBuild2(const size_t aHashDirSize, const uint32_t aHtLog2ChunkSizeMain, const uint32_t aHtLog2ChunkSizeSub)
      : AlgOpBase2(algop_et::kAlgOpJoin3dRpBuild),
        _hashtable(aHashDirSize, aHtLog2ChunkSizeMain, aHtLog2ChunkSizeSub) {}
    AlgOpJoinHash3dRpBuild2(const size_t aHashDirSize, const uint32_t aHtLog2ChunkSize)
      : AlgOpJoinHash3dRpBuild2(aHashDirSize, aHtLog2ChunkSize, aHtLog2ChunkSize) {}
  public:
    static constexpr uint RP_BUFFER_SIZE = 12;
    // one entry in the ring buffer
    struct rp_t {
      const tuple_t*       _tuple = nullptr;  // the tuple ptr to insert into the hashtable
            entry_main_t*  _ht_entry;         // the hash table diretory entry where _tuple should be inserted
            hashval_t      _hashval;          // the tuple's hash value, needed for insert
      inline void clear() { _tuple = nullptr; }
      inline bool is_empty() const { return (_tuple == nullptr); }
      inline bool is_occupied() const { return !is_empty(); }
    };
  public:
    inline bool mem_alloc() { return _hashtable.mem_alloc(); }
    inline bool mem_init()  { return _hashtable.mem_init(); }
    inline void init();
           void step(const input_t* aInput);
    inline void fin();
    inline void clear() { _hashtable.clear(); }
    inline void mem_free() { _hashtable.mem_free(); }
  public:
    inline const hashtable_t& hashtable() const { return _hashtable; }
    inline const consumer_t* consumer() const { return nullptr; }
    inline const AlgOpBase2* consumer_poly() const { return nullptr; }
  public:
    inline size_t mem_consumption_ht() const { return _hashtable.mem_consumption(); }
    inline size_t ht_num_main_nodes() const { return _hashtable.rsv_main_card(); }
    inline size_t ht_num_sub_nodes() const { return _hashtable.rsv_sub_card(); }
  private:
    // increment _curr in range [0, AMAC_BUFFER_SIZE - 1]
    inline void rp_advance() { _curr = (++_curr >= RP_BUFFER_SIZE ? 0 : _curr); }
    // process (= probe) a single buffer entry, NO sanity
    inline void process_buffer_entry(rp_t& aRp);
  private:
    hashtable_t _hashtable;
    rp_t        _rp_buffer[RP_BUFFER_SIZE];  // the RP ring buffer
    uint        _curr;                       // current index into the ring buffer, values [0, RP_BUFFER_SIZE - 1]
};

template <alg_hashfun_c Thashfun, alg_binary_predicate_c Tequalfun, uint Tnoptrmain, uint Tnoptrsub>
void
AlgOpJoinHash3dRpBuild2<Thashfun, Tequalfun, Tnoptrmain, Tnoptrsub>
::init() {
  reset();
  for (uint i = 0; i < RP_BUFFER_SIZE; ++i) { _rp_buffer[i].clear(); }
  _curr = 0;
}

template <alg_hashfun_c Thashfun, alg_binary_predicate_c Tequalfun, uint Tnoptrmain, uint Tnoptrsub>
void
AlgOpJoinHash3dRpBuild2<Thashfun, Tequalfun, Tnoptrmain, Tnoptrsub>
::step(const input_t* aTuple) {
  constexpr bool lTrace = false;
  constexpr bool lInsertByOperator = true;

  if constexpr (lTrace) { std::cout << "  AlgOpJoinHash3dRpBuild2::step() -> build tuple #" << count() << std::endl; }
  inc();

  const hashval_t lHashval = hashfun_t::eval(aTuple);
  entry_main_t* lDirEntry = _hashtable.get_entry(lHashval);

  if constexpr (lInsertByOperator) {
    if constexpr (lTrace) {
      //std::cout << __PRETTY_FUNCTION__ << ": " << "Insert by operator." << std::endl;
      std::cout << "    " << "aTuple = " << *aTuple << std::endl;
      std::cout << "    " << "_curr = " << _curr << std::endl;
      std::cout << "    " << "_rp_buffer[_curr].is_empty() = " << _rp_buffer[_curr].is_empty() << std::endl;
    }
    rp_t& lRp = _rp_buffer[_curr];
    if (lRp.is_occupied()) {
      process_buffer_entry(lRp);
      lRp.clear();
    }

    // insert the new build tuple aTuple into the newly cleared buffer entry
    if constexpr (lTrace) {
      std::cout << "    " << "inserting tuple " << *aTuple << " into the RP buffer." << std::endl;
    }
    assert(lRp.is_empty());
    lRp._tuple = aTuple;
    lRp._hashval = lHashval;
    lRp._ht_entry = lDirEntry;
    PREFETCH_BUILD(lDirEntry);
    rp_advance();

  } else {  // lInsertByOperator == false
    // XXX testing only, insert should be handled by operator, not hashtable
    if constexpr (lTrace) { std::cout << __PRETTY_FUNCTION__ << ": " << "Insert by hashtable." << std::endl; }
    _hashtable.template insert<eqfun_t>(lHashval, aTuple);
  }
}

template <alg_hashfun_c Thashfun, alg_binary_predicate_c Tequalfun, uint Tnoptrmain, uint Tnoptrsub>
void
AlgOpJoinHash3dRpBuild2<Thashfun, Tequalfun, Tnoptrmain, Tnoptrsub>
::fin() {
  constexpr bool lTrace = false;
  if constexpr (lTrace) { std::cout << "  AlgOpJoinHash3dRpBuild2::fin()" << std::endl; }

  // move body of for-loop into function and use both here and in step(), code is identical
  for (_curr = 0; _curr < RP_BUFFER_SIZE; ++_curr) {
    rp_t& lRp = _rp_buffer[_curr];
    if (lRp.is_empty()) { continue; }
    if constexpr (lTrace) {
      std::cout << "    " << "processing the RP buffer: " << "_curr = " << _curr << ", tuple " << *lRp._tuple << std::endl;
    }
    process_buffer_entry(lRp);
    lRp.clear();
  }
}

template <alg_hashfun_c Thashfun, alg_binary_predicate_c Tequalfun, uint Tnoptrmain, uint Tnoptrsub>
void
AlgOpJoinHash3dRpBuild2<Thashfun, Tequalfun, Tnoptrmain, Tnoptrsub>
::process_buffer_entry(rp_t& aRp) {
  assert(aRp.is_occupied());

  /*
   * Note: This is the same code as in hashtable_t::insert(...).
   *
   * To insert lHashval and aTuple...
   * - go to directory entry d at lHashval % dir_size()
   * - check if d's content is empty.
   *   - empty: add aTuple to d's content and set d's hashvalue to lHashval. end.
   *   - not empty: Walk the main collision chain nodes e until you either
   *                (a) find a node with matching hashvalue and eqfun_t::eval(...) returns true, or
   *                (b) reach the end (e == nullptr).
   *     - (a): Check if matching node e is full
   *       - e is full: Check if sub chain node s exists and if it still has space.
   *         - true: Insert into s.
   *         - false: Create a new subchain node s and insert it at the head of e's sub chain.
   *       - e is not full: Insert into e.
   *     - (b): Create a new main node e2.
   *            Insert lHashval and aTuple into e2 and insert at head of collision chain: d -> e2 -> ...
   */


  entry_main_t* d = aRp._ht_entry;
  [[maybe_unused]] bool b = false;  // result of insert, only for NDEBUG

  if (d->content().is_empty()) {
    // directory entry is empty, can just insert here
    assert(0 == d->content().size());
    b = d->content().do_insert(aRp._tuple);
    d->content().hashval(aRp._hashval);
    assert(b);
    assert(1 == d->content().size());
  } else {
    // walk the main collision chain
    entry_main_t* e = d;
    while ((nullptr != e) &&
           (e->content().hashval() != aRp._hashval) &&                        // compare hash values
           (!(eqfun_t::eval(e->content().get_first_tuple(), aRp._tuple)))) {  // evaluate join predicate
      e = e->next();
    }
    if (nullptr == e) {
      // didn't find a main node with matching hashvalue and join key,
      // need new main node.
      // situation before: d -> e0 -> e1 -> ...
      // situation after:  d -> e_new -> e0 -> e1 -> ...
      entry_main_t* e_new = _hashtable.new_main_entry();
      assert(0 == e_new->content().size());
      b = e_new->content().do_insert(aRp._tuple);
      e_new->content().hashval(aRp._hashval);
      assert(b);
      assert(1 == e_new->content().size());
      e_new->next(d->next());
      d->next(e_new);
    } else {
      // found a node e with matching hashvalue and join key --> tuple belongs here
      if (e->content().is_full()) {
        // main node full, must insert into sub node
        if constexpr (1 == Tnoptrsub) {
          // only one tuple per sub node, so no need to look at the current state of the sub chain,
          // just create a new one and insert
          entry_sub_t* s = _hashtable.new_sub_entry();
          assert(s->content().is_empty());
          b = s->content().do_insert(aRp._tuple);
          assert(b);
          assert(1 == s->content().size());
          s->next(e->sub());
          e->sub(s);
          assert(s->content().size() == 1);
        } else {
          entry_sub_t* s = e->sub();
          if (nullptr == s || s->content().is_full()) {
            // no sub chain yet or first subchain node is full. create a new node.
            s = _hashtable.new_sub_entry();
            assert(s->content().is_empty());
            b = s->content().do_insert(aRp._tuple);
            assert(b);
            assert(1 == s->content().size());
            s->next(e->sub());
            e->sub(s);
            assert(s->content().size() == 1);
          } else {
            // s exists and is not yet full. can safely insert.
            assert(!s->content().is_full());
            b = s->content().do_insert(aRp._tuple);
            assert(b);
            assert(s->content().size() >= 2);
          }
        }
      } else {
        // matching main node e is not yet full, so insert there
        b = e->content().do_insert(aRp._tuple);
        assert(b);
        assert(!e->content().is_empty());
        assert(aRp._hashval == e->content().hashval());
      }
    }
  }
}

/*
 * 3D hash join probe operator with unnest
 * with Rolling Prefetching
 *
 * Probe into the nested hashtable provided by the respective AlgOpJoinHash3dRpBuild2 operator.
 *
 * Note:
 * AlgOpBase2::_count counts the number of output tuples going from this operator to its consumer.
 * AlgOpJoinHash3dAmacProbe2::_c_mid counts the number of distinct probe tuples finding at least one match,
 * i.e., the size of Prb \lsjoin Bld.
 *
 * Template parameters:
 * - Tbuild: the build operator whose hashtable is probed by this operator
 * - Thashfun: probe hash function mapping tuples to hash values (usually uint64_t) using the eval() function
 * - Tjoinpred: join predicate (equality comparison) function that compares a probe and a build tuple
 *   for join attribute equality using its eval() function
 * - Tconcatfun: concatenation function that concatenates a matching probe and build tuple into an output_t tuple.
 */
template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun>
class AlgOpJoinHash3dRpProbe2 : public AlgOpBase2 {
  public:
    using consumer_t   = Tconsumer;
    using build_t      = Tbuild;
    using hashfun_t    = Thashfun;
    using input_t      = typename hashfun_t::input_t;
    using output_t     = typename consumer_t::input_t;
    using joinpred_t   = Tjoinpred;
    using concatfun_t  = Tconcatfun;
    using hashval_t    = typename hashfun_t::output_t;
    using hashtable_t  = typename build_t::hashtable_t;
    using entry_main_t = typename hashtable_t::entry_main_t;
    using entry_sub_t  = typename hashtable_t::entry_sub_t;
    using tuple_bld_t  = typename build_t::input_t;
    using tuple_prb_t  = input_t;
    using join_res_t   = output_t;
  public:
    AlgOpJoinHash3dRpProbe2(consumer_t* aConsumer, build_t* aBuildOperator)
      : AlgOpBase2(algop_et::kAlgOpJoin3dRpProbeWithUnnest),
        _consumer(aConsumer), _build(aBuildOperator), _output(), _numCmps1(0), _numCmps2(0), _c_mid(0),
        _rp_buffer(), _curr(0) {}
  public:
    static constexpr uint RP_BUFFER_SIZE = 12;
    // one entry in the ring buffer
    struct rp_t {
      const tuple_prb_t*   _tuple_prb = nullptr;  // the probe tuple ptr
      const entry_main_t*  _ht_entry;             // the hash table diretory corresponding to _tuple_prb
            hashval_t      _hashval;              // the hash value of _tuple_prb
      inline void clear() { _tuple_prb = nullptr; }
      inline bool is_empty() const { return (_tuple_prb == nullptr); }
      inline bool is_occupied() const { return !is_empty(); }
    };
  public:
    inline bool mem_alloc() { return _consumer->mem_alloc(); }
    inline bool mem_init() { return _consumer->mem_init(); }
    inline void init();
    inline void step(const input_t* aProbeTuple);
    inline void fin();
    inline void clear() { _consumer->clear(); }
    inline void mem_free() { _consumer->mem_free(); }
  public:
    inline const consumer_t* consumer() const { return _consumer; }
    inline const AlgOpBase2* consumer_poly() const { return _consumer; }
    inline       uint64_t    numCmps1()  const { return _numCmps1; }
    inline       uint64_t    numCmps2()  const { return _numCmps2; }
  private:
    // increment _curr in range [0, AMAC_BUFFER_SIZE - 1]
    inline void rp_advance() { _curr = (++_curr >= RP_BUFFER_SIZE ? 0 : _curr); }
    // process (= probe) a single buffer entry, NO sanity
    inline void process_buffer_entry(rp_t& aRp);
  private:
    consumer_t* _consumer;
    build_t*    _build;     // AlgOpJoinHash3dRpBuild2, has getter method hashtable()
    output_t    _output;    // unnested output tuple
    uint64_t    _numCmps1;  // total number of collision chain *hashval* comparisons (requires entry_main_t access)
    uint64_t    _numCmps2;  // total number of collision chain *joinpred* evaluations (requires tuple access)
    uint64_t    _c_mid;     // number of *nested* result tuples before unnesting
                            // (= number of qualifying probe tuples = |Prb lsjoin Bld|).
                            // Let R := key relation, S := foreign key relation.
                            // For build unique (S join R): c_mid = |S lsjoin R| = |S join R| = c_res ==> unnesting: 1 -> 1
                            // For build non-unique (R join S): c_mid = |R lsjoin S| <= |R join S| = c_res ==> unnesting 1 -> n
    rp_t        _rp_buffer[RP_BUFFER_SIZE];  // the RP ring buffer
    uint        _curr;                       // current index into the ring buffer, values [0, RP_BUFFER_SIZE - 1]
};

template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun>
void
AlgOpJoinHash3dRpProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun>
::init() {
  reset();
  _numCmps1 = 0;
  _numCmps2 = 0;
  _c_mid = 0;
  for (uint i = 0; i < RP_BUFFER_SIZE; ++i) { _rp_buffer[i].clear(); }
  _curr = 0;
  _consumer->init();
}

template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun>
void
AlgOpJoinHash3dRpProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun>
::step(const input_t* aProbeTuple) {
  constexpr bool lTrace = false;
  if constexpr (lTrace) { std::cout << "  AlgOpJoinHash3dRpProbe2::step() -> probe tuple: " << *aProbeTuple << std::endl; }

  using hashval_bld_t = typename build_t::hashval_t;
  static_assert(std::is_same<hashval_t, hashval_bld_t>::value);

  const hashtable_t& lHt = _build->hashtable();
  const hashval_t lArgTupleHash = hashfun_t::eval(aProbeTuple);
  const entry_main_t* lArgTupleDirEntry = lHt.get_entry(lArgTupleHash);

  // check and process current buffer entry
  rp_t& lRp = _rp_buffer[_curr];
  if (lRp.is_occupied()) {
    process_buffer_entry(lRp);
    lRp.clear();
  }

  // insert the argument tuple aProbeTuple into the newly cleared buffer entry
  if constexpr (lTrace) {
    std::cout << "    " << "inserting tuple " << *aProbeTuple << " into the RP buffer." << std::endl;
  }
  assert(lRp.is_empty());
  lRp._tuple_prb = aProbeTuple;
  lRp._hashval = lArgTupleHash;
  lRp._ht_entry = lArgTupleDirEntry;
  PREFETCH_PROBE(lArgTupleDirEntry);
  rp_advance();
}
template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun>
void
AlgOpJoinHash3dRpProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun>
::fin() {
  constexpr bool lTrace = false;
  if constexpr (lTrace) { std::cout << "  AlgOpJoinHash3dRpProbe2::fin()" << std::endl; }

  for (_curr = 0; _curr < RP_BUFFER_SIZE; ++_curr) {
    rp_t& lRp = _rp_buffer[_curr];
    if (lRp.is_empty()) { continue; }
    if constexpr (lTrace) {
      std::cout << "    " << "processing the RP buffer: " << "_curr = " << _curr << ", tuple " << *lRp._tuple_prb << std::endl;
    }
    process_buffer_entry(lRp);
    lRp.clear();
  }
  _consumer->fin();
}

// Processes a single buffer element to completion = probe.
// Does NOT check if buffer element aRp is occupied.
// Does NOT clear the buffer element aRp.
// Does NOT advance _curr.
template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun>
void
AlgOpJoinHash3dRpProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun>
::process_buffer_entry(rp_t& aRp) {
  constexpr bool lTrace = false;
  assert(aRp.is_occupied());

  /*
   * Probe:
   * Find a matching build tuples given a probe tuple using a comparison operator (join predicate).
   *
   * Let t := aProbeTuple.
   * 1. Calculate hashvalue h from t using hashfun_t.
   * 2. Get entry_main_t d in the hash diretory using h.
   * 3. Check if d contains at least one tuple.
   *    - If not, no match. Terminate.
   *    - Otherwise, walk the main collision chain (entry_main_t e) and compare each node by hashvalue and join predicate (joinpred_t).
   * 4. If matching node e is found
   *    - concat() and step() for all tuples stored in e.
   *    - Walk e's sub chain (nodes s of type entry_sub_t) if it exists: concat() and step() for all tuples stored in s.
   */

  const entry_main_t* e = aRp._ht_entry;

  // must catch the case where directory entry is empty
  if (e->content().is_empty()) {  // access *e
    if constexpr (lTrace) {
      std::cout << "    " << "empty directory entry" << std::endl;
    }
  } else {
    uint64_t lNumComparisons1 = 0;  // number of collision chain *hashval* comparisons (-> access e)
    uint64_t lNumComparisons2 = 0;  // number of collision chain *joinpred* evaluations (-> access tuple)

    while ((nullptr != e)) {
      ++lNumComparisons1;
      if (e->content().hashval() == aRp._hashval) {
        ++lNumComparisons2;
        if (joinpred_t::eval(aRp._tuple_prb, e->content().get_first_tuple())) {
          break;
        }
      }
      e = e->next();  // access *e
    }

    _numCmps1 += lNumComparisons1;
    _numCmps2 += lNumComparisons2;

    if (nullptr != e) {
      // e is a match
      ++_c_mid;  // found a matching main node = a matching distinct join key
      if constexpr (lTrace) {
        std::cout << "    " << "found matching main node: " << e << std::endl;
      }
      // UNNEST
      // step 1: unnest e
      if constexpr (lTrace) {
        std::cout << "    " << "unnest main node" << std::endl;
      }
      concatfun_t::eval_prb(_output, aRp._tuple_prb);
      for (uint i = e->content().begin(); i < e->content().end(); ++i) {
        concatfun_t::eval_bld(_output, e->content().tuple(i));
        if constexpr (lTrace) {
          std::cout << "      " << "rs[ " << _output << " ]" << std::endl;
        }
        inc();
        _consumer->step(&_output);
      }
      // step 2: unnest subchain
      if constexpr (lTrace) {
        std::cout << "    " << "unnest sub nodes" << std::endl;
      }
      const entry_sub_t* s = e->sub();
      while (nullptr != s) {
        // access *s
        for (uint i = s->content().begin(); i < s->content().end(); ++i) {
          concatfun_t::eval_bld(_output, s->content().tuple(i));
          if constexpr (lTrace) {
            std::cout << "      " << "rs[ " << _output << " ]" << std::endl;
          }
          inc();
          _consumer->step(&_output);
        }
        s = s->next();
      }
    } else {
      if constexpr (lTrace) {
        std::cout << "    " << "no match found." << std::endl;
      }
    }
  }
}
