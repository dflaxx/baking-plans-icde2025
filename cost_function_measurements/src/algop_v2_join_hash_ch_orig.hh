#pragma once

#include "dfinfra/standard_includes.hh"

#include "algop_v2_base.hh"
#include "hashtable_chained_v2.hh"
#include "hj_htnode_content.hh"
#include "concepts.hh"

#include <concepts>


/*
 * Regular hash join build operator using the chaining hash table v2
 *
 * In the chaining hash table (HashtableChained2),
 * both directory entries and collision chain entries are of type entry_t.
 * An entry_t has a next pointer and a payload (Tcontent).
 * The payload is a hj_ch_content_c that contains an array of tuple pointers 
 * and an equally-sized array of hash values.
 * The size is directed by this class's template parameter Tnoptr.
 *
 * Template parameters
 * - Thashfun: hash function mapping tuples to uint64_t (hash values) using eval() function
 * - Tnoptr: number of tuple pointers in the hash table's content's pointer array (default: 1)
 */
template <alg_hashfun_c Thashfun,
          uint Tnoptr = 1>
class AlgOpJoinHashChOrigBuild2 : public AlgOpBase2 {
  public:
    using hashfun_t    = Thashfun;
    using hashval_t    = typename hashfun_t::output_t;
    using input_t      = typename hashfun_t::input_t;
    using output_t     = void;
    using consumer_t   = void;
    using ht_content_t = ht_content_ch_tt<input_t, hashval_t, Tnoptr>;
    using hashtable_t  = HashtableChained2<ht_content_t>;
  public:
    AlgOpJoinHashChOrigBuild2(const size_t aHashDirSize, const uint32_t aHtLog2ChunkSize)
      : AlgOpBase2(algop_et::kAlgOpJoinChOrigBuild),
        _hashtable(aHashDirSize, aHtLog2ChunkSize) {}
  public:
    inline bool mem_alloc() { return _hashtable.mem_alloc(); }
    inline bool mem_init()  { return _hashtable.mem_init(); }
    inline void init()      { reset(); }
           void step(const input_t* aTuple);
    inline void fin() {}
    inline void clear() { _hashtable.clear(); }
    inline void mem_free() { _hashtable.mem_free(); }
  public:
    inline const hashtable_t& hashtable() const { return _hashtable; }
    inline const consumer_t* consumer() const { return nullptr; }
    inline const AlgOpBase2* consumer_poly() const { return nullptr; }
  public:
    inline size_t mem_consumption_ht() const { return _hashtable.mem_consumption(); }
    inline size_t ht_num_nodes() const { return _hashtable.rsv_card(); }
  private:
    hashtable_t _hashtable;
};

template <alg_hashfun_c Thashfun, uint Tnoptr>
void
AlgOpJoinHashChOrigBuild2<Thashfun, Tnoptr>
::step(const input_t* aTuple) {
  inc();
  const hashval_t lHashval = hashfun_t::eval(aTuple);

  constexpr bool lTrace = false;
  constexpr bool lInsertByOperator = true;  // insert by operator is the default case
  static_assert(lInsertByOperator);  // make sure we use the operator for inserts (for measurements)

  if constexpr (lInsertByOperator) {
    if constexpr (lTrace) { std::cout << __PRETTY_FUNCTION__ << ": " << "Insert by operator." << std::endl; }
    /*
     * Note: This is the same code as in hashtable_t::insert(...).
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
    using entry_t = typename hashtable_t::entry_t;
    entry_t* e0 = _hashtable.get_entry(lHashval);  // dir entry
    [[maybe_unused]] bool b = false;  // only for assertions. unused for NDEBUG.
    if (!(e0->content().is_full())) {
      b = e0->content().do_insert(aTuple, lHashval);
      assert(b);
    } else {
      // dir entry e0 is full.
      // situation before insert: e0 -> e1 (e1 may be nullptr)
      entry_t* e1 = e0->next();
      if constexpr (1 == Tnoptr) {
        // if there is only one tuple per node, we don't need to look at the state of e1
        entry_t* e2 = _hashtable.new_entry();
        assert(e2->content().is_empty());
        b = e2->content().do_insert(aTuple, lHashval);
        assert(b);
        e2->next(e1);
        e0->next(e2);
      } else {
        if (nullptr == e1 || e1->content().is_full()) {
          // cannot insert into e1, so we need a new entry e2.
          // situation after insert:  e0 -> e2 -> e1 (e1 may be nullptr)
          entry_t* e2 = _hashtable.new_entry();
          assert(e2->content().is_empty());
          b = e2->content().do_insert(aTuple, lHashval);
          assert(b);
          e2->next(e1);
          e0->next(e2);
        } else {
          // can insert into e1
          b = e1->content().do_insert(aTuple, lHashval);
          assert(b);
        }
      }
    }
  } else {
    // testing only, insert should be handled by operator, not by hashtable
    if constexpr (lTrace) { std::cout << __PRETTY_FUNCTION__ << ": " << "Insert by hashtable." << std::endl; }
    _hashtable.insert(lHashval, aTuple);
  }
}


/*
 * Regular hash join probe operator
 *
 * Probe into the chaining hashtable provided by the respective AlgOpJoinHashChOrigBuild operator.
 *
 * Note:
 * AlgBase::_count counts the number of matches when probing,
 * i.e., the number of output tuples.
 *
 * Template parameters:
 * - alg2_operator_build_c Tbuild
 *   type of the build operator that owns the hashtable (AlgOpJoinHashChOrigBuild<...>)
 * - alg_hashfun_c Thashfun
 *   type of the probe hash function, has eval(...) for probe tuples
 * - alg_binary_predicate_c Tjoinpred
 *   type of the equality join predicate, has binary eval(...) for a pair of probe/build tuples
 * - alg2_concatfun_c Tconcatfun
 *   type of the concat function, has eval_prb(...), eval_bld(...)
 * - bool IsBuildKeyUnique
 *   If the hash table was built on a unique (key) attribute, probe can terminate after the first match.
 */
template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun,
          bool IsBuildKeyUnique = false>
class AlgOpJoinHashChOrigProbe2 : public AlgOpBase2 {
  public:
    using consumer_t  = Tconsumer;
    using build_t     = Tbuild;
    using hashfun_t   = Thashfun;
    using input_t     = typename hashfun_t::input_t;
    using output_t    = typename consumer_t::input_t;
    using hashval_t   = typename hashfun_t::output_t;
    using joinpred_t  = Tjoinpred;
    using concatfun_t = Tconcatfun;
    using hashtable_t = typename build_t::hashtable_t;
    using entry_t     = typename build_t::hashtable_t::entry_t;
    using tuple_bld_t = typename build_t::input_t;
    using tuple_prb_t = input_t;
    using join_res_t  = output_t;
  public:
      inline AlgOpJoinHashChOrigProbe2(consumer_t* aConsumer, build_t* aBuildOperator)
        : AlgOpBase2(algop_et::kAlgOpJoinChOrigProbe),
          _consumer(aConsumer), _build(aBuildOperator), _output(), _numCmps(0) {}
  public:
    inline bool mem_alloc() { return _consumer->mem_alloc(); }
    inline bool mem_init() { return _consumer->mem_init(); }
    inline void init();
    inline void step(const input_t* aProbeTuple);
    inline void fin() { _consumer->fin(); }
    inline void clear() { _consumer->clear(); }
    inline void mem_free() { _consumer->mem_free(); }
  public:
    inline const consumer_t* consumer() const { return _consumer; }
    inline const AlgOpBase2* consumer_poly() const { return _consumer; }
    inline       uint64_t    numCmps()  const { return _numCmps; }
  private:
    consumer_t* _consumer;  
    build_t*    _build;     // AlgOpJoinHashChOrigBuild, has getter method hashtable()
    output_t    _output;    // concat'd output tuple: <tuple_left_t*, tuple_right_t*>
    uint64_t    _numCmps;   // total number of collision chain comparisons
};

template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun,
          bool IsBuildKeyUnique>
void
AlgOpJoinHashChOrigProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun, IsBuildKeyUnique>
::init() {
  //std::cout << "  AlgOpJoinHashChOrigProbe2::init()" << std::endl;
  reset();
  _numCmps = 0;
  _consumer->init();
}

template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun,
          bool IsBuildKeyUnique>
void
AlgOpJoinHashChOrigProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun, IsBuildKeyUnique>
::step(const input_t* aProbeTuple) {
  constexpr bool lTrace = false;
  if constexpr (lTrace) { std::cout << "  AlgOpJoinHashChOrigProbe2::step() -> probe tuple: " << *aProbeTuple << "\n"; }

  /*
   * To probe:
   * Let t := aProbeTuple.
   * 1. Compute hashvalue h from t using Thashfun.
   * 2. Get entry_t e in the hash directory using h.
   * 3. Repeat:
   *    - Compare h to all of e's hash values. e's content might hold more than one hash value if Tnoptr > 1.
   *      If hash values match, evaluate join predicate on t and the matching tuple s from e.
   *      If join predicate is true, concatenate t and s using Tconcatfun and push to consumer.
   *      If IsBuildKeyUnique, terminate.
   *    - e := e->next
   *    ... while e != nullptr.
   */

  using hashval_bld_t = typename build_t::hashval_t;
  static_assert(std::is_same<hashval_t, hashval_bld_t>::value);

  const hashtable_t& lHt = _build->hashtable();

  hashval_t lProbeHash = hashfun_t::eval(aProbeTuple);

  concatfun_t::eval_prb(_output, aProbeTuple);

  const entry_t* e = lHt.get_entry(lProbeHash);  // ptr calc, no mem access
  bool lMore = true;  // helper variable to break outer while loop from inner for loop

  // walk the collision chain
  while ((nullptr != e) && lMore) {
    // check all hash values and tuple pointers in the ptr array
    for (uint i = e->content().begin(); i < e->content().end(); ++i) {
      assert(!e->content().is_empty());
      const tuple_bld_t* s = e->content().tuple(i);          // get build tuple s
      const hashval_t lBuildHash = e->content().hashval(i);  // and its hash value
      ++_numCmps;
      if ((lProbeHash == lBuildHash) && (joinpred_t::eval(aProbeTuple, s))) {  // access *s
        concatfun_t::eval_bld(_output, s);
        _consumer->step(&_output);
        inc();
        if constexpr (IsBuildKeyUnique) {
          // this is the only join match
          if constexpr (lTrace) {
            std::cout << "    break from collision chain iteration because IsBuildKeyUnique = " << IsBuildKeyUnique << std::endl;
          }
          lMore = false;  // = break outer while
          break;  // break for
        }
      }
    }
    e = e->next();
  }
}
