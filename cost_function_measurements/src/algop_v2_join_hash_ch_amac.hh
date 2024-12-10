#pragma once

#include "dfinfra/standard_includes.hh"

#include "algop_v2_base.hh"
#include "hashtable_chained_v2.hh"
#include "hj_htnode_content.hh"
#include "concepts.hh"
#include "prefetch.hh"

#include "gminfra/bit_subsets.hh"

#include <bitset>
#include <concepts>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <format>


/*
 * Regular hash join build operator using the chaining hash table v2
 * with AMAC
 *
 * In the chaining hash table (HashtableChained2),
 * both directory entries and collision chain entries are of type entry_t.
 * An entry_t has a next pointer and a payload (Tcontent).
 * The payload is a hj_ch_content_c that contains an array of tuple pointers 
 * and an equally-sized array of hash values.
 * The size if directed by this class's template parameter Tnoptr.
 *
 * Template parameters
 * - Thashfun: hash function mapping tuples to uint64_t (hash values) using eval() function
 * - Tnoptr: number of tuple pointers in the hash table's content's pointer array (default: 1)
 */
template <alg_hashfun_c Thashfun,
          uint Tnoptr = 1>
class AlgOpJoinHashChAmacBuild2 : public AlgOpBase2 {
  public:
    using hashfun_t    = Thashfun;
    using hashval_t    = typename hashfun_t::output_t;
    using input_t      = typename hashfun_t::input_t;
    using output_t     = void;
    using tuple_t      = input_t;
    using consumer_t   = void;
    using ht_content_t = ht_content_ch_tt<input_t, hashval_t, Tnoptr>;
    using hashtable_t  = HashtableChained2<ht_content_t>;
    using entry_t      = typename hashtable_t::entry_t;
  public:
    AlgOpJoinHashChAmacBuild2(const size_t aHashDirSize, const uint32_t aHtLog2ChunkSize)
      : AlgOpBase2(algop_et::kAlgOpJoinChAmacBuild),
        _hashtable(aHashDirSize, aHtLog2ChunkSize),
        _amac(),
        _curr(0) {}
  public:
    static constexpr uint AMAC_BUFFER_SIZE = 12;  // size of the AMAC ring buffer
    // the different AMAC stages
    enum class amac_todo_et {
      kAmacEmpty          = 0,
      kAmacInsertDirEntry = 1,
      kAmacInsertChain    = 2
    };
    // one entry in the AMAC ring buffer
    struct amac_t {
      const tuple_t*     _tuple;     // the tuple ptr to insert into the hashtable
            entry_t*     _ht_entry;  // the hash table directory entry where _tuple should be inserted
            hashval_t    _hashval;   // the _tuple's hash value
            amac_todo_et _todo;      // the AMAC stage
      inline void clear() { _todo = amac_todo_et::kAmacEmpty; }
      inline bool is_empty() const { return (_todo == amac_todo_et::kAmacEmpty); }
      inline bool is_occupied() const { return (_todo != amac_todo_et::kAmacEmpty); }
      // get _todo as an integer for easier output
      inline std::underlying_type_t<amac_todo_et> todo() const { return static_cast<std::underlying_type_t<amac_todo_et>>(_todo); }
    };
  public:
    inline bool mem_alloc() { return _hashtable.mem_alloc(); }
    inline bool mem_init()  { return _hashtable.mem_init(); }
    inline void init();
           void step(const input_t* aTuple);
           void fin();
    inline void clear();
    inline void mem_free() { _hashtable.mem_free(); }
  public:
    inline const hashtable_t& hashtable() const { return _hashtable; }
    inline const consumer_t* consumer() const { return nullptr; }
    inline const AlgOpBase2* consumer_poly() const { return nullptr; }
  public:
    inline size_t mem_consumption_ht() const { return _hashtable.mem_consumption(); }
    inline size_t ht_num_nodes() const { return _hashtable.rsv_card(); }  // collision chain nodes in rsv
  private:
    // increment _curr in range [0, AMAC_BUFFER_SIZE - 1]
    inline void amac_advance() { _curr = (++_curr >= AMAC_BUFFER_SIZE ? 0 : _curr); }
  private:
    hashtable_t _hashtable;
    amac_t      _amac[AMAC_BUFFER_SIZE];  // the AMAC ring buffer
    uint        _curr;                    // current index into the AMAC ring buffer, value range [0, AMAC_BUFFER_SIZE - 1]
};

template <alg_hashfun_c Thashfun, uint Tnoptr>
void
AlgOpJoinHashChAmacBuild2<Thashfun, Tnoptr>
::init() {
  reset();
  for (uint i = 0; i < AMAC_BUFFER_SIZE; ++i) { _amac[i].clear(); }
  _curr = 0;
}

template <alg_hashfun_c Thashfun, uint Tnoptr>
void
AlgOpJoinHashChAmacBuild2<Thashfun, Tnoptr>
::clear() {
  _hashtable.clear();
  for (uint i = 0; i < AMAC_BUFFER_SIZE; ++i) { _amac[i].clear(); }
}

template <alg_hashfun_c Thashfun, uint Tnoptr>
void
AlgOpJoinHashChAmacBuild2<Thashfun, Tnoptr>
::step(const input_t* aTuple) {
  constexpr bool lTrace = false;
  constexpr bool lUseAmac = true;
  static_assert(lUseAmac);  // make sure we use AMAC for the measurements

  if constexpr (lTrace) { std::cout << "  AlgOpJoinHashChAmacBuild2::step() -> build tuple #" << count() << std::endl; }
  inc();

  if constexpr (lUseAmac) {
    const hashval_t lHashval = hashfun_t::eval(aTuple);   // access *aTuple
    entry_t* lDirEntry = _hashtable.get_entry(lHashval);  // no memory access, only pointer computation

    if constexpr (lTrace) {
      std::cout << "    " << "aTuple = " << *aTuple << std::endl;
      std::cout << "    " << "_curr = " << _curr << std::endl;
      std::cout << "    " << "_amac[_curr].is_empty() = " << _amac[_curr].is_empty() << std::endl;
    }

    while (_amac[_curr].is_occupied()) {
      amac_t& lAmac = _amac[_curr];
      if constexpr (lTrace) {
        std::cout << "    " << "processing the AMAC buffer: "
          << "_curr = " << _curr << ", _todo = " << lAmac.todo() << ", tuple " << *lAmac._tuple << std::endl;
      }
      // process tuples in the AMAC ring buffer until a free spot is found for the new tuple aTuple
      entry_t* e_dir = lAmac._ht_entry;  // directory entry
      
      // stage 1: directory entry d has been prefetched, test it for insert
      if (amac_todo_et::kAmacInsertDirEntry == lAmac._todo) {
        if (e_dir->content().do_insert(lAmac._tuple, lAmac._hashval)) {
          // Insert was successful. Done.
          lAmac.clear();
          if constexpr (lTrace) { std::cout << "      " << "todo: 1 -> 0. done." << std::endl; }
          break;  // from the while loop. Don't advance _curr, so next tuple is inserted in this free spot.
        } else {
          // Insert not successful (e_dir is full). Must insert into the collision chain.
          assert(e_dir->content().is_full());
          lAmac._todo = amac_todo_et::kAmacInsertChain;
          PREFETCH_BUILD(e_dir->next());
          if constexpr (lTrace) { std::cout << "      " << "todo: 1 -> 2." << std::endl; }
        }

      // stage 2: collision chain entry e_nxt has been prefetched in stage 1. Test if its full and insert.
      } else if (amac_todo_et::kAmacInsertChain == lAmac._todo) {
        entry_t* e_nxt = e_dir->next();
        // check if e_nxt is a) present and b) has space
        if ((nullptr == e_nxt) || !(e_nxt->content().do_insert(lAmac._tuple, lAmac._hashval))) {
          // e_nxt does not yet exist or is full; need to create new entry in the chain
          entry_t* e_new = _hashtable.new_entry();
          [[maybe_unused]] bool b = e_new->content().do_insert(lAmac._tuple, lAmac._hashval);
          assert(b);
          e_dir->_next = e_new;
          e_new->_next = e_nxt;
          if constexpr (lTrace) { std::cout << "      " << "todo: 2 -> 0 [by insert in e_new]. done." << std::endl; }
        } else {
          if constexpr (lTrace) { std::cout << "      " << "todo: 2 -> 0 [by insert in e_nxt]. done." << std::endl; }
        }
        lAmac.clear();
        break;

      } else {
        throw std::logic_error("Unexpected AMAC stage in AlgOpJoinHashChAmacBuild2::step(): " + std::to_string(_amac[_curr].todo()));
      }
      amac_advance();
    }  // end while: amac occupied

    // insert aTuple into the free spot in the ring buffer and issue a prefetch
    if constexpr (lTrace) {
      std::cout << "    " << "inserting tuple " << *aTuple << " into the AMAC buffer." << std::endl;
    }
    assert(_amac[_curr].is_empty());
    _amac[_curr]._tuple = aTuple;
    _amac[_curr]._ht_entry = lDirEntry;
    _amac[_curr]._hashval = lHashval;
    _amac[_curr]._todo = amac_todo_et::kAmacInsertDirEntry;
    PREFETCH_BUILD(lDirEntry);
    amac_advance();

  } else {
    //std::cout << "  Warning: AlgOpJoinHashChAmacBuild2::step() does not use AMAC!" << std::endl;
    const hashval_t lHashval = hashfun_t::eval(aTuple);

    constexpr bool lInsertByOperator = true;  // insert by operator is the default case

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
      entry_t* e0 = _hashtable.get_entry(lHashval);
      [[maybe_unused]] bool b = false;  // for asserts only, else unused
      if (!(e0->content().is_full())) {
        // can insert into e0
        b = e0->content().do_insert(aTuple, lHashval);
        assert(b);
      } else {
        // cannot insert into e0
        // situation before insert: e0 -> e1 (e1 may be nullptr)
        entry_t* e1 = e0->next();
        if (nullptr == e1 || e1->content().is_full()) {
          // cannot insert into e1, must create new entry e2.
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
    } else {
      // testing only, insert should be handled by operator, not by hashtable
      if constexpr (lTrace) { std::cout << __PRETTY_FUNCTION__ << ": " << "Insert by hashtable." << std::endl; }
      _hashtable.insert(lHashval, aTuple);
    }
  }
}

template <alg_hashfun_c Thashfun, uint Tnoptr>
void
AlgOpJoinHashChAmacBuild2<Thashfun, Tnoptr>
::fin() {
  // Process the remaining non-empty entries in the AMAC ring buffer.
  // Note: No more AMAC-style processing.
  // Each buffer element is considered once and finished to completion.
  constexpr bool lTrace = false;
  if constexpr (lTrace) { std::cout << "  AlgOpJoinHashChAmacBuild2::fin()" << std::endl; }

  // check each buffer element once and finish it to completion
  for (uint _curr = 0; _curr < AMAC_BUFFER_SIZE; ++_curr) {
    amac_t& lAmac = _amac[_curr];
    if constexpr (lTrace) {
      std::cout << "    " << "_curr = " << _curr << std::endl;
      std::cout << "    " << "_todo = " << lAmac.todo() << std::endl;
    }

    if (amac_todo_et::kAmacEmpty == lAmac._todo) { continue; }

    if constexpr (lTrace) {
      std::cout << "      " << "_tuple = " << *lAmac._tuple << std::endl;
    }

    entry_t* e_dir = lAmac._ht_entry; 

    if (amac_todo_et::kAmacInsertDirEntry == lAmac._todo) {
      // stage 1: directory entry d has been prefetched, test it for insert
      if (e_dir->content().do_insert(lAmac._tuple, lAmac._hashval)) {
        // Insert was successful. Done.
        lAmac.clear();
        if constexpr (lTrace) { std::cout << "      " << "todo: 1 -> 0. done." << std::endl; }
        continue;  // next buffer entry in next for loop iteration
      } else {
        // Insert not successful (e_dir is full). Must insert into the collision chain.
        assert(e_dir->content().is_full());
        lAmac._todo = amac_todo_et::kAmacInsertChain;
        PREFETCH_BUILD(e_dir->next());
        if constexpr (lTrace) { std::cout << "      " << "todo: 1 -> 2." << std::endl; }
      }
    }
    // note: no 'else if' -> fall through from stage 1 to 2 directly
    if (amac_todo_et::kAmacInsertChain == lAmac._todo) {
      // stage 2: collision chain entry e_nxt has been prefetched in stage 1. Test if its full and insert.
      entry_t* e_nxt = e_dir->next();
      if ((nullptr == e_nxt) || !(e_nxt->content().do_insert(lAmac._tuple, lAmac._hashval))) {
        // e_nxt does not yet exist or is full; need to create new entry in the chain
        entry_t* e_new = _hashtable.new_entry();
        [[maybe_unused]] bool b = e_new->content().do_insert(lAmac._tuple, lAmac._hashval);
        assert(b);
        e_dir->_next = e_new;
        e_new->_next = e_nxt;
        if constexpr (lTrace) { std::cout << "      " << "todo: 2 -> 0 [by insert in e_new]. done." << std::endl; }
      } else {
        if constexpr (lTrace) { std::cout << "      " << "todo: 2 -> 0 [by insert in e_nxt]. done." << std::endl; }
      }
      lAmac.clear();  // done with stage 2 and with this buffer entry
    } else {
      throw std::logic_error("Unexpected AMAC stage in AlgOpJoinHashChAmacBuild2::fin(): " + std::to_string(_amac[_curr].todo()));
    }
  }
}

/*
 * Regular hash join probe operator
 * with AMAC
 *
 * Probe into the chaining hashtable provided by the respective AlgOpJoinHashChAmacBuild2 operator.
 *
 * Note:
 * AlgBase::_count counts the number of matches when probing,
 * i.e., the number of output tuples.
 *
 * Template parameters:
 * - bool IsBuildKeyUnique: If the hash table was built on a unique (key) attribute,
 *                          probe can terminate after the first match.
 */
template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun,
          bool IsBuildKeyUnique = false>
class AlgOpJoinHashChAmacProbe2 : public AlgOpBase2 {
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
    inline AlgOpJoinHashChAmacProbe2(consumer_t* aConsumer, build_t* aBuildOperator)
      : AlgOpBase2(algop_et::kAlgOpJoinChOrigProbe),
        _consumer(aConsumer), _build(aBuildOperator), _output(), _numCmps(),
        _amac(), _curr(0) {}
  public:
    static constexpr uint AMAC_BUFFER_SIZE = 12;  // size of the AMAC ring buffer
    // the different AMAC stages for this operator
    enum class amac_todo_et {
      kAmacEmpty        = 0,
      kAmacEvalHashval  = 1,
      kAmacEvalJoinpred = 2
    };
    // one entry in the AMAC ring buffer
    struct amac_t {
      const tuple_prb_t* _tuple_prb;    // probe tuple from step()
      const tuple_bld_t* _tuple_bld;    // build tuple from hash table
      const entry_t*     _ht_entry;     // build hash table entry
            hashval_t    _hashval_prb;  // hashvalue of the probe tuple
            amac_todo_et _todo;         // AMAC stage
            uint         _idx;          // index into the content of the current HT entry_t
      inline void clear() { _todo = amac_todo_et::kAmacEmpty; }
      inline bool is_empty() const { return (amac_todo_et::kAmacEmpty == _todo); }
      inline bool is_occupied() const { return (amac_todo_et::kAmacEmpty != _todo); }
      // get _todo as an integer for easier output
      inline std::underlying_type_t<amac_todo_et> todo() const { return static_cast<std::underlying_type_t<amac_todo_et>>(_todo); }
    };

  public:
    inline bool mem_alloc() { return _consumer->mem_alloc(); }
    inline bool mem_init() { return _consumer->mem_init(); }
    inline void init();
    inline void step(const input_t* aProbeTuple);
    inline void fin();
    inline void clear();
    inline void mem_free() { _consumer->mem_free(); }
  public:
    inline const consumer_t* consumer() const { return _consumer; }
    inline const AlgOpBase2* consumer_poly() const { return _consumer; }
    inline       uint64_t    numCmps()  const {
      throw "Function not implemented.";
      return _numCmps;
    }
  private:
    // increment _curr in range [0, AMAC_BUFFER_SIZE - 1]
    inline void amac_advance() { _curr = (++_curr >= AMAC_BUFFER_SIZE ? 0 : _curr); }
           void amac_print_buffer(const size_t aIndent);
  private:
    consumer_t* _consumer;  
    build_t*    _build;     // AlgOpJoinHashChAmacBuild, has getter method hashtable()
    output_t    _output;    // concat'd output tuple: <tuple_left_t*, tuple_right_t*>
    uint64_t    _numCmps;   // total number of collision chain comparisons
    amac_t      _amac[AMAC_BUFFER_SIZE];  // the AMAC ring buffer
    uint        _curr;                    // current index into the AMAC ring buffer, value range [0, AMAC_BUFFER_SIZE - 1]
};

template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun,
          bool IsBuildKeyUnique>
void
AlgOpJoinHashChAmacProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun, IsBuildKeyUnique>
::init() {
  //std::cout << "  AlgOpJoinHashChAmacProbe2::init()" << std::endl;
  reset();
  _numCmps = 0;
  for (uint i = 0; i < AMAC_BUFFER_SIZE; ++i) { _amac[i].clear(); }
  _curr = 0;
  _consumer->init();
}

template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun,
          bool IsBuildKeyUnique>
void
AlgOpJoinHashChAmacProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun, IsBuildKeyUnique>
::clear() {
  for (uint i = 0; i < AMAC_BUFFER_SIZE; ++i) { _amac[i].clear(); }
  _consumer->clear();
}

template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun,
          bool IsBuildKeyUnique>
void
AlgOpJoinHashChAmacProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun, IsBuildKeyUnique>
::step(const input_t* aProbeTuple) {
  constexpr bool lTrace = false;
  constexpr bool lUseAmac = true;
  if constexpr (lTrace) { std::cout << "  AlgOpJoinHashChAmacProbe2::step() -> probe tuple: " << *aProbeTuple << std::endl; }
  if constexpr (lTrace) { amac_print_buffer(4); }

  using hashval_bld_t = typename build_t::hashval_t;
  static_assert(std::is_same<hashval_t, hashval_bld_t>::value);

  if constexpr (lUseAmac) {
    // argument = this function call's aProbeTuple; differentiate from other hash values and pointers from the AMAC buffer
    const hashval_t lProbeHashArgument = hashfun_t::eval(aProbeTuple);
    const entry_t* lProbeHtEntryArgument = _build->hashtable().get_entry(lProbeHashArgument);

    while (_amac[_curr].is_occupied()) {
      // process tuples in the AMAC buffer until a free spot is found to insert aProbeTuple
      amac_t& lAmac = _amac[_curr];
      if constexpr (lTrace) {
        std::cout << "    " << "processing the AMAC buffer: "
          << "_curr = " << _curr << ", _todo = " << lAmac.todo()
          << ", probe tuple " << *lAmac._tuple_prb << std::endl;
      }
      if constexpr (lTrace) { amac_print_buffer(4); }
      const entry_t* e = lAmac._ht_entry;

      // AMAC stage 2
      if (amac_todo_et::kAmacEvalJoinpred == lAmac._todo) {
        if constexpr (lTrace) { std::cout << "      " << "STAGE: " << lAmac.todo() << std::endl; }
        if constexpr (lTrace) { std::cout << "      " << "checking build tuple #" << lAmac._idx << " in entry " << lAmac._ht_entry << ": " << *lAmac._tuple_bld << std::endl; }
        // stage 2: tuple lAmac._tuple_bld has been prefetched, evaluate the join pred
        if (joinpred_t::eval(lAmac._tuple_prb, lAmac._tuple_bld)) {
          //if constexpr (lTrace) { std::cout << "      " << "join match with build tuple " << *lAmac._tuple_bld << std::endl; }
          if constexpr (lTrace) { std::cout << "      " << "concat: " << (*lAmac._tuple_prb) << " + " << (*lAmac._tuple_bld) << std::endl; }
          concatfun_t::eval_prb(_output, lAmac._tuple_prb);
          concatfun_t::eval_bld(_output, lAmac._tuple_bld);
          inc();
          _consumer->step(&_output);
          if constexpr (IsBuildKeyUnique) {
            // this is the one and only matching build tuple
            lAmac.clear();
            if constexpr (lTrace) { std::cout << "      " << "todo: 2 -> 0. done." << std::endl; }
            break;  // from while loop, don't advance _curr
          }
        }
        if constexpr (lTrace) { std::cout << "      " << "checking remaining hash values in entry " << lAmac._ht_entry << std::endl;; }
        for (uint i = lAmac._idx + 1; i < e->content().end(); ++i) {
          // check the rest of the hash values in e->content()
          if (lAmac._hashval_prb == e->content().hashval(i)) {
            lAmac._tuple_bld = e->content().tuple(i);
            PREFETCH_PROBE(lAmac._tuple_bld);
            lAmac._idx = i;
            lAmac._todo = amac_todo_et::kAmacEvalJoinpred;
            if constexpr (lTrace) { std::cout << "      " << "todo: 1 -> 2." << std::endl; }
            goto PROBE_STEP_WHILE_END;  // advance _curr and go to next while iteration
          }
        }
        // if we arrive here, the for loop ran through and didn't find a matching hash value.
        // continue with next collision chain entry, if it exists.
        if (nullptr != e->next()) {
          PREFETCH_PROBE(e->next());
          lAmac._ht_entry = e->next();
          lAmac._todo = amac_todo_et::kAmacEvalHashval;
        } else {
          // if there is no next collision chain entry, we're done
          lAmac.clear();
          break;  // don't advane _curr
        }

      // AMAC stage 1
      } else if (amac_todo_et::kAmacEvalHashval == lAmac._todo) {
        if constexpr (lTrace) { std::cout << "      " << "STAGE: " << lAmac.todo() << std::endl; }
        if constexpr (lTrace) { std::cout << "        " << "entry: " << lAmac._ht_entry << std::endl; }
        // stage 1: hash table entry e has been prefetched, check its hash values
        for (uint i = e->content().begin(); i < e->content().end(); ++i) {
          // loop through e's content and check the hash values
          if (lAmac._hashval_prb == e->content().hashval(i)) {
            lAmac._tuple_bld = e->content().tuple(i);
            PREFETCH_PROBE(lAmac._tuple_bld);
            lAmac._idx = i;
            lAmac._todo = amac_todo_et::kAmacEvalJoinpred;
            if constexpr (lTrace) { std::cout << "      " << "found a hash value match. todo: 1 -> 2." << std::endl; }
            goto PROBE_STEP_WHILE_END;  // advance _curr and go to next while iteration
          }
        }  // end for e->content()...
        
        // if we arrive here, the for loop ran through and didn't find a matching hash value.
        // continue with next collision chain entry, if it exists.
        if (nullptr != e->next()) {
          assert(amac_todo_et::kAmacEvalHashval == lAmac._todo);
          PREFETCH_PROBE(e->next());
          lAmac._ht_entry = e->next();
          lAmac._todo = amac_todo_et::kAmacEvalHashval;
          if constexpr (lTrace) { std::cout << "        " << "next entry: " << lAmac._ht_entry << std::endl; }
        } else {
          // if there is no next collision chain entry, we're done
          lAmac.clear();
          if constexpr (lTrace) { std::cout << "      " << "todo: 1 -> 0. done." << std::endl; }
          break;  // don't advane _curr
        }

      } else {
        throw std::logic_error("Unexpected AMAC stage in AlgOpJoinHashChAmacProbe2::step(): " + std::to_string(_amac[_curr].todo()));
      }
      PROBE_STEP_WHILE_END:
      amac_advance();
    }
    
    // found a free spot in AMAC ring buffer, insert aProbeTuple there
    if constexpr (lTrace) {
      std::cout << "    " << "inserting probe tuple " << *aProbeTuple << " into the AMAC buffer at index " << _curr << "." << std::endl;
    }
    _amac[_curr]._tuple_prb = aProbeTuple;
    _amac[_curr]._tuple_bld = nullptr;
    _amac[_curr]._hashval_prb = lProbeHashArgument;
    _amac[_curr]._ht_entry = lProbeHtEntryArgument;
    _amac[_curr]._todo = amac_todo_et::kAmacEvalHashval;
    PREFETCH_PROBE(lProbeHtEntryArgument);
    amac_advance();

  } else {
    //std::cout << "  Warning: AlgOpJoinHashChAmacProbe2::step() does not use AMAC!" << std::endl;
    /*
     * To probe:
     * Let t := aProbeTuple.
     * 1. Compute hashvalue h from t using Thashfun.
     * 2. Get entry_t e in the hash directory using h.
     * 3. Repeat:
     *    - Compare h to all of e's hash values.
     *      If match, evaluate join predicate on t and the matching tuple s from e.
     *      If match, contact t and s using Tconcatfun and push to consumer.
     *      If IsBuildKeyUnique, terminate.
     *    - e := e->next
     *    ... while e != nullptr.
     */

    using hashval_bld_t = typename build_t::hashval_t;
    static_assert(std::is_same<hashval_t, hashval_bld_t>::value);

    const hashtable_t& lHt = _build->hashtable();

    const hashval_t lProbeHash = hashfun_t::eval(aProbeTuple);

    join_res_t lJoinRes;
    concatfun_t::eval_prb(lJoinRes, aProbeTuple);

    const entry_t* e = lHt.get_entry(lProbeHash);  // ptr calc, no mem access
    bool lMore = true;

    // walk the collision chain
    while ((nullptr != e) && lMore) {
      // check all hash values and tuple pointers in the ptr array
      for (uint i = e->content().begin(); i < e->content().end(); ++i) {
        assert(!e->content().is_empty());
        const tuple_bld_t* s = e->content().tuple(i);            // get build tuple s
        const hashval_t lBuildHash = e->content().hashval(i);  // and its hash value hb
        if ((lProbeHash == lBuildHash) && (joinpred_t::eval(aProbeTuple, s))) {  // access *s
          concatfun_t::eval_bld(lJoinRes, s);
          inc();
          _consumer->step(&lJoinRes);
          if constexpr (IsBuildKeyUnique) {
            // this is the only join match
            if constexpr (lTrace) {
              std::cout << "    break from collision chain iteration because IsBuildKeyUnique = " << IsBuildKeyUnique << std::endl;
            }
            lMore = false;
            break;
          }
        }
      }
      e = e->next();
    }
  }
}

template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun,
          bool IsBuildKeyUnique>
void
AlgOpJoinHashChAmacProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun, IsBuildKeyUnique>
::fin() {
  static constexpr bool lTrace = false;
  constexpr bool lUseAmac = true;
  if constexpr (lTrace) { std::cout << "  AlgOpJoinHashChAmacProbe2::fin()" << std::endl; }
  
  if constexpr (lUseAmac) {
    // build a bitmap that indicates which AMAC buffer elements are not empty
    uint64_t lBitmapAmacTodo = 0;  // bit i is set iff _amac[i].is_occupied()
    for (uint i = 0; i < AMAC_BUFFER_SIZE; ++i) {
      lBitmapAmacTodo |= (_amac[i].is_occupied() << i);
    }

    // loop until bitmap is empty, i.e., all AMAC buffer elements have been completed
    while (0 != lBitmapAmacTodo) {
      if constexpr (lTrace) {
        std::cout << "    " << "AMAC buffer | bitmap = " << std::bitset<AMAC_BUFFER_SIZE>(lBitmapAmacTodo) << std::endl;
        amac_print_buffer(4);
      }
      // iterate all '1's in the bitmap in ascending index order
      for (BvMemberIdxAsc lIter(lBitmapAmacTodo); lIter.isValid(); ++lIter) {
        _curr = *lIter;  // *lIter return the index of the bit set
        amac_t& lAmac = _amac[_curr];
        if constexpr (lTrace) { std::cout << "    " << "_curr = " << _curr << ", todo = " << lAmac.todo() << std::endl; }
        assert(lAmac.is_occupied());

        const entry_t* e = lAmac._ht_entry;

        // AMAC stage 2
        if (amac_todo_et::kAmacEvalJoinpred == lAmac._todo) {
          if constexpr (lTrace) { std::cout << "      " << "STAGE: " << lAmac.todo() << std::endl; }
          if constexpr (lTrace) { std::cout << "      " << "checking build tuple #" << lAmac._idx << " in entry " << lAmac._ht_entry << ": " << *lAmac._tuple_bld << std::endl; }
          // stage 2: tuple lAmac._tuple_bld has been prefetched, evaluate the join pred
          if (joinpred_t::eval(lAmac._tuple_prb, lAmac._tuple_bld)) {
            //if constexpr (lTrace) { std::cout << "      " << "join match with build tuple " << *lAmac._tuple_bld << std::endl; }
            if constexpr (lTrace) { std::cout << "      " << "concat: " << (*lAmac._tuple_prb) << " + " << (*lAmac._tuple_bld) << std::endl; }
            concatfun_t::eval_prb(_output, lAmac._tuple_prb);
            concatfun_t::eval_bld(_output, lAmac._tuple_bld);
            inc();
            _consumer->step(&_output);
            if constexpr (IsBuildKeyUnique) {
              // this is the one and only matching build tuple
              lAmac.clear();
              lBitmapAmacTodo &= ~(1 << _curr);
              if constexpr (lTrace) { std::cout << "      " << "todo: 2 -> 0. done." << std::endl; }
              break;  // from while loop, don't advance _curr
            }
          }
          if constexpr (lTrace) { std::cout << "      " << "checking remaining hash values in entry " << lAmac._ht_entry << std::endl;; }
          for (uint i = lAmac._idx + 1; i < e->content().end(); ++i) {
            // check the rest of the hash values in e->content()
            if (lAmac._hashval_prb == e->content().hashval(i)) {
              lAmac._tuple_bld = e->content().tuple(i);
              PREFETCH_PROBE(lAmac._tuple_bld);
              lAmac._idx = i;
              lAmac._todo = amac_todo_et::kAmacEvalJoinpred;
              if constexpr (lTrace) { std::cout << "      " << "todo: 2 -> 2." << std::endl; }
              goto PROBE_FIN_WHILE_END;  // advance _curr and go to next while iteration
            }
          }
          // if we arrive here, the for loop ran through and didn't find a matching hash value.
          // continue with next collision chain entry, if it exists.
          if (nullptr != e->next()) {
            PREFETCH_PROBE(e->next());
            lAmac._ht_entry = e->next();
            lAmac._todo = amac_todo_et::kAmacEvalHashval;
          } else {
            // if there is no next collision chain entry, we're done
            lAmac.clear();
            lBitmapAmacTodo &= ~(1 << _curr);
            break;  // don't advane _curr
          }

        // AMAC stage 1
        } else if (amac_todo_et::kAmacEvalHashval == lAmac._todo) {
          if constexpr (lTrace) { std::cout << "      " << "STAGE: " << lAmac.todo() << std::endl; }
          if constexpr (lTrace) { std::cout << "        " << "entry: " << lAmac._ht_entry << std::endl; }
          // stage 1: hash table entry e has been prefetched, check its hash values
          for (uint i = e->content().begin(); i < e->content().end(); ++i) {
            // loop through e's content and check the hash values
            if (lAmac._hashval_prb == e->content().hashval(i)) {
              lAmac._tuple_bld = e->content().tuple(i);
              PREFETCH_PROBE(lAmac._tuple_bld);
              lAmac._idx = i;
              lAmac._todo = amac_todo_et::kAmacEvalJoinpred;
              if constexpr (lTrace) { std::cout << "      " << "found a hash value match. todo: 1 -> 2." << std::endl; }
              goto PROBE_FIN_WHILE_END;  // advance _curr and go to next while iteration
            }
          }  // end for e->content()...
          
          // if we arrive here, the for loop ran through.
          // continue with next collision chain entry, if it exists.
          if (nullptr != e->next()) {
            assert(amac_todo_et::kAmacEvalHashval == lAmac._todo);
            PREFETCH_PROBE(e->next());
            lAmac._ht_entry = e->next();
            lAmac._todo = amac_todo_et::kAmacEvalHashval;
            if constexpr (lTrace) { std::cout << "        " << "next entry: " << lAmac._ht_entry << std::endl; }
          } else {
            // if there is no next collision chain entry, we're done
            lAmac.clear();
            lBitmapAmacTodo &= ~(1 << _curr);
            if constexpr (lTrace) { std::cout << "      " << "todo: 1 -> 0. done." << std::endl; }
            break;  // don't advane _curr
          }

        } else {
          throw std::logic_error("Unexpected AMAC stage in AlgOpJoinHashChAmacProbe2::step(): " + std::to_string(_amac[_curr].todo()));
        }
      }  // end for every bitmap member
      PROBE_FIN_WHILE_END: ;
    }  // end while bitmap != 0

  } else {
    // no AMAC/prefetching in fin: process the remaining AMAC entries to completion, i.e., consider each AMAC entry once
    for (_curr = 0; _curr < AMAC_BUFFER_SIZE; ++_curr) {
      amac_t& lAmac = _amac[_curr];
      if constexpr (lTrace) { std::cout << "    " << "_curr = " << _curr << ", _todo = " << lAmac.todo() << std::endl; }

      if (lAmac.is_empty()) { continue; }  // nothing to do

      bool lMore = true;  // indicate if there can be more matches in subsequent collision chain nodes

      // each AMAC entry is at least in stage 1
      const entry_t* e = lAmac._ht_entry;
      const tuple_prb_t* lProbeTuple = lAmac._tuple_prb;
      const hashval_t lProbeHash = lAmac._hashval_prb;
      concatfun_t::eval_prb(_output, lProbeTuple);  // prepare output tuple
      if constexpr (lTrace) { std::cout << "      checking entry " << e << std::endl; }

      // if we're already in stage 2 (kAmacEvalJoinpred),
      // then we might have already looked at some of the hash values and tuples in e,
      // so we need to continue where we left of (lAmac._idx).
      const uint lIdxBegin = (amac_todo_et::kAmacEvalJoinpred == lAmac._todo) ? lAmac._idx
                                                                              : e->content().begin();
      for (uint i = lIdxBegin; i < e->content().end(); ++i) {
        const tuple_bld_t* lBuildTuple = e->content().tuple(i);  // potential join partner
        const hashval_t lBuildHash = e->content().hashval(i);
        // no more processing in AMAC stages, so just test everything in one go
        if ((lProbeHash == lBuildHash) && joinpred_t::eval(lProbeTuple, lBuildTuple)) {
          // join partner found!
          concatfun_t::eval_bld(_output, lBuildTuple);
          if constexpr (lTrace) { std::cout << "      " << "join partner found: " << (*lProbeTuple) << " + " << *(lBuildTuple) << std::endl; }
          inc();
          _consumer->step(&_output);
          if constexpr (IsBuildKeyUnique) {
            // this is the one and only match
            lAmac.clear();
            lMore = false;
            break;
          }
        }
      }  // end for e->content()

      // there might be more matches in e->next(), e->next()->next()...
      if (lMore) {
        e = e->next();
        while ((nullptr != e) && lMore) {
          if constexpr (lTrace) { std::cout << "      checking entry " << e << std::endl; }
          for (uint i = e->content().begin(); i < e->content().end(); ++i) {
            const tuple_bld_t* lBuildTuple = e->content().tuple(i);
            const hashval_t lBuildHash = e->content().hashval(i);
            if ((lProbeHash == lBuildHash) && joinpred_t::eval(lProbeTuple, lBuildTuple)) {
              // join partner found!
              concatfun_t::eval_bld(_output, lBuildTuple);
              if constexpr (lTrace) { std::cout << "      " << "join partner found: " << (*lProbeTuple) << " + " << *(lBuildTuple) << std::endl; }
              inc();
              _consumer->step(&_output);
              if constexpr (IsBuildKeyUnique) {
                // this is the one and only match
                lAmac.clear();
                lMore = false;
                break;
              }
            }
          }
          e = e->next();
        }  // end while e
      }  // end if more
      lAmac.clear();  // here, we're done with this buffer entry for sure
    }  // end loop over AMAC buffer
  }  // end else lUseAmac

  if constexpr (lTrace) { amac_print_buffer(4); }
#ifndef NDEBUG
  // make sure that AMAC buffer is empty
  for (size_t i = 0; i < AMAC_BUFFER_SIZE; ++i) {
    assert(_amac[i].is_empty() && "AMAC ring buffer has non-empty elements at the end of fin.");
  }
#endif

  _consumer->fin();
}

template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun,
          bool IsBuildKeyUnique>
void
AlgOpJoinHashChAmacProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun, IsBuildKeyUnique>
::amac_print_buffer(const size_t aIndent) {
  const std::string lIndent(aIndent, ' ');
  std::cout << lIndent << std::string(AMAC_BUFFER_SIZE * 3 + 1, '-') << std::endl;
  std::cout << lIndent;
  for (uint i = 0; i < AMAC_BUFFER_SIZE; ++i) {
    std::cout << "|" << std::setw(2) << i;
    if (i == AMAC_BUFFER_SIZE - 1) {
      std::cout << "|" << std::endl;
    }
  }
  std::cout << lIndent;
  for (uint i = 0; i < AMAC_BUFFER_SIZE; ++i) {
    std::cout << "|" << std::setw(2) << _amac[i].todo();
    if (i == AMAC_BUFFER_SIZE - 1) {
      std::cout << "|" << std::endl;
    }
  }
  std::cout << lIndent << std::string(AMAC_BUFFER_SIZE * 3 + 1, '-') << std::endl;
}

