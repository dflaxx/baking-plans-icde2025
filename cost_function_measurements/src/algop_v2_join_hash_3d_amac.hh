#pragma once

#include "algop_v2_join_hash_ch_amac.hh"
#include "dfinfra/standard_includes.hh"
#include "dfinfra/exceptions.hh"

#include "algop_v2_base.hh"
#include "hashtable_nested_v2.hh"
#include "hj_htnode_content.hh"
#include "concepts.hh"

/*
 * 3D hash join build operator using the nested hash table v2
 * with AMAC
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
class AlgOpJoinHash3dAmacBuild2 : public AlgOpBase2 {
  public:
    using hashfun_t         = Thashfun;
    using hashval_t         = typename hashfun_t::output_t;
    using input_t           = typename hashfun_t::input_t;
    using output_t          = void;
    using tuple_t           = input_t;
    using consumer_t        = void;
    using eqfun_t           = Tequalfun;
    using ht_content_sub_t  = ht_content_3d_sub_tt<input_t, Tnoptrsub>;
    using ht_content_main_t = ht_content_3d_main_tt<hashval_t, ht_content_sub_t, Tnoptrmain>;
    using hashtable_t       = HashtableNested2<ht_content_main_t, ht_content_sub_t>;
    using entry_main_t      = typename hashtable_t::entry_main_t;
    using entry_sub_t       = typename hashtable_t::entry_sub_t;
  public:
    AlgOpJoinHash3dAmacBuild2(const size_t aHashDirSize,
                     const uint32_t aHtLog2ChunkSizeMain,
                     const uint32_t aHtLog2ChunkSizeSub)
      : AlgOpBase2(algop_et::kAlgOpJoin3dAmacBuild),
        _hashtable(aHashDirSize, aHtLog2ChunkSizeMain, aHtLog2ChunkSizeSub),
        _amac(), _curr(0) {}
    AlgOpJoinHash3dAmacBuild2(const size_t aHashDirSize, const uint32_t aHtLog2ChunkSize)
      : AlgOpJoinHash3dAmacBuild2(aHashDirSize, aHtLog2ChunkSize, aHtLog2ChunkSize) {}
  public:
    /* Here comes the AMAC stuff */
    static constexpr uint AMAC_BUFFER_SIZE = 12;  // size of the AMAC ring buffer
    // the different AMAC stages for 3D build
    enum class amac_todo_et {
      kAmacEmpty          = 0,
      kAmacInsertDirEntry = 1,  // access the (already prefetched) directory entry and check if it's empty              // if empty, insert
      kAmacFindMainByHash = 2,  // access a main entry's hashvalue (after prefetching the node) and compare it          // find main node
      kAmacVerifyMainByEq = 3,  // access a main entry's tuple (after prefetching it) and eval the equality predicate   // try insert main node
      kAmacInsertSub      = 4   // access matching (but full) main entry's sub node (after prefetching it)              // try insert sub node
    };
    // one entry in the AMAC ring buffer
    struct amac_t {
      const tuple_t* _tuple;               // build tuple to be inserted
            entry_main_t* _ht_entry_main;  // hash table main entry
            hashval_t     _hashval;        // hash value of build tuple
            amac_todo_et  _todo;           // AMAC stage
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
           void step(const input_t* aInput);
           void fin();
    inline void clear() { _hashtable.clear(); }
    inline void mem_free() { _hashtable.mem_free(); }
  public:
    inline const hashtable_t& hashtable() const { return _hashtable; }
    inline const consumer_t* consumer() const { return nullptr; }
    inline const AlgOpBase2* consumer_poly() const { return nullptr; }
  public:
    inline size_t mem_consumption_ht() const { return _hashtable.mem_consumption(); }
    inline size_t ht_num_main_nodes() const { return _hashtable.rsv_main_card(); }  // main collision chain nodes in rsv
    inline size_t ht_num_sub_nodes() const { return _hashtable.rsv_sub_card(); }
  private:
    // increment _curr in range [0, AMAC_BUFFER_SIZE - 1]
    inline void amac_advance() { _curr = (++_curr >= AMAC_BUFFER_SIZE ? 0 : _curr); }
           void amac_print_buffer(const size_t aIndent);
  private:
    hashtable_t _hashtable;
    amac_t      _amac[AMAC_BUFFER_SIZE];  // the AMAC ring buffer
    uint        _curr;                    // current index into the AMAC ring buffer, value range [0, AMAC_BUFFER_SIZE - 1]
};

template <alg_hashfun_c Thashfun, alg_binary_predicate_c Tequalfun, uint Tnoptrmain, uint Tnoptrsub>
void
AlgOpJoinHash3dAmacBuild2<Thashfun, Tequalfun, Tnoptrmain, Tnoptrsub>
::init() {
  reset();
  for (uint i = 0; i < AMAC_BUFFER_SIZE; ++i) { _amac[i].clear(); }
  _curr = 0;
}

template <alg_hashfun_c Thashfun, alg_binary_predicate_c Tequalfun, uint Tnoptrmain, uint Tnoptrsub>
void
AlgOpJoinHash3dAmacBuild2<Thashfun, Tequalfun, Tnoptrmain, Tnoptrsub>
::step(const input_t* aTuple) {
  constexpr bool lTrace = false;
  constexpr bool lUseAmac = true;

  if constexpr (lTrace) { std::cout << "  AlgOpJoinHash3dAmacBuild2::step() -> build tuple #" << count() << std::endl; }
  inc();

  if constexpr (lUseAmac) {
    const hashval_t lHashvalArgument = hashfun_t::eval(aTuple);   // access *aTuple
    entry_main_t* lMainEntryArgument = _hashtable.get_entry(lHashvalArgument);  // no memory access, only pointer computation

    if constexpr (lTrace) {
      std::cout << "    " << "aTuple = " << *aTuple << std::endl;
      std::cout << "    " << "_curr = " << _curr << std::endl;
      std::cout << "    " << "_amac[_curr].is_empty() = " << _amac[_curr].is_empty() << std::endl;
    }
    
    // process tuples in the AMAC ring buffer until a free spot if found for the new tuple aTuple
    while (_amac[_curr].is_occupied()) {
      amac_t& lAmac = _amac[_curr];
      entry_main_t* e = lAmac._ht_entry_main;
      if constexpr (lTrace) {
        std::cout << "    " << "processing the AMAC buffer: "
          << "_curr = " << _curr << ", _todo = " << lAmac.todo() << ", tuple " << *lAmac._tuple << std::endl;
      }

      /*
       * Process the possible AMAC stages:
       * - kAmacInsertDirEntry = 1,  // access the (already prefetched) directory entry and check if it's empty
       * - kAmacFindMainByHash = 2,  // access a main entry's hashvalue (after prefetching the node) and compare it
       * - kAmacVerifyMainByEq = 3,  // access a main entry's tuple (after prefetching it) and eval the equality predicate
       * - kAmacInsertSub      = 4   // access matching (but full) main entry's sub node (after prefetching it)
       */

      [[maybe_unused]] bool b = false;  // result of insert operations

      // stage 1: directory entry has been prefetched, test it for insert
      if (amac_todo_et::kAmacInsertDirEntry == lAmac._todo) {
        if (e->content().is_empty()) {
          // nothing in directory entry yet, just insert here
          assert(0 == e->content().size());
          b = e->content().do_insert(lAmac._tuple);
          assert(b);
          e->content().hashval(lAmac._hashval);
          lAmac.clear();
          if constexpr (lTrace) { std::cout << "      " << "todo: 1 -> 0. done." << std::endl; }
          break;  // from the while loop. Don't advance _curr, so next tuple is inserted in this free spot.
        } else {
          // directory entry is not empty, must find matching main node for insert
          lAmac._todo = amac_todo_et::kAmacFindMainByHash;
          if constexpr (lTrace) { std::cout << "      " << "todo: 1 -> 2." << std::endl; }
        }

      // stage 2: current main entry has already been prefetched, compare its hashvalue
      } else if (amac_todo_et::kAmacFindMainByHash == lAmac._todo) {
        if (e->content().hashval() == lAmac._hashval) {
          // hashval matches, next check equality predicate on tuple
          PREFETCH_BUILD(e->content().get_first_tuple());
          lAmac._todo = amac_todo_et::kAmacVerifyMainByEq;
          if constexpr (lTrace) { std::cout << "      " << "todo: 2 -> 3." << std::endl; }
        } else {
          // hashval does not match, check if there are more main entries
          if (nullptr != e->next()) {
            // more main entries. prefetch next.
            lAmac._ht_entry_main = e->next();
            PREFETCH_BUILD(e->next());
            if constexpr (lTrace) { std::cout << "      " << "todo: 2 -> 2." << std::endl; }
          } else {
            // no more main entries but none matched so far -> create a new one and insert.
            entry_main_t* e_new = _hashtable.new_main_entry();
            assert(e_new->content().is_empty());
            b = e_new->content().do_insert(lAmac._tuple);
            e_new->content().hashval(lAmac._hashval);
            assert(b);
            e->_next = e_new;
            lAmac.clear();  // done with this entry
            if constexpr (lTrace) { std::cout << "      " << "todo: 2 -> 0. done." << std::endl; }
            break;  // from the while loop. Don't advance _curr, so next tuple is inserted in this free spot.
          }
        }

      // stage 3: tuple of current main entry has already been prefetched, evaluate equality predicate now
      } else if (amac_todo_et::kAmacVerifyMainByEq == lAmac._todo) {
        assert(!e->content().is_empty());
        assert(e->content().hashval() == lAmac._hashval);
        if (eqfun_t::eval(e->content().get_first_tuple(), lAmac._tuple)) {  // evaluate equality predicate
          // this main entry matches -> insert here
          if (e->content().do_insert(lAmac._tuple)) {
            // insert succeeded, done.
            lAmac.clear();
            if constexpr (lTrace) { std::cout << "      " << "todo: 3 -> 0. done." << std::endl; }
            break;  // from the while loop. Don't advance _curr, so next tuple is inserted in this free spot.
          } else {
            // insert did not succeed because main entry is full, check first sub node
            if (nullptr != e->sub()) {
              // sub node exists, prefetch it
              PREFETCH_BUILD(e->sub());
              lAmac._todo = amac_todo_et::kAmacInsertSub;
              if constexpr (lTrace) { std::cout << "      " << "todo: 3 -> 4. done." << std::endl; }
            } else {
              // no sub nodes exists so far; create it and insert
              entry_sub_t* s_new = _hashtable.new_sub_entry();
              assert(s_new->content().is_empty());
              b = s_new->content().do_insert(lAmac._tuple);
              assert(b);
              assert(nullptr == e->sub());
              e->_sub = s_new;
              lAmac.clear();
              if constexpr (lTrace) { std::cout << "      " << "todo: 3 -> 0. done." << std::endl; }
              break;  // from the while loop. Don't advance _curr, so next tuple is inserted in this free spot.
            }
          }
        } else {
          // this main entry does not match, check next
          if (nullptr != e->next()) {
            // more main entries. prefetch next.
            lAmac._ht_entry_main = e->next();
            PREFETCH_BUILD(e->next());
            lAmac._todo = amac_todo_et::kAmacFindMainByHash;
            if constexpr (lTrace) { std::cout << "      " << "todo: 3 -> 2." << std::endl; }
          } else {
            // no more main entries but none matched so far. create a new one and insert.
            entry_main_t* e_new = _hashtable.new_main_entry();
            assert(e_new->content().is_empty());
            b = e_new->content().do_insert(lAmac._tuple);
            e_new->content().hashval(lAmac._hashval);
            assert(b);
            e->_next = e_new;
            lAmac.clear();  // done with this entry
            if constexpr (lTrace) { std::cout << "      " << "todo: 3 -> 0. done." << std::endl; }
            break;  // from the while loop. Don't advance _curr, so next tuple is inserted in this free spot.
          }
        }

      // stage 4: sub entry of current entry has already been prefetched. try to insert here.
      } else if (amac_todo_et::kAmacInsertSub == lAmac._todo) {
        // situation before: main entry e -(sub)-> s -> s' -> ...
        entry_sub_t* s = e->sub();  // assume (hope) that e is still in the cache
        assert(!e->content().is_empty());
        assert(e->content().hashval() == lAmac._hashval);
        assert(nullptr != e->sub());
        if (!(s->content().do_insert(lAmac._tuple))) {
          // insert did not succeed, need new sub entry
          if constexpr (lTrace) { std::cout << "        " << "need new sub entry" << std::endl; }
          assert(s->content().is_full());
          entry_sub_t* s_new = _hashtable.new_sub_entry();
          assert(s_new->content().is_empty());
          b = s_new->content().do_insert(lAmac._tuple);
          assert(b);
          s_new->_next = s;
          e->_sub = s_new;
          // situation now: e -> s_new -> s -> s' -> ...
        }
        // insert into (potentially new) sub entry succeeded, done.
        lAmac.clear();
        if constexpr (lTrace) { std::cout << "      " << "todo: 4 -> 0. done." << std::endl; }
        break;  // from the while loop. Don't advance _curr, so next tuple is inserted in this free spot.

      } else {
        throw std::logic_error("Unexpected AMAC stage in AlgOpJoinHash3dAmacBuild2::step(): " + std::to_string(_amac[_curr].todo()));
      }

      amac_advance();
    }  // end while: amac occupied

    // if we arrive here, we found a free spot in the AMAC ring buffer.
    // --> insert aTuple into the free spot in the ring buffer and issue a prefetch
    if constexpr (lTrace) {
      std::cout << "    " << "inserting tuple " << *aTuple << " into the AMAC buffer. todo: 0 -> 1." << std::endl;
    }
    assert(_amac[_curr].is_empty());
    _amac[_curr]._tuple = aTuple;
    _amac[_curr]._ht_entry_main = lMainEntryArgument;
    _amac[_curr]._hashval = lHashvalArgument;
    _amac[_curr]._todo = amac_todo_et::kAmacInsertDirEntry;
    PREFETCH_BUILD(lMainEntryArgument);
    amac_advance();

  } else {
    const hashval_t lHashval = hashfun_t::eval(aTuple);
    constexpr bool lInsertByOperator = true;

    if constexpr (lInsertByOperator) {
      if constexpr (lTrace) { std::cout << __PRETTY_FUNCTION__ << ": " << "Insert by operator." << std::endl; }
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

      using entry_main_t = typename hashtable_t::entry_main_t;
      using entry_sub_t = typename hashtable_t::entry_sub_t;
      entry_main_t* d = _hashtable.get_entry(lHashval);
      [[maybe_unused]] bool b = false;

      if (d->content().is_empty()) {
        // directory entry is empty, can just insert here
        assert(0 == d->content().size());
        b = d->content().do_insert(aTuple);
        d->content().hashval(lHashval);
        assert(b);
        assert(1 == d->content().size());
      } else {
        // walk the main collision chain
        entry_main_t* e = d;
        while ((nullptr != e) &&
               (e->content().hashval() != lHashval) &&                        // compare hash values
               (!(eqfun_t::eval(e->content().get_first_tuple(), aTuple)))) {  // evaluate join predicate
          e = e->next();
        }
        if (nullptr == e) {
          // didn't find a main node with matching hashvalue and join key,
          // need new main node.
          // situation before: d -> e0 -> e1 -> ...
          // situation after:  d -> e_new -> e0 -> e1 -> ...
          entry_main_t* e_new = _hashtable.new_main_entry();
          assert(0 == e_new->content().size());
          b = e_new->content().do_insert(aTuple);
          e_new->content().hashval(lHashval);
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
              s = _hashtable.new_sub_entry();
              assert(s->content().is_empty());
              b = s->content().do_insert(aTuple);
              assert(b);
              assert(1 == s->content().size());
              s->next(e->sub());
              e->sub(s);
              assert(s->content().size() == 1);
            } else {
              // s exists and is not yet full. can safely insert.
              assert(!s->content().is_full());
              b = s->content().do_insert(aTuple);
              assert(b);
              assert(s->content().size() >= 2);
            }
          } else {
            // matching main node e is not yet full, so insert there
            b = e->content().do_insert(aTuple);
            assert(b);
            assert(!e->content().is_empty());
            assert(lHashval == e->content().hashval());
          }
        }
      }
    } else {  // lInsertByOperator == false
      // XXX testing only, insert should be handled by operator, not hashtable
      if constexpr (lTrace) { std::cout << __PRETTY_FUNCTION__ << ": " << "Insert by hashtable." << std::endl; }
      _hashtable.template insert<eqfun_t>(lHashval, aTuple);
    }
  }  // end else lUseAmac
}

template <alg_hashfun_c Thashfun, alg_binary_predicate_c Tequalfun, uint Tnoptrmain, uint Tnoptrsub>
void
AlgOpJoinHash3dAmacBuild2<Thashfun, Tequalfun, Tnoptrmain, Tnoptrsub>
::fin() {
  // Process the remaining non-empty entries in the AMAC ring buffer.
  // Note: No more AMAC-style round-robin processing. Consider each buffer element once, finish it to completion.

  constexpr bool lTrace = false;
  constexpr bool lAmacInStep = true;  // AMAC has been used in step() and we must finish the buffer
  if constexpr (lTrace) { std::cout << "  AlgOpJoinHash3dAmacBuild2::fin()" << std::endl; }

  if constexpr (lAmacInStep) {
    // check each buffer element once and finish it to completion. no prefetching any more.
    for (_curr = 0; _curr < AMAC_BUFFER_SIZE; ++_curr) {
      if constexpr (lTrace) { amac_print_buffer(4); }
      amac_t& lAmac = _amac[_curr];
      if constexpr (lTrace) {
        std::cout << "    " << "_curr = " << _curr << std::endl;
        std::cout << "    " << " todo = " << lAmac.todo() << std::endl;
      }

      if (amac_todo_et::kAmacEmpty == lAmac._todo) { continue; }  // nothing to do for this entry

      assert(nullptr != lAmac._tuple);
      if constexpr (lTrace) {
        std::cout << "    " << " tuple = " << *(lAmac._tuple) << std::endl;
      }

      entry_main_t* e = lAmac._ht_entry_main;
      [[maybe_unused]] bool b = false;

      // stage 1
      if (amac_todo_et::kAmacInsertDirEntry == lAmac._todo) {
        if (e->content().is_empty()) {
          assert(0 == e->content().size());
          b = e->content().do_insert(lAmac._tuple);
          e->content().hashval(lAmac._hashval);
          assert(b);
          lAmac.clear();
          if constexpr (lTrace) { std::cout << "      " << "todo: 1 -> 0. done." << std::endl; }
          continue;
        } else {
          lAmac._todo = amac_todo_et::kAmacFindMainByHash;
          if constexpr (lTrace) { std::cout << "      " << "todo: 1 -> 2." << std::endl; }
        }
      }

      // stages 2 & 3
      if ((amac_todo_et::kAmacFindMainByHash == lAmac._todo) || amac_todo_et::kAmacVerifyMainByEq == lAmac._todo) {
        while (nullptr != e) {
          if ((e->content().hashval() == lAmac._hashval) &&
              (eqfun_t::eval(e->content().get_first_tuple(), lAmac._tuple))) {
            // found the correct (matching) main entry
            if (e->content().do_insert(lAmac._tuple)) {
              // done
              if constexpr (lTrace) { std::cout << "      " << "todo: 2/3 -> 0. done." << std::endl; }
              lAmac.clear();
              break;
            } else {
              // insert not successful due to no space in e
              lAmac._ht_entry_main = e;
              lAmac._todo = amac_todo_et::kAmacInsertSub;
              if constexpr (lTrace) { std::cout << "      " << "todo: 2/3 -> 4." << std::endl; }
              break;
            }
          }
          lAmac._ht_entry_main = e;
          e = e->next();
        }
        if (nullptr == e) {
          // we're at the end of the main chain, need new main node
          entry_main_t* e_new = _hashtable.new_main_entry();
          assert(e_new->content().is_empty());
          b = e_new->content().do_insert(lAmac._tuple);
          assert(b);
          e_new->content().hashval(lAmac._hashval);
          lAmac._ht_entry_main->_next = e_new;
          lAmac.clear();
          if constexpr (lTrace) { std::cout << "      " << "todo: 2/3 -> 0. done." << std::endl; }
          continue;
        }
      }

      // stage 4
      if (amac_todo_et::kAmacInsertSub == lAmac._todo) {
        e = lAmac._ht_entry_main;
        assert ((e->content().hashval() == lAmac._hashval) && (eqfun_t::eval(e->content().get_first_tuple(), lAmac._tuple)));
        entry_sub_t* s = e->sub();
        if ((nullptr != s) && (s->content().do_insert(lAmac._tuple))) {
          // done
          if constexpr (lTrace) { std::cout << "      " << "insert into existing sub entry" << std::endl; }
          if constexpr (lTrace) { std::cout << "      " << "todo: 4 -> 0. done." << std::endl; }
        } else {
          // need new sub node
          entry_sub_t* s_new = _hashtable.new_sub_entry();
          assert(s_new->content().is_empty());
          b = s_new->content().do_insert(lAmac._tuple);
          assert(b);
          s_new->_next = s;
          e->_sub = s_new;
          if constexpr (lTrace) { std::cout << "      " << "insert into newly created sub entry" << std::endl; }
          if constexpr (lTrace) { std::cout << "      " << "todo: 4 -> 0. done." << std::endl; }
        }
        lAmac.clear();
        continue;
      }
    }  // end for each AMAC buffer element
    if constexpr (lTrace) { amac_print_buffer(4); }
#ifndef NDEBUG
    // make sure that all AMAC buffer is empty
    for (size_t i = 0; i < AMAC_BUFFER_SIZE; ++i) {
      assert(_amac[i].is_empty() && "AMAC ring buffer has non-empty elements at the end of fin.");
    }
#endif
  }
  // for lAmacInStep == false, there is nothing to do, so fin() is empty.
}

template <alg_hashfun_c Thashfun, alg_binary_predicate_c Tequalfun, uint Tnoptrmain, uint Tnoptrsub>
void
AlgOpJoinHash3dAmacBuild2<Thashfun, Tequalfun, Tnoptrmain, Tnoptrsub>
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


/*
 * 3D hash join probe operator with unnest using the nested hash table v2
 * with AMAC
 *
 * Probe into the nested hashtable provided by the respective AlgOpJoinHash3dAmacBuild2 operator.
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
class AlgOpJoinHash3dAmacProbe2 : public AlgOpBase2 {
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
    AlgOpJoinHash3dAmacProbe2(consumer_t* aConsumer, build_t* aBuildOperator)
      : AlgOpBase2(algop_et::kAlgOpJoin3dAmacProbeWithUnnest),
        _consumer(aConsumer), _build(aBuildOperator), _output(),
        _amac(), _curr(0),
        _numCmps1(0), _numCmps2(0), _c_mid(0) {}
  public:
    /* Here comes the AMAC stuff */
    static constexpr uint AMAC_BUFFER_SIZE = 12;  // size of the AMAC ring buffer
    // the different AMAC stages for 3D probe
    enum class amac_todo_et {
      kAmacEmpty          = 0,
      kAmacCheckDirectory = 1,
      kAmacFindMainByHash = 2,
      kAmacVerifyMainByEq = 3,
      kAmacUnnestSub      = 4,
    };
    // one entry in the AMAC ring buffer
    struct amac_t {
      const tuple_prb_t*   _tuple_prb;     // pointer to probe tuple from step()
      const entry_main_t*  _ht_entry_main; // pointer to build hash table main entry
      const entry_sub_t*   _ht_entry_sub;  // pointer to build hash table sub entry
            hashval_t      _hashval_prb;   // hashvalue of the probe tuple
            amac_todo_et   _todo;          // AMAC stage
      inline void clear() { _todo = amac_todo_et::kAmacEmpty; }
      inline bool is_empty() const { return (_todo == amac_todo_et::kAmacEmpty); }
      inline bool is_occupied() const { return (_todo != amac_todo_et::kAmacEmpty); }
      // get _todo as an integer for easier output
      inline std::underlying_type_t<amac_todo_et> todo() const { return static_cast<std::underlying_type_t<amac_todo_et>>(_todo); }
    };
  public:
    inline bool mem_alloc() { return _consumer->mem_alloc(); }
    inline bool mem_init() { return _consumer->mem_init(); }
    inline void init();
           void step(const input_t* aProbeTuple);
           void fin();
    inline void clear();
    inline void mem_free() { _consumer->mem_free(); }
  public:
    inline const consumer_t* consumer() const { return _consumer; }
    inline const AlgOpBase2* consumer_poly() const { return _consumer; }
    inline       uint64_t    numCmps1()  const { return _numCmps1; }
    inline       uint64_t    numCmps2()  const { return _numCmps2; }
  private:
    // increment _curr in range [0, AMAC_BUFFER_SIZE - 1]
    inline void amac_advance() { _curr = (++_curr >= AMAC_BUFFER_SIZE ? 0 : _curr); }
           void amac_print_buffer(const size_t aIndent);
  private:
    consumer_t* _consumer;
    build_t*    _build;                   // AlgOpJoinHash3dAmacBuild2, has getter method hashtable()
    output_t    _output;                  // unnested output tuple (_consumer->step(&_output))
    amac_t      _amac[AMAC_BUFFER_SIZE];  // the AMAC ring buffer
    uint        _curr;                    // current index into the AMAC ring buffer, value range [0, AMAC_BUFFER_SIZE - 1]
    uint64_t    _numCmps1;                // total number of collision chain *hashval* comparisons (requires entry_main_t access)
    uint64_t    _numCmps2;                // total number of collision chain *joinpred* evaluations (requires tuple access)
    uint64_t    _c_mid;                   // number of *nested* result tuples before unnesting
                                          // (= number of qualifying probe tuples = |Prb lsjoin Bld|).
};

template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun>
void
AlgOpJoinHash3dAmacProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun>
::init() {
  reset();
  _numCmps1 = 0;
  _numCmps2 = 0;
  for (uint i = 0; i < AMAC_BUFFER_SIZE; ++i) { _amac[i].clear(); }
  _curr = 0;
  _consumer->init();
}

template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun>
void
AlgOpJoinHash3dAmacProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun>
::clear() {
  for (uint i = 0; i < AMAC_BUFFER_SIZE; ++i) { _amac[i].clear(); }
  _consumer->clear();
}

template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun>
void
AlgOpJoinHash3dAmacProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun>
::step(const input_t* aProbeTuple) {
  constexpr bool lTrace = false;
  constexpr bool lUseAmac = true;

  if constexpr (lTrace) { std::cout << "  AlgOpJoinHash3dAmacProbe2::step() -> probe tuple: " << *aProbeTuple << std::endl; }

  using hashval_bld_t = typename build_t::hashval_t;
  static_assert(std::is_same<hashval_t, hashval_bld_t>::value);

  if constexpr (lUseAmac) {
    const hashval_t lHashvalArgument = hashfun_t::eval(aProbeTuple);
    const entry_main_t* lMainEntryArgument = _build->hashtable().get_entry(lHashvalArgument);

    // process tuples in the AMAC ring buffer until a free spot is found to insert aProbeTuple
    while (_amac[_curr].is_occupied()) {
      amac_t& lAmac = _amac[_curr];
      assert(!lAmac.is_empty());

      if constexpr (lTrace) {
        std::cout << "    " << "processing the AMAC buffer: "
          << "_curr = " << _curr << ", _todo = " << lAmac.todo()
          << ", probe tuple " << *lAmac._tuple_prb << std::endl;
      }

      const entry_main_t* e = lAmac._ht_entry_main;

      /* Process the different AMAC stages:
       * - kAmacCheckDirectory = 1,
       * - kAmacFindMainByHash = 2,
       * - kAmacVerifyMainByEq = 3,
       * - kAmacUnnestSub      = 4,
       */

      // stage 1: directory entry has been prefetched, check if it's empty
      if (amac_todo_et::kAmacCheckDirectory == lAmac._todo) {
        if (e->content().is_empty()) {
          if constexpr (lTrace) { std::cout << "    " << "directory entry is empty. 1 -> 0. done." << std::endl; }
          lAmac.clear();
          break;
        } else {
          if constexpr (lTrace) { std::cout << "    " << "1 -> 2." << std::endl; }
          lAmac._todo = amac_todo_et::kAmacFindMainByHash;
        }

      // stage 2: main node has been prefetched, check hash value equality
      } else if (amac_todo_et::kAmacFindMainByHash == lAmac._todo) {
        assert((nullptr != e) && (!e->content().is_empty()));
        ++_numCmps1;
        if (e->content().hashval() == lAmac._hashval_prb) {
          // found a potentially matching main node, check join predicate next
          if constexpr (lTrace) { std::cout << "    " << "hash values match, check join predicate next. 2 -> 3." << std::endl; }
          PREFETCH_PROBE(e->content().get_first_tuple());
          lAmac._todo = amac_todo_et::kAmacVerifyMainByEq;
        } else {
          // hash value does not match, check the next node
          if constexpr (lTrace) { std::cout << "    " << "hash values differ." << std::endl; }
          if (nullptr == e->next()) {
            // no more main nodes to examine, no more matches, done
            if constexpr (lTrace) { std::cout << "    " << "no more main nodes to examine. 2 -> 0. done." << std::endl; }
            lAmac.clear();
            break;
          } else {
            // more main nodes to examine
            PREFETCH_PROBE(e->next());
            lAmac._ht_entry_main = e->next();
            if constexpr (lTrace) { std::cout << "    " << "more main nodes to examine. 2 -> 2." << std::endl; }
          }
        }

      // stage 3: main node (incl content) and first tuple have been prefetched, check join predicate
      } else if (amac_todo_et::kAmacVerifyMainByEq == lAmac._todo) {
        const tuple_prb_t* lTuplePrb = lAmac._tuple_prb;
        ++_numCmps2;
        if (joinpred_t::eval(lTuplePrb, e->content().get_first_tuple())) {
          // join predicate evaluates to true, so all tuples in this entry match
          ++_c_mid;
          if constexpr (lTrace) { std::cout << "    " << "join predicate returns true." << std::endl; }
          // produce join results
          concatfun_t::eval_prb(_output, lTuplePrb);
          for (uint i = e->content().begin(); i < e->content().end(); ++i) {
            assert(joinpred_t::eval(lTuplePrb, e->content().tuple(i)));
            concatfun_t::eval_bld(_output, e->content().tuple(i));
            if constexpr(lTrace) {
              std::cout << "      " << "rs[" << _output._r->k << ',' << _output._s->k << ']' << std::endl;
            }
            inc();
            _consumer->step(&_output);
          }
          // if there is a sub chain, every sub chain entry also matches
          if (nullptr != e->sub()) {
            lAmac._ht_entry_sub = e->sub();
            lAmac._todo = amac_todo_et::kAmacUnnestSub;
            PREFETCH_PROBE(lAmac._ht_entry_sub);
            if constexpr (lTrace) { std::cout << "    " << "sub chain exists. 3 -> 4." << std::endl; }
          } else {
            // this is the only matching main node. there are no sub nodes, so we're done
            lAmac.clear();
            if constexpr (lTrace) { std::cout << "    " << "sub chain does not exist. 3 -> 0. done." << std::endl; }
            break;

          }
        } else {
          // join predicate evaluates to false. no tuple in this entry matches. continue with the next main entry.
          if constexpr (lTrace) { std::cout << "    " << "join predicate returns false." << std::endl; }
          if (nullptr != e->next()) {
            lAmac._ht_entry_main = e->next();
            PREFETCH_PROBE(lAmac._ht_entry_main);
            if constexpr (lTrace) { std::cout << "    " << "more main nodes to examine. 3 -> 2." << std::endl; }
          } else {
            // no more main entries to examine
            if constexpr (lTrace) { std::cout << "    " << "no more main nodes to examine. 3 -> 0. done." << std::endl; }
            lAmac.clear();
            break;
          }
        }

      // stage 4: walk the sub chain (= unnest) and produce output tuples
      } else if (amac_todo_et::kAmacUnnestSub == lAmac._todo) {
        const tuple_prb_t* lTuplePrb = lAmac._tuple_prb;
        const entry_sub_t* s = lAmac._ht_entry_sub;
        concatfun_t::eval_prb(_output, lTuplePrb);

        assert(nullptr != s);
        for (uint i = s->content().begin(); i < s->content().end(); ++i) {
          concatfun_t::eval_bld(_output, s->content().tuple(i));
          if constexpr(lTrace) {
            std::cout << "      " << "rs[" << _output._r->k << ',' << _output._s->k << ']' << std::endl;
          }
          inc();
          _consumer->step(&_output);
        }

        if (nullptr != s->next()) {
          lAmac._ht_entry_sub = s->next();
          PREFETCH_PROBE(lAmac._ht_entry_sub);
          if constexpr (lTrace) { std::cout << "    " << "more sub nodes to unnest. 4 -> 4." << std::endl; }
        } else {
          if constexpr (lTrace) { std::cout << "    " << "no more sub nodes to unnest. 4 -> 0. done." << std::endl; }
          lAmac.clear();
          break;
        }

      } else {
        throw std::logic_error("Unexpected AMAC stage in AlgOpJoinHash3dAmacProbe2::step(): " + std::to_string(_amac[_curr].todo()));
      }

      amac_advance();
    }

    // if we arrived here, _amac[_curr] is a free slot in the ring buffer and we can insert aProbeTuple here
    if constexpr (lTrace) {
      std::cout
        << "    " << "inserting probe tuple " << *aProbeTuple
        << " into the AMAC buffer at index " << _curr << "." << " 0 -> 1."
        << std::endl;
    }
    _amac[_curr]._tuple_prb = aProbeTuple;
    _amac[_curr]._ht_entry_main = lMainEntryArgument;
    _amac[_curr]._ht_entry_sub = nullptr;
    _amac[_curr]._hashval_prb = lHashvalArgument;
    _amac[_curr]._todo = amac_todo_et::kAmacCheckDirectory;
    PREFETCH_PROBE(lMainEntryArgument);
    amac_advance();

  }
  else {
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
    const hashtable_t& lHt = _build->hashtable();
    const hashval_t lProbeHash = hashfun_t::eval(aProbeTuple);

    const entry_main_t* e = lHt.get_entry(lProbeHash);  // pointer calc, no mem access

    // must catch the case where directory entry is empty
    if (e->content().is_empty()) {  // access *e
      if constexpr (lTrace) {
        std::cout << "    " << "empty directory entry" << std::endl;
      }
    } else {
      uint64_t lNumComparisons1 = 0;  // number of collision chain *hashval* comparisons (-> access e)
      uint64_t lNumComparisons2 = 0;  // number of collision chain *joinpred* evaluations (-> access tuple)

      while (nullptr != e) {
        ++lNumComparisons1;
        if (e->content().hashval() == lProbeHash) {
          ++lNumComparisons2;
          if (joinpred_t::eval(aProbeTuple, e->content().get_first_tuple())) {
            break;
          }
        }
        e = e->next();  // access *e
      }

      _numCmps1 += lNumComparisons1;
      _numCmps2 += lNumComparisons2;

      if (nullptr != e) {
        // e is a match
        if constexpr (lTrace) {
          std::cout << "    " << "found matching main node: " << e << std::endl;
        }
        // UNNEST
        // step 1: unnest e
        if constexpr (lTrace) {
          std::cout << "    " << "unnest main node" << std::endl;
        }
        concatfun_t::eval_prb(_output, aProbeTuple);
        for (uint i = e->content().begin(); i < e->content().end(); ++i) {
          concatfun_t::eval_bld(_output, e->content().tuple(i));
          if constexpr (lTrace) {
              std::cout << "    " << "rs[" << _output._r->k << ',' << _output._s->k << ']' << std::endl;
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
              std::cout << "    " << "rs[" << _output._r->k << ',' << _output._s->k << ']' << std::endl;
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
}

template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun>
void
AlgOpJoinHash3dAmacProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun>
::fin() {
  constexpr bool lTrace = false;
  constexpr bool lAmacInStep = true;  // AMAC was used in step, must finish the buffer
  constexpr bool lSimple = false;  // true -> just finish the AMAC buffer to completion without further prefetching

  if constexpr (lTrace) {
    std::cout
      << std::format("  AlgOpJoinHash3dAmacProbe2::fin() [amac_in_step={}, simple_amac={}]", lAmacInStep, lSimple)
      << std::endl;
  }

  if constexpr (lAmacInStep) {
    if constexpr (!lSimple) {
      // continue AMAC-style processing with prefetching.
      // build a bitmap that indicates which AMAC buffer elements are not empty
      uint64_t lBitmapAmacTodo = 0;  // bit i is set iff _amac[i].is_occupied()
      for (uint i = 0; i < AMAC_BUFFER_SIZE; ++i) {
        lBitmapAmacTodo |= (_amac[i].is_occupied() << i);
      }

      // while there are occupied AMAC buffer elements, keep processing them like in step() one stage after the other, with prefetching
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
          const entry_main_t* e = lAmac._ht_entry_main;
          // stage 0
          if (amac_todo_et::kAmacEmpty == lAmac._todo) { continue; }

          // note that stage 1 falls through to stage 2,
          // but stages 2-4 are if/else if constructs.

          // stage 1
          if (amac_todo_et::kAmacCheckDirectory == lAmac._todo) {
            if (e->content().is_empty()) {
              if constexpr (lTrace) { std::cout << "      " << "directory entry is empty. 1 -> 0. done." << std::endl; }
              lAmac.clear();
              lBitmapAmacTodo &= ~(1 << _curr);
              continue;  // for: next buffer entry
            } else {
              if constexpr (lTrace) { std::cout << "      " << "1 -> 2." << std::endl; }
              lAmac._todo = amac_todo_et::kAmacFindMainByHash;
            }
          }

          // stage 2: main node has been prefetched, check hash value equality
          if (amac_todo_et::kAmacFindMainByHash == lAmac._todo) {
            assert((nullptr != e) && (!e->content().is_empty()));
            ++_numCmps1;
            if (e->content().hashval() == lAmac._hashval_prb) {
              // found a potentially matching main node, check join predicate next
              if constexpr (lTrace) { std::cout << "      " << "hash values match, check join predicate next. 2 -> 3." << std::endl; }
              PREFETCH_PROBE(e->content().get_first_tuple());
              lAmac._todo = amac_todo_et::kAmacVerifyMainByEq;
            } else {
              // hash value does not match, check the next node
              if constexpr (lTrace) { std::cout << "      " << "hash values differ." << std::endl; }
              if (nullptr == e->next()) {
                // no more main nodes to examine, no more matches, done
                if constexpr (lTrace) { std::cout << "      " << "no more main nodes to examine. 2 -> 0. done." << std::endl; }
                lAmac.clear();
                lBitmapAmacTodo &= ~(1 << _curr);
                continue;
              } else {
                // more main nodes to examine
                PREFETCH_PROBE(e->next());
                lAmac._ht_entry_main = e->next();
                if constexpr (lTrace) { std::cout << "      " << "more main nodes to examine. 2 -> 2." << std::endl; }
              }
            }
          }
          // stage 3: main node (incl content) and first tuple have been prefetched, check join predicate
          else if (amac_todo_et::kAmacVerifyMainByEq == lAmac._todo) {
            const tuple_prb_t* lTuplePrb = lAmac._tuple_prb;
            ++_numCmps2;
            if (joinpred_t::eval(lTuplePrb, e->content().get_first_tuple())) {
              ++_c_mid;
              // join predicate evaluates to true, so all tuples in this entry match
              if constexpr (lTrace) { std::cout << "      " << "join predicate returns true." << std::endl; }
              // produce join results
              concatfun_t::eval_prb(_output, lTuplePrb);
              for (uint i = e->content().begin(); i < e->content().end(); ++i) {
                assert(joinpred_t::eval(lTuplePrb, e->content().tuple(i)));
                concatfun_t::eval_bld(_output, e->content().tuple(i));
                if constexpr(lTrace) {
                  std::cout << "        " << "rs[" << _output._r->k << ',' << _output._s->k << ']' << std::endl;
                }
                inc();
                _consumer->step(&_output);
              }
              // if there is a sub chain, every sub chain entry also matches
              if (nullptr != e->sub()) {
                lAmac._ht_entry_sub = e->sub();
                lAmac._todo = amac_todo_et::kAmacUnnestSub;
                PREFETCH_PROBE(lAmac._ht_entry_sub);
                if constexpr (lTrace) { std::cout << "      " << "sub chain exists. 3 -> 4." << std::endl; }
              } else {
                // this is the only matching main node. there are no sub nodes, so we're done
                lAmac.clear();
                lBitmapAmacTodo &= ~(1 << _curr);
                if constexpr (lTrace) { std::cout << "      " << "sub chain does not exist. 3 -> 0. done." << std::endl; }
                continue;

              }
            } else {
              // join predicate evaluates to false. no tuple in this entry matches. continue with the next main entry.
              if constexpr (lTrace) { std::cout << "      " << "join predicate returns false." << std::endl; }
              if (nullptr != e->next()) {
                lAmac._ht_entry_main = e->next();
                PREFETCH_PROBE(lAmac._ht_entry_main);
                if constexpr (lTrace) { std::cout << "      " << "more main nodes to examine. 3 -> 2." << std::endl; }
              } else {
                // no more main entries to examine
                if constexpr (lTrace) { std::cout << "      " << "no more main nodes to examine. 3 -> 0. done." << std::endl; }
                lAmac.clear();
                lBitmapAmacTodo &= ~(1 << _curr);
                continue;
              }
            }
          }
          // stage 4
          else if (amac_todo_et::kAmacUnnestSub == lAmac._todo) {
            const tuple_prb_t* lTuplePrb = lAmac._tuple_prb;
            const entry_sub_t* s = lAmac._ht_entry_sub;
            concatfun_t::eval_prb(_output, lTuplePrb);

            assert(nullptr != s);
            for (uint i = s->content().begin(); i < s->content().end(); ++i) {
              concatfun_t::eval_bld(_output, s->content().tuple(i));
              if constexpr(lTrace) {
                std::cout << "      " << "rs[" << _output._r->k << ',' << _output._s->k << ']' << std::endl;
              }
              inc();
              _consumer->step(&_output);
            }

            if (nullptr != s->next()) {
              lAmac._ht_entry_sub = s->next();
              PREFETCH_PROBE(lAmac._ht_entry_sub);
              if constexpr (lTrace) { std::cout << "    " << "more sub nodes to unnest. 4 -> 4." << std::endl; }
            } else {
              if constexpr (lTrace) { std::cout << "    " << "no more sub nodes to unnest. 4 -> 0. done." << std::endl; }
              lAmac.clear();
              lBitmapAmacTodo &= ~(1 << _curr);
              continue;
            }
          }
          else {
            throw std::logic_error("Unexpected AMAC stage in AlgOpJoinHash3dAmacProbe2::fin(): " + std::to_string(_amac[_curr].todo()));
          }
        }  // end for
      }  // end while
      assert(lBitmapAmacTodo == 0);
    }
    else {
      // simple approach with AMAC:
      // go through each AMAC buffer element and finish it to completion, no prefetching
      for (_curr = 0; _curr < AMAC_BUFFER_SIZE; ++_curr) {
        amac_t& lAmac = _amac[_curr];
        if constexpr (lTrace) {
          std::cout << "    " << "_curr = " << _curr <<  " | " << "todo = " << lAmac.todo() << std::endl;
        }

        if (lAmac.is_empty()) { continue; }  // nothing to do for this entry

        if constexpr (lTrace) {
          std::cout << "    " << "probe tuple = " << *lAmac._tuple_prb << std::endl;
        }

        const entry_main_t* e = lAmac._ht_entry_main;

        // stage 1
        if (amac_todo_et::kAmacCheckDirectory == lAmac._todo) {
          if constexpr (lTrace) { std::cout << "    " << "todo = " << lAmac.todo() << std::endl; }
          if (e->content().is_empty()) {
            if constexpr (lTrace) { std::cout << "      " << "empty directory. 1 -> 0. done." << std::endl; }
            lAmac.clear();
            continue;
          }
          if constexpr (lTrace) { std::cout << "      " << "1 -> 2." << std::endl; }
          lAmac._todo = amac_todo_et::kAmacFindMainByHash;
        }

        // stage 2 & 3
        if ((amac_todo_et::kAmacFindMainByHash == lAmac._todo) ||
            (amac_todo_et::kAmacVerifyMainByEq == lAmac._todo)) {
          if constexpr (lTrace) { std::cout << "    " << "todo = " << lAmac.todo() << std::endl; }

          uint64_t lNumComparisons1 = 0;  // number of collision chain *hashval* comparisons (-> access e)
          uint64_t lNumComparisons2 = 0;  // number of collision chain *joinpred* evaluations (-> access tuple)

          const tuple_prb_t* lProbeTuple = lAmac._tuple_prb;

          while (nullptr != e) {
            ++lNumComparisons1;
            if (e->content().hashval() == lAmac._hashval_prb) {
              ++lNumComparisons2;
              if (joinpred_t::eval(lProbeTuple, e->content().get_first_tuple())) {
                break;  // while
              }
            }
            e = e->next();  // access *e
          }

          _numCmps1 += lNumComparisons1;
          _numCmps2 += lNumComparisons2;

          if (nullptr != e) {
            // e is a match
            if constexpr (lTrace) { std::cout << "      " << "found matching main node: " << e << std::endl; }
            // UNNEST
            // step 1: unnest e
            if constexpr (lTrace) { std::cout << "      " << "unnest main node" << std::endl; }
            concatfun_t::eval_prb(_output, lProbeTuple);
            for (uint i = e->content().begin(); i < e->content().end(); ++i) {
              concatfun_t::eval_bld(_output, e->content().tuple(i));
              if constexpr (lTrace) {
                std::cout << "      " << "rs[" << _output._r->k << ',' << _output._s->k << ']' << std::endl;
              }
              inc();
              _consumer->step(&_output);
            }
            if (nullptr != e->sub()) {
              if constexpr (lTrace) { std::cout << "      " << "check sub nodes next. 2/3 -> 4." << std::endl; }
              lAmac._todo = amac_todo_et::kAmacUnnestSub;
              lAmac._ht_entry_sub = e->sub();
            } else {
              if constexpr (lTrace) { std::cout << "      " << "no sub nodes. 2/3 -> 0. done." << std::endl; }
              lAmac.clear();
              continue;
            }
          } else {
            if constexpr (lTrace) { std::cout << "      " << "no match found. 2/3 -> 0. done." << std::endl; }
            lAmac.clear();
            continue;
          }
        }

        // stage 4
        if (amac_todo_et::kAmacUnnestSub == lAmac._todo) {
          if constexpr (lTrace) { std::cout << "    " << "todo = " << lAmac.todo() << std::endl; }
          // step 2: unnest subchain
          if constexpr (lTrace) {
            std::cout << "      " << "unnest sub nodes" << std::endl;
          }
          concatfun_t::eval_prb(_output, lAmac._tuple_prb);
          const entry_sub_t* s = lAmac._ht_entry_sub;
          while (nullptr != s) {
            // access *s
            for (uint i = s->content().begin(); i < s->content().end(); ++i) {
              concatfun_t::eval_bld(_output, s->content().tuple(i));
              if constexpr (lTrace) {
                std::cout << "      " << "rs[" << _output._r->k << ',' << _output._s->k << ']' << std::endl;
              }
              inc();
              _consumer->step(&_output);
            }
            s = s->next();
          }
          if constexpr (lTrace) { std::cout << "      " << "4 -> 0. done." << std::endl; }
          lAmac.clear();
        }
      }


    }  // end lSimple == true
    if constexpr (lTrace) { amac_print_buffer(4); }
#ifndef NDEBUG
    // make sure that all AMAC buffer elements are empty
    for (size_t i = 0; i < AMAC_BUFFER_SIZE; ++i) {
      assert(_amac[i].is_empty() && "AMAC ring buffer has non-empty elements at the end of fin.");
    }
#endif
  }  // end if lAmacInStep
  _consumer->fin();
}

template <alg2_consumer_c Tconsumer, alg2_operator_build_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg2_concatfun_c Tconcatfun>
void
AlgOpJoinHash3dAmacProbe2<Tconsumer, Tbuild, Thashfun, Tjoinpred, Tconcatfun>
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


