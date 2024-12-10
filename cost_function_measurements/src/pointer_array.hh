#ifndef POINTER_ARRAY_FIXED_SIZE_TT_HH
#define POINTER_ARRAY_FIXED_SIZE_TT_HH
#pragma once

#include "dfinfra/standard_includes.hh"
#include "cnt_ptr_ut.hh"

/*
 * tiny fixed size array/container for pointers
 * first element of _ptr holds idx to the next free element in _ptr
 *
 * states
 * - _ptr[0] = 0/nullptr : indicates empty array: see is_empty()
 * - _ptr[0] & 0x1       : indicate counter = index of first empty array slot
 * - _ptr[0] & 0x1 == 0  : pointer, array is full
 *
 * _ptr is filled from end, see insert()
 *
 * begin() also returns the index of the last successful insert
 * end() always returns n (one index past the last valid index) -> cheap for-loop boundary check
 *
 * do not insert nullptr
 * pointers must have the least significant bit equal to zero
 *
 * Author: Guido Moerkotte
 *
 * Details
 *
 *   0            n-1
 *  -----------------
 * | x |   | ... |   |
 *  -----------------
 *   ^--contains counter
 *    <-- insert --
 *
 * Example: n = 4, indexes 0 ... 3
 *   n   #inserts   begin()   end()   _ptr[0]   --(>> 1)-->  count_unchecked
 *   4   0          4         4       0000 = 0               0000 = 0
 *   4   1          3         4       0101 = 5               0010 = 2 (= next_empty)
 *   4   2          2         4       0011 = 3               0001 = 1 (= next_empty)
 *   4   3          1         4       0001 = 1               0000 = 0 (= next_empty)
 *   4   4          0         4       xxx0                   (not meaningful)
 */

template<typename Tptr, uint n>
struct pointer_array_tt {
  using ptr_t = Tptr;

  static constexpr uint cnt_idx = 0;
  static constexpr uint capacity = n;

  inline __attribute__((always_inline)) bool is_empty() const { return (nullptr == _ptr[cnt_idx]); }
  inline __attribute__((always_inline)) bool is_full()  const { return ((nullptr != _ptr[cnt_idx]) && (0 == (((luint_t)(_ptr[cnt_idx])) & 0x1))); }
  // next_empty: only valid if not full
  inline __attribute__((always_inline)) luint_t next_empty() const {
    //return ((luint_t)(_ptr[cnt_idx]) >> 1);  // Guido
    assert(!is_full());
    return (is_empty() ? n-1 : (is_full() ? 0 : (cnt_ptr_ut<ptr_t>(_ptr[cnt_idx]).count_unchecked()))); 
  }
  inline __attribute__((always_inline)) luint_t size() const { return end() - begin(); }

  inline __attribute__((always_inline)) uint begin() const {
    // return (is_empty() ? n : (is_full() ? 0 : (cnt_ptr_ut(_ptr[cnt_idx]).count_unchecked()) + 1)); 
    // old gcc on ampere requires
    return (is_empty() ? n : (is_full() ? 0 : (cnt_ptr_ut<ptr_t>(_ptr[cnt_idx]).count_unchecked()) + 1)); 
  }
  inline constexpr uint end() const { return n; }
 
  // unchecked, make sure you have i in a valid range 
  inline __attribute__((always_inline)) ptr_t ptr(const uint i) const { return _ptr[i]; }
  inline __attribute__((always_inline)) ptr_t operator[](const uint i) const { return _ptr[i]; }

  // return the first inserted tuple
  // make sure that it is not empty!
  inline __attribute__((always_inline)) ptr_t get_first_ptr() const { return _ptr[n - 1]; }
  // inserts into this, returns true if not full, returns false if full
  inline __attribute__((always_inline)) bool insert(const ptr_t aPtr);

  inline __attribute__((always_inline)) void clear() {
    /*
    for(uint i = 0; i < n; ++i) {
      _ptr[i] = nullptr;
    }
    */
    _ptr[cnt_idx] = nullptr; // sic!, wuerde reichen
  }

  ptr_t _ptr[n];
};

template<typename T, uint n>
bool
pointer_array_tt<T,n>::insert(const ptr_t aPtr) {
  assert(nullptr != aPtr);  // don't insert nullptr
  assert(0 == (reinterpret_cast<luint_t>(aPtr) & 0x1));  // all pointers must have LSB == 0

  cnt_ptr_ut lPtr0 = _ptr[cnt_idx];  // get counter/pointer at element _ptr[0]

  // first insert (empty array)
  if (is_empty()) {
    _ptr[n - 1] = aPtr;  // insert at last position
    if constexpr (2 <= n) {
      // if array can hold 2 or more pointers, update the counter
      lPtr0 = luint_t{n} - 2;
      _ptr[cnt_idx] = lPtr0.ptr_unchecked();
    }
    return true;
  }
  
  // full (insert not possible)
  if (lPtr0.is_ptr()) { assert(is_full()); return false; }

  // other inserts (not full, not empty)
  const luint_t lCnt = lPtr0.count_unchecked();  // _ptr[0]'s counter holds index of next free element

  if (0 != lCnt) {
    // this is not the last possible insert, so we must update the counter
    lPtr0.dec();
    _ptr[cnt_idx] = lPtr0.ptr_unchecked();
  }
  _ptr[lCnt] = aPtr;
  return true;
}

#endif
