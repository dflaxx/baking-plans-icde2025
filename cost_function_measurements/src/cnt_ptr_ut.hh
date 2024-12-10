#ifndef CNT_PTR_UT_HH
#define CNT_PTR_UT_HH
#pragma once

#include "dfinfra/standard_includes.hh"
#include <limits>

/*
 * helper struct containing 
 * the union of a counter and a pointer
 *
 * Author: Guido Moerkotte
 */

using luint_t = long unsigned int;

template<class Tptr>
union cnt_ptr_ut {
   using ptr_t = Tptr;
   luint_t _cnt;
   ptr_t   _ptr;
   cnt_ptr_ut() : _cnt(0x1) {}
   cnt_ptr_ut(const luint_t c) : _cnt((c << 1) | 0x1) {}
   cnt_ptr_ut(ptr_t  aPtr) : _ptr(aPtr) {}
   inline bool     is_cnt() const { return (0 != (_cnt & 0x1)); }
   inline bool     is_ptr() const { return (0 == (_cnt & 0x1)); }
   inline luint_t  count_unchecked() const { return (_cnt >> 1); }
   inline luint_t  count() const { return (is_ptr() ? std::numeric_limits<luint_t>::max() : (_cnt >> 1)); }
   inline void     inc() { _cnt += 2; }
   inline void     dec() { _cnt -= 2; }
   inline ptr_t    ptr() { return (is_ptr() ? _ptr : nullptr); }
   inline ptr_t    ptr_unchecked() { return _ptr; }
   inline void     init() { _cnt = 0x1; }
   inline cnt_ptr_ut& operator=(ptr_t x) { _ptr = x; return (*this); }
   inline cnt_ptr_ut& operator=(const luint_t c) { _cnt = (c << 1) | 0x1; return (*this); }
};

#endif
