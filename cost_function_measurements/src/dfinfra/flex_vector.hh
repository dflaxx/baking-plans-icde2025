#ifndef FLEX_VECTOR_HH
#define FLEX_VECTOR_HH
#pragma once

#include "standard_includes.hh"

#include <cstddef>
#include <cstring>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>


/*
 * A container with a maximum fixed size
 * that stores its elements contiguously on the heap.
 *
 * It's somewhat in the middle between a std::vector and a std::array:
 * The size does not need to be know at compile time (in contrast to a std::array),
 * but the capacity cannot be changed dynamically (in contrast to a std::vector).
 * In constrast to FixedVectorNC, the size can be adjusted at runtime as long 
 * as it is smaller than the capacity.
 * A resize does not touch the contents of the vector,
 * no matter whether it increases or decreases the current size.
 * Values at positions size <= i < capacity are merely "hidden".
 *
 * Upon memory allocation, size == capacity.
 *
 * Memory must be explicitly allocated before use and deallocated after use (no RAII).
 * Content of the array is not initialized, hence the "NC" in the class name.
 *
 * Template parameter BoomOnMissingDealloc (BOMD) directs if missing deallocs triggers a warning.
 *
 * Internals
 * - uses std::aligned_alloc and std::free to alloc and dealloc memory.
 */

namespace df::infra {

template<typename T, bool BoomOnMissingDealloc = false>
class FlexVectorNC {
  static_assert(std::is_copy_assignable_v<T> && std::is_move_assignable_v<T>);

  public:
    using value_type = T;
    using pointer_type = T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = std::size_t;

    //static constexpr std::size_t DEFAULT_ALIGNMENT = __STDCPP_DEFAULT_NEW_ALIGNMENT__;
    static constexpr std::size_t DEFAULT_ALIGNMENT = 64;

    template <bool IsConst>
    class Iterator;
    using iterator_type = Iterator<false>;
    using const_iterator_type = Iterator<true>;

    static_assert(std::contiguous_iterator<iterator_type>);
    static_assert(std::contiguous_iterator<const_iterator_type>);

  public:
    FlexVectorNC() : _data(nullptr), _size(0), _capacity(0) {}

    FlexVectorNC(const size_type aCapacity, const bool aWithAlloc = false, const size_type aAlignment = DEFAULT_ALIGNMENT);

    // delete copy constructor
    FlexVectorNC(const FlexVectorNC<T, BoomOnMissingDealloc>&) = delete;

    // delete copy assignment operator
    const FlexVectorNC<T, BoomOnMissingDealloc>&
    operator=(const FlexVectorNC<T, BoomOnMissingDealloc>&) = delete;

    ~FlexVectorNC();

  public:
    // aSize: number of elements (not bytes)
    bool mem_alloc(const size_type aSize, const size_type aAlignment = DEFAULT_ALIGNMENT);
    bool mem_init();
    bool mem_free();

    inline reference operator[](const size_type aPos) { return _data[aPos]; }
    inline const_reference operator[](const size_type aPos) const { return _data[aPos]; }
    inline reference at(const size_type aPos) { assert_pos_capa(aPos); return _data[aPos]; }
    inline const_reference at(const size_type aPos) const { assert_pos_capa(aPos); return _data[aPos]; }

    inline reference front() { return _data[0]; }
    inline const_reference front() const { return _data[0]; }
    inline reference back() { return _data[size() - 1]; }
    inline const_reference back() const { return _data[size() - 1]; }

    // note: range-based for-loop never calls cbegin()/cend(), even for const objects!
    // must provide const-overloads for begin()/end().
    inline iterator_type begin() { return iterator_type(&front()); }
    inline iterator_type end() { return iterator_type(&back() + 1); }
    inline const_iterator_type begin() const { return const_iterator_type(&front()); }
    inline const_iterator_type end() const { return const_iterator_type(&back() + 1); }
    inline const_iterator_type cbegin() const { return const_iterator_type(&front()); }
    inline const_iterator_type cend() const { return const_iterator_type(&back() + 1); }

  public:
    // getters
    inline bool is_allocated() const { return (size() > 0) && (_data != nullptr); }
    inline size_type size() const { return _size; }
    inline size_type capacity() const { return _capacity; }
    inline size_type byte_size() const { return (capacity() * sizeof(T)); }  // XXX returns capacity, not size!
    // setters
    inline void resize(const size_type aSize) { assert(aSize <= capacity()); _size = aSize; }

  private:
    inline void assert_pos_size(const size_type pos);
    inline void assert_pos_capa(const size_type pos);

  private:
    T*        _data;
    size_type _size;
    size_type _capacity;

  public:
    template <bool IsConst>
    class Iterator {
      public:
        using iterator_category [[maybe_unused]] = std::contiguous_iterator_tag;
        using iterator_concept [[maybe_unused]] = std::contiguous_iterator_tag;
        using self_type = Iterator<IsConst>;
        // value_type should be cv-unqualified, element_type should be cv_qualified.
        // See https://stackoverflow.com/a/66050521.
        using value_type = T;
        using element_type = typename std::conditional<IsConst, const T, T>::type;
        using difference_type   = std::ptrdiff_t;
        using pointer   = typename std::conditional<IsConst, const element_type*, element_type*>::type;
        using reference = typename std::conditional<IsConst, const element_type&, element_type&>::type;

      public:
        inline Iterator() : _current(nullptr) {}
        inline Iterator(pointer p) : _current(p) {}

        // Always permit conversion (copy constructor) from iterator to const iterator
        inline Iterator(const Iterator<false>& other) : _current(other._current) {}

        // Always permit copy-assignment from iterator to const iterator.
        inline self_type& operator=(const Iterator<false>& other) {
          // must compare the raw addresses when mixing Iterator<true> and Iterator<false>.
          // const void* for const-correctness (cannot cast const Iterator<true>* to void*.)
          if (static_cast<const void*>(this) == static_cast<const void*>(&other)) {
            return *this;
          }
          _current = other._current;
          return *this;
        }

      public:
        reference operator*() const;
        pointer operator->() const;
        reference operator[](const size_t) const;

        self_type& operator++();
        self_type operator++(int);
        self_type& operator--();
        self_type operator--(int);

        self_type& operator+=(const difference_type);
        self_type& operator-=(const difference_type);
        
        // covers all comparisons: <, <=, =, !=, >=, >
        friend auto operator<=>(self_type, self_type) = default;

        // free arithmetic functions
        // XXX declared here because not easy to declare out-of-line.

        // it + 7
        friend self_type operator+(const self_type& lhs, const difference_type rhs) {
          return self_type(lhs._current + rhs);
        }
        // it - 7
        friend self_type operator-(const self_type& lhs, const difference_type rhs) {
          return self_type(lhs._current - rhs);
        }
        // 7 + it
        friend self_type operator+(const difference_type lhs, const self_type& rhs) {
          return self_type(lhs + rhs._current);
        }
        friend difference_type operator-(const self_type& lhs, const self_type& rhs) {
          return lhs._current - rhs._current;
        }

        // allow const Iterator to access private members of non-const Iterator (needed for copy constructor)
        friend Iterator<true>;

      private:
        pointer _current;
    };
};


/* FlexVectorNC */

template <typename T, bool BOMD>
FlexVectorNC<T, BOMD>::FlexVectorNC(const size_type aCapacity, const bool aWithAlloc, const size_type aAlignment)
    : _data(nullptr), _size(0), _capacity(aCapacity) {
  if (aWithAlloc) {
    mem_alloc(aCapacity, aAlignment);
  }
}

template <typename T, bool BOMD>
FlexVectorNC<T, BOMD>::~FlexVectorNC() {
  if constexpr (BOMD) {
    if (is_allocated()) {
      std::cerr << __PRETTY_FUNCTION__ << ": " << "Memory has not been freed before dtor call." << std::endl;
      //throw std::logic_error("Memory has not been freed before dtor call.");
      // XXX destructors should not throw exeptions
    }
  }
  mem_free();
}

template <typename T, bool BOMD>
bool
FlexVectorNC<T, BOMD>::mem_alloc(const size_type aCapacity, const size_type aAlignment) {
  assert(nullptr == _data);
  _capacity = aCapacity;
  _size     = aCapacity;
  void* p = std::aligned_alloc(aAlignment, capacity() * sizeof(T));
  _data = static_cast<pointer_type>(p);
  assert(nullptr != _data);
  return (nullptr != _data);
}

template <typename T, bool BOMD>
bool
FlexVectorNC<T, BOMD>::mem_init() {
  if(nullptr == _data) { return false; }
  std::memset(_data, 0, byte_size());
  return true;
}

template <typename T, bool BOMD>
bool
FlexVectorNC<T, BOMD>::mem_free() {
  // Note from https://en.cppreference.com/w/cpp/memory/c/free:
  // The function accepts (and does nothing with) the null pointer to reduce the amount of special-casing.
  // Whether allocation succeeds or not, the pointer returned by an allocation function can be passed to std::free. 
  std::free(_data);
  _data = nullptr;
  _capacity = 0;
  _size     = 0;
  return true;
}

template <typename T, bool BOMD>
void
FlexVectorNC<T, BOMD>::assert_pos_size(const size_type aPos) {
  if (aPos >= size()) {
    throw std::out_of_range("FlexVectorNC: pos (" + std::to_string(aPos) + ") >= size (" + std::to_string(size()) + ")");
  }
}

template <typename T, bool BOMD>
void
FlexVectorNC<T, BOMD>::assert_pos_capa(const size_type aPos) {
  if (aPos >= capacity()) {
    throw std::out_of_range("FlexVectorNC: pos (" + std::to_string(aPos) + ") >= size (" + std::to_string(size()) + ")");
  }
}


/* Iterator */

template <typename T, bool BOMD>
template <bool IsConst>
typename FlexVectorNC<T, BOMD>::template Iterator<IsConst>::reference
FlexVectorNC<T, BOMD>::Iterator<IsConst>::operator*() const {
  return *_current; 
}

template <typename T, bool BOMD>
template <bool IsConst>
typename FlexVectorNC<T, BOMD>::template Iterator<IsConst>::pointer
FlexVectorNC<T, BOMD>::Iterator<IsConst>::operator->() const {
  return _current;
}

template <typename T, bool BOMD>
template <bool IsConst>
typename FlexVectorNC<T, BOMD>::template Iterator<IsConst>::reference
FlexVectorNC<T, BOMD>::Iterator<IsConst>::operator[](const size_t n) const {
  return *(_current + n); 
}

template <typename T, bool BOMD>
template <bool IsConst>
typename FlexVectorNC<T, BOMD>::template Iterator<IsConst>::self_type&
FlexVectorNC<T, BOMD>::Iterator<IsConst>::operator++() {
  _current++;
  return *this;
}

template <typename T, bool BOMD>
template <bool IsConst>
typename FlexVectorNC<T, BOMD>::template Iterator<IsConst>::self_type
FlexVectorNC<T, BOMD>::Iterator<IsConst>::operator++(int) {
  self_type tmp = *this;
  ++(*this);
  return tmp;
}

template <typename T, bool BOMD>
template <bool IsConst>
typename FlexVectorNC<T, BOMD>::template Iterator<IsConst>::self_type&
FlexVectorNC<T, BOMD>::Iterator<IsConst>::operator--() {
  _current--;
  return *this;
}

template <typename T, bool BOMD>
template <bool IsConst>
typename FlexVectorNC<T, BOMD>::template Iterator<IsConst>::self_type
FlexVectorNC<T, BOMD>::Iterator<IsConst>::operator--(int) {
  self_type tmp = *this;
  ++(*this);
  return tmp;
}

template <typename T, bool BOMD>
template <bool IsConst>
typename FlexVectorNC<T, BOMD>::template Iterator<IsConst>::self_type&
FlexVectorNC<T, BOMD>::Iterator<IsConst>::operator+=(const difference_type aOffset) {
  _current += aOffset;
  return *this;
}

template <typename T, bool BOMD>
template <bool IsConst>
typename FlexVectorNC<T, BOMD>::template Iterator<IsConst>::self_type&
FlexVectorNC<T, BOMD>::Iterator<IsConst>::operator-=(const difference_type aOffset) {
  _current -= aOffset;
  return *this;
}


}

#endif

