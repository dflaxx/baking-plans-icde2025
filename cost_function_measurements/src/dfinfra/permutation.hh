#pragma once

#include "standard_includes.hh"
#include <sstream>
#include <type_traits>

namespace df::infra {

/*
 * Permutation with n elements from 0, ..., n-1
 * Elements are of type Tuint (unsigned int type).
 *
 * Class has functions to
 * - invert
 * - compose
 * - rank (permutation -> uint) and unrank (uint -> permutation)
 *   (algorithms by Myrvold/Ruskey)
 * - determine if permutation is reducible
 *   or has common subset prefix (CSP) with another permutation
 * 
 * TODO:
 * - comparison operators <, <=, >=, > (lexicographic ordering)
 * - assignment operator
 * - isCyclic(), toCylceNotation()
 */
template <typename Tuint>
class Permutation {
public:
  using self_t   = Permutation<Tuint>;
  using value_t  = Tuint;
  using value_vt = std::vector<value_t>;

public:
  enum class special_perm { id, rev };

public:
  /* constructors */
  Permutation()                   : _data()      {}
  explicit Permutation(const size_t aSize) : _data(aSize) {
    init(aSize, special_perm::id);
  }
  explicit Permutation(const size_t aSize, const special_perm aSpecial) : _data(aSize) {
    init(aSize, aSpecial);
  }
  Permutation(const value_vt& aData) : _data(aData)                   {}
  Permutation(value_vt&& aData)      : _data(std::move(aData))        {}
  Permutation(const self_t& aOther)  : _data(aOther._data)            {}
  Permutation(self_t&& aOther)       : _data(std::move(aOther._data)) {}

  Permutation(const size_t aSize, const uint64_t aRank) : _data(aSize) {
    unrank(aSize, aRank);
  }
public:
  inline self_t& operator=(const self_t& aOther) {  // copy-assignment
    if (this != &aOther) {
      _data = aOther._data;
    }
    return *this;
  }

  inline self_t& operator=(self_t&& aOther) {  // move-assignment
    if (this != &aOther) {
      _data = std::move(aOther._data);
    }
    return *this;
  }

public:
  void init(const size_t aSize, const special_perm aSpecial) {
    _data.resize(aSize);
    switch(aSpecial) {
      case special_perm::id:
        for (size_t i = 0; i < aSize; ++i) { _data[i] = i; };
        break;
      case special_perm::rev:
        for (size_t i = 0; i < aSize; ++i) { _data[i] = aSize - 1 - i; };
        break;
    }
  }

  // returns a permutation that has a single GIP at index aGipIdx
  static self_t constructUniqGip(const size_t aSize, const size_t aGipIdx) {
    // allow aGipIdx = aSize-1 although by def. GIP can only occurr at indices [0, aSize-2].
    self_t lRes(aSize);
    lRes.assertIdx(aGipIdx);

    value_t lValue = 0;
    // fill upper half
    for (size_t i = aGipIdx+1; i < lRes.size(); ++i) {
      lRes._data.at(i) = lValue;
      ++lValue;
    }
    // fill lower half
    for (size_t i = 0; i < aGipIdx+1; ++i) {
      lRes._data.at(i) = lValue;
      ++lValue;
    }
    assert(lRes.hasGipAt(aGipIdx));
    assert((lRes.getNoGip() == 1) || (aGipIdx == aSize-1));
    return lRes;
  }
 
  // TODO: This method does not work.
  //       It's not possible to construct such permutations for arbitrary
  //       aSize, aGipLb and aGipUb!
  // returns a permutation that has GIPs at all indices [aGipLb, aGipUb]
  static self_t constructGips(const size_t aSize, const size_t aGipLb, const size_t aGipUb) {
    // allow aGipIdx = aSize-1 although by def. GIP can only occurr at indices [0, aSize-2].
    assert(aGipUb >= aGipLb);
    const size_t lNumGips = aGipUb - aGipLb + 1;
    assert(lNumGips <= aSize);

    if (lNumGips == 1) { return constructUniqGip(aSize, aGipLb); }

    self_t lRes(aSize);
    lRes.assertIdx(aGipLb);
    lRes.assertIdx(aGipUb);

    // fill GIP interval with values decreasing from left to right, starting with the maximum
    value_t lValue = aSize - 1;
    for (size_t i = aGipLb; i <= aGipUb; ++i) {
      lRes._data.at(i) = lValue;
      --lValue;
    }
    // fill part left of the GIP interval with values decreasing from right to left
    if (0 < aGipLb) {
      for (int i = aGipLb-1; i >= 0; --i) {
        lRes._data.at(i) = lValue;
        --lValue;
      }
    }
    // fill part rigth of the GIP interval with values decreasing from right to left
    if (aGipUb < aSize - 1) {
      for (int i = aSize-1; i > aGipUb; --i) {
        lRes._data.at(i) = lValue;
        --lValue;
      }
    }
    assert(lValue == static_cast<value_t>(-1));
    
    return lRes;
  }

public:
  void invert() {
    value_vt lInv(_data.size());
    for (size_t i = 0; i < _data.size(); ++i) {
      lInv[_data[i]] = i;
    }
    _data = std::move(lInv);
  }

  self_t getInverse() const {
    self_t lInvPerm(*this);
    lInvPerm.invert();
    return lInvPerm;
  }

  // return f(this())
  self_t getComposition(const self_t f) const {
    assert(size() == f.size());
    self_t lResult(size());
    for (size_t i = 0; i < size(); ++i) {
      lResult[i] = f[_data[i]];
    }
    return lResult;
  }

  // is identity permutation
  bool isId() const {
    for (size_t i = 0; i < size(); ++i) {
      if (_data[i] != i) { return false; }
    }
    return true;
  }

  // is reverse identity permutation
  bool isRev() const {
    for (size_t i = 0; i < size(); ++i) {
      if (_data[i] != (size() - 1 - i)) { return false; }
    }
    return true;
  }

  // has common subset prefix (CSP) with another permutation.
  // permutations p1, p2 of length n have csp
  // :<=> there is an index i, 0 <= i < n-1,
  //      such that {p1(0), ..., p1(i)} = {p2(0), ..., p2(i)}.
  bool hasCsp(const self_t& aOtherPerm) const {
    uint64_t lRelsetSelf = 0;
    uint64_t lRelsetOther = 0;
    assert(size() == aOtherPerm.size());
    const size_t l = size() - 1;
    //std::cout << "  : size=" << size() << ", l=" << l << std::endl;
    for(size_t i = 0; i < l; ++i) {
      lRelsetSelf |= (1 << get(i));
      lRelsetOther |= (1 << aOtherPerm(i));
      //std::cout << "  :   " << i << " : " << lRelsetSelf << " - " << lRelsetOther << "\n";
      if(lRelsetSelf == lRelsetOther) {
        return true;
      }
    }
    return false;
  }

  // has global inversion point
  bool hasGip() const {
    for (size_t j = 0; j < size(); ++j) {
      if (hasGipAt(j)) {
        return true;
      }
    }
    return false;
  }

  // has a global inversion point (GIP) at aPos
  bool hasGipAt(const size_t aPos, const bool aTrace = false) const {
    if (aTrace) { std::cout << "    " << __PRETTY_FUNCTION__ << std::endl; }
    if (aTrace) { std::cout << "    aPos = " << aPos << std::endl; }

    const Tuint lMinBound = (size()-1) - (aPos+1) + 1;
    if (aTrace) { std::cout << "    lMinBound = " << lMinBound << std::endl; }

    for (size_t i = 0; i <= aPos; ++i) {
      if (get(i) >= lMinBound) {
        continue;
      } else {
        return false;
      }
    }
    return true;
  }

  // number of GIPs
  size_t getNoGip() const {
    size_t lRes = 0;
    for (size_t i = 0; i < size()-1; ++i) {
      lRes += hasGipAt(i);
    }
    return lRes;
  }

  bool isReducible() const {
    const self_t lId(size(), special_perm::id);  // id perm
    return hasCsp(lId);
  }

  bool isIrreducible() const { return !isReducible(); }

  // inflection point ip as defined by DSQ/GM:
  // Let pi \in Perm(z) be a permutation; let x \in bv(z) be a bitvector.
  // ip(x, pi) := max{ g | 0 <= g < z, pi{g} \subseteq x }
  // => Largest "overlap" of permutation and bitvector.
  size_t inflectionPoint(const uint64_t aBv) const {
    int64_t lRes = -1;
    uint64_t lSet = 0;
    for (size_t i = 0; i < size(); ++i) {
      lSet |= (1 << get(i));
      if (!(aBv == (lSet | aBv))) {  // !(lSet \subseteq aBv)
        lRes = i;
        break;
      }
    }
    if (0 > lRes) {
      lRes = size() - 1;
    }

    return lRes;
  }



public:
        size_t   size()                 const { return _data.size();  }
  const value_t& get (const size_t idx) const { return _data.at(idx); }
        value_t& get (const size_t idx)       { return _data.at(idx); }

public:
  const value_t& operator[](const size_t idx) const { return get(idx); }
        value_t& operator[](const size_t idx)       { return get(idx); }
  const value_t& operator()(const size_t idx) const { return get(idx); }
        value_t& operator()(const size_t idx)       { return get(idx); }

public:
  // conversion (cast) to std::vector<Tuint> a.k.a. value_vt
  operator value_vt() const { return _data; }

  // conversion to std::vector<Tnum> of different type
  template <typename Tnum> 
  std::vector<Tnum> toVec() const;
  //template <typename Tnum> 
  //inline std::vector<Tnum> toVec() const {
  //  //if (std::is_same<value_t, Tnum>::value) {  // true for uint64_t and uint?
  //  //  return _data;
  //  //} else {
  //    return std::vector<Tnum>(_data.cbegin(), _data.cend());
  //  //}
  //}


public:
  // Myrvold/Ruskey algorithm to rank a permutation
  // implementation taken from Guido's src/infra/permutation.hh
  // non-destructive version, i.e. no change to *this
  uint64_t rank() const {
    self_t lSelf(*this);
    self_t lInv = getInverse();
    uint64_t lRes = 0;
    uint64_t f = 1;
    const size_t& n = size();
    for (uint64_t i = n - 1; i > 0; --i) {
      const uint64_t s = lSelf(i);

      const uint64_t x = lSelf(i); // swap
      lSelf(n-1) = lSelf(lInv(i));
      lSelf(lInv(i)) = x;

      const uint64_t y = lInv(s); // swap
      lInv(s) = lInv(i);
      lInv(i) = y;

      lRes += (s * f);
      f *= (i+1);
    }
    return lRes;
  }

public:
  // to bitvector (uint)
  uint64_t toBv(const size_t aNumElements) const {
    assert(64 > size());
    uint64_t lRes = 0;
    for (size_t i = 0; i < aNumElements; ++i) {
      lRes |= (1 << get(i));
    }
    return lRes;
  }

  uint64_t toBv() const {
    return toBv(size());
  }


  std::string toString() const {
    std::stringstream lStream;
    lStream << *this;
    return lStream.str();
  }


private:
  void swapElements(const size_t i, const size_t j) {
    assertIdx(i);
    assertIdx(j);
    value_t tmp = _data[i];
    _data[i] = _data[j];
    _data[j] = tmp;
  }

  // Myrvold/Ruskey algorithm to unrank a permutation
  // implementation taken from Guido's src/infra/permutation.hh
  void unrank(const size_t aSize, uint64_t aRank) {
    init(aSize, special_perm::id);
    for(size_t i = aSize; i > 0; --i) {
      Tuint x = _data[i-1];
      _data[i-1] = _data[aRank % i];
      _data[aRank % i] = x;
      aRank = aRank / i;
    }
  }

private:
  bool isValidIdx(const size_t idx) const { return (0 <= idx) && (idx < _data.size()); }
  void assertIdx (const size_t idx) const { assert(isValidIdx(idx)); }


private:
  std::vector<Tuint> _data;

};


/* Operators */

template <typename Tuint>
inline std::ostream&
operator<<(std::ostream& os, const Permutation<Tuint>& aPermutation) {
  os << "(";
  if (aPermutation.size() > 0) {
    for (size_t i = 0; i < aPermutation.size()-1; ++i) {
      os << aPermutation[i] << " ";
    }
    os << aPermutation[aPermutation.size()-1];
  }
  os << ")";
  return os;
}


template <typename Tuint>
inline bool
operator==(const Permutation<Tuint>& aP1, const Permutation<Tuint>& aP2) {
  if (aP1.size() != aP2.size()) { return false; }
  for (size_t i = 0; i < aP1.size(); ++i) {
    if (aP1[i] != aP2[i]) { return false; }
  }
  return true;
}

template <typename Tuint>
inline bool
operator!=(const Permutation<Tuint>& aP1, const Permutation<Tuint>& aP2) {
  return !(aP1 == aP2);
}


// conversion to std::vector<Tnum> of different type
template <typename Tuint>
template <typename Tnum> 
std::vector<Tnum>
Permutation<Tuint>::toVec() const {
  return std::vector<Tnum>(_data.cbegin(), _data.cend());
}

// specialization of member does not work without specialization of the class
//template <typename Tuint>
//template <> 
//std::vector<Tuint>
//Permutation<Tuint>::toVec() const {
//  return _data;
//}

}
