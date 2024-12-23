#ifndef INFRA_VECTOR_OPERATORS_HH
#define INFRA_VECTOR_OPERATORS_HH

#include <iostream>
#include <iomanip>
#include <vector>
#include <assert.h>
#include "float_ops.hh"


template<class T>
std::vector<T>&
operator+=(std::vector<T>& x, const std::vector<T>& y) {
  x.insert(x.end(), y.begin(), y.end());
  return x;
}

template<class T>
void 
concat(std::vector<T>& x, const std::vector<T>& y, const std::vector<T>& z) {
  x.resize(y.size() + z.size());
  x = y;
  x += z;
}

template<typename T>
std::ostream&
operator<<(std::ostream& os, const std::vector<T>& x) {
  const int lWidth = os.width();
  // os << std::setw(1) << '[';
  for(size_t i = 0; i < x.size(); ++i) {
    os << std::setw(lWidth) << x[i];
    if((i + 1) < x.size()) {
      os << ' ';
    }
  }
  // os << std::setw(1) << ']';
  return os;
}

template<class T>
std::istream&
operator>>(std::istream& aIs, std::vector<T>& x) {
  for(size_t i = 0; i < x.size(); ++i) {
    aIs >> x[i];
  }
  return aIs;
}

template<class T>
bool
is_prefix(const std::vector<T>& aPrefix, const std::vector<T>& aBig) {
  if(aPrefix.size() > aBig.size()) {
    return false;
  }
  for(size_t i = 0; i < aPrefix.size(); ++i) {
    if(aPrefix[i] != aBig[i]) {
      return false;
    }
  }
  return true;
}

template<class T>
void
subvector_idx(std::vector<T>& aResult, const std::vector<T> aVec, const std::vector<unsigned int>& aIdx) {
  aResult.resize(aIdx.size());
  for(size_t i = 0; i < aIdx.size(); ++i) {
    aResult[i] = aVec[aIdx[i]];
  }
}

// intersection of sorted vectors
template<class T>
void
vec_isec_srt(std::vector<T>& aResult, const std::vector<T> aVec1, const std::vector<T> aVec2) {
  typename std::vector<T>::const_iterator i1 = aVec1.begin();
  typename std::vector<T>::const_iterator i2 = aVec2.begin();
  aResult.clear();
  while((aVec1.end() != i1) && (aVec2.end() != i2)) {
    if( ((*i1) < (*i2)) ) {
      ++i1;
    } else
    if( ((*i1) > (*i2)) ) {
      ++i2;
    } else {
      aResult.push_back(*i1);
      ++i1;
      ++i2;
    }
  }
}

template<class T>
bool
is_sorted_uniq(const std::vector<T>& v) {
  for(size_t i = 1; i < v.size(); ++i) {
    if(v[i-1] >= v[i]) {
      return false;
    }
  }
  return true;
}

template<class T>
void
elim_one(std::vector<T>& aResult, const T& aElem, const std::vector<T>& v) {
  aResult.clear();
  for(const auto& x : v) {
    if(x != aElem) {
      aResult.push_back(x);
    }
  }
}

template<class T>
T
vec_min(const std::vector<T>& aIn) {
  assert(0 < aIn.size());
  T lMin = aIn[0];
  for(size_t i = 1; i < aIn.size(); ++i) {
    if(aIn[i] < lMin) {
      lMin = aIn[i];
    }
  }
  return lMin;
}

template<class T>
T
vec_max(const std::vector<T>& aIn) {
  assert(0 < aIn.size());
  T lMax = aIn[0];
  for(size_t i = 1; i < aIn.size(); ++i) {
    if(aIn[i] > lMax) {
      lMax = aIn[i];
    }
  }
  return lMax;
}

template<class T>
double
vec_avg(const std::vector<T>& aIn) {
  assert(0 < aIn.size());
  double lSum = 0;
  for(size_t i = 0; i < aIn.size(); ++i) {
    lSum += (double) aIn[i];
  }
  return lSum / ((double) aIn.size());
}

template<class T>
bool
vec_all_equal(const std::vector<T>& aIn, const double aValue, const double aEps) {
  for(uint i = 0; i < aIn.size(); ++i) {
    if(!is_equal(aIn[i], aValue, aEps)) {
      return false;
    }
  }
  return true;
}

template<class T>
bool
vec_is_subset(const std::vector<T>& aSub, const std::vector<T>& aSup) {
  if(aSub.size() > aSup.size()) { return false; }
  for(const auto& x : aSub) {
    if(aSup.end() == std::find(aSup.begin(), aSup.end(), x)) {
      return false;}
  }
  return true;
}



#endif
