#include "meas_eval.hh"
#include <algorithm>

void
meas_eval_t::fin() {
  _first = _meas[0];
  std::sort(_meas.begin(), _meas.end());
}

uint64_t
meas_eval_t::sum() const {
  const uint n = _meas.size();
  double lRes = 0;
  for(uint i = 0; i < n; ++i) {
    lRes += _meas[i];
  }
  return lRes;
}

double
meas_eval_t::avg() const {
  const double n = _meas.size();
  const double lSum = sum();
  return (lSum / n);
}

double
meas_eval_t::median() const {
  const uint n = _meas.size();
  if(0 == (n & 0x1)) {
    return (((double) _meas[(n >> 1) - 1] + (double) _meas[n >> 1]) / double{2});
  }
  return _meas[n >> 1];
}

std::ostream& 
meas_eval_t::print(std::ostream& os) const {
  os << '[';
  for(uint i = 0; i < _meas.size(); ++i) {
    os << ' ' << _meas[i];
  }
  os << " ; " << min() << ',' << max() << ',' << avg() << ',' << median() << ',' << max2() << ',' << first() ;
  os << " ]";
  return os;
}

std::ostream& 
operator<<(std::ostream& os, const meas_eval_t& x) {
 return x.print(os);
}

