#ifndef MEASURE_EVAL_HH
#define MEASURE_EVAL_HH

#include "dfinfra/standard_includes.hh"

/*
 * Simple container of time measurements
 * 
 * Initialize to size using ctor.
 * Prepare with init(). Add values with step(). Postprocessing with fin().
 *
 * Statistical functions: min, max, sum, avg, median, 2nd largest value
 *
 * Author: Guido Moerkotte
 */
struct meas_eval_t {
  using uint64_vt = std::vector<uint64_t>;

  uint64_vt _meas;   // the values, added using step(...)
  uint64_t  _first;  // the (chronologially) first value, _meas[0] before sorting
  uint      _curr;   // next free element in _meas

  // ctor: n = vector size
  meas_eval_t(const uint n) : _meas(n), _first(0), _curr(0) { assert(2 <= _meas.size()); };

  inline void init() { _curr = 0; }
  // inline void step(const uint64_t x) { assert(_meas.size() > _curr); _meas[_curr++] = x; } // inline failed
  inline void step(const uint64_t x) { _meas[_curr++] = x; }
         void fin();

  inline uint64_t min()  const { return _meas[0]; }
  inline uint64_t max()  const { return _meas[_meas.size() - 1]; }
         uint64_t sum()  const;
         double   avg()  const;
         double   median() const;
  inline uint64_t max2()  const { return _meas[_meas.size() - 2]; }
  inline uint64_t first() const { return _first; }
  inline void     resize(const uint n) { _meas.resize(n); } // no assert here
  std::ostream& print(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream&, const meas_eval_t&);

#endif
