#include "standard_includes.hh"

#include <algorithm>
#include <concepts>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <ranges>
#include <stdexcept>

template <typename T>
concept Numeric = std::integral<T> || std::floating_point<T>;



/*
 * A simple statistics accumulator.
 * Can add numeric values and return measures like count, average, min, max.
 * Can optionally materialize the values to compute quantiles.
 */
template <Numeric T>
class SimpleStat {
  public:
    using value_t = T;

  public:
    inline SimpleStat() = default;
    inline SimpleStat(const bool aMaterialize) : _materialize(aMaterialize) {}

  public:
    inline void add(const T x);
           void add(const std::initializer_list<T>);

    // always available
    inline T min() const { return _min; }
    inline T max() const { return _max; }

    inline size_t count() const { return _count; }
    inline T sum() const { return _sum; }
    double avg() const;

    // requires _materialize = true
    inline double median() { return quantile_lazy(0.5); }
    double quantile(const double p);  // enforces sort()
    double quantile_lazy(const double p);  // does not enforce total order
    // p <= 0.0 -> min(). p >= 1.0 -> max().

  public:
    inline bool maintain_sorting() const { return _maintain_sorting; }
    inline void maintain_sorting(const bool b) { _maintain_sorting = b; }

    bool sorted() const { return _sorted; }
    bool check_sorted() const { return std::is_sorted(_values.cbegin(), _values.cend()); }

    void print_values(std::ostream& = std::cout) const;

  private:
    void sort();

    size_t get_quantile_idx_upper(const double p) const;
    size_t get_quantile_idx_upper(const size_t num, const size_t denom) const;

    bool is_quantile_from_two(const double p) const;  // is the requested quantile a single value
                                                      // or the mean of two values?

    static inline bool
    is_integer(const double d) { double dummy; return (std::modf(d, &dummy) == 0.0); }

  private:
    // aggregate values
    T _min {std::numeric_limits<T>::max()};
    T _max {std::numeric_limits<T>::min()};
    T _sum {0};
    size_t _count {0};

    // materialization of values if requested by _materialize
    std::vector<T> _values {};

    const bool _materialize {false}; // indicate if values should be materialized
    bool _maintain_sorting {false};  // if _values are sorted, maintain this sorting when adding
    bool _sorted {false};            // internal flag to indicate if _values are sorted currently
};

/* public member function */

template <Numeric T>
void
SimpleStat<T>::add(const T x) {
  _count++;
  _sum += x;
  if (x < _min) { _min = x; }
  if (x > _max) { _max = x; }

  if (_materialize) {
    if (count() == 1 && _maintain_sorting) {
      _sorted = true;
    }
    if (_sorted && _maintain_sorting) {
      _values.insert(std::upper_bound(_values.begin(), _values.end(), x), x);
    } else {
      _values.push_back(x);
      _sorted = false;
    }
  }
}

template <Numeric T>
void
SimpleStat<T>::add(const std::initializer_list<T> aValues) {
  for (const auto& elem : aValues) {
    add(elem);
  }
}

template <Numeric T>
double
SimpleStat<T>::avg() const {
  assert(count() > 0);
  if (count() == 0) {
    throw std::runtime_error("SimpleStat::avg(): Error: Division by 0. count() is 0.");
  }
  return (sum() / static_cast<double>(count()));
}

template <Numeric T>
double
SimpleStat<T>::quantile(const double p) {
  assert(_materialize && count() > 0);
  if (!_materialize || count() == 0) {
    throw std::runtime_error("SimpleStat::quantile(): _materialize is false or count() is 0.");
  }

  if (p <= 0.0) { return min(); }
  if (p >= 1.0) { return max(); }

  if (!_sorted) { sort(); }
  assert(check_sorted());

  size_t lIdx = get_quantile_idx_upper(p);
  if (is_quantile_from_two(p)) {
    return (_values[lIdx - 1] + _values[lIdx]) / 2.0;
  } else {
    return _values[lIdx];
  }
}

template <Numeric T>
double
SimpleStat<T>::quantile_lazy(const double p) {
  assert(_materialize && count() > 0);
  if (!_materialize || count() == 0) {
    throw std::runtime_error("SimpleStat::quantile(): _materialize is false or count() is 0.");
  }

  if (p <= 0.0) { return min(); }
  if (p >= 1.0) { return max(); }

  size_t lIdx0 = get_quantile_idx_upper(p);
  std::nth_element(_values.begin(), _values.begin() + lIdx0, _values.end());
  if (is_quantile_from_two(p)) {
    size_t lIdx1 = lIdx0 - 1;
    std::nth_element(_values.begin(), _values.begin() + lIdx1, _values.end());
    return (_values[lIdx0] + _values[lIdx1]) / 2.0;
  } else {
    return _values[lIdx0];
  }
}

template <Numeric T>
void
SimpleStat<T>::print_values(std::ostream& os) const {
  os << "[";
    if (_values.size() > 0) {
    for (const auto& elem : _values | std::views::take(_values.size() - 1)) {
      os << elem << ", ";
    }
    os << _values.back();
  }
  os << "]" << std::endl;;
}


/* private member function */

template <Numeric T>
void
SimpleStat<T>::sort() {
  std::sort(_values.begin(), _values.end());
  _sorted = true;
}

template <Numeric T>
size_t
SimpleStat<T>::get_quantile_idx_upper(const double p) const {
  return count() * p;
}

template <Numeric T>
size_t
SimpleStat<T>::get_quantile_idx_upper(const size_t num, const size_t denom) const {
  return (count() / static_cast<double>(denom)) * num;
}

template <Numeric T>
bool
SimpleStat<T>::is_quantile_from_two(const double p) const {
  // If number_of_values * p is an integer,
  // then the quantile is the mean of two values.
  // Otherwise, its a single value.
  return is_integer(count() * p);
}
