#pragma once
#include "standard_includes.hh"

#include <cmath>

#include <concepts>
#include <limits>

/*
 * Some math functions
 */

namespace df::infra {

// factorial of an unsigned integer
template <std::unsigned_integral Tuint>
Tuint factorial(const Tuint n) {
  Tuint lRes = 1;
  for (Tuint i = 2; i <= n; ++i) {
    Tuint lNewRes = lRes * i;
    if (lRes != 0 && lNewRes / lRes != i) {
      std::cout << __PRETTY_FUNCTION__ << ": Warning. Overflow." << std::endl;
    }
    lRes = lNewRes;
  }
  return lRes;
}

// binomial coefficient "n over k" for unsigned integers n, k
// Source: https://www.codewhoop.com/numbers/binomial-coefficient.html
template <std::unsigned_integral Tuint>
Tuint binomial(const Tuint n, Tuint k) {
  assert(n >= k);
  Tuint lRes = 1;  // binom(n, 0) = binom(n, n) = 1
  if (k > (n - k)) {
    k = n - k;     // binom(n, k) = binom(n, n-k)  (symmetry)
  }
  for (size_t i = 0; i < k; ++i) {
    lRes *= (n - i);  // numerator: n * (n-1) * ... * (n-k+1)
    lRes /= (i + 1);  // denominator: 1 * 2 * ... * k
  }
  return lRes;
}

/*
 * given insigned integer n >= 0 in base 10, return number of digits needed to represent in base b.
 * - approach with while and division: https://stackoverflow.com/a/1489873
 * - approach with log: https://www.geeksforgeeks.org/given-number-n-decimal-base-find-number-digits-base-base-b/
 */
template <std::unsigned_integral Tuint>
size_t number_of_digits(Tuint n, const Tuint b = 10) {
  if (0 == n) { return 1; }
  // log has precision issues for large n, which result in wrong number of digits
  //return static_cast<size_t>(std::floor(std::log(n) / std::log(b)) + 1);
  size_t digits = 0;
  while (n) {
    n /= b;
    ++digits;
  }
  return digits;
}

// given an unsigned integer n >= 0, check if n is a power of b
// https://codereview.stackexchange.com/a/117204
template <std::unsigned_integral Tuint>
bool is_power_of(Tuint n, const Tuint b = 10) {
  while (n >= b && n % b == 0) {
    n /= b;
  }
  return (n == 1);
}

inline double
qerr(const double x_true, const double x_est) {
  return std::max(x_true / x_est, x_est / x_true);
}

}
