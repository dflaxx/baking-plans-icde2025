#include "standard_includes.hh"

#include <chrono>
#include <functional>

#include "chrono_helpers.hh"

namespace df::infra {

/*
 * Repeat a function aFunc until its cumulative execution time is at least aMinTime.
 * Return the cumulative execution time and the number of repetitions.
 * aTeardown is a std::function that is executed after each repetition
 * to prepare for the next, e.g. to clear data structures.
 *
 * If aMinTime is not reached or exceeded after n repetitions,
 * the measurement is restarted with 2*n repetitions.
 *
 * aTeardownBlock is used to teardown after one round of repetitions,
 * not after each individual repetition.
 *
 * XXX aTeardown() is part of the measurement!
 */
inline
std::pair<std::chrono::nanoseconds, size_t>
repeat_mintime(const std::chrono::nanoseconds aMinTime,
               const std::function<void(void)> aFunc,
               const std::function<void(void)> aTeardown = [](){},
               const bool aTeardownAfterLast = false,
               const size_t aMinRepeat = 1,
               const std::function<void(void)> aTeardownBlock = [](){}) {
  using namespace std::chrono_literals;
  using clock_t = std::chrono::high_resolution_clock;

  constexpr bool lTrace = false;

  size_t n = aMinRepeat;
  std::chrono::time_point<clock_t> t0, t1;
  std::chrono::nanoseconds lTotalTime = 0s;
  
  do {
    t0 = clock_t::now();
    for (size_t i = 0; i < n; ++i) {
      aFunc();
      if (i != n-1) { aTeardown(); }  // don't tear down last iteration
    }
    t1 = clock_t::now();
    lTotalTime = t1 - t0;
    
    // mintime not yet reached, must continue
    if (lTotalTime < aMinTime) {
      if constexpr (lTrace) {
        std::cout << __FUNCTION__ << ": mintime not yet reached, " << n << " iterations not sufficient." << std::endl;
      }
      n *= 2;
      aTeardown();  // teardown last iteration of for-loop
      aTeardownBlock();
    }
  } while (lTotalTime < aMinTime);

  if (aTeardownAfterLast) {
    aTeardown();
  }

  return {lTotalTime, n};
}

};
