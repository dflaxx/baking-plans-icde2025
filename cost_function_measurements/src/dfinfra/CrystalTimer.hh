#pragma once
/*
 * Wrapper for measuring the passage of time in terms of processor clock ticks.
 * Based on the rdtsc instruction.
 * Clock ticks can be converted to time using the (measured) frequency of the crystal.
 *
 * - x86_64: __rdtsc() intrinsic
 * - ARMv7:  no intrinsic, not implemented
 * - ARMv8:  no intrinsic, using inline asm
 *
 */

#include "standard_includes.hh"
#include <chrono>

#ifdef __x86_64__
  #include <x86intrin.h>
#elif (__ARM_ARCH == 8)
  // supported
#else 
  #error CrystalTimer: Architecture not supported.
#endif 

namespace df::infra {

class CrystalTimer {
  public:
    CrystalTimer() = default;
    void init();  // force measurement of crystal frequency
  private:
    using time_point = uint64_t;
  public:
    static inline time_point now();
    static inline time_point now_serial();
           inline void       start()        { _beg = now(); }
           inline void       stop()         { _end = now(); }
           inline void       start_serial() { _beg = now_serial(); }
           inline void       stop_serial()  { _end = now_serial(); }
  public:
    inline time_point begin() const { return _beg; }
    inline time_point end()   const { return _end; }
    inline bool       wraps() const { return end() < begin(); }
  public:
    inline uint64_t duration_cycles() const;  // num of cycles between begin and end, allows one wraparound
    inline uint64_t cycles() const { return duration_cycles(); }  // compatibility with Guido

    // returns duration as std::chrono::duration, defaults to nanoseconds
    template<class Rep = int64_t, class Period = std::nano>
    std::chrono::duration<Rep, Period> duration() const;

    double duration_ns() const;
    double duration_us() const;
    double duration_ms() const;
    double duration_s() const;
    
  public:
     static double crystal_frequency();
  private:
    inline static uint64_t cycles(const time_point aBegin, const time_point aEnd);

    // always uses now_serial() for frequency measurement
           static double   measure_frequency(const size_t aNumIterations);
    inline static double   measure_frequency() { return measure_frequency(1); }
  private:
    time_point    _beg, _end;
    static double _crystal_frequency;
    
};


// from Guido
inline CrystalTimer::time_point
CrystalTimer::now() {
  uint64_t lCyc = 0;
  #ifdef __x86_64__
    lCyc =  __rdtsc();
  #elif (__ARM_ARCH == 8)
    asm volatile("mrs %0, cntvct_el0" : "=r" (lCyc));
  #endif
  //#elif (__ARM_ARCH == 7)
  //  lCyc = pm_get_cycle_counter(); // ACHTUNG: hat nur 32 bit
  //#endif
  return lCyc; 
}

inline CrystalTimer::time_point
CrystalTimer::now_serial() {
  uint64_t lCyc = 0;
  #ifdef __x86_64__
    _mm_lfence();
    lCyc =  __rdtsc();
    _mm_lfence();
  #endif
  return lCyc; 
}

inline uint64_t
CrystalTimer::duration_cycles() const {
  return cycles(_beg, _end);
}

template<class Rep, class Period>
std::chrono::duration<Rep, Period>
CrystalTimer::duration() const {
  const uint64_t lNumCycle = duration_cycles();
  const double lFreq = crystal_frequency();
  const double lSeconds = lNumCycle / lFreq;
  const auto lDuration = std::chrono::duration<double, std::ratio<1>>(lSeconds);
  return std::chrono::duration_cast<std::chrono::duration<Rep, Period>>(lDuration);
}

inline uint64_t
CrystalTimer::cycles(const time_point aBegin, const time_point aEnd) {
  if (aBegin > aEnd) {
    return std::numeric_limits<time_point>::max() - aBegin + aEnd;
  }
  return aEnd - aBegin;
}

}
