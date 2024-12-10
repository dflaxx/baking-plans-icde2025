#include "CrystalTimer.hh"
#include <limits>

extern "C" {
  #include "cmeasure.h"
}

namespace df::infra {

double CrystalTimer::_crystal_frequency = -1;  // init value, freq not set

void
CrystalTimer::init() {
  _crystal_frequency = measure_frequency();
}

double
CrystalTimer::duration_ns() const {
  return duration_s() * static_cast<double>(1'000'000'000);
}
double
CrystalTimer::duration_us() const {
  return duration_s() * static_cast<double>(1'000'000);
}
double
CrystalTimer::duration_ms() const {
  return duration_s() * static_cast<double>(1'000);
}
double
CrystalTimer::duration_s() const {
  return static_cast<double>(duration_cycles()) / crystal_frequency();
}

// static methods

double
CrystalTimer::crystal_frequency() {
  if (0 > _crystal_frequency) {
    _crystal_frequency = measure_frequency();
  }
  return _crystal_frequency;
}


// from Guido, with changes
double
CrystalTimer::measure_frequency(size_t aNumIterations) {
  aNumIterations = std::max(aNumIterations, 1UL);  // at least one iteration
  double lFreq = 0;
  constexpr bool lTrace = false;
  constexpr double lMeasureTime = 1.0; // seconds

  for (size_t i = 0; i < aNumIterations; ++i) {
    cmeasure_t lMeas;
    cmeasure_start(&lMeas);
    const uint64_t lBegin = now_serial();
    do {
      cmeasure_stop(&lMeas);
    } while(lMeasureTime > cmeasure_total_s(&lMeas));
    const uint64_t lEnd = now_serial();
    const uint64_t lNoCyc = cycles(lBegin, lEnd);
    lFreq += ((static_cast<double>(lNoCyc)) / cmeasure_total_s(&lMeas));
    if(lTrace) {
      std::cout << "CrystalClock::measure_frequency:" << std::endl;
      std::cout << "   time " << cmeasure_total_s(&lMeas) << " [s]" << std::endl;
      if(lEnd < lBegin) {
        std::cout << "   wrap arround." << std::endl;
      }
      std::cout << "   cycles    :" << lNoCyc << std::endl;
      std::cout << "   frequency :" << lFreq   << std::endl;
    }
  }

  return lFreq / aNumIterations;  // mean of n iterations
}


}
