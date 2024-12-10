#pragma once

#include "dfinfra/standard_includes.hh"

#include <random>

enum class algop_et {
  kAlgOpInvalid                    =  0,
  kAlgOpScan                       =  1,
  kAlgOpScanMat                    =  2,
  kAlgOpMat                        =  3,
  kAlgOpTop                        =  4,
  kAlgOpJoinChOrigBuild            =  5,
  kAlgOpJoinChOrigProbe            =  6,
  kAlgOpJoinChRpBuild              =  7,
  kAlgOpJoinChRpProbe              =  8,
  kAlgOpJoinChAmacBuild            =  9,
  kAlgOpJoinChAmacProbe            = 10,
  //
  kAlgOpJoin3dOrigBuild            = 11,
  kAlgOpJoin3dOrigProbeWithUnnest  = 12,
  kAlgOpJoin3dOrigProbeNoUnnest    = 13,
  kAlgOpJoin3dOrigUnnest           = 14,
  kAlgOpJoin3dRpBuild              = 15,
  kAlgOpJoin3dRpProbeWithUnnest    = 16,
  kAlgOpJoin3dRpProbeNoUnnest      = 17,
  kAlgOpJoin3dRpUnnest             = 18,
  kAlgOpJoin3dAmacBuild            = 19,
  kAlgOpJoin3dAmacProbeWithUnnest  = 20,
  kAlgOpJoin3dAmacProbeNohUnnest   = 21,
  kAlgOpJoin3dAmacUnnest           = 22,
  //
  kNoAlgOp                         = 23
};

// ----

enum meas_todo_et {
  kMeasureBuildBUN = 0,
  kMeasureBuildBNU = 1,
  kMeasureProbeBUN = 2,
  kMeasureProbeBNU = 3,
  kNoMeasTodo      = 4
};

// prefetch kind
enum prefetch_et {
  kPrefetchNone       = 0,
  kRollingPrefetching = 1,
  kAmac               = 2,
  kNumPrefetch        = 3
};

// TODO: fix and unify naming scheme
// TODO: align with header naming in line.hh
enum plan_et {
  kPlanChOrigRS  =  0,
  kPlanChOrigSR  =  1,
  kPlanChImprRS  =  2,
  kPlanChImprSR  =  3,

  kPlanChRpUpkRS =  4, // neu
  kPlanChRpUpkSR =  5, // neu
  kPlanChRpPkdRS =  6, // neu
  kPlanChRpPkdSR =  7, // neu

  kPlanChAmacRS  =  8,
  kPlanChAmacSR  =  9,
  kPlanChImacRS  = 10,
  kPlanChImacSR  = 11,
  //
  kPlan3dOrigRS  = 12,
  kPlan3dOrigSR  = 13,
  kPlan3dImprRS  = 14,
  kPlan3dImprSR  = 15,

  kPlan3dRpUpkRS = 16, // neu
  kPlan3dRpUpkSR = 17, // neu
  kPlan3dRpPkdRS = 18, // neu
  kPlan3dRpPkdSR = 19, // neu

  kPlan3dAmacRS  = 20,
  kPlan3dAmacSR  = 21,
  kPlan3dImacRS  = 22,
  kPlan3dImacSR  = 23,
  kNoPlan        = 24
};

inline std::string
to_string(const plan_et e) {
  switch (e) {
    case kPlanChOrigRS:  return "kPlanChOrigRS";  break;
    case kPlanChOrigSR:  return "kPlanChOrigSR";  break;
    case kPlanChImprRS:  return "kPlanChImprRS";  break;
    case kPlanChImprSR:  return "kPlanChImprSR";  break;
    case kPlanChRpUpkRS: return "kPlanChRpUpkRS"; break;
    case kPlanChRpUpkSR: return "kPlanChRpUpkSR"; break;
    case kPlanChRpPkdRS: return "kPlanChRpPkdRS"; break;
    case kPlanChRpPkdSR: return "kPlanChRpPkdSR"; break;
    case kPlanChAmacRS:  return "kPlanChAmacRS";  break;
    case kPlanChAmacSR:  return "kPlanChAmacSR";  break;
    case kPlanChImacRS:  return "kPlanChImacRS";  break;
    case kPlanChImacSR:  return "kPlanChImacSR";  break;
    case kPlan3dOrigRS:  return "kPlan3dOrigRS";  break;
    case kPlan3dOrigSR:  return "kPlan3dOrigSR";  break;
    case kPlan3dImprRS:  return "kPlan3dImprRS";  break;
    case kPlan3dImprSR:  return "kPlan3dImprSR";  break;
    case kPlan3dRpUpkRS: return "kPlan3dRpUpkRS"; break;
    case kPlan3dRpUpkSR: return "kPlan3dRpUpkSR"; break;
    case kPlan3dRpPkdRS: return "kPlan3dRpPkdRS"; break;
    case kPlan3dRpPkdSR: return "kPlan3dRpPkdSR"; break;
    case kPlan3dAmacRS:  return "kPlan3dAmacRS";  break;
    case kPlan3dAmacSR:  return "kPlan3dAmacSR";  break;
    case kPlan3dImacRS:  return "kPlan3dImacRS";  break;
    case kPlan3dImacSR:  return "kPlan3dImacSR";  break;
    case kNoPlan:        return "kNoPlan";        break;
    default: assert(0 ==1); return "to_string(const plan_et e) -- this should not happen";
  }
}

using planset_t = uint64_t;
inline bool
rs_contains(const uint aPlanNo, const planset_t P) {
  return (0 != ((P >> aPlanNo) & 0x1));
}


enum join_order_et {
  kJoinOrderRS = 0, // R probe S build
  kJoinOrderSR = 1, // S probe R build
  kNoJoinOrder = 2
};

enum join_impl_et {
  kJoinImplInvalid = 0,
  kJoinImplChOrig  = 1,
  kJoinImplChAmac  = 2,
  kJoinImpl3dOrig  = 3,
  kJoinImpl3dAmac  = 4,
  kNoJoinImpl      = 5
};

using rng32_t = std::mt19937;
using rng64_t = std::mt19937_64;
          
template<typename Tuint> 
struct rng_tt;
          
template<>      
struct rng_tt<uint32_t> {
  using result_type = uint32_t;
  constexpr static auto min() { return rng32_t::min(); }
  constexpr static auto max() { return rng32_t::max(); }
  rng_tt() : _rng() {}
  uint32_t operator()() { return _rng(); }
  rng32_t _rng;
};

template<>
struct rng_tt<uint64_t> {
  using result_type = uint64_t;
  constexpr static auto min() { return rng64_t::min(); }
  constexpr static auto max() { return rng64_t::max(); }
  rng_tt() : _rng() {}
  uint64_t operator()() { return _rng(); }
  rng64_t _rng;
};
