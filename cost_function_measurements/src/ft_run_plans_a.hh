#ifndef HASHJOIN2_FT_RUN_PLANS_A_HH
#define HASHJOIN2_FT_RUN_PLANS_A_HH
#pragma once

#include "dfinfra/standard_includes.hh"
#include "gminfra/prime_table_1_1.hh"

extern "C" {
  #include "dfinfra/cbind_to_hw_thread.h"
  #include "dfinfra/cmeasure.h"
}

#include "cb_glob.hh"
#include "line.hh"
#include "RelationRS.hh"
#include "plan.hh"

// more includes, thus include it but last in main

template<typename T>
void
ft_check_result_a(const T& X, const T& Y, const std::string& u, const std::string& v, const bool aPrint) {
  const bool lPrint = aPrint || (X != Y);
  if(lPrint) {
    std::cout << u << " vs " << v << std::endl;
    std::cout << u << ':' << std::endl;
    for(const auto& x : X) {
      std::cout << "  " << x._key_r << ',' << x._key_s << std::endl;
    }
    std::cout << v << ':' << std::endl;
    for(const auto& y : Y) {
      std::cout << "  " << y._key_r << ',' << y._key_s << std::endl;
    }
  }
  assert(X == Y && "Sets of result tuples of the two plans differ.");
}

template<typename T_rel_R, typename T_rel_S, bool Ttest, bool Ttrace>
void
ft_run_plans_a(const T_rel_R& R, 
               const T_rel_S& S,
               const size_t   aHashDirSizeR,
               const size_t   aHashDirSizeS,
               const uint     aLog2ChunkSize,
                     line_t&  aLine,
               const CbGlob&  aCbGlob) {

  static constexpr bool lManyPtrNo  = false;
  static constexpr bool lManyPtrYes = true;
  static constexpr bool lCh = false;
  static constexpr bool l3d = true;
  static constexpr prefetch_et lPrefetchOff  = prefetch_et::kPrefetchNone;
  static constexpr prefetch_et lPrefetchRp = prefetch_et::kRollingPrefetching;
  static constexpr prefetch_et lPrefetchAmac = prefetch_et::kAmac;

  using tuple_r_t = T_rel_R::tuple_t;
  using attr_t    = tuple_r_t::attr_t;
  using bun_t     = bun_tt<attr_t>;
  using hashfun_t = ht::HashMurmur<attr_t>;
  using R_t       = T_rel_R;
  using S_t       = T_rel_S;


  if constexpr (Ttrace) {
    std::cout << "=== run_measurement(" << sizeof(attr_t) << ") ===" << std::endl;
    print_RS_wh<bun_t, hashfun_t>(std::cout, R, S, R.card(), aHashDirSizeS) << std::endl;
  }

  using plan_ch_orig_rs_t = PlanJoin<R_t, S_t, kJoinOrderRS, lCh, lPrefetchOff,  lManyPtrNo, Ttrace, Ttest>;
  using plan_ch_orig_sr_t = PlanJoin<R_t, S_t, kJoinOrderSR, lCh, lPrefetchOff,  lManyPtrNo, Ttrace, Ttest>;
  using plan_ch_impr_rs_t = PlanJoin<R_t, S_t, kJoinOrderRS, lCh, lPrefetchOff,  lManyPtrYes, Ttrace, Ttest>;
  using plan_ch_impr_sr_t = PlanJoin<R_t, S_t, kJoinOrderSR, lCh, lPrefetchOff,  lManyPtrYes, Ttrace, Ttest>;

  using plan_3d_orig_rs_t = PlanJoin<R_t, S_t, kJoinOrderRS, l3d, lPrefetchOff,  lManyPtrNo, Ttrace, Ttest>;
  using plan_3d_orig_sr_t = PlanJoin<R_t, S_t, kJoinOrderSR, l3d, lPrefetchOff,  lManyPtrNo, Ttrace, Ttest>;
  using plan_3d_impr_rs_t = PlanJoin<R_t, S_t, kJoinOrderRS, l3d, lPrefetchOff,  lManyPtrYes, Ttrace, Ttest>;
  using plan_3d_impr_sr_t = PlanJoin<R_t, S_t, kJoinOrderSR, l3d, lPrefetchOff,  lManyPtrYes, Ttrace, Ttest>;

  using plan_ch_amac_rs_t = PlanJoin<R_t, S_t, kJoinOrderRS, lCh, lPrefetchAmac, lManyPtrNo, Ttrace, Ttest>;
  using plan_ch_amac_sr_t = PlanJoin<R_t, S_t, kJoinOrderSR, lCh, lPrefetchAmac, lManyPtrNo, Ttrace, Ttest>;
  using plan_ch_imac_rs_t = PlanJoin<R_t, S_t, kJoinOrderRS, lCh, lPrefetchAmac, lManyPtrYes, Ttrace, Ttest>;
  using plan_ch_imac_sr_t = PlanJoin<R_t, S_t, kJoinOrderSR, lCh, lPrefetchAmac, lManyPtrYes, Ttrace, Ttest>;

  using plan_3d_amac_rs_t = PlanJoin<R_t, S_t, kJoinOrderRS, l3d, lPrefetchAmac, lManyPtrNo, Ttrace, Ttest>;
  using plan_3d_amac_sr_t = PlanJoin<R_t, S_t, kJoinOrderSR, l3d, lPrefetchAmac, lManyPtrNo, Ttrace, Ttest>;
  using plan_3d_imac_rs_t = PlanJoin<R_t, S_t, kJoinOrderRS, l3d, lPrefetchAmac, lManyPtrYes, Ttrace, Ttest>;
  using plan_3d_imac_sr_t = PlanJoin<R_t, S_t, kJoinOrderSR, l3d, lPrefetchAmac, lManyPtrYes, Ttrace, Ttest>;


  plan_ch_orig_rs_t  lPlanChOrigRS;
  plan_ch_orig_sr_t  lPlanChOrigSR;
  plan_ch_impr_rs_t  lPlanChImprRS;
  plan_ch_impr_sr_t  lPlanChImprSR;
    
  plan_3d_orig_rs_t  lPlan3dOrigRS;
  plan_3d_orig_sr_t  lPlan3dOrigSR;
  plan_3d_impr_rs_t  lPlan3dImprRS;
  plan_3d_impr_sr_t  lPlan3dImprSR;
    
  plan_ch_amac_rs_t  lPlanChAmacRS;
  plan_ch_amac_sr_t  lPlanChAmacSR;
  plan_ch_imac_rs_t  lPlanChImacRS;
  plan_ch_imac_sr_t  lPlanChImacSR;
    
  plan_3d_amac_rs_t  lPlan3dAmacRS;
  plan_3d_amac_sr_t  lPlan3dAmacSR;
  plan_3d_imac_rs_t  lPlan3dImacRS;
  plan_3d_imac_sr_t  lPlan3dImacSR;


  static_assert(kPlanChOrigRS == lPlanChOrigRS.PLAN_ID);
  static_assert(kPlanChOrigSR == lPlanChOrigSR.PLAN_ID);
  static_assert(kPlanChImprRS == lPlanChImprRS.PLAN_ID);
  static_assert(kPlanChImprSR == lPlanChImprSR.PLAN_ID);
  static_assert(kPlanChAmacRS == lPlanChAmacRS.PLAN_ID);
  static_assert(kPlanChAmacSR == lPlanChAmacSR.PLAN_ID);
  static_assert(kPlanChImacRS == lPlanChImacRS.PLAN_ID);
  static_assert(kPlanChImacSR == lPlanChImacSR.PLAN_ID);

  static_assert(kPlan3dOrigRS == lPlan3dOrigRS.PLAN_ID);
  static_assert(kPlan3dOrigSR == lPlan3dOrigSR.PLAN_ID);
  static_assert(kPlan3dImprRS == lPlan3dImprRS.PLAN_ID);
  static_assert(kPlan3dImprSR == lPlan3dImprSR.PLAN_ID);
  static_assert(kPlan3dAmacRS == lPlan3dAmacRS.PLAN_ID);
  static_assert(kPlan3dAmacSR == lPlan3dAmacSR.PLAN_ID);
  static_assert(kPlan3dImacRS == lPlan3dImacRS.PLAN_ID);
  static_assert(kPlan3dImacSR == lPlan3dImacSR.PLAN_ID);


  lPlanChOrigRS.build_plan(aHashDirSizeS, aLog2ChunkSize);
  lPlanChOrigSR.build_plan(aHashDirSizeR, aLog2ChunkSize);
  lPlanChImprRS.build_plan(aHashDirSizeS, aLog2ChunkSize);
  lPlanChImprSR.build_plan(aHashDirSizeR, aLog2ChunkSize);

  lPlan3dOrigRS.build_plan(aHashDirSizeS, aLog2ChunkSize);
  lPlan3dOrigSR.build_plan(aHashDirSizeR, aLog2ChunkSize);
  lPlan3dImprRS.build_plan(aHashDirSizeS, aLog2ChunkSize);
  lPlan3dImprSR.build_plan(aHashDirSizeR, aLog2ChunkSize);

  lPlanChAmacRS.build_plan(aHashDirSizeS, aLog2ChunkSize);
  lPlanChAmacSR.build_plan(aHashDirSizeR, aLog2ChunkSize);
  lPlanChImacRS.build_plan(aHashDirSizeS, aLog2ChunkSize);
  lPlanChImacSR.build_plan(aHashDirSizeR, aLog2ChunkSize);

  lPlan3dAmacRS.build_plan(aHashDirSizeS, aLog2ChunkSize);
  lPlan3dAmacSR.build_plan(aHashDirSizeR, aLog2ChunkSize);
  lPlan3dImacRS.build_plan(aHashDirSizeS, aLog2ChunkSize);
  lPlan3dImacSR.build_plan(aHashDirSizeR, aLog2ChunkSize);

  if constexpr (Ttrace) {
    std::cout << "done building plans" << std::endl;
  }

  if constexpr (Ttest) {
    lPlanChOrigRS.run_once(R, S);
    lPlanChOrigSR.run_once(R, S);
    lPlanChImprRS.run_once(R, S);
    lPlanChImprSR.run_once(R, S);
    lPlan3dOrigRS.run_once(R, S);
    lPlan3dOrigSR.run_once(R, S);
    lPlan3dImprRS.run_once(R, S);
    lPlan3dImprSR.run_once(R, S);
    lPlanChAmacRS.run_once(R, S);
    lPlanChAmacSR.run_once(R, S);
    lPlanChImacRS.run_once(R, S);
    lPlanChImacSR.run_once(R, S);
    lPlan3dAmacRS.run_once(R, S);
    lPlan3dAmacSR.run_once(R, S);
    lPlan3dImacRS.run_once(R, S);
    lPlan3dImacSR.run_once(R, S);
    aLine._card_result = lPlan3dOrigSR.top_count();
  } else {
    lPlanChOrigRS.run_many(R, S, aLine);
    lPlanChOrigSR.run_many(R, S, aLine);
    lPlanChImprRS.run_many(R, S, aLine);
    lPlanChImprSR.run_many(R, S, aLine);
    lPlan3dOrigRS.run_many(R, S, aLine);
    lPlan3dOrigSR.run_many(R, S, aLine);
    lPlan3dImprRS.run_many(R, S, aLine);
    lPlan3dImprSR.run_many(R, S, aLine);
    lPlanChAmacRS.run_many(R, S, aLine);
    lPlanChAmacSR.run_many(R, S, aLine);
    lPlanChImacRS.run_many(R, S, aLine);
    lPlanChImacSR.run_many(R, S, aLine);
    lPlan3dAmacRS.run_many(R, S, aLine);
    lPlan3dAmacSR.run_many(R, S, aLine);
    lPlan3dImacRS.run_many(R, S, aLine);
    lPlan3dImacSR.run_many(R, S, aLine);
    aLine._card_result = lPlan3dOrigSR.top_count();
  }

  if constexpr (Ttest) {
    std::cout << "|R join^ch.o S| = " << lPlanChOrigRS.result().size() << std::endl;
    std::cout << "|S join^ch.o R| = " << lPlanChOrigSR.result().size() << std::endl;
    std::cout << "|R join^ch.i S| = " << lPlanChImprRS.result().size() << std::endl;
    std::cout << "|S join^ch.i R| = " << lPlanChImprSR.result().size() << std::endl;

    std::cout << "|R join^3d.o S| = " << lPlan3dOrigRS.result().size() << std::endl;
    std::cout << "|S join^3d.o R| = " << lPlan3dOrigSR.result().size() << std::endl;
    std::cout << "|R join^3d.i S| = " << lPlan3dImprRS.result().size() << std::endl;
    std::cout << "|S join^3d.i R| = " << lPlan3dImprSR.result().size() << std::endl;

    std::cout << "|R join^ch.a S| = " << lPlanChAmacRS.result().size() << std::endl;
    std::cout << "|S join^ch.a R| = " << lPlanChAmacSR.result().size() << std::endl;
    std::cout << "|R join^ch.A S| = " << lPlanChImacRS.result().size() << std::endl;
    std::cout << "|S join^ch.A R| = " << lPlanChImacSR.result().size() << std::endl;

    std::cout << "|R join^3d.a S| = " << lPlan3dAmacRS.result().size() << std::endl;
    std::cout << "|S join^3d.a R| = " << lPlan3dAmacSR.result().size() << std::endl;
    std::cout << "|R join^3d.A S| = " << lPlan3dImacRS.result().size() << std::endl;
    std::cout << "|S join^3d.A R| = " << lPlan3dImacSR.result().size() << std::endl;

    ft_check_result_a(lPlanChOrigRS.result(), lPlanChOrigSR.result(), "lPlanChOrigRS", "lPlanChOrigSR", false);
    ft_check_result_a(lPlanChOrigRS.result(), lPlanChImprRS.result(), "lPlanChOrigRS", "lPlanChImprRS", false);
    ft_check_result_a(lPlanChOrigRS.result(), lPlanChImprSR.result(), "lPlanChOrigRS", "lPlanChImprSR", false);

    ft_check_result_a(lPlanChOrigRS.result(), lPlanChAmacRS.result(), "lPlanChOrigRS", "lPlanChAmacRS", false);
    ft_check_result_a(lPlanChOrigRS.result(), lPlanChAmacSR.result(), "lPlanChOrigRS", "lPlanChAmacSR", false);
    ft_check_result_a(lPlanChOrigRS.result(), lPlanChImacRS.result(), "lPlanChOrigRS", "lPlanChImacRS", false);
    ft_check_result_a(lPlanChOrigRS.result(), lPlanChImacSR.result(), "lPlanChOrigRS", "lPlanChImacSR", false);

    ft_check_result_a(lPlanChOrigRS.result(), lPlan3dOrigRS.result(), "lPlanChOrigRS", "lPlan3dOrigRS", false);
    ft_check_result_a(lPlanChOrigRS.result(), lPlan3dOrigSR.result(), "lPlanChOrigRS", "lPlan3dOrigSR", false);
    ft_check_result_a(lPlanChOrigRS.result(), lPlan3dImprRS.result(), "lPlanChOrigRS", "lPlan3dImprRS", false);
    ft_check_result_a(lPlanChOrigRS.result(), lPlan3dImprSR.result(), "lPlanChOrigRS", "lPlan3dImprSR", false);

    ft_check_result_a(lPlanChOrigRS.result(), lPlan3dAmacRS.result(), "lPlanChOrigRS", "lPlan3dAmacRS", false);
    ft_check_result_a(lPlanChOrigRS.result(), lPlan3dAmacSR.result(), "lPlanChOrigRS", "lPlan3dAmacSR", false);
    ft_check_result_a(lPlanChOrigRS.result(), lPlan3dImacRS.result(), "lPlanChOrigRS", "lPlan3dImacRS", false);
    ft_check_result_a(lPlanChOrigRS.result(), lPlan3dImacSR.result(), "lPlanChOrigRS", "lPlan3dImacSR", false);
  }

  assert(lPlanChOrigRS.top_count() == lPlanChOrigSR.top_count());
  assert(lPlanChOrigRS.top_count() == lPlanChImprRS.top_count());
  assert(lPlanChOrigRS.top_count() == lPlanChImprSR.top_count());
  assert(lPlanChOrigRS.top_count() == lPlanChAmacRS.top_count());
  assert(lPlanChOrigRS.top_count() == lPlanChAmacSR.top_count());
  assert(lPlanChOrigRS.top_count() == lPlanChImacRS.top_count());
  assert(lPlanChOrigRS.top_count() == lPlanChImacSR.top_count());
  assert(lPlanChOrigRS.top_count() == lPlan3dOrigRS.top_count());
  assert(lPlanChOrigRS.top_count() == lPlan3dOrigSR.top_count());
  assert(lPlanChOrigRS.top_count() == lPlan3dImprRS.top_count());
  assert(lPlanChOrigRS.top_count() == lPlan3dImprSR.top_count());
  assert(lPlanChOrigRS.top_count() == lPlan3dAmacRS.top_count());
  assert(lPlanChOrigRS.top_count() == lPlan3dAmacSR.top_count());
  assert(lPlanChOrigRS.top_count() == lPlan3dImacRS.top_count());
  assert(lPlanChOrigRS.top_count() == lPlan3dImacSR.top_count());
}

#endif
