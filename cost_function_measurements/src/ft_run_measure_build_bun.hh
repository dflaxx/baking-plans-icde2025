#pragma once

#include "dfinfra/standard_includes.hh"

#include "dfinfra/CrystalTimer.hh"

#include "gminfra/prime_table_1_1.hh"
#include "gminfra/tmath.hh"
#include "gminfra/bit_subsets.hh"

extern "C" {
  #include "dfinfra/cbind_to_hw_thread.h"
  #include "dfinfra/cmeasure.h"
}

#include "RelationRS.hh"
#include "cb_glob.hh"
#include "line.hh"
#include "plan_container.hh"

// more includes
// thus include last in main

/*
 * Measure hash table build for unique build side.
 *
 * Template parameters:
 * - Tattr: attribute type, works for uint64_t or uint32_t
 * - Ttest: boolean that directs if the join result is tested for correctness
 * - Ttrace: boolean, print tracing output
 */
template<typename Tattr, bool Ttest, bool Ttrace>
void
ft_run_measure_build_bun(const CbGlob& aCbGlob) {
  using attr_t = Tattr;
  using rng_t  = rng_tt<attr_t>;
  using bun_t  = bun_tt<attr_t>;  // 2-tuple: [k, a]
  using R_t    = RelationRS2<bun_t>;  // key relation = build relation
  using S_t    = RelationRS2<bun_t>;

  // constexpr bool    lTraceInner    = true;
  constexpr uint    lLog2ChunkSize = 12;

  meas_eval_t lMeasBuild(17);  // store several measurements

  // plans to be executed/measured (note: only "SR" plans because "build unique" (R))
  PlanContainer<R_t, S_t, Ttest, Ttrace> lPlanContainer;
  planset_t lPlanSet = (1 << kPlanChOrigSR)  |
                       (1 << kPlanChImprSR)  |
                       (1 << kPlanChRpUpkSR) |
                       (1 << kPlanChRpPkdSR) |
                       (1 << kPlanChAmacSR)  |
                       (1 << kPlanChImacSR)  |
                       (1 << kPlan3dOrigSR)  |
                       (1 << kPlan3dImprSR)  |
                       (1 << kPlan3dRpUpkSR) |
                       (1 << kPlan3dRpPkdSR) |
                       (1 << kPlan3dAmacSR)  |
                       (1 << kPlan3dImacSR);

  constexpr uint64_t lMinHtDirSize  = (Ttest ? 17 : 191);  // mininum hash directory size

  // calculate hash table sizes
  df::infra::CrystalTimer lClock;
  lClock.start();
  uint64_vt lHashTableSizesR(aCbGlob.max_log_card_key() + 1);  // one entry for each size 0 <= n <= max_log_card_key()
  for (uint i = 0; i <= aCbGlob.max_log_card_key(); ++i) {
     const uint64_t lCardKey = std::ceil((uint64_t{1} << i) * aCbGlob.card_key_factor());
     lHashTableSizesR[i] = std::max(lMinHtDirSize, get_next_prime_1_1(lCardKey));
  }
  lClock.stop();
  std::cout << "# calc hashtable sizes for R: " << lClock.duration_ns() << " [ns]" << std::endl;

  line_t lLine;  // measurement output
  rng_t  lRng;

  cmeasure_t lMeasure;
  cmeasure_start(&lMeasure);
  R_t R;
  // init R.k? to truely allocate physical memory (YES)
  const uint64_t lCapR = std::ceil((uint64_t{1} << aCbGlob.max_log_card_key()) * aCbGlob.card_key_factor());
  //std::cout << "@@ max: " << (uint64_t{1} << aCbGlob.max_log_card_key()) << std::endl;
  //std::cout << "@@ max w/ factor: " << std::ceil((uint64_t{1} << aCbGlob.max_log_card_key()) * aCbGlob.card_key_factor()) << std::endl;
  R.mem_alloc(lCapR);
  build_rel_key<bun_t>(R, lCapR);
  //std::cout << "@@ |R| = " << R.size() << std::endl;
  cmeasure_stop(&lMeasure);
  std::cout << "# time to build R: " << cmeasure_total_s(&lMeasure) << " [s]" << std::endl;

  if(nullptr != aCbGlob.os_ptr()) {
    line_t::print_header_line(aCbGlob.os());
  } else {
    // line_t::print_header_line(std::cout);  // XXX
  }

  lLine._attr_size = sizeof(attr_t);           // (4)

  cmeasure_t lMeasureLoop;
  cmeasure_start(&lMeasureLoop);
  for (uint lLog2CardKey =  aCbGlob.min_log_card_key();
            lLog2CardKey <= aCbGlob.max_log_card_key();
          ++lLog2CardKey) {
    const uint64_t lCardKey = std::ceil((uint64_t{1} << lLog2CardKey) * aCbGlob.card_key_factor());
    const double   lLog2CardKeyScaled = std::log2(lCardKey);
    const uint     lNoRep   = aCbGlob.get_no_rep(lLog2CardKey);
    const uint64_t lHashtableSizeR = lHashTableSizesR[lLog2CardKey];

    //std::cout << "@@  lLog2CardKey      : " << lLog2CardKey << std::endl
    //          << "@@  lCardKey          : " << lCardKey << std::endl
    //          << "@@  lLog2CardKeyScaled: "<< lLog2CardKeyScaled << std::endl;

    lMeasBuild.resize(lNoRep);

    R.resize(lCardKey);  // resize only sets size, no change to stored values

    // cmeasure_start(&lMeasure); // lohnt nicht
    for (BvMemberIdxDesc<planset_t> lIter(lPlanSet); lIter.isValid(); ++lIter) {
      const plan_et lPlanId = (plan_et) (*lIter);
      assert(lPlanId < kNoPlan);

      //std::cout << "++ " << __FUNCTION__ << ": " << lPlanId << " = " << to_string(lPlanId) << std::endl;

      lPlanContainer.build_plan(lPlanId, lHashtableSizeR, lLog2ChunkSize);
      lMeasBuild.init();
      lPlanContainer.run_build_part_a(lPlanId, R, lNoRep, lMeasBuild);
      const size_t lMemConsumptionHt = lPlanContainer.mem_consumption_ht(lPlanId);
      const size_t lHtNumNodes = lPlanContainer.ht_num_nodes(lPlanId);
      const size_t lHtNumNodesSub = lPlanContainer.ht_num_nodes_sub(lPlanId);
      const double lHtDirFillFactor = lPlanContainer.ht_dir_fill_factor(lPlanId);
      const size_t lCntBld = lPlanContainer.count_build(lPlanId);
      const size_t lCntPrb = lPlanContainer.count_probe(lPlanId);
      lPlanContainer.run_build_part_b(lPlanId, lMeasBuild);
      lMeasBuild.fin();

      assert(R.size() == lCntBld);
      assert(0 == lCntPrb);

      lLine._rt_stat_build[lPlanId]._cc_min = lMeasBuild.min();
      lLine._rt_stat_build[lPlanId]._cc_max = lMeasBuild.max();
      lLine._rt_stat_build[lPlanId]._cc_med = lMeasBuild.median();
      lLine._rt_stat_build[lPlanId]._cc_avg = lMeasBuild.avg();
      lLine._rt_stat_build[lPlanId]._cc_mx2 = lMeasBuild.max2();
      lLine._rt_stat_build[lPlanId]._cc_fst = lMeasBuild.first();

      // ht statistics: memory consumption, number of collision chain nodes, fill degree of dir
      lLine._mem_stat_ht.at(lPlanId) = lMemConsumptionHt;
      lLine._ht_num_nodes_main.at(lPlanId) = lHtNumNodes;
      lLine._ht_num_nodes_sub.at(lPlanId) = lHtNumNodesSub;
      lLine._ht_dir_fill.at(lPlanId) = lHtDirFillFactor;

      // build and probe counts
      lLine._cnt_build.at(lPlanId) = lCntBld;
      lLine._cnt_top.at(lPlanId) = lCntPrb;

      //std::cout
      //  << "@@@ " << to_string(lPlanId) << std::endl
      //  << "@@@   mem stat:    " << lLine._mem_stat_ht.at(lPlanId) << " = " << lMemConsumptionHt << std::endl
      //  << "@@@   ht main:     " << lLine._ht_num_nodes_main.at(lPlanId) << " = " << lHtNumNodes << std::endl
      //  << "@@@   ht sub:      " << lLine._ht_num_nodes_sub.at(lPlanId) << " = " << lHtNumNodesSub << std::endl
      //  << "@@@   ht dir fill: " << lLine._ht_dir_fill.at(lPlanId) << " = " << lHtDirFillFactor << std::endl;

      lLine._id  = ((((lLog2CardKey * 100 + 0) * 100 + 0) * 10 + 0) * 10 + (aCbGlob.card_key_factor() != 1.0)) * 10 + 0;  // (1)
      lLine._no_rep_key       = lNoRep;                // (3)
      lLine._log_card_key     = lLog2CardKeyScaled;    // (11)
      lLine._card_key         = lCardKey;              // (12)
      lLine._hashdir_size_key = lHashtableSizeR;       // (15)
    } // end loop plans
    // cmeasure_stop(&lMeasure);
    if(nullptr != aCbGlob.os_ptr()) {
      aCbGlob.os() << lLine << std::endl;
    } else {
      std::cout << lLine << std::endl;
    }
  } // end loop log2_card_key
  cmeasure_stop(&lMeasureLoop);
  std::cout << "# runtime loop: " 
            << cmeasure_total_s(&lMeasureLoop) << " [s], "
            << cmeasure_total_h(&lMeasureLoop) << " [h]"
            << std::endl;
}
