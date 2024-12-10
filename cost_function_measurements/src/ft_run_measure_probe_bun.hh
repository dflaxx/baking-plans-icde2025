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

#include "HllErtl.hh"
#include "RelationRS.hh"
#include "cb_glob.hh"
#include "line.hh"
#include "plan_container.hh"

// more includes
// thus include last in main

template<typename Tattr, bool Ttest, bool Ttrace>
void
ft_run_measure_probe_bun(const CbGlob& aCbGlob) {
  using attr_t = Tattr;
  using rng_t  = rng_tt<attr_t>;
  using hll_t  = HllErtl<attr_t>;
  using bun_t  = bun_tt<attr_t>;
  using R_t    = RelationRS2<bun_t>;
  using S_t    = RelationRS2<bun_t>;

  // constexpr bool    lTraceInner    = true;
  constexpr uint    lLog2ChunkSize = 12;
  // constexpr double  lHllFactor     = 1.03;

  meas_eval_t lMeasBuild(17);  // required for run_build_part_a/b, measured values not used
  meas_eval_t lMeasProbe(17);

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

  constexpr uint64_t lMinHtDirSize  = (Ttest ? 17 : 191);

  // calculate hash table sizes
  df::infra::CrystalTimer lClock;
  lClock.start();
  uint64_vt lHashTableSizesR(aCbGlob.max_log_card_key() + 1);
  for(uint i = 0; i <= aCbGlob.max_log_card_key(); ++i) {
     const uint64_t lCardKey = std::ceil((uint64_t{1} << i) * aCbGlob.card_key_factor());
     lHashTableSizesR[i] = std::max(lMinHtDirSize, get_next_prime_1_1(lCardKey));
  }
  lClock.stop();
  std::cout << "# calc hashtable sizes for R: " << lClock.duration_ns() << " [ns]" << std::endl;

  line_t lLine[aCbGlob.max_log_card_fk() + 1]; // max log2 card S = 30
  rng_t  lRng;


  /* base relations: memory alloc & data gen */

  cmeasure_t lMeasure;
  cmeasure_start(&lMeasure);
  R_t R;
  // init R.k? to truely allocate physical memory (YES)
  const uint64_t lCapR = std::ceil((uint64_t{1} << aCbGlob.max_log_card_key()) * aCbGlob.card_key_factor());  // max size of R
  //std::cout << "@@ max: " << (uint64_t{1} << aCbGlob.max_log_card_key()) << std::endl;
  //std::cout << "@@ max w/ factor: " << lCapR << std::endl;
  R.mem_alloc(lCapR);
  build_rel_key<bun_t>(R, lCapR);
  //std::cout << "@@ |R| = " << R.size() << std::endl;
  cmeasure_stop(&lMeasure);
  std::cout << "# time to build R: " << cmeasure_total_s(&lMeasure) << " [s]" << std::endl;
  assert(R.size() == lCapR);

  // init S.k? to truely allocate physical memory (YES)
  cmeasure_start(&lMeasure);
  S_t S;
  const uint64_t lCapS = std::ceil((uint64_t{1} << aCbGlob.max_log_card_fk()) * aCbGlob.card_fk_factor());  // max size of S
  //std::cout << "@@ FK max: " << (uint64_t{1} << aCbGlob.max_log_card_fk()) << std::endl;
  //std::cout << "@@ FK max w/ factor: " << lCapS << std::endl;
  S.mem_alloc(lCapS);
  S.mem_init();
  cmeasure_stop(&lMeasure);
  //std::cout << "@@ |S| = " << S.size() << std::endl;
  std::cout << "# time to alloc S: " << cmeasure_total_s(&lMeasure) << " [s]" << std::endl;
  assert(S.size() == lCapS);


  if(nullptr != aCbGlob.os_ptr()) {
    line_t::print_header_line(aCbGlob.os());
  } else {
    // line_t::print_header_line(std::cout);  // XXX
  }

  for(uint i = 0; i <= aCbGlob.max_log_card_fk(); ++i) {
    lLine[i]._attr_size = sizeof(attr_t);           // (4)
  }

  uint64_t lResCard = 0;
  assert(aCbGlob.max_log_dom_fk_reduction() <= aCbGlob.max_log_card_fk());

  // constants for the log2 of some measures of the foreign key relation S
  const uint lLog2CardMin      = aCbGlob.min_log_card_fk();
  const uint lLog2CardMax      = aCbGlob.max_log_card_fk();
  const uint lLog2DomRedMin    = aCbGlob.min_log_dom_fk_reduction();
  const uint lLog2DomRedMax    = aCbGlob.max_log_dom_fk_reduction();

  // compute lower and upper limit for effective domain size of FK relation S
  const uint lLog2DomSizeBegin = std::max(
                                   aCbGlob.min_log_domsize_fk(),
                                   (lLog2CardMin < lLog2DomRedMax) ?
                                    lLog2CardMin : (lLog2CardMin - lLog2DomRedMax));
  const uint lLog2DomSizeLimit = std::min(
                                   aCbGlob.max_log_domsize_fk(),
                                  (lLog2CardMax < lLog2DomRedMin) ?
                                   lLog2CardMax : (lLog2CardMax - lLog2DomRedMin));

  std::cout << "# card  fk = " << lLog2CardMin      << " .. " << lLog2CardMax << std::endl;
  std::cout << "# dom  red = " << lLog2DomRedMin    << " .. " << lLog2DomRedMax << std::endl;
  std::cout << "# dom size = " << lLog2DomSizeBegin << " .. " << lLog2DomSizeLimit << std::endl;

  // @loop_order
  // 1. skew: for skew in {0, 1}
  // 2. FK dom size: for fk_dom_size between fk_dom_size_begin, fk_dom_size_limit
  // 3. key card: for key_card between key_card_min, key_card_max
  // 4. plans: for each plan in plan_set
  // 5. FK card: for fk_card between fk_card_begin, fk_card_limit

  uint lCount = 0;
  uint lCountHeavyS = 0;
  cmeasure_t lMeasureLoop; 
  cmeasure_start(&lMeasureLoop);
  cmeasure_t lMeasureBodySkew;
  cmeasure_t lMeasureLoopCardKey;
  cmeasure_t lMeasureLoopPlan;

  hll_t lHll(14);
  // loop 1.
  for(uint lSkew = aCbGlob.min_skew(); lSkew <= aCbGlob.max_skew(); ++lSkew) {
    cmeasure_start(&lMeasureBodySkew);
    std::cout << "#  skew = " << lSkew << std::endl;
    // loop 2.
    for(uint lLog2DomSizeFk  = lLog2DomSizeBegin;
             lLog2DomSizeFk <= lLog2DomSizeLimit;
           ++lLog2DomSizeFk) {
      const uint     lLog2CardFkBegin = std::max<uint>(lLog2CardMin, lLog2DomSizeFk + lLog2DomRedMin);
      const uint     lLog2CardFkLimit = std::min<uint>(lLog2CardMax, lLog2DomSizeFk + lLog2DomRedMax);
      const uint64_t lCardLimit = std::ceil((uint64_t{1} << lLog2CardFkLimit) * aCbGlob.card_fk_factor());
      const uint64_t lDomSizeFk = std::ceil((uint64_t{1} << lLog2DomSizeFk) * aCbGlob.card_fk_factor());

      const double lLog2DomSizeFkScaled = std::log2(lDomSizeFk);

      // std::cout << "#    lLog2DomSizeFk = " <<   lLog2DomSizeFk << std::endl;
      // std::cout << "#    lLog2Card      = " <<   lLog2CardFkBegin << " .. "
      //                                       <<   lLog2CardFkLimit << std::endl;

      cmeasure_start(&lMeasure);
      build_rel_fk<bun_t,rng_t>(S, lCardLimit, lDomSizeFk, lSkew, lRng); // sic!
      cmeasure_stop(&lMeasure);
      std::cout << "#   time build S(" << lSkew << ',' 
                << lLog2DomSizeFkScaled << "): " 
                << cmeasure_total_s(&lMeasure) << " [s]" 
                << std::endl;
      ++lCountHeavyS;

      cmeasure_start(&lMeasureLoopCardKey);
      // loop 3.
      for(uint lLog2CardKey =  aCbGlob.min_log_card_key();
               lLog2CardKey <= aCbGlob.max_log_card_key();
               ++lLog2CardKey) {
        const uint64_t lCardKey = std::ceil((uint64_t{1} << lLog2CardKey) * aCbGlob.card_key_factor());
        const uint64_t lHashtableSizeR = lHashTableSizesR[lLog2CardKey];  // lookup of precomputed HT sizes

        const double lLog2CardKeyScaled = std::log2(lCardKey);

        R.resize(lCardKey);  // sic!

        cmeasure_start(&lMeasureLoopPlan);
        // loop 4.
        for(BvMemberIdxDesc<planset_t> lIter(lPlanSet); lIter.isValid(); ++lIter) {
          const plan_et lPlanId = (plan_et) (*lIter);
          assert(lPlanId < kNoPlan);

          //std::cout << "++ " << __FUNCTION__ << ": " << lPlanId << " = " << to_string(lPlanId) << std::endl;

          lPlanContainer.build_plan(lPlanId, lHashtableSizeR, lLog2ChunkSize);
          lMeasBuild.init();
          lPlanContainer.run_build_part_a(lPlanId, R, 1, lMeasBuild);

          // HT stats, BLD count
          size_t lMemConsumptionHt = lPlanContainer.mem_consumption_ht(lPlanId);
          const size_t lHtNumNodes = lPlanContainer.ht_num_nodes(lPlanId);
          const size_t lHtNumNodesSub = lPlanContainer.ht_num_nodes_sub(lPlanId);
          const double lHtDirFillFactor = lPlanContainer.ht_dir_fill_factor(lPlanId);
          const size_t lCntBld = lPlanContainer.count_build(lPlanId);

          for(uint lLog2CardFk = lLog2CardFkBegin; lLog2CardFk <= lLog2CardFkLimit; ++lLog2CardFk) {
            lLine[lLog2CardFk]._time = 0;
          }

          // loop 5.
          for(uint lLog2CardFk = lLog2CardFkBegin; lLog2CardFk <= lLog2CardFkLimit; ++lLog2CardFk) {
            const uint64_t lCardFk       = std::ceil((uint64_t{1} << lLog2CardFk) * aCbGlob.card_key_factor());
            const uint     lLog2DomRedFk = lLog2CardFk - lLog2DomSizeFk;
            const uint     lNoRep        = aCbGlob.get_no_rep(lLog2CardFk);

            const double lLog2CardFkScaled = std::log2(lCardFk);

            S.resize(lCardFk);
            // determine number of distinct values
            lHll.clear();
            for(uint i = 0; i < S.card(); ++i) {
              lHll.insert(ht::murmur_hash<attr_t>(S[i].a));
            }
            const double   lEstNoDvDouble = lHll.estimate();
            const uint64_t lEstNoDvInt    = (uint64_t) std::ceil(lEstNoDvDouble);

            assert(lLog2DomSizeFk >= aCbGlob.min_log_card_fk());
            assert(lLog2DomSizeFk <= aCbGlob.max_log_card_fk());
            assert(lLog2DomRedFk  >= aCbGlob.min_log_dom_fk_reduction());
            assert(lLog2DomRedFk  <= aCbGlob.max_log_dom_fk_reduction());
            assert((lLog2DomSizeFk + lLog2DomRedFk) == lLog2CardFk);
            assert(S.card() == lCardFk);

            ++lCount;

            lMeasProbe.resize(lNoRep);
            lMeasProbe.init();
            lResCard = lPlanContainer.run_probe(lPlanId,
                                                S,
                                                lNoRep,
                                                lMeasProbe);
            lMeasProbe.fin();

            lLine[lLog2CardFk]._time += lMeasProbe.sum();

            lLine[lLog2CardFk]._rt_stat_probe[lPlanId]._cc_min = lMeasProbe.min();
            lLine[lLog2CardFk]._rt_stat_probe[lPlanId]._cc_max = lMeasProbe.max();
            lLine[lLog2CardFk]._rt_stat_probe[lPlanId]._cc_med = lMeasProbe.median();
            lLine[lLog2CardFk]._rt_stat_probe[lPlanId]._cc_avg = lMeasProbe.avg();
            lLine[lLog2CardFk]._rt_stat_probe[lPlanId]._cc_mx2 = lMeasProbe.max2();
            lLine[lLog2CardFk]._rt_stat_probe[lPlanId]._cc_fst = lMeasProbe.first();

            // ht statistics: memory consumption, number of collision chain nodes, fill degree of dir
            lLine[lLog2CardFk]._mem_stat_ht.at(lPlanId) = lMemConsumptionHt;
            lLine[lLog2CardFk]._ht_num_nodes_main.at(lPlanId) = lHtNumNodes;
            lLine[lLog2CardFk]._ht_num_nodes_sub.at(lPlanId) = lHtNumNodesSub;
            lLine[lLog2CardFk]._ht_dir_fill.at(lPlanId) = lHtDirFillFactor;

            // build and top (~probe) counts
            lLine[lLog2CardFk]._cnt_build.at(lPlanId) = lCntBld;
            lLine[lLog2CardFk]._cnt_top.at(lPlanId) = lResCard;

            // id format: KKFFDDSkf (key, fk, domain, skew, scale key, scale fk)
            lLine[lLog2CardFk]._id =  // (1)
              ( ( ((lLog2CardKey * 100 + lLog2CardFk) * 100 + lLog2DomSizeFk) * 10 + lSkew
                ) * 10 + (aCbGlob.card_key_factor() != 1.0)
              ) * 10 + (aCbGlob.card_fk_factor() != 1.0);

            lLine[lLog2CardFk]._no_rep_fk        = lNoRep;                // (3)
            lLine[lLog2CardFk]._attr_size        = sizeof(attr_t);        // (4)
            lLine[lLog2CardFk]._log_card_fk      = lLog2CardFkScaled;     // (5)
            lLine[lLog2CardFk]._card_fk          = lCardFk;               // (6)
            lLine[lLog2CardFk]._log_dom_red_fk   = lLog2DomRedFk;         // (7)
            lLine[lLog2CardFk]._log_dom_fk       = lLog2DomSizeFkScaled;  // (8)
            //std::cout << "@@@ lLog2DomSizeFkScaled = " << lLog2DomSizeFkScaled << std::endl;
            //std::cout << "@@@ lLine[lLog2CardFk]._log_dom_fk  = " << lLine[lLog2CardFk]._log_dom_fk << std::endl;
            lLine[lLog2CardFk]._dom_fk           = lDomSizeFk;            // (9)
            lLine[lLog2CardFk]._skew_fk          = lSkew;                 // (10)
            lLine[lLog2CardFk]._log_card_key     = lLog2CardKeyScaled;    // (11)
            lLine[lLog2CardFk]._card_key         = lCardKey;              // (12)
            lLine[lLog2CardFk]._nodv_fk_hll      = lEstNoDvInt;           // (13)
            lLine[lLog2CardFk]._hashdir_size_key = lHashtableSizeR;       // (15)
            lLine[lLog2CardFk]._card_result      = lResCard;              // (17)
            if(kPlan3dOrigSR == lPlanId) {
              lLine[lLog2CardFk]._card_sj_sr = lPlanContainer.card_sj_sr(); // (16)
            }
          } // end loop 5: log2_card_fk

          lPlanContainer.run_build_part_b(lPlanId, lMeasBuild);

        } // end loop 4: plans
        cmeasure_stop(&lMeasureLoopPlan);
        std::cout << "#     runtime loop plan (" << lLog2CardKey << "): " << cmeasure_total_s(&lMeasureLoopPlan) << " [s]" << std::endl;

        for(uint lLog2CardFk = lLog2CardFkBegin; lLog2CardFk <= lLog2CardFkLimit; ++lLog2CardFk) {
          lLine[lLog2CardFk]._time /= 2600000000.0;
        }

        if(nullptr != aCbGlob.os_ptr()) {
          for(uint i = lLog2CardFkBegin; i <= lLog2CardFkLimit; ++i) {
            aCbGlob.os() << lLine[i] << std::endl;
          }
        } else {
          for(uint i = lLog2CardFkBegin; i <= lLog2CardFkLimit; ++i) {
            std::cout << lLine[i] << std::endl;
          }
        }

      } // end loop 3: log2_card_key
      cmeasure_stop(&lMeasureLoopCardKey);
      std::cout << "#   runtime loop card_key(" << lLog2DomSizeFk
                << "): " << cmeasure_total_h(&lMeasureLoopCardKey) << " [h]"
                << std::endl;
    } // all dom sizes for S considered
    cmeasure_stop(&lMeasureBodySkew);
    std::cout << "#   runtime body skew(" << lSkew << "): " 
              << cmeasure_total_s(&lMeasureBodySkew) << " [s]" << std::endl;
  } // end loop skew
  cmeasure_stop(&lMeasureLoop);
  std::cout << "# no run     = " << lCount << std::endl;
  std::cout << "# no heavy s = " << lCountHeavyS << std::endl;
  std::cout << "# runtime loop: " 
            << cmeasure_total_s(&lMeasureLoop) << " [s], "
            << cmeasure_total_h(&lMeasureLoop) << " [h]"
            << std::endl;
}
