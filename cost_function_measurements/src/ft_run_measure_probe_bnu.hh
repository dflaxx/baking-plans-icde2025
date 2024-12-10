#pragma once

#include "dfinfra/standard_includes.hh"
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

/*
 * run measurements for probe with build on none-unique
 *
 * loop order of parameter iteration, data gen and measurement,
 * see below at @loop_order
 */
template<typename Tattr, bool Ttest, bool Ttrace>
void
ft_run_measure_probe_bnu(const CbGlob& aCbGlob) {
  using attr_t = Tattr;
  using rng_t  = rng_tt<attr_t>;
  using hll_t  = HllErtl<attr_t>;
  using bun_t  = bun_tt<attr_t>;
  using R_t    = RelationRS2<bun_t>;
  using S_t    = RelationRS2<bun_t>;

  constexpr bool    lTraceInner    = true;
  constexpr uint    lLog2ChunkSize = 12;
  constexpr double  lHllFactor     = 1.03;  // HLL security margin: use HLL estimate x lHllFactor

  // set up measurement containers
  meas_eval_t lMeasBuild(17);
  meas_eval_t lMeasProbe(17);

  // plan container, plans to run
  PlanContainer<R_t, S_t, Ttest, Ttrace> lPlanContainer;
  planset_t lPlanSet = (1 << kPlanChOrigRS)  |
                       (1 << kPlanChImprRS)  |
                       (1 << kPlanChRpUpkRS) |
                       (1 << kPlanChRpPkdRS) |
                       (1 << kPlanChAmacRS)  |
                       (1 << kPlanChImacRS)  |
                       (1 << kPlan3dOrigRS)  |
                       (1 << kPlan3dRpUpkRS) |
                       (1 << kPlan3dRpPkdRS) |
                       (1 << kPlan3dImprRS)  |
                       (1 << kPlan3dAmacRS)  |
                       (1 << kPlan3dImacRS);

  // minimum hash table direcroty size
  constexpr uint64_t lMinHtDirSize  = (Ttest ? 17 : 191);

  line_t lLine[aCbGlob.max_log_card_key() + 1]; // max log2 card R = 30
  rng_t  lRng;


  /* base relations: memory alloc & data gen */

  cmeasure_t lMeasure;
  cmeasure_start(&lMeasure);
  R_t R;  // key relation
  // init R.k? to truely allocate physical memory (YES)
  const uint64_t lCapR = std::ceil((uint64_t{1} << aCbGlob.max_log_card_key()) * aCbGlob.card_key_factor());  // max size of R
  //std::cout << "@@ max: " << (uint64_t{1} << aCbGlob.max_log_card_key()) << std::endl;
  //std::cout << "@@ max w/ factor: " << std::ceil((uint64_t{1} << aCbGlob.max_log_card_key()) * aCbGlob.card_key_factor()) << std::endl;
  R.mem_alloc(lCapR);
  build_rel_key<bun_t>(R, lCapR);
  //std::cout << "@@ |R| = " << R.size() << std::endl;
  cmeasure_stop(&lMeasure);
  std::cout << "# time to build R: " << cmeasure_total_s(&lMeasure) << " [s]" << std::endl;

  // init S.k? to truely allocate physical memory (YES)
  cmeasure_start(&lMeasure);
  S_t S;  // FK relation
  const uint64_t lCapS = std::ceil((uint64_t{1} << aCbGlob.max_log_card_fk()) * aCbGlob.card_fk_factor());  // max size of S
  //std::cout << "@@ FK max: " << (uint64_t{1} << aCbGlob.max_log_card_fk()) << std::endl;
  //std::cout << "@@ FK max w/ factor: " << std::ceil((uint64_t{1} << aCbGlob.max_log_card_fk()) * aCbGlob.card_fk_factor()) << std::endl;
  S.mem_alloc(lCapS);
  S.mem_init();
  cmeasure_stop(&lMeasure);
  //std::cout << "@@ |S| = " << S.size() << std::endl;
  std::cout << "# time to alloc S: " << cmeasure_total_s(&lMeasure) << " [s]" << std::endl;


  if(nullptr != aCbGlob.os_ptr()) {
    line_t::print_header_line(aCbGlob.os());
  } else {
    // line_t::print_header_line(std::cout);  // XXX
  }

  for(uint i = 0; i <= aCbGlob.max_log_card_key(); ++i) {
    lLine[i]._attr_size = sizeof(attr_t);           // (4)
  }

  uint64_t lResCard = 0;
  assert(aCbGlob.max_log_dom_fk_reduction() <= aCbGlob.max_log_card_fk());

  // constants for the log2 of some measures of the foreign key relation S
  const uint lLog2CardMin      = aCbGlob.min_log_card_fk();
  const uint lLog2CardMax      = aCbGlob.max_log_card_fk();
  const uint lLog2DomRedMin    = aCbGlob.min_log_dom_fk_reduction();
  const uint lLog2DomRedMax    = aCbGlob.max_log_dom_fk_reduction();

  // compute lower and upper limit for effective domain size of FK relation
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
  // 3. FK card: for fk_card between fk_card_begin, fk_card_limit
  // 4. plans: for each plan in plan_set
  // 5. key card: for key_card between key_card_min, key_card_max

  uint lCount = 0;
  uint lCountHeavyS = 0;
  cmeasure_t lMeasureLoop;
  cmeasure_start(&lMeasureLoop);
  cmeasure_t lMeasureBodySkew;  // loop 1.
  cmeasure_t lMeasureBodyDomSizeFk;  // loop 2.
  cmeasure_t lMeasureLoopPlan;  // loop 4.
  hll_t lHll(14);

  // loop 1.
  for(uint lSkew = aCbGlob.min_skew(); lSkew <= aCbGlob.max_skew(); ++lSkew) {
    cmeasure_start(&lMeasureBodySkew);
    std::cout << "#  skew = " << lSkew << std::endl;

    // loop 2.
    for(uint lLog2DomSizeFk  = lLog2DomSizeBegin;
             lLog2DomSizeFk <= lLog2DomSizeLimit;
           ++lLog2DomSizeFk) {
      // compute lower and upper limit for effetive cardinality of foreign key relation
      const uint     lLog2CardFkBegin = std::max<uint>(lLog2CardMin, lLog2DomSizeFk + lLog2DomRedMin);
      const uint     lLog2CardFkLimit = std::min<uint>(lLog2CardMax, lLog2DomSizeFk + lLog2DomRedMax);
      const uint64_t lCardLimit = std::ceil((uint64_t{1} << lLog2CardFkLimit) * aCbGlob.card_fk_factor());
      const uint64_t lDomSizeFk = std::ceil((uint64_t{1} << lLog2DomSizeFk) * aCbGlob.card_fk_factor());

      const double lLog2DomSizeFkScaled = std::log2(lDomSizeFk);

      cmeasure_start(&lMeasureBodyDomSizeFk);

      // std::cout << "# lLog2DomSizeFk = " <<   lLog2DomSizeFk << std::endl;
      // std::cout << "# lLog2Card      = " <<   lLog2CardFkBegin << " .. "
      //                                       <<   lLog2CardFkLimit << std::endl;

      // build FK relation S
      cmeasure_start(&lMeasure);
      build_rel_fk<bun_t,rng_t>(S, lCardLimit, lDomSizeFk, lSkew, lRng); // sic!
      cmeasure_stop(&lMeasure);
      std::cout << "# time build S   : " << cmeasure_total_s(&lMeasure) << " [s]" << std::endl;
      ++lCountHeavyS;


      // loop 3.
      for(uint lLog2CardFk = lLog2CardFkBegin; lLog2CardFk <= lLog2CardFkLimit; ++lLog2CardFk) {
        const uint64_t lCardFk       = std::ceil((uint64_t{1} << lLog2CardFk) * aCbGlob.card_fk_factor());
        const uint     lLog2DomRedFk = lLog2CardFk - lLog2DomSizeFk;

        const double lLog2CardFkScaled = std::log2(lCardFk);

        if constexpr (lTraceInner) {
          std::cout << "#      lLog2CardFk = " << lLog2CardFk   << std::endl;
          std::cout << "#      lCardFk     = " << lCardFk       << std::endl;
          std::cout << "##     lLog2DomRed = " << lLog2DomRedFk << std::endl;
        }

        S.resize(lCardFk);

        for(uint i = aCbGlob.min_log_card_key(); i <= aCbGlob.max_log_card_key(); ++i) {
          lLine[i]._no_rep_fk = aCbGlob.get_no_rep(lLog2CardFk); // (3)
          lLine[i]._log_card_fk = lLog2CardFkScaled;             // (5)
          lLine[i]._card_fk = lCardFk;                           // (6)
          lLine[i]._log_dom_red_fk = lLog2DomRedFk;              // (7)
          lLine[i]._log_dom_fk = lLog2DomSizeFkScaled;           // (8)
          lLine[i]._dom_fk = lDomSizeFk;                         // (9)
          lLine[i]._skew_fk = lSkew;                             // (10)
          lLine[i]._time = 0;
        }

        assert(lLog2DomSizeFk >= aCbGlob.min_log_card_fk());
        assert(lLog2DomSizeFk <= aCbGlob.max_log_card_fk());
        assert(lLog2DomRedFk  >= aCbGlob.min_log_dom_fk_reduction());
        assert(lLog2DomRedFk  <= aCbGlob.max_log_dom_fk_reduction());
        assert((lLog2DomSizeFk + lLog2DomRedFk) == lLog2CardFk);
        assert(S.card() == lCardFk);

        // determine number of distinct values via HLL
        lHll.clear();
        for(uint i = 0; i < S.card(); ++i) {
          lHll.insert(ht::murmur_hash<attr_t>(S[i].a));
        }
        const double   lEstNoDvDouble = lHll.estimate();
        const uint64_t lEstNoDvInt    = (uint64_t) std::ceil(lEstNoDvDouble);
        const uint64_t lHashDirSzMinS = (uint64_t) std::ceil(lEstNoDvDouble * lHllFactor);
        const uint64_t lHashtableSizeS = get_next_prime_1_1(std::max(lMinHtDirSize, lHashDirSzMinS));

        if constexpr (lTraceInner) {
          std::cout << "##     hll : " << mt::roundXt(lEstNoDvDouble)
                    << ' ' << lHashDirSzMinS
                    << ' ' << lHashtableSizeS
                    << ' ' << (((double) lCardFk) / lEstNoDvDouble)
                    << std::endl;
        }

        for(uint i = aCbGlob.min_log_card_key(); i <= aCbGlob.max_log_card_key(); ++i) {
          lLine[i]._nodv_fk_hll     = lEstNoDvInt;     // (13)
          lLine[i]._hashdir_size_fk = lHashtableSizeS; // (14)
        }

        cmeasure_start(&lMeasureLoopPlan); 
        // loop 4.
        for(BvMemberIdxDesc<planset_t> lIter(lPlanSet); lIter.isValid(); ++lIter) {
          const plan_et lPlanId = (plan_et) (*lIter);
          assert(lPlanId < kNoPlan);

          //std::cout << "++ " << __FUNCTION__ << ": " << lPlanId << " = " << to_string(lPlanId) << std::endl;

          lPlanContainer.build_plan(lPlanId, lHashtableSizeS, lLog2ChunkSize);
          lMeasBuild.init();
          lPlanContainer.run_build_part_a(lPlanId, S, 1, lMeasBuild);  // build operator: build HT

          // HT stats, BLD count
          size_t lMemConsumptionHt = lPlanContainer.mem_consumption_ht(lPlanId);
          const size_t lHtNumNodes = lPlanContainer.ht_num_nodes(lPlanId);
          const size_t lHtNumNodesSub = lPlanContainer.ht_num_nodes_sub(lPlanId);
          const double lHtDirFillFactor = lPlanContainer.ht_dir_fill_factor(lPlanId);
          const size_t lCntBld = lPlanContainer.count_build(lPlanId);

          cmeasure_start(&lMeasure);
          // loop 5.
          for(uint lLog2CardKey =  aCbGlob.min_log_card_key();
                   lLog2CardKey <= aCbGlob.max_log_card_key();
                 ++lLog2CardKey) {
            const uint64_t lCardKey           = std::ceil((uint64_t{1} << lLog2CardKey) * aCbGlob.card_key_factor());
            const double   lLog2CardKeyScaled = std::log2(lCardKey);
            // const uint64_t lHashtableSizeR = lHashTableSizesR[lLog2CardKey];
            // const uint64_t lHashtableSizeR = get_next_prime_1_1(std::max(lMinHtDirSize, lCardKey));
            // lLine[lLog2CardKey]._hashdir_size_key = lHashtableSizeR;                  // (15)

            // id format: KKFFDDSkf (key, fk, domain, skew, scale key, scale fk)
            lLine[lLog2CardKey]._id =
              ( ( ((lLog2CardKey * 100 + lLog2CardFk) * 100 + lLog2DomSizeFk) * 10 + lSkew
                ) * 10 + (aCbGlob.card_key_factor() != 1.0)
              ) * 10 + (aCbGlob.card_fk_factor() != 1.0); // (1)
            lLine[lLog2CardKey]._no_rep_key       = aCbGlob.get_no_rep(lLog2CardKey); // (2)
            lLine[lLog2CardKey]._log_card_key     = lLog2CardKeyScaled;               // (11)
            lLine[lLog2CardKey]._card_key         = lCardKey;                         // (12)

            R.resize(lCardKey);  // sic!
                                 // resize only sets _size, does not touch data or dealloc memory!
                                 // cheap

            ++lCount;

            lMeasProbe.resize(lLine[lLog2CardKey]._no_rep_key);

            // the actual measurement: probe operator
            // PlanContainer::run_probe runs the plan multiple times and updates lMeasProbe
            lMeasProbe.init();
            lResCard = lPlanContainer.run_probe(lPlanId, 
                                                R, 
                                                lLine[lLog2CardKey]._no_rep_key, 
                                                lMeasProbe);
            lMeasProbe.fin();
            lLine[lLog2CardKey]._time += lMeasProbe.sum();

            lLine[lLog2CardKey]._rt_stat_probe[lPlanId]._cc_min = lMeasProbe.min();
            lLine[lLog2CardKey]._rt_stat_probe[lPlanId]._cc_max = lMeasProbe.max();
            lLine[lLog2CardKey]._rt_stat_probe[lPlanId]._cc_med = lMeasProbe.median();
            lLine[lLog2CardKey]._rt_stat_probe[lPlanId]._cc_avg = lMeasProbe.avg();
            lLine[lLog2CardKey]._rt_stat_probe[lPlanId]._cc_mx2 = lMeasProbe.max2();
            lLine[lLog2CardKey]._rt_stat_probe[lPlanId]._cc_fst = lMeasProbe.first();

            // ht statistics: memory consumption, number of collision chain nodes, fill degree of dir
            lLine[lLog2CardKey]._mem_stat_ht.at(lPlanId) = lMemConsumptionHt;
            lLine[lLog2CardKey]._ht_num_nodes_main.at(lPlanId) = lHtNumNodes;
            lLine[lLog2CardKey]._ht_num_nodes_sub.at(lPlanId) = lHtNumNodesSub;
            lLine[lLog2CardKey]._ht_dir_fill.at(lPlanId) = lHtDirFillFactor;

            // build and top (~probe) counts
            lLine[lLog2CardKey]._cnt_build.at(lPlanId) = lCntBld;
            lLine[lLog2CardKey]._cnt_top.at(lPlanId) = lResCard;

            if(kPlan3dOrigRS == lPlanId) {
              lLine[lLog2CardKey]._card_result = lResCard;
              lLine[lLog2CardKey]._card_sj_sr  = lPlanContainer.card_sj_sr();
            }
          } // done probing one plan for all cardinalities of R (loop 5.)
          cmeasure_stop(&lMeasure);
          std::cout << "#     runtime probe plan_" << lPlanId << ": " 
                    << cmeasure_total_s(&lMeasure) << " [s] = " 
                    << cmeasure_total_h(&lMeasure) << " [h]" 
                    << std::endl;

          lPlanContainer.run_build_part_b(lPlanId, lMeasBuild);  // build operator: cleanup & dealloc memory
          lMeasBuild.fin();

        } // done probing all plans (loop 4.)
        cmeasure_stop(&lMeasureLoopPlan); 
        std::cout << "#   runtime measure loop plan(" 
                  << lLog2CardFk << ',' << lLog2DomRedFk << "): "
                  << cmeasure_total_s(&lMeasureLoopPlan) << " [s]" << std::endl;

        // should call Line::print_probe_rs
        if(nullptr != aCbGlob.os_ptr()) {
          for(uint i = aCbGlob.min_log_card_key(); i <= aCbGlob.max_log_card_key(); ++i) {
            aCbGlob.os() << lLine[i] << std::endl;
          }
        } else {
          // std::cout << lLine << std::endl;  // XXX
          for(uint i = aCbGlob.min_log_card_key(); i <= aCbGlob.max_log_card_key(); ++i) {
            std::cout << lLine[i] << std::endl;
          }
        }
      } // end loop lLog2CardFk (loop 3.)
      cmeasure_stop(&lMeasureBodyDomSizeFk);
      std::cout << "#   runtime body log2DomSizeFk(" 
                << lSkew << ',' << lLog2DomSizeFk << "): " 
                << cmeasure_total_s(&lMeasureBodyDomSizeFk) << " [s]" << std::endl;
    } // done: all dom sizes for S considered (loop 2.)
    cmeasure_stop(&lMeasureBodySkew);
    std::cout << "  # runtime body skew(" << lSkew << "): " 
              << cmeasure_total_s(&lMeasureBodySkew) << " [s]" << std::endl;
  } // end loop skew (loop 1.)

  cmeasure_stop(&lMeasureLoop);
  std::cout << "# no run     = " << lCount << std::endl;
  std::cout << "# no heavy s = " << lCountHeavyS << std::endl;
  std::cout << "# runtime loop: " 
            << cmeasure_total_s(&lMeasureLoop) << " [s], "
            << cmeasure_total_h(&lMeasureLoop) << " [h]"
            << std::endl;
}
