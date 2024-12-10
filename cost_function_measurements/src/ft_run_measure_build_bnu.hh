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

template<typename Tattr, bool Ttest, bool Ttrace>
void
ft_run_measure_build_bnu(const CbGlob& aCbGlob) {
  using attr_t = Tattr;
  using rng_t  = rng_tt<attr_t>;
  using hll_t  = HllErtl<attr_t>;
  using bun_t  = bun_tt<attr_t>;
  using R_t    = RelationRS2<bun_t>;
  using S_t    = RelationRS2<bun_t>;

  constexpr bool    lTraceInner    = false;
  constexpr uint    lLog2ChunkSize = 12;
  constexpr double  lHllFactor     = 1.03;

  meas_eval_t lMeasBuild(17);

  // plans to be executed/measured (note: only "RS" plans because "build non-unique" (S))
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

  constexpr uint64_t lMinHtDirSize  = (Ttest ? 17 : 191);

  line_t lLine;
  rng_t  lRng;

  cmeasure_t lMeasure;
  cmeasure_start(&lMeasure);

  // init S.k? to truely allocate physical memory (YES)
  cmeasure_start(&lMeasure);
  S_t S;
  const uint64_t lCapS = std::ceil((uint64_t{1} << aCbGlob.max_log_card_fk()) * aCbGlob.card_fk_factor());
  //std::cout << "@@ max: " << (uint64_t{1} << aCbGlob.max_log_card_fk()) << std::endl;
  //std::cout << "@@ max w/ factor: " << std::ceil((uint64_t{1} << aCbGlob.max_log_card_fk()) * aCbGlob.card_fk_factor()) << std::endl;
  S.mem_alloc(lCapS);
  S.mem_init();
  cmeasure_stop(&lMeasure);
  //std::cout << "@@ |S| = " << S.size() << std::endl;
  std::cout << "# time to alloc S: " << cmeasure_total_s(&lMeasure) << " [s]" << std::endl;

  if(nullptr != aCbGlob.os_ptr()) {
    line_t::print_header_line(aCbGlob.os());
  } else {
    line_t::print_header_line(std::cout);
  }

  lLine._attr_size = sizeof(attr_t);           // (4)

  // constants for the log2 of some measures of the foreign key relation S
  const uint lLog2CardMin      = aCbGlob.min_log_card_fk();
  const uint lLog2CardMax      = aCbGlob.max_log_card_fk();
  const uint lLog2DomRedMin    = aCbGlob.min_log_dom_fk_reduction();
  const uint lLog2DomRedMax    = aCbGlob.max_log_dom_fk_reduction();
  const uint lLog2DomSizeBegin = (lLog2CardMin < lLog2DomRedMax) ? 
                                  lLog2CardMin : (lLog2CardMin - lLog2DomRedMax);
  const uint lLog2DomSizeLimit = (lLog2CardMax < lLog2DomRedMin) ?
                                  lLog2CardMax : (lLog2CardMax - lLog2DomRedMin);

  std::cout << "# card  fk = " << lLog2CardMin      << " .. " << lLog2CardMax << std::endl;
  std::cout << "# dom  red = " << lLog2DomRedMin    << " .. " << lLog2DomRedMax << std::endl;
  std::cout << "# dom size = " << lLog2DomSizeBegin << " .. " << lLog2DomSizeLimit << std::endl;

  uint lCount = 0;
  uint lCountHeavyS = 0;
  cmeasure_t lMeasureLoop;
  cmeasure_t lMeasureBodySkew;
  cmeasure_start(&lMeasureLoop);
  hll_t lHll(14);
  for(uint lSkew = aCbGlob.min_skew(); lSkew <= aCbGlob.max_skew(); ++lSkew) {
    cmeasure_start(&lMeasureBodySkew);
    std::cout << "# skew = " << lSkew << std::endl;
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
                << lLog2DomSizeFk << "): " 
                << cmeasure_total_s(&lMeasure) << " [s]"
                << std::endl;
      ++lCountHeavyS;

      for(uint lLog2CardFk = lLog2CardFkBegin; lLog2CardFk <= lLog2CardFkLimit; ++lLog2CardFk) {
        const uint64_t lCardFk       = std::ceil((uint64_t{1} << lLog2CardFk) * aCbGlob.card_fk_factor());
        const uint     lLog2DomRedFk = lLog2CardFk - lLog2DomSizeFk;

        const double lLog2CardFkScaled = std::log2(lCardFk);

        if constexpr (lTraceInner) {
          std::cout << "#      lLog2CardFk = " << lLog2CardFk   << std::endl;
          std::cout << "#      lCardFk     = " << lCardFk       << std::endl;
          std::cout << "##     lLog2DomRed = " << lLog2DomRedFk << std::endl;
        } else {
          std::cout << "#     lLog2CardFk = " << lLog2CardFk   << std::endl;
        }

        S.resize(lCardFk);

        const uint lNoRep      = aCbGlob.get_no_rep(lLog2CardFk);
        lLine._no_rep_fk       = lNoRep;               // (3)
        lLine._log_card_fk     = lLog2CardFkScaled;    // (5)  <= scaled
        lLine._card_fk         = lCardFk;              // (6)
        lLine._log_dom_red_fk  = lLog2DomRedFk;        // (7)
        lLine._log_dom_fk      = lLog2DomSizeFkScaled; // (8)  <= scaled
        lLine._dom_fk          = lDomSizeFk ;          // (9)
        lLine._skew_fk         = lSkew;                // (10)

        assert(lLog2DomSizeFk >= aCbGlob.min_log_card_fk());
        assert(lLog2DomSizeFk <= aCbGlob.max_log_card_fk());
        assert(lLog2DomRedFk  >= aCbGlob.min_log_dom_fk_reduction());
        assert(lLog2DomRedFk  <= aCbGlob.max_log_dom_fk_reduction());
        assert((lLog2DomSizeFk + lLog2DomRedFk) == lLog2CardFk);
        assert(S.card() == lCardFk);

        lMeasBuild.resize(lNoRep);

        // determine number of distinct values
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

        lLine._nodv_fk_hll     = lEstNoDvInt;     // (13)
        lLine._hashdir_size_fk = lHashtableSizeS; // (14)

        for(BvMemberIdxDesc<planset_t> lIter(lPlanSet); lIter.isValid(); ++lIter) {
          const plan_et lPlanId = (plan_et) (*lIter);
          assert(lPlanId < kNoPlan);

          //std::cout << "++ " << __FUNCTION__ << ": " << lPlanId << " = " << to_string(lPlanId) << std::endl;

          cmeasure_start(&lMeasure);

          lPlanContainer.build_plan(lPlanId, lHashtableSizeS, lLog2ChunkSize);

          lMeasBuild.init();
          lPlanContainer.run_build_part_a(lPlanId, S, lNoRep, lMeasBuild);
          size_t lMemConsumptionHt = lPlanContainer.mem_consumption_ht(lPlanId);
          const size_t lHtNumNodes = lPlanContainer.ht_num_nodes(lPlanId);
          const size_t lHtNumNodesSub = lPlanContainer.ht_num_nodes_sub(lPlanId);
          const double lHtDirFillFactor = lPlanContainer.ht_dir_fill_factor(lPlanId);
          const size_t lCntBld = lPlanContainer.count_build(lPlanId);
          const size_t lCntPrb = lPlanContainer.count_probe(lPlanId);
          lPlanContainer.run_build_part_b(lPlanId, lMeasBuild);
          lMeasBuild.fin();

          assert(S.size() == lCntBld);
          assert(0 == lCntPrb);

          // id format: KKFFDDSkf (key, fk, domain, skew, scale key, scale fk)
          lLine._id = ((((0 * 100 + lLog2CardFk) * 100 + lLog2DomSizeFk) * 10 + lSkew) * 10 + 0) * 10 + (aCbGlob.card_fk_factor() != 1.0); // (1)
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

          cmeasure_stop(&lMeasure);
          std::cout << "#       runtime build plan_" << lPlanId << "(" << to_string(lPlanId) << "): " << cmeasure_total_s(&lMeasure) << " [s]" << std::endl;
        } // done build all plans

        // should call Line::print_build_bnu
        if(nullptr != aCbGlob.os_ptr()) {
          aCbGlob.os() << lLine << std::endl;
        } else {
          // std::cout << lLine << std::endl;  // XXX
          std::cout << lLine << std::endl;
        }
      } // log2CardFk
    } // all dom sizes for S considered
    cmeasure_stop(&lMeasureBodySkew);
    std::cout << "  # runtime body skew(" << lSkew << "): " 
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
