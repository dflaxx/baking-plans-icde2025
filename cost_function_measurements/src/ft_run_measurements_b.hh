#pragma once

#include "dfinfra/standard_includes.hh"

#include "dfinfra/CrystalTimer.hh"

#include "gminfra/prime_table_1_1.hh"
#include "gminfra/tmath.hh"

extern "C" {
  #include "dfinfra/cbind_to_hw_thread.h"
  #include "dfinfra/cmeasure.h"
}

#include "cb_glob.hh"
#include "line.hh"

#include "HllErtl.hh"

#include "RelationRS.hh"
#include "plan_container.hh"

// more includes
// thus include last in main

template<typename Tattr, bool Ttest, bool Ttrace>
void
ft_run_measurements_b(const CbGlob& aCbGlob) {
  using attr_t = Tattr;
  using rng_t  = rng_tt<attr_t>;
  using hll_t  = HllErtl<attr_t>;
  using bun_t  = bun_tt<attr_t>;
  using R_t    = RelationRS2<bun_t>;
  using S_t    = RelationRS2<bun_t>;

  constexpr bool     lTraceInner    = false;
  constexpr uint     lLog2ChunkSize = 12;
  constexpr double   lHllFactor     = 1.1;


  df::infra::CrystalTimer lClock;
  lClock.start();
  constexpr uint64_t lMinHtDirSize  = (Ttest ? 17 : 191);
  uint64_vt lHashTableSizesR(aCbGlob.max_log_card_key() + 1);
  for(uint i = 0; i <= aCbGlob.max_log_card_key(); ++i) {
     const uint64_t lCardR = uint64_t{1} << i;
     lHashTableSizesR[i] = get_next_prime_1_1(std::max(lMinHtDirSize, lCardR));
  }
  lClock.stop();
  std::cout << "# calc hashtable sizes for R: " << lClock.duration_ns() << " [ns]" << std::endl;


  line_t lLine;
  rng_t  lRng;


  cmeasure_t lMeasure; // measure total time to produce a line
  cmeasure_start(&lMeasure);
  R_t R;
  // init R.k? to truely allocate physical memory
  const uint64_t lCapR = uint64_t{1} << aCbGlob.max_log_card_key();
  R.mem_alloc(lCapR);
  build_rel_key<bun_t>(R, lCapR);
  cmeasure_stop(&lMeasure);
  std::cout << "# time to build R: " << cmeasure_total_s(&lMeasure) << " [s]" << std::endl;

  // init S.k? to truely allocate physical memory
  cmeasure_start(&lMeasure);
  S_t S;
  const uint64_t lCapS = uint64_t{1} << aCbGlob.max_log_card_fk();
  S.mem_alloc(lCapS);
  S.mem_init();
  cmeasure_stop(&lMeasure);
  std::cout << "# time to alloc S: " << cmeasure_total_s(&lMeasure) << " [s]" << std::endl;


  if(nullptr != aCbGlob.os_ptr()) {
    line_t::print_header_line(aCbGlob.os());
  } else {
    // line_t::print_header_line(std::cout);  // XXX
  }

  lLine._attr_size = sizeof(attr_t);           // (4)

  assert(aCbGlob.max_log_dom_fk_reduction() <= aCbGlob.max_log_card_fk());
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
  cmeasure_start(&lMeasureLoop);
  cmeasure_t lMeasureBodySkew;
  for(uint lSkew = aCbGlob.min_skew(); lSkew <= aCbGlob.max_skew(); ++lSkew) {
    cmeasure_start(&lMeasureBodySkew);
    std::cout << "#  skew = " << lSkew << std::endl;
    for(uint lLog2DomSizeFk  = lLog2DomSizeBegin;
             lLog2DomSizeFk <= lLog2DomSizeLimit;
           ++lLog2DomSizeFk) {
      const uint     lLog2CardBegin = std::max<uint>(lLog2CardMin, lLog2DomSizeFk + lLog2DomRedMin);
      const uint     lLog2CardLimit = std::min<uint>(lLog2CardMax, lLog2DomSizeFk + lLog2DomRedMax);
      const uint64_t lCardLimit = uint64_t{1} << lLog2CardLimit;
      const uint64_t lDomSizeFk = uint64_t{1} << lLog2DomSizeFk;

      // std::cout << "#    lLog2DomSizeFk = " <<   lLog2DomSizeFk << std::endl;
      // std::cout << "#    lLog2Card      = " <<   lLog2CardBegin << " .. "
      //                                       <<   lLog2CardLimit << std::endl;

      cmeasure_start(&lMeasure);
      build_rel_fk<bun_t,rng_t>(S, lCardLimit, lDomSizeFk, lSkew, lRng); // sic!
      cmeasure_stop(&lMeasure);
      // std::cout << "#    time build S(" << lLog2DomSizeFk << "): " 
      //           << cmeasure_total_s(&lMeasure) 
      //           << " [s]" << std::endl;
      ++lCountHeavyS;

      for(uint lLog2CardFk = lLog2CardBegin; lLog2CardFk <= lLog2CardLimit; ++lLog2CardFk) {
        const uint64_t lCardFk       = uint64_t{1} << lLog2CardFk;
        const uint     lLog2DomRedFk = lLog2CardFk - lLog2DomSizeFk;

        if constexpr (lTraceInner) {
          std::cout << "#      lLog2CardFk = " << lLog2CardFk   << std::endl;
          std::cout << "#      lCardFk     = " << lCardFk       << std::endl;
          std::cout << "##     lLog2DomRed = " << lLog2DomRedFk << std::endl;
        }

        S.resize(lCardFk);

        lLine._no_rep_fk       = aCbGlob.get_no_rep(lLog2CardFk); // (3)
        lLine._log_card_fk     = lLog2CardFk;    // (5)
        lLine._card_fk         = lCardFk;        // (6)
        lLine._log_dom_red_fk  = lLog2DomRedFk;  // (7)
        lLine._log_dom_fk      = lLog2DomSizeFk; // (8)
        lLine._dom_fk          = lDomSizeFk ;    // (9)
        lLine._skew_fk         = lSkew;          // (10)

        assert(lLog2DomSizeFk >= aCbGlob.min_log_card_fk());
        assert(lLog2DomSizeFk <= aCbGlob.max_log_card_fk());
        assert(lLog2DomRedFk  >= aCbGlob.min_log_dom_fk_reduction());
        assert(lLog2DomRedFk  <= aCbGlob.max_log_dom_fk_reduction());
        assert((lLog2DomSizeFk + lLog2DomRedFk) == lLog2CardFk);
        assert(S.card() == lCardFk);

        hll_t lHll(14);
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

        lLine._nodv_fk_hll = lEstNoDvInt;         // (13)
        lLine._hashdir_size_fk = lHashtableSizeS; // (14)
        
        for(uint lLog2CardKey =  aCbGlob.min_log_card_key();
                 lLog2CardKey <= aCbGlob.max_log_card_key();
               ++lLog2CardKey) {
          const uint64_t lCardR = uint64_t{1} << lLog2CardKey;
          const uint64_t lHashtableSizeR = lHashTableSizesR[lLog2CardKey];
          // const uint64_t lHashtableSizeR = get_next_prime_1_1(std::max(lMinHtDirSize, lCardR));

          lLine._id = ((lLog2CardKey * 100 + lLog2CardFk) * 100 + lLog2DomSizeFk) * 10 + lSkew; // (1)
          lLine._no_rep_key = aCbGlob.get_no_rep(lLog2CardKey); // (2)
          lLine._log_card_key     = lLog2CardKey;               // (11)
          lLine._card_key         = lCardR;                     // (12)
          lLine._hashdir_size_key = lHashtableSizeR;            // (15)

          R.resize(lLog2CardKey);  // sic!

          ++lCount;
          cmeasure_start(&lMeasure);
          // ft_run_plans_b<R_t,S_t,Ttest,Ttrace>(R, S, lHashtableSizeR, lHashtableSizeS, lLog2ChunkSize, lLine, aCbGlob);
          cmeasure_stop(&lMeasure);
          lLine._time = cmeasure_total_s(&lMeasure); // (18)
          if(nullptr != aCbGlob.os_ptr()) {
             aCbGlob.os() << lLine << std::endl;
          } else {
            // std::cout << lLine << std::endl;  // XXX
          }
        }
      }
    }
    cmeasure_stop(&lMeasureBodySkew);
    std::cout << "  # runtime body skew(" << lSkew << "): " 
              << cmeasure_total_s(&lMeasureBodySkew) << " [s]" << std::endl;
  }
  cmeasure_stop(&lMeasureLoop);
  std::cout << "# no run     = " << lCount << std::endl;
  std::cout << "# no heavy s = " << lCountHeavyS << std::endl;
  std::cout << "# runtime loop: " << cmeasure_total_s(&lMeasureLoop) << " [s]" << std::endl;
}
