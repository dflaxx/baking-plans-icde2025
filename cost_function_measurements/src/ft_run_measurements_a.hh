#ifndef HASHJOIN2_FT_RUN_MEASUREMENTS_A_HH
#define HASHJOIN2_FT_RUN_MEASUREMENTS_A_HH
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
#include "ft_run_plans_a.hh"

#include "HllErtl.hh"

// more includes
// thus include last in main

template<typename Tattr, bool Ttest, bool Ttrace>
void
ft_run_measurements_a(const CbGlob& aCbGlob) {
  using attr_t = Tattr;
  using rng_t  = rng_tt<attr_t>;
  using hll_t  = HllErtl<attr_t>;
  using bun_t  = bun_tt<attr_t>;
  using R_t    = RelationRS2<bun_t>;
  using S_t    = RelationRS2<bun_t>;

  constexpr uint     lLog2ChunkSize = 12;
  constexpr uint64_t lMinHtDirSize  = (Ttest ? 17 : 191);
  constexpr double   lHllFactor     = 1.1;  // protect against HLL underestimates
  line_t lLine;
  rng_t  lRng;
  R_t    R;
  S_t    S;
  const uint64_t lCapR = uint64_t{1} << aCbGlob.max_log_card_key();
  const uint64_t lCapS = uint64_t{1} << aCbGlob.max_log_card_fk();
  R.mem_alloc(lCapR);
  for(attr_t i = 0; i < lCapR; ++i) {
    R[i] = {i, i};
  }
  // init R.k? to truely allocate physical memory
  S.mem_alloc(lCapS);
  // init S.k? to truely allocate physical memory

  cmeasure_t lMeasure; // measure total time to produce a line

  if(nullptr != aCbGlob.os_ptr()) {
    line_t::print_header_line(aCbGlob.os());
  } else {
    line_t::print_header_line(std::cout);
  }

  lLine._attr_size = sizeof(attr_t);           // (4)

  uint   lCount = 0;
  // loop 1: foreign key card
  for(uint lLogCardFk =  aCbGlob.min_log_card_fk();
           lLogCardFk <= aCbGlob.max_log_card_fk();
         ++lLogCardFk) {
    // std::cout << "# lLogCardFk = " << lLogCardFk << std::endl;
    lLine._no_rep_fk  = aCbGlob.get_no_rep(lLogCardFk); // (3)
    const uint64_t lCardFk = uint64_t{1} << lLogCardFk;
    const uint lMaxLog2DomReduction = std::min(lLogCardFk, aCbGlob.max_log_dom_fk_reduction());
    lLine._log_card_fk = lLogCardFk;  // (5)
    lLine._card_fk     = lCardFk;     // (6)
    // loop 2: domain reduction
    for(uint lLogDomRedFk = aCbGlob.min_log_dom_fk_reduction();
             lLogDomRedFk <= lMaxLog2DomReduction;
           ++lLogDomRedFk) {
      const uint lLogDomSizeFk = lLogCardFk - lLogDomRedFk;
      assert(lLogDomSizeFk <= aCbGlob.max_log_card_fk());
      lLine._log_dom_red_fk  = lLogDomRedFk;                 // (7)
      lLine._log_dom_fk      = lLogDomSizeFk;                // (8)
      lLine._dom_fk          = uint64_t{1} << lLogDomSizeFk; // (9)
      // loop 3: skew
      for(uint lSkew = aCbGlob.min_skew(); lSkew <= aCbGlob.max_skew(); ++lSkew) {
        lLine._skew_fk = lSkew;  // (10)
        // data generation: foreign key relation
        build_rel_fk<bun_t,rng_t>(S, lLine.card_fk(), lLine.dom_fk(), lLine.skew_fk(), lRng);
        // use HLL sketch to estimate number of distinct values in foreign key column
        hll_t lHll(14);  // hll_t(uint n: number of registers in log2)
        for(const auto& s : S.tuples()) {
          lHll.insert(ht::murmur_hash<attr_t>(s.a));  // XXX why insert hash(x) and not x?
        }
        const double   lEstNoDvDouble = lHll.estimate();
        const uint64_t lEstNoDvInt    = (uint64_t) std::ceil(lEstNoDvDouble);
        const uint64_t lHashDirSzMinS = (uint64_t) std::ceil(lEstNoDvDouble * lHllFactor);
        const uint64_t lHashtableSizeS = get_next_prime_1_1(std::max(lMinHtDirSize, lHashDirSzMinS));  // get prime number that is larger than a lower bound
        lLine._nodv_fk_hll = lEstNoDvInt;         // (13)
        lLine._hashdir_size_fk = lHashtableSizeS; // (14)
        
        // loop 4: key card
        for(uint lLogCardKey =  aCbGlob.min_log_card_key();
                 lLogCardKey <= aCbGlob.max_log_card_key();
               ++lLogCardKey) {
          const uint64_t lCardR = uint64_t{1} << lLogCardKey;
          const uint64_t lHashtableSizeR = get_next_prime_1_1(std::max(lMinHtDirSize, lCardR));
          lLine._no_rep_key = aCbGlob.get_no_rep(lLogCardKey); // (2)
          lLine._id = ((lLogCardKey * 100 + lLogCardFk) * 100 + lLogDomSizeFk) * 10 + lSkew; // (1)
          lLine._log_card_key     = lLogCardKey;     // (11)
          lLine._card_key         = lCardR;          // (12)
          lLine._hashdir_size_key = lHashtableSizeR; // (15)
          // data generation: key relation
          build_rel_key<bun_t>(R, lLine._card_key);
          cmeasure_start(&lMeasure);
          ft_run_plans_a<R_t,S_t,Ttest,Ttrace>(R, S, lHashtableSizeR, lHashtableSizeS, lLog2ChunkSize, lLine, aCbGlob);
          cmeasure_stop(&lMeasure);
          lLine._time = cmeasure_total_s(&lMeasure); // (18)
          if(nullptr != aCbGlob.os_ptr()) {
             aCbGlob.os() << lLine << std::endl;
          } else {
            std::cout << lLine << std::endl;
          }
        }
      }
    }
  }
  std::cout << "# no run = " << lCount << std::endl;
}

#endif

