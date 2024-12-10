/*
 * Original file name: /home/moer/src/HashJoin2/main_measure.cpp
 * Author: Guido Moerkotte
 *
 * Measure runtime (build and probe separately) of a single key/foreign key join 
 * using two hash join algorithms (CH and 3D)
 * in different implementations variants (unpacked/packed, no prefetching/AMAC).
 *
 * Distinguish four measurement configurations:
 * - bld+bun: build on unique key side
 * - bld+bnu: build on non-unique foreign key side
 * - prb+bun: probe with fk side on hash table built on unique key side
 * - prb+bnu: probe with key side on hash table built on non-unique fk side
 *
 * Input parameters:
 * All inputs prefixed with "l" are in log2.
 * - key side cardinality min/max
 * - fk side cardinality min/max
 * - fk domain and domain reduction min/max:
 *   fk_domain_size = fk_card / fk_domain_reduction (in normal)
 *   fk_domain_size = fk_card - fk_domain_reduciton (in log2)
 * - skew: 0, 1 (uniform or zipf(1.0) distribution
 * - key and fk scale factor: scale (non-log) cardinalities up or down by this factor
 */
#include "dfinfra/standard_includes.hh"
#include "gminfra/prime_table_1_1.hh"

extern "C" {
  #include "dfinfra/cbind_to_hw_thread.h"
  #include "dfinfra/cmeasure.h"
}

#include "cb_glob.hh"
#include "line.hh"
//#include "arg.hh"

#include "RelationRS.hh"
#include "hj_htnode_content.hh"
#include "hashtable_chained_v2.hh"
#include "hashtable_nested_v2.hh"
#include "meas_eval.hh"
#include "hasht.hh"
#include "HllErtl.hh"

#include "algebra_v2.hh"
#include "plan.hh"

// #include "ft_run_plans_b.hh"
#include "ft_run_measurements_b.hh"

#include "ft_run_measure_build_bun.hh"
#include "ft_run_measure_build_bnu.hh"
#include "ft_run_measure_probe_bnu.hh"
#include "ft_run_measure_probe_bun.hh"

#include <sys/resource.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Winline"
#include "lib/CLI11.hpp"

void cli_option_setup(CLI::App& aApp, CbGlob& aCb) {
  aApp.add_option("--l_key_min" , aCb.min_log_card_key()        , "minimum log2 card key")                  ->default_val(0)  ->check(CLI::Range(0, 30));
  aApp.add_option("--l_key_max" , aCb.max_log_card_key()        , "maximum log2 card key")                  ->default_val(30) ->check(CLI::Range(0, 30));
  aApp.add_option("--l_fk_min"  , aCb.min_log_card_fk()         , "minimum log2 card fk")                   ->default_val(0)  ->check(CLI::Range(0, 30));
  aApp.add_option("--l_fk_max"  , aCb.max_log_card_fk()         , "maximum log2 card fk")                   ->default_val(30) ->check(CLI::Range(0, 30));
  aApp.add_option("--l_red_min" , aCb.min_log_dom_fk_reduction(), "minimum log2 fk dom reduction")          ->default_val(0)  ->check(CLI::Range(0, 10));
  aApp.add_option("--l_red_max" , aCb.max_log_dom_fk_reduction(), "maximum log2 fk dom reduction")          ->default_val(10) ->check(CLI::Range(0, 10));
  aApp.add_option("--l_dom_min" , aCb.min_log_domsize_fk()      , "minimum log2 fk dom size")               ->default_val(0)  ->check(CLI::Range(0, 30));
  aApp.add_option("--l_dom_max" , aCb.max_log_domsize_fk()      , "maximum log2 fk dom size")               ->default_val(10) ->check(CLI::Range(0, 30));
  aApp.add_option("--skew_min"  , aCb.min_skew()                , "minimum skew")                           ->default_val(0)  ->check(CLI::Range(0, 1));
  aApp.add_option("--skew_max"  , aCb.max_skew()                , "maximum skew")                           ->default_val(1)  ->check(CLI::Range(0, 1));
  aApp.add_option("--hwtno"     , aCb.hwtno()                   , "hardware thread number")                 ->default_val(0);
  aApp.add_option("--od"        , aCb.directory()               , "output directory")                       ->default_val("./")->check(CLI::ExistingDirectory);
  aApp.add_option("--of"        , aCb.filename()                , "output file name prefix")                ->default_val("zOut");
  aApp.add_option("-m,--measure", aCb.measure()                 , "what to measure <build|probe>_<bun|bnu>")->required();
  //
  aApp.add_option("--key_scale_factor", aCb.card_key_factor()   , "factor by which key card is scaled up or down: card = 2^(l_key) * key_scale_factor")->default_val(1.0);
  aApp.add_option("--fk_scale_factor" , aCb.card_fk_factor()    , "factor by which fk card is scaled up or down")->default_val(1.0);
}

int
main(const int argc, const char* argv[]) {
  CbGlob lCbGlob;  // global control block; stores experiment parameters
  CLI::App lApp;
  cli_option_setup(lApp, lCbGlob);
  CLI11_PARSE(lApp, argc, argv);

  cbind_to_hw_thread(lCbGlob.hwtno(), 0);

  std::cout << "control block (" << lCbGlob.to_string() << "):" << std::endl;
  lCbGlob.print(std::cout);
  std::cout << "--- end cb ---" << std::endl;

  if (lCbGlob.min_log_card_key() > lCbGlob.max_log_card_key()) {
    std::cout << "Error: min_log_card_key > max_log_card_key" << std::endl;
    return -1;
  }
  if (lCbGlob.min_log_card_fk() > lCbGlob.max_log_card_fk()) {
    std::cout << "Error: min_log_card_fk > max_log_card_fk" << std::endl;
    return -1;
  }
  if (lCbGlob.min_log_dom_fk_reduction() > lCbGlob.max_log_dom_fk_reduction()) {
    std::cout << "Error: min_log_dom_fk_reduction > max_log_dom_fk_reduction" << std::endl;
    return -1;
  }
  if (lCbGlob.min_log_domsize_fk() > lCbGlob.max_log_domsize_fk()) {
    std::cout << "Error: min_log_domsize_fk > max_log_domsize_fk" << std::endl;
    return -1;
  }
  if (30 < lCbGlob.max_log_card_key()) {
    std::cout << "max_log_card_key = " << lCbGlob.max_log_card_key() << " too high." << std::endl;
    return -1;
  }
  if (30 < lCbGlob.max_log_card_fk()) {
    std::cout << "max_log_card_fk = " << lCbGlob.max_log_card_fk() << " too high." << std::endl;
    return -1;
  }
  if ((1 < lCbGlob.min_skew()) || (1 < lCbGlob.max_skew()) || lCbGlob.min_skew() > lCbGlob.max_skew()) {
    std::cout << "bad skew definition: " 
              << lCbGlob.min_skew() 
              << " <= " 
              << lCbGlob.max_skew() 
              << " <= 1 " 
              << std::endl;
    return -1;
  }
  if ((0 >= lCbGlob.card_key_factor()) || (0 >= lCbGlob.card_fk_factor())) {
    std::cout << "bad card factor definition: "
      << "card_key_factor = " << lCbGlob.card_key_factor() << ", card_fk_factor = " << lCbGlob.card_fk_factor()
      << std::endl;
    return -1;
  }

  if ((0 < lCbGlob.filename().size()) && (!lCbGlob.alloc_os())) {
    std::cout << "Problems allocating output stream" << std::endl;
    return -1;
  }


  // constexpr bool lRunMeasurement = false;
  // if(lRunMeasurement) {
  //   ft_run_measurements_b<uint64_t, false, false>(lCbGlob);
  // }

  /*
   * Signature:
   * template<typename Tattr, bool Ttest, bool Ttrace>
   * void ft_run_measure_XXX(const CbGlob& aCbGlob);
   */
  if ("build_bun" == lCbGlob.measure()) {
    ft_run_measure_build_bun<uint64_t, false, false>(lCbGlob);
  } else
  if ("build_bnu" == lCbGlob.measure()) {
    ft_run_measure_build_bnu<uint64_t, false, false>(lCbGlob);
  } else
  if ("probe_bnu" == lCbGlob.measure()) {
    ft_run_measure_probe_bnu<uint64_t, false, false>(lCbGlob);
  } else
  if ("probe_bun" == lCbGlob.measure()) {
    ft_run_measure_probe_bun<uint64_t, false, false>(lCbGlob);
  } else
  if ("gen_heavy" == lCbGlob.measure()) {
     // make sure ft_run_plans_b is commented out in ft_run_measurements_b
     ft_run_measurements_b<uint64_t, false, false>(lCbGlob);
  } else {
    std::cout << "unknown measure: '" << lCbGlob.measure() << "'." << std::endl;
    return -1;
  }

  struct rusage lRusage;
  int lRc = getrusage(RUSAGE_SELF, &lRusage);
  if (lRc) {
     std::cout << "failed to call getrusage" << std::endl;
  }
  std::cout << "max resident memory [GB]: " << ((double) lRusage.ru_maxrss / ((double) (1024*1024))) << std::endl;
  std::cout << "num of major page faults: " << lRusage.ru_majflt << std::endl;
  std::cout << "num of swaps:             " << lRusage.ru_nswap << std::endl;
  /*
   * https://www.gnu.org/software/libc/manual/html_node/Resource-Usage.html
   *
   * struct rusage
   *
   * long int ru_maxrss
   * The maximum resident set size used, in kilobytes. That is, the maximum number of kilobytes of physical memory that processes used simultaneously.
   *
   * long int ru_majflt
   * The number of page faults which were serviced by doing I/O.
   *
   * long int ru_nswap
   * The number of times processes was swapped entirely out of main memory.
   */

  return 0;
}

#pragma GCC diagnostic pop
