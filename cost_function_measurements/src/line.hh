#ifndef LINE_HH
#define LINE_HH
#pragma once

#include "dfinfra/standard_includes.hh"
//#include "hj_types.hh"
#include "algebra_types_v2.hh"

struct rt_stat_t {
  uint64_t _cc_min;
  uint64_t _cc_max;
  double   _cc_med;
  double   _cc_avg;
  uint64_t _cc_mx2;
  uint64_t _cc_fst;

  std::ostream& print(std::ostream& os, const char aSep) const;
};

std::ostream& operator<<(std::ostream& os, const rt_stat_t& aRtStat);

/*
 * line contains major information
 * additional information may be written for the different
 * join algorithms to a separate file * using _id for a potential join
 * with other result files
 */

struct line_t {
  uint      _id;
  uint      _no_rep_key; // for R as build/probe
  uint      _no_rep_fk;  // for S as build/probe
  uint      _attr_size;
  double    _log_card_fk;
  uint64_t  _card_fk;  
  uint      _log_dom_red_fk;
  double    _log_dom_fk;
  uint64_t  _dom_fk;
  uint      _skew_fk;
  double    _log_card_key;
  uint64_t  _card_key;
  uint64_t  _nodv_fk_hll;
  uint64_t  _hashdir_size_fk;  // hash table size (number of buckets) on FK side
  uint64_t  _hashdir_size_key;  // hash table size (number of buckets) on key side
  uint64_t  _card_sj_sr;  // semi-join size (R lsjoin S), R key side, S foreign key side
                          // -> how many unique probe keys from R find at least one join partner
                          //    in non-unique build side S.
                          // = number of nested output tuples created: [r, {s0, s1, ...}]
  uint64_t  _card_result; // join size (S join R), R key side, S foreign key side
  double    _time;        // runtime [s] to produce this line

  rt_stat_t _rt_stat_build[kNoPlan]; // build times statistics per plan
  rt_stat_t _rt_stat_probe[kNoPlan]; // probe times statistics per plan

  uint64_vt _mem_stat_ht;  // memory consumption of hash table (directory + reservoir) per plan

  uint64_vt _ht_num_nodes_main;  // number of (main) collision chain nodes of chaining and 3d hash table
  uint64_vt _ht_num_nodes_sub;   // number of sub chain nodes; 3d hash table only, set to 0 for chaining HT

  double_vt _ht_dir_fill;  // fraction of non-empty hash table directory entries

  uint64_vt _cnt_build;  // _count from algebraic build operator
  uint64_vt _cnt_top;    //                       top

  line_t() : _id(0),
             _no_rep_key(0),
             _no_rep_fk(0),
             _attr_size(0),
             _log_card_fk(0),
             _card_fk(0),
             _log_dom_red_fk(0),
             _log_dom_fk(0),
             _dom_fk(0),
             _skew_fk(0),
             _log_card_key(0.0),
             _card_key(0),
             _nodv_fk_hll(0),
             _hashdir_size_fk(0),
             _hashdir_size_key(0),
             _card_sj_sr(0),
             _card_result(0),
             _time(0),
             _rt_stat_build(),
             _rt_stat_probe(),
             _mem_stat_ht(kNoPlan, 0),
             _ht_num_nodes_main(kNoPlan, 0),
             _ht_num_nodes_sub(kNoPlan, 0),
             _ht_dir_fill(kNoPlan, 0.0),
             _cnt_build(kNoPlan, 0),
             _cnt_top(kNoPlan, 0) {}

  inline uint     id()             const { return _id;       }
  inline uint     no_rep_key()     const { return _no_rep_key; }
  inline uint     no_rep_fk()      const { return _no_rep_fk; }
  inline uint     attr_size()      const { return _attr_size; }
  inline double   log_card_fk()    const { return _log_card_fk; }
  inline uint64_t card_fk()        const { return _card_fk;  }
  inline uint     log_dom_red_fk() const { return _log_dom_red_fk; }
  inline double   log_dom_fk()     const { return _log_dom_fk; }
  inline uint64_t dom_fk()         const { return _dom_fk;   }
  inline uint     skew_fk()        const { return _skew_fk;  }
  inline double   log_card_key()   const { return _log_card_key; }
  inline uint64_t card_key()       const { return _card_key; }
  inline uint64_t nodv_fk_hll()    const { return _nodv_fk_hll; }

  inline uint64_t hashdir_size_fk()  const { return _hashdir_size_fk;  }
  inline uint64_t hashdir_size_key() const { return _hashdir_size_key; }

  inline uint64_t card_sj_sr()     const { return _card_sj_sr; }
  inline uint64_t card_result()    const { return _card_result; }

  inline double   time()           const { return _time; }

  std::ostream& print(std::ostream& os, const char aSep) const;
  std::ostream& print_build_rs(std::ostream& os, const char* aSep) const;
  std::ostream& print_build_sr(std::ostream& os, const char* aSep) const;
  std::ostream& print_probe_rs(std::ostream& os, const char* aSep) const;
  std::ostream& print_probe_sr(std::ostream& os, const char* aSep) const;
  // TODO: need multiple header lines for each case
  static std::ostream& print_header_line(std::ostream& os);
};

std::ostream& operator<<(std::ostream& os, const line_t& aLine);

#endif
