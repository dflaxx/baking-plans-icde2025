#ifndef HASHJOIN_PLAN_HH
#define HASHJOIN_PLAN_HH
#pragma once

#include "dfinfra/CrystalTimer.hh"
#include "algebra_v2.hh"
#include "RelationRS.hh"
#include "line.hh"
#include "meas_eval.hh"

/*
 * plans for joins
 */

// template<...,bool Timpr>
template<typename Tr, typename Ts,
         join_order_et aJoinOrder,
         bool          NESTED,
         prefetch_et   PREFETCH,
         bool          aManyPtr,
         bool          Ttrace,
         bool          Ttest>
class PlanJoin {
  static_assert(prefetch_et::kPrefetchNone <= PREFETCH && PREFETCH < prefetch_et::kNumPrefetch);
  public:
    using R_t            = Tr;
    using S_t            = Ts;
    using tuple_r_t      = typename Tr::tuple_t;
    using tuple_s_t      = typename Ts::tuple_t;
    using attr_t         = typename tuple_r_t::attr_t;
    using bun_t          = bun_tt<attr_t>;
    using hash_r_t       = hash_R_tt<bun_t>;
    using hash_s_t       = hash_S_tt<bun_t>;
    using eq_r_t         = eq_join_attr_R<bun_t>;
    using eq_s_t         = eq_join_attr_S<bun_t>;
    using join_pred_rs_t = join_pred_RS<bun_t,bun_t>;
    using join_pred_sr_t = join_pred_SR<bun_t,bun_t>;
    using concat_rs_t    = concat_RS<bun_t>;
    using concat_sr_t    = concat_SR<bun_t>;
    using join_res_t     = join_res_tt<bun_t>;

    static constexpr uint NO_PTR_CH      = aManyPtr ? 3 : 1;
    static constexpr uint NO_PTR_3d_MAIN = aManyPtr ? 5 : 1;
    static constexpr uint NO_PTR_3d_SUB  = aManyPtr ? 3 : 1;
    //static constexpr plan_et PLAN_ID = (plan_et) ((NESTED * 8) + (AMAC * 4) + (aManyPtr * 2) + aJoinOrder);
    // kNumPrefetch = 3
    static constexpr plan_et PLAN_ID = (plan_et) (((((NESTED * 3 + PREFETCH) << 1) + aManyPtr) << 1)  + aJoinOrder);
    /*
     * N P M J
     * 0 0 0 0 -> 0   ch nopf upk RS
     * 0 0 0 1 -> 1   ch nopf upk SR
     * 0 0 1 0 -> 2   ch nopf pkd RS
     * 0 0 1 1 -> 3   ch nopf pkd SR
     * 0 1 0 0 -> 4   ch roll upk RS
     * 0 1 0 1 -> 5   ch roll upk SR
     * 0 1 1 0 -> 6   ch roll pkd RS
     * 0 1 1 1 -> 7   ch roll pkd SR
     * 0 2 0 0 -> 8   ch amac upk RS
     * 0 2 0 1 -> 9   ch amac upk SR
     * 0 2 1 0 -> 10  ch amac pkd RS
     * 0 2 1 1 -> 11  ch amac pkd SR
     * 1 0 0 0 -> 12  3d nopf upk RS
     * 1 0 0 1 -> 13  3d nopf upk SR
     * 1 0 1 0 -> 14  3d nopf pkd RS
     * 1 0 1 1 -> 15  3d nopf pkd SR
     * 1 1 0 0 -> 16  3d roll upk RS
     * 1 1 0 1 -> 17  3d roll upk SR
     * 1 1 1 0 -> 18  3d roll pkd RS
     * 1 1 1 1 -> 19  3d roll pkd SR
     * 1 2 0 0 -> 20  3d amac upk RS
     * 1 2 0 1 -> 21  3d amac upk SR
     * 1 2 1 0 -> 22  3d amac pkd RS
     * 1 2 1 1 -> 23  3d amac pkd SR
    */


    // join R and S: build on S, probe with R
    struct bundle_rs_t {
      using rel_build_t = Ts;
      using rel_probe_t = Tr;
      using tuple_bld_t = tuple_s_t;
      using tuple_prb_t = tuple_r_t;
      using hash_bld    = hash_s_t;
      using hash_prb    = hash_r_t;
      using eq_bld      = eq_s_t;
      using join_pred   = join_pred_rs_t;
      using concat      = concat_rs_t;
      using join_res    = join_res_t;
      static constexpr bool build_uniq = false;

      using join_ch_build_orig_t = AlgOpJoinHashChOrigBuild2<hash_bld, NO_PTR_CH>;
      template<typename Tconsumer> 
      using join_ch_probe_orig_t = AlgOpJoinHashChOrigProbe2<Tconsumer, 
                                                             join_ch_build_orig_t,
                                                             hash_prb,
                                                             join_pred,
                                                             concat, 
                                                             build_uniq>;

      using join_ch_build_rp_t = AlgOpJoinHashChRpBuild2<hash_bld, NO_PTR_CH>;
      template<typename Tconsumer> 
      using join_ch_probe_rp_t = AlgOpJoinHashChRpProbe2<    Tconsumer, 
                                                             join_ch_build_rp_t,
                                                             hash_prb,
                                                             join_pred,
                                                             concat, 
                                                             build_uniq>;

      using join_ch_build_amac_t = AlgOpJoinHashChAmacBuild2<hash_bld, NO_PTR_CH>;
      template<typename Tconsumer> 
      using join_ch_probe_amac_t = AlgOpJoinHashChAmacProbe2<Tconsumer, 
                                                             join_ch_build_amac_t,
                                                             hash_prb,
                                                             join_pred,
                                                             concat, 
                                                             build_uniq>;

      using join_3d_build_orig_t = AlgOpJoinHash3dOrigBuild2<hash_bld, eq_bld, NO_PTR_3d_MAIN, NO_PTR_3d_SUB>;
      template<typename Tconsumer>
      using join_3d_probe_orig_t = AlgOpJoinHash3dOrigProbe2<Tconsumer,
                                                             join_3d_build_orig_t,
                                                             hash_prb,
                                                             join_pred,
                                                             concat>;

      using join_3d_build_rp_t = AlgOpJoinHash3dRpBuild2<hash_bld, eq_bld, NO_PTR_3d_MAIN, NO_PTR_3d_SUB>;
      template<typename Tconsumer>
      using join_3d_probe_rp_t = AlgOpJoinHash3dRpProbe2<    Tconsumer,
                                                             join_3d_build_rp_t,
                                                             hash_prb,
                                                             join_pred,
                                                             concat>;

      using join_3d_build_amac_t = AlgOpJoinHash3dAmacBuild2<hash_bld, eq_bld, NO_PTR_3d_MAIN, NO_PTR_3d_SUB>;
      template<typename Tconsumer>
      using join_3d_probe_amac_t = AlgOpJoinHash3dAmacProbe2<Tconsumer,
                                                             join_3d_build_amac_t,
                                                             hash_prb,
                                                             join_pred,
                                                             concat>;

      static inline const S_t& build([[maybe_unused]] const R_t& R, const S_t&S) { return S; };
      static inline const R_t& probe(const R_t& R, [[maybe_unused]] const S_t&S) { return R; };
      static inline       uint no_rep_build([[maybe_unused]] const uint aNoRepR, const uint aNoRepS) { return aNoRepS; }
      static inline       uint no_rep_probe(const uint aNoRepR, [[maybe_unused]] const uint aNoRepS) { return aNoRepR; }
    };

    // join S and R: build on R, probe with S
    struct bundle_sr_t {
      using rel_build_t = Tr;
      using rel_probe_t = Ts;
      using tuple_bld_t = tuple_r_t;
      using tuple_prb_t = tuple_s_t;
      using hash_bld    = hash_r_t;
      using hash_prb    = hash_s_t;
      using eq_bld      = eq_r_t;
      using join_pred   = join_pred_sr_t;
      using concat      = concat_sr_t;
      using join_res    = join_res_t;
      static constexpr bool build_uniq = true;

      using join_ch_build_orig_t = AlgOpJoinHashChOrigBuild2<hash_bld, NO_PTR_CH>;
      template<typename Tconsumer> 
      using join_ch_probe_orig_t = AlgOpJoinHashChOrigProbe2<Tconsumer, 
                                                             join_ch_build_orig_t, 
                                                             hash_prb, 
                                                             join_pred, 
                                                             concat,
                                                             build_uniq>;

      using join_ch_build_rp_t = AlgOpJoinHashChRpBuild2<hash_bld, NO_PTR_CH>;
      template<typename Tconsumer> 
      using join_ch_probe_rp_t = AlgOpJoinHashChRpProbe2<    Tconsumer, 
                                                             join_ch_build_rp_t, 
                                                             hash_prb, 
                                                             join_pred, 
                                                             concat,
                                                             build_uniq>;

      using join_ch_build_amac_t = AlgOpJoinHashChAmacBuild2<hash_bld, NO_PTR_CH>;
      template<typename Tconsumer> 
      using join_ch_probe_amac_t = AlgOpJoinHashChAmacProbe2<Tconsumer, 
                                                             join_ch_build_amac_t, 
                                                             hash_prb, 
                                                             join_pred, 
                                                             concat,
                                                             build_uniq>;

      using join_3d_build_orig_t = AlgOpJoinHash3dOrigBuild2<hash_bld, eq_bld, NO_PTR_3d_MAIN, NO_PTR_3d_SUB>;
      template<typename Tconsumer>
      using join_3d_probe_orig_t = AlgOpJoinHash3dOrigProbe2<Tconsumer,
                                                             join_3d_build_orig_t,
                                                             hash_prb,
                                                             join_pred,
                                                             concat>;

      using join_3d_build_rp_t = AlgOpJoinHash3dRpBuild2<hash_bld, eq_bld, NO_PTR_3d_MAIN, NO_PTR_3d_SUB>;
      template<typename Tconsumer>
      using join_3d_probe_rp_t = AlgOpJoinHash3dRpProbe2<    Tconsumer,
                                                             join_3d_build_rp_t,
                                                             hash_prb,
                                                             join_pred,
                                                             concat>;

      using join_3d_build_amac_t = AlgOpJoinHash3dAmacBuild2<hash_bld, eq_bld, NO_PTR_3d_MAIN, NO_PTR_3d_SUB>;
      template<typename Tconsumer>
      using join_3d_probe_amac_t = AlgOpJoinHash3dAmacProbe2<Tconsumer,
                                                             join_3d_build_amac_t,
                                                             hash_prb,
                                                             join_pred,
                                                             concat>;

      static inline const R_t& build(const R_t& R, [[maybe_unused]] const S_t&S) { return R; };
      static inline const S_t& probe([[maybe_unused]] const R_t& R, const S_t&S) { return S; };
      static inline       uint no_rep_build(const uint aNoRepR, [[maybe_unused]] const uint aNoRepS) { return aNoRepR; }
      static inline       uint no_rep_probe([[maybe_unused]] const uint aNoRepR, const uint aNoRepS) { return aNoRepS; }
    };
    using bundle_t = typename std::conditional<kJoinOrderRS == aJoinOrder, bundle_rs_t, bundle_sr_t>::type;

    using rel_build_t = typename bundle_t::rel_build_t;
    using rel_probe_t = typename bundle_t::rel_probe_t;

    using top_res = AlgOpTop2<join_res_t, Ttest, Ttrace>;
  
    using join_build = typename std::conditional<NESTED,
                         typename std::conditional<prefetch_et::kPrefetchNone == PREFETCH,
                                   typename bundle_t::join_3d_build_orig_t,
                                   typename std::conditional<prefetch_et::kRollingPrefetching == PREFETCH,
                                            typename bundle_t::join_3d_build_rp_t,
                                            typename bundle_t::join_3d_build_amac_t>::type >::type,
                         typename std::conditional<prefetch_et::kPrefetchNone == PREFETCH,
                                  typename bundle_t::join_ch_build_orig_t,
                                  typename std::conditional<prefetch_et::kRollingPrefetching == PREFETCH,
                                           typename bundle_t::join_ch_build_rp_t,
                                           typename bundle_t::join_ch_build_amac_t>::type >::type >::type;

    using join_probe = typename std::conditional<NESTED,
                         typename std::conditional<prefetch_et::kPrefetchNone == PREFETCH,
                                  typename bundle_t::template join_3d_probe_orig_t<top_res>,
                                  typename std::conditional<prefetch_et::kRollingPrefetching == PREFETCH,
                                           typename bundle_t::template join_3d_probe_rp_t<top_res>,
                                           typename bundle_t::template join_3d_probe_amac_t<top_res> >::type >::type,
                         typename std::conditional<prefetch_et::kPrefetchNone == PREFETCH,
                                  typename bundle_t::template join_ch_probe_orig_t<top_res>, 
                                  typename std::conditional<prefetch_et::kRollingPrefetching == PREFETCH,
                                           typename bundle_t::template join_ch_probe_rp_t<top_res>,
                                           typename bundle_t::template join_ch_probe_amac_t<top_res> >::type >::type >::type;

    using scan_build = AlgOpScan2<join_build>;
    using scan_probe = AlgOpScan2<join_probe>;

  private:
    PlanJoin(const PlanJoin&) = delete;
    PlanJoin& operator=(const PlanJoin&) = delete;
  public:
    PlanJoin() : _top(nullptr), 
                 _build(nullptr), 
                 _probe(nullptr), 
                 _scan_build(nullptr), 
                 _scan_probe(nullptr),
                 _3d_mat_probe(nullptr),
                 _3d_scan_mat(nullptr),
                 _3d_unnest(nullptr),
                 _meas_build(7),
                 _meas_probe(7),
                 _meas_unnest(7),
                 _card_sj_sr(0) {}
    ~PlanJoin();
  public:
    void build_plan(const size_t aHashDirSize, const uint32_t aLog2ChunkSize);
    void run_once(const R_t& R, const S_t& S);
    void run_many(const R_t& R, const S_t& S, line_t& aLine);
  public:
    // run plan in parts:
    // 1. build part a: build & clear the hash table aNoRep times,
    //                  then build it again and leave it filled for subsequent probe.
    // 2. probe: probe hash table from build part a aNoRep times
    // 3. build part b: clear the hash table from build part a,
    //                  adding to the total time for build
    void     run_build_part_a(const rel_build_t& R_bld,
                              const uint aNoRep, meas_eval_t& aMeasBuild);
    uint64_t run_probe(const rel_probe_t& R_prb,
                       const uint aNoRep, meas_eval_t& aMeasProbe);
    void     run_build_part_b(meas_eval_t& aMeasBuild);
  public: 
    inline auto     result()      const { return _top->result(); }
    inline uint64_t top_count()   const { return _top->count(); }

    inline uint64_t meas_build()  const { return _meas_build; }
    inline uint64_t meas_probe()  const { return _meas_probe; }
    inline uint64_t meas_unnest() const { return _meas_unnest; }

    inline uint64_t c_mid()       const { return _card_sj_sr; }
    inline uint64_t card_sj_sr()  const { return _card_sj_sr; }
  public:
    inline const join_build* build() const { return _build; }
    inline const join_probe* probe() const { return _probe; }
    inline size_t mem_consumption_ht() const { return _build->mem_consumption_ht(); }
  public:
    void clear();
  private:
    using CrystalTimer = df::infra::CrystalTimer;
  private:
    top_res*    _top;
    join_build* _build;
    join_probe* _probe;
    scan_build* _scan_build;
    scan_probe* _scan_probe;
    void*       _3d_mat_probe; // only for measuring unnest after materializing of nested probe result
    void*       _3d_scan_mat;
    void*       _3d_unnest;
    uint64_t    _meas_build;
    uint64_t    _meas_probe;
    uint64_t    _meas_unnest;
    uint64_t    _card_sj_sr;
};

template<typename Tr, typename Ts,
         join_order_et aJoinOrder,
         bool          NESTED,
         prefetch_et   PREFETCH,
         bool          aManyPtr,
         bool          Ttrace,
         bool          Ttest>
PlanJoin<Tr,Ts,aJoinOrder,NESTED,PREFETCH,aManyPtr,Ttrace,Ttest>::~PlanJoin() {
  if(nullptr != _top)   { delete _top; }
  if(nullptr != _build) { delete _build; }
  if(nullptr != _probe) { delete _probe; }
  if(nullptr != _scan_build) { delete _scan_build; }
  if(nullptr != _scan_probe) { delete _scan_probe; }
}

template<typename Tr, typename Ts,
         join_order_et aJoinOrder,
         bool          NESTED,
         prefetch_et   PREFETCH,
         bool          aManyPtr,
         bool          Ttrace,
         bool          Ttest>
void
PlanJoin<Tr,Ts,aJoinOrder,NESTED,PREFETCH,aManyPtr,Ttrace,Ttest>
::build_plan(const size_t aHashDirSize, const uint32_t aLog2ChunkSize) {
   clear();
   _build = new join_build(aHashDirSize, aLog2ChunkSize);
   _scan_build = new scan_build(_build);
   _top = new top_res(std::cout);
   _probe = new join_probe(_top, _build);
   _scan_probe = new scan_probe(_probe);
}

template<typename Tr, typename Ts,
         join_order_et aJoinOrder,
         bool          NESTED,
         prefetch_et   PREFETCH,
         bool          aManyPtr,
         bool          Ttrace,
         bool          Ttest>
void
PlanJoin<Tr,Ts,aJoinOrder,NESTED,PREFETCH,aManyPtr,Ttrace,Ttest>
::run_once(const R_t& R, const S_t& S) {
  if constexpr (Ttrace) {
    std::cout << "run_once ("
              << (kJoinOrderRS == aJoinOrder ? "R join S" : "S join R")
              << ") [NESTED = " << NESTED
              << ", IMPR = " << aManyPtr
              << ", PREFETCH = " << PREFETCH
              << "]:" << std::endl;
  }
  CrystalTimer lClock;
  lClock.start();
  _scan_build->mem_alloc();
  _scan_build->mem_init();
  _scan_build->run(bundle_t::build(R,S));
  lClock.stop();
  _meas_build = lClock.cycles();

  lClock.start();
  _scan_probe->mem_alloc();
  _scan_probe->mem_init();
  _scan_probe->run(bundle_t::probe(R,S));
  _scan_probe->clear();
  _scan_probe->mem_free();
  lClock.stop();
  _meas_probe = lClock.cycles();

  lClock.start();
  _scan_build->clear();
  _scan_build->mem_free();
  lClock.stop();
  _meas_build += lClock.cycles();

  if constexpr (Ttrace) {
    std::cout << "  count:" << std::endl;
    std::cout << "  scan build: " << _scan_build->count() << std::endl;
    std::cout << "  scan probe: " << _scan_probe->count() << std::endl;
    std::cout << "  top: count: " << _top->count() << std::endl;
    std::cout << "  cycles build: " << _meas_build << std::endl;
    std::cout << "  cycles probe: " << _meas_build << std::endl;
    std::cout << "  hashtable:" << std::endl;
    _build->hashtable().print(std::cout) << std::endl;
  }
}

template<typename Tr, typename Ts,
         join_order_et aJoinOrder,
         bool          NESTED,
         prefetch_et   PREFETCH,
         bool          aManyPtr,
         bool          Ttrace,
         bool          Ttest>
void
PlanJoin<Tr,Ts,aJoinOrder,NESTED,PREFETCH,aManyPtr,Ttrace,Ttest>
::run_many(const R_t& R, const S_t& S, line_t& aLine) {
  if constexpr (Ttrace) {
    std::cout << "run_many ("
              << (kJoinOrderRS == aJoinOrder ? "R join S" : "S join R")
              << ") [NESTED = " << NESTED
              << ", IMPR = " << aManyPtr
              << ", PREFETCH = " << PREFETCH
              << "]: plan_id: " << PLAN_ID << std::endl;
  }

  const uint lNoRepBuild = bundle_t::no_rep_build(aLine.no_rep_key(), aLine.no_rep_fk());
  const uint lNoRepProbe = bundle_t::no_rep_probe(aLine.no_rep_key(), aLine.no_rep_fk());
  meas_eval_t lMeasBuild(lNoRepBuild);
  meas_eval_t lMeasProbe(lNoRepProbe);

  lMeasBuild.init();
  lMeasProbe.init();

  CrystalTimer lClock;
  lClock.start();
  _scan_build->mem_alloc();
  _scan_build->mem_init();
  _scan_build->run(bundle_t::build(R,S));
  for(uint i = 1; i < lNoRepBuild; ++i) {
    _scan_build->clear();
    _scan_build->mem_free();
    lClock.stop();
    lMeasBuild.step(lClock.cycles());
    lClock.start();
    _scan_build->mem_alloc();
    _scan_build->mem_init();
    _scan_build->run(bundle_t::build(R,S));
  }
  lClock.stop();
  _meas_build = lClock.cycles();

  for(uint i = 0; i < lNoRepProbe; ++i) {
    lClock.start();
    _scan_probe->mem_alloc();
    _scan_probe->mem_init();
    _scan_probe->run(bundle_t::probe(R,S));
    _scan_probe->clear();
    _scan_probe->mem_free();
    lClock.stop();
    lMeasProbe.step(lClock.cycles());
  }

  lClock.start();
  _scan_build->clear();
  _scan_build->mem_free();
  lClock.stop();
  _meas_build += lClock.cycles();
  lMeasBuild.step(meas_build());
  
  lMeasBuild.fin();
  lMeasProbe.fin();

  if constexpr (Ttrace) {
    std::cout << "measurements build: " << lMeasBuild << std::endl;
    std::cout << "measurements probe: " << lMeasProbe << std::endl;
  }

  aLine._rt_stat_build[PLAN_ID]._cc_min = lMeasBuild.min();
  aLine._rt_stat_build[PLAN_ID]._cc_max = lMeasBuild.max();
  aLine._rt_stat_build[PLAN_ID]._cc_med = lMeasBuild.median();
  aLine._rt_stat_build[PLAN_ID]._cc_avg = lMeasBuild.avg();
  aLine._rt_stat_build[PLAN_ID]._cc_mx2 = lMeasBuild.max2();
  aLine._rt_stat_build[PLAN_ID]._cc_fst = lMeasBuild.first();

  aLine._rt_stat_probe[PLAN_ID]._cc_min = lMeasProbe.min();
  aLine._rt_stat_probe[PLAN_ID]._cc_max = lMeasProbe.max();
  aLine._rt_stat_probe[PLAN_ID]._cc_med = lMeasProbe.median();
  aLine._rt_stat_probe[PLAN_ID]._cc_avg = lMeasProbe.avg();
  aLine._rt_stat_probe[PLAN_ID]._cc_mx2 = lMeasProbe.max2();
  aLine._rt_stat_probe[PLAN_ID]._cc_fst = lMeasProbe.first();

  // XXX
  // constexpr selects plan (R join3d S),                            <== consider this correct
  // indicates how many unique probe keys form R
  // find at least one join partner in non-unique build side S.
  // However, variable name refer to size of (S lsjoin R)            <== consider naming wrong
  // TODO: fix name
  if constexpr (kPlan3dOrigRS == PLAN_ID) {
    aLine._card_sj_sr = _probe->c_mid();  // set semijoin size of R lsjoin S.
  }
}



template<typename Tr, typename Ts,
         join_order_et aJoinOrder,
         bool          NESTED,
         prefetch_et   PREFETCH,
         bool          aManyPtr, 
         bool          Ttrace,
         bool          Ttest>
void
PlanJoin<Tr,Ts,aJoinOrder,NESTED,PREFETCH,aManyPtr,Ttrace,Ttest>
::run_build_part_a(const rel_build_t& R_bld, 
                   const uint   aNoRepBuild, 
                   meas_eval_t& aMeasBuild) {
  CrystalTimer lClock;
  lClock.start();
  _scan_build->mem_alloc();
  _scan_build->mem_init();
  _scan_build->run(R_bld);
  for(uint i = 1; i < aNoRepBuild; ++i) {
    _scan_build->clear();
    _scan_build->mem_free();
    lClock.stop();
    aMeasBuild.step(lClock.cycles());
    lClock.start();
    _scan_build->mem_alloc();
    _scan_build->mem_init();
    _scan_build->run(R_bld);
  }
  lClock.stop();
  _meas_build = lClock.cycles(); // remember the first part for adding to the rest in part_b
}

template<typename Tr, typename Ts,
         join_order_et aJoinOrder,
         bool          NESTED,
         prefetch_et   PREFETCH,
         bool          aManyPtr, 
         bool          Ttrace,
         bool          Ttest>
uint64_t
PlanJoin<Tr,Ts,aJoinOrder,NESTED,PREFETCH,aManyPtr,Ttrace,Ttest>
::run_probe(const rel_probe_t& R_prb, 
            const uint         aNoRepProbe, 
                  meas_eval_t& aMeasProbe) {
  CrystalTimer lClock;
  for(uint i = 0; i < aNoRepProbe; ++i) {
    lClock.start();
    _scan_probe->mem_alloc();
    _scan_probe->mem_init();
    _scan_probe->run(R_prb);
    _scan_probe->clear();
    _scan_probe->mem_free();
    lClock.stop();
    aMeasProbe.step(lClock.cycles());
  }
  if constexpr (kPlan3dOrigRS == PLAN_ID) {
    _card_sj_sr = _probe->c_mid();
  }
  return top_count();
}

template<typename Tr, typename Ts,
         join_order_et aJoinOrder,
         bool          NESTED,
         prefetch_et   PREFETCH,
         bool          aManyPtr, 
         bool          Ttrace,
         bool          Ttest>
void
PlanJoin<Tr,Ts,aJoinOrder,NESTED,PREFETCH,aManyPtr,Ttrace,Ttest>
::run_build_part_b(meas_eval_t& aMeasBuild) {
  CrystalTimer lClock;
  lClock.start();
  _scan_build->clear();
  _scan_build->mem_free();
  lClock.stop();
  _meas_build += lClock.cycles(); // contains runtime from part_a
  aMeasBuild.step(meas_build());
}

template<typename Tr, typename Ts,
         join_order_et aJoinOrder,
         bool          NESTED,
         prefetch_et   PREFETCH,
         bool          aManyPtr,
         bool          Ttrace,
         bool          Ttest>
void
PlanJoin<Tr,Ts,aJoinOrder,NESTED,PREFETCH,aManyPtr,Ttrace,Ttest>
::clear() {
  if(nullptr != _top)   { delete _top; }
  if(nullptr != _build) { delete _build; }
  if(nullptr != _probe) { delete _probe; }
  if(nullptr != _scan_build) { delete _scan_build; }
  if(nullptr != _scan_probe) { delete _scan_probe; }
  _top   = nullptr;
  _build = nullptr;
  _probe = nullptr;
  _scan_build = nullptr;
  _scan_probe = nullptr;
}


#endif
