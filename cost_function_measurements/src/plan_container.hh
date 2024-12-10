#ifndef HASHJOIN2_PLAN_CONTAINER_HH
#define HASHJOIN2_PLAN_CONTAINER_HH
#pragma once

#include "gminfra/bit_subsets.hh"

#include "plan.hh"

template<typename T_rel_R, typename T_rel_S, bool Ttest, bool Ttrace>
class PlanContainer {
  public:
    static constexpr bool lManyPtrNo  = false;
    static constexpr bool lManyPtrYes = true;
    static constexpr bool lCh = false;
    static constexpr bool l3d = true;
    static constexpr prefetch_et lPrefetchOff  = prefetch_et::kPrefetchNone;
    static constexpr prefetch_et lPrefetchRp = prefetch_et::kRollingPrefetching;
    static constexpr prefetch_et lPrefetchAmac = prefetch_et::kAmac;
  public:
    using R_t = T_rel_R;
    using S_t = T_rel_S;
  public:
    PlanContainer() : _planChOrigRS(),
                      _planChOrigSR(),
                      _planChImprRS(),
                      _planChImprSR(),
                      _planChUpkRpRS(),
                      _planChUpkRpSR(),
                      _planChPkdRpRS(),
                      _planChPkdRpSR(),
                      _planChAmacRS(),
                      _planChAmacSR(),
                      _planChImacRS(),
                      _planChImacSR(),
                      _plan3dOrigRS(),
                      _plan3dOrigSR(),
                      _plan3dImprRS(),
                      _plan3dImprSR(),
                      _plan3dUpkRpRS(),
                      _plan3dUpkRpSR(),
                      _plan3dPkdRpRS(),
                      _plan3dPkdRpSR(),
                      _plan3dAmacRS(),
                      _plan3dAmacSR(),
                      _plan3dImacRS(),
                      _plan3dImacSR(),
                      _card_sj_sr(0) {}
    ~PlanContainer() {}
  public:
    using plan_ch_orig_rs_t   = PlanJoin<R_t, S_t, kJoinOrderRS, lCh, lPrefetchOff,  lManyPtrNo,  Ttrace, Ttest>;
    using plan_ch_orig_sr_t   = PlanJoin<R_t, S_t, kJoinOrderSR, lCh, lPrefetchOff,  lManyPtrNo,  Ttrace, Ttest>;
    using plan_ch_impr_rs_t   = PlanJoin<R_t, S_t, kJoinOrderRS, lCh, lPrefetchOff,  lManyPtrYes, Ttrace, Ttest>;
    using plan_ch_impr_sr_t   = PlanJoin<R_t, S_t, kJoinOrderSR, lCh, lPrefetchOff,  lManyPtrYes, Ttrace, Ttest>;
    using plan_ch_upk_rp_rs_t = PlanJoin<R_t, S_t, kJoinOrderRS, lCh, lPrefetchRp,   lManyPtrNo,  Ttrace, Ttest>;
    using plan_ch_upk_rp_sr_t = PlanJoin<R_t, S_t, kJoinOrderSR, lCh, lPrefetchRp,   lManyPtrNo,  Ttrace, Ttest>;
    using plan_ch_pkd_rp_rs_t = PlanJoin<R_t, S_t, kJoinOrderRS, lCh, lPrefetchRp,   lManyPtrYes, Ttrace, Ttest>;
    using plan_ch_pkd_rp_sr_t = PlanJoin<R_t, S_t, kJoinOrderSR, lCh, lPrefetchRp,   lManyPtrYes, Ttrace, Ttest>;
    using plan_ch_amac_rs_t   = PlanJoin<R_t, S_t, kJoinOrderRS, lCh, lPrefetchAmac, lManyPtrNo,  Ttrace, Ttest>;
    using plan_ch_amac_sr_t   = PlanJoin<R_t, S_t, kJoinOrderSR, lCh, lPrefetchAmac, lManyPtrNo,  Ttrace, Ttest>;
    using plan_ch_imac_rs_t   = PlanJoin<R_t, S_t, kJoinOrderRS, lCh, lPrefetchAmac, lManyPtrYes, Ttrace, Ttest>;
    using plan_ch_imac_sr_t   = PlanJoin<R_t, S_t, kJoinOrderSR, lCh, lPrefetchAmac, lManyPtrYes, Ttrace, Ttest>;
  
    using plan_3d_orig_rs_t   = PlanJoin<R_t, S_t, kJoinOrderRS, l3d, lPrefetchOff,  lManyPtrNo,  Ttrace, Ttest>;
    using plan_3d_orig_sr_t   = PlanJoin<R_t, S_t, kJoinOrderSR, l3d, lPrefetchOff,  lManyPtrNo,  Ttrace, Ttest>;
    using plan_3d_impr_rs_t   = PlanJoin<R_t, S_t, kJoinOrderRS, l3d, lPrefetchOff,  lManyPtrYes, Ttrace, Ttest>;
    using plan_3d_impr_sr_t   = PlanJoin<R_t, S_t, kJoinOrderSR, l3d, lPrefetchOff,  lManyPtrYes, Ttrace, Ttest>;
    using plan_3d_upk_rp_rs_t = PlanJoin<R_t, S_t, kJoinOrderRS, l3d, lPrefetchRp,   lManyPtrNo,  Ttrace, Ttest>;
    using plan_3d_upk_rp_sr_t = PlanJoin<R_t, S_t, kJoinOrderSR, l3d, lPrefetchRp,   lManyPtrNo,  Ttrace, Ttest>;
    using plan_3d_pkd_rp_rs_t = PlanJoin<R_t, S_t, kJoinOrderRS, l3d, lPrefetchRp,   lManyPtrYes, Ttrace, Ttest>;
    using plan_3d_pkd_rp_sr_t = PlanJoin<R_t, S_t, kJoinOrderSR, l3d, lPrefetchRp,   lManyPtrYes, Ttrace, Ttest>;
    using plan_3d_amac_rs_t   = PlanJoin<R_t, S_t, kJoinOrderRS, l3d, lPrefetchAmac, lManyPtrNo,  Ttrace, Ttest>;
    using plan_3d_amac_sr_t   = PlanJoin<R_t, S_t, kJoinOrderSR, l3d, lPrefetchAmac, lManyPtrNo,  Ttrace, Ttest>;
    using plan_3d_imac_rs_t   = PlanJoin<R_t, S_t, kJoinOrderRS, l3d, lPrefetchAmac, lManyPtrYes, Ttrace, Ttest>;
    using plan_3d_imac_sr_t   = PlanJoin<R_t, S_t, kJoinOrderSR, l3d, lPrefetchAmac, lManyPtrYes, Ttrace, Ttest>;
  public:
    static_assert(kPlanChOrigRS  == plan_ch_orig_rs_t::PLAN_ID);
    static_assert(kPlanChOrigSR  == plan_ch_orig_sr_t::PLAN_ID);
    static_assert(kPlanChImprRS  == plan_ch_impr_rs_t::PLAN_ID);
    static_assert(kPlanChImprSR  == plan_ch_impr_sr_t::PLAN_ID);
    static_assert(kPlanChRpUpkRS == plan_ch_upk_rp_rs_t::PLAN_ID); 
    static_assert(kPlanChRpUpkSR == plan_ch_upk_rp_sr_t::PLAN_ID); 
    static_assert(kPlanChRpPkdRS == plan_ch_pkd_rp_rs_t::PLAN_ID); 
    static_assert(kPlanChRpPkdSR == plan_ch_pkd_rp_sr_t::PLAN_ID); 
    static_assert(kPlanChAmacRS  == plan_ch_amac_rs_t::PLAN_ID);
    static_assert(kPlanChAmacSR  == plan_ch_amac_sr_t::PLAN_ID);
    static_assert(kPlanChImacRS  == plan_ch_imac_rs_t::PLAN_ID);
    static_assert(kPlanChImacSR  == plan_ch_imac_sr_t::PLAN_ID);

    static_assert(kPlan3dOrigRS  == plan_3d_orig_rs_t::PLAN_ID);
    static_assert(kPlan3dOrigSR  == plan_3d_orig_sr_t::PLAN_ID);
    static_assert(kPlan3dImprRS  == plan_3d_impr_rs_t::PLAN_ID);
    static_assert(kPlan3dImprSR  == plan_3d_impr_sr_t::PLAN_ID);
    static_assert(kPlan3dRpUpkRS == plan_3d_upk_rp_rs_t::PLAN_ID);
    static_assert(kPlan3dRpUpkSR == plan_3d_upk_rp_sr_t::PLAN_ID);
    static_assert(kPlan3dRpPkdRS == plan_3d_pkd_rp_rs_t::PLAN_ID);
    static_assert(kPlan3dRpPkdSR == plan_3d_pkd_rp_sr_t::PLAN_ID);
    static_assert(kPlan3dAmacRS  == plan_3d_amac_rs_t::PLAN_ID);
    static_assert(kPlan3dAmacSR  == plan_3d_amac_sr_t::PLAN_ID);
    static_assert(kPlan3dImacRS  == plan_3d_imac_rs_t::PLAN_ID);
    static_assert(kPlan3dImacSR  == plan_3d_imac_sr_t::PLAN_ID);
  public:
    void build_plan(const plan_et aPlanId, 
                    const size_t  aHashDirSize, 
                    const uint    aLog2ChunkSize);
    void build_plan_set(const planset_t aPlanSet,
                        const size_t    aHashDirSize,
                        const uint      aLog2ChunkSize);
  public:
    template<typename Trel>
    void run_build_part_a(const plan_et aPlanId,
                          const Trel&   aRel,
                          const uint    aNoRep,
                          meas_eval_t&  aMeasBuild);
    // returns result size
    template<typename Trel>
    uint64_t run_probe(const plan_et aPlanId,
                   const Trel&   aRel,
                   const uint    aNoRep,
                   meas_eval_t&  aMeasProbe);

    void run_build_part_b(const plan_et aPlanId,
                          meas_eval_t&  aMeasBuild);
  public:
    inline uint64_t c_mid()      const { return _card_sj_sr; }
    inline uint64_t card_sj_sr() const { return _card_sj_sr; }
           size_t   mem_consumption_ht(const plan_et aPlanId) const;
           size_t   ht_num_nodes(const plan_et aPlanId) const;  // ch + 3d
           size_t   ht_num_nodes_sub(const plan_et aPlanId) const;  // 3d only, ch = 0
           double   ht_dir_fill_factor(const plan_et aPlanId) const;
           uint64_t count_build(const plan_et aPlanId) const;  // _count from build operator
           uint64_t count_probe(const plan_et aPlanId) const;  //             probe
  public:
    plan_ch_orig_rs_t   _planChOrigRS;
    plan_ch_orig_sr_t   _planChOrigSR;
    plan_ch_impr_rs_t   _planChImprRS;
    plan_ch_impr_sr_t   _planChImprSR;
    plan_ch_upk_rp_rs_t _planChUpkRpRS;
    plan_ch_upk_rp_sr_t _planChUpkRpSR;
    plan_ch_pkd_rp_rs_t _planChPkdRpRS;
    plan_ch_pkd_rp_sr_t _planChPkdRpSR;
    plan_ch_amac_rs_t   _planChAmacRS;
    plan_ch_amac_sr_t   _planChAmacSR;
    plan_ch_imac_rs_t   _planChImacRS;
    plan_ch_imac_sr_t   _planChImacSR;

    plan_3d_orig_rs_t  _plan3dOrigRS;
    plan_3d_orig_sr_t  _plan3dOrigSR;
    plan_3d_impr_rs_t  _plan3dImprRS;
    plan_3d_impr_sr_t  _plan3dImprSR;
    plan_3d_upk_rp_rs_t _plan3dUpkRpRS;
    plan_3d_upk_rp_sr_t _plan3dUpkRpSR;
    plan_3d_pkd_rp_rs_t _plan3dPkdRpRS;
    plan_3d_pkd_rp_sr_t _plan3dPkdRpSR;
    plan_3d_amac_rs_t  _plan3dAmacRS;
    plan_3d_amac_sr_t  _plan3dAmacSR;
    plan_3d_imac_rs_t  _plan3dImacRS;
    plan_3d_imac_sr_t  _plan3dImacSR;

    uint64_t _card_sj_sr; // set by kPlan3dOrigRS
};

template<typename T_rel_R, typename T_rel_S, bool Ttest, bool Ttrace>
void
PlanContainer<T_rel_R,T_rel_S,Ttest,Ttrace>::build_plan_set(const planset_t aPlanSet,
                                                            const size_t    aHashDirSize,
                                                            const uint      aLog2ChunkSize) {
  for(BvMemberIdxDesc<planset_t> lIter(aPlanSet); lIter.isValid(); ++lIter) {
    build_plan((plan_et) (*lIter), aHashDirSize, aLog2ChunkSize);
  }
}

template<typename T_rel_R, typename T_rel_S, bool Ttest, bool Ttrace>
void
PlanContainer<T_rel_R,T_rel_S,Ttest,Ttrace>::build_plan(const plan_et aPlanId, 
                                                        const size_t  aHashDirSize, 
                                                        const uint    aLog2ChunkSize) {
  switch(aPlanId) {
    case kPlanChOrigRS:  _planChOrigRS.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlanChOrigSR:  _planChOrigSR.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlanChImprRS:  _planChImprRS.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlanChImprSR:  _planChImprSR.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlanChRpUpkRS: _planChUpkRpRS.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlanChRpUpkSR: _planChUpkRpSR.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlanChRpPkdRS: _planChPkdRpRS.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlanChRpPkdSR: _planChPkdRpSR.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlanChAmacRS:  _planChAmacRS.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlanChAmacSR:  _planChAmacSR.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlanChImacRS:  _planChImacRS.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlanChImacSR:  _planChImacSR.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlan3dOrigRS:  _plan3dOrigRS.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlan3dOrigSR:  _plan3dOrigSR.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlan3dImprRS:  _plan3dImprRS.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlan3dImprSR:  _plan3dImprSR.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlan3dRpUpkRS: _plan3dUpkRpRS.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlan3dRpUpkSR: _plan3dUpkRpSR.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlan3dRpPkdRS: _plan3dPkdRpRS.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlan3dRpPkdSR: _plan3dPkdRpSR.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlan3dAmacRS:  _plan3dAmacRS.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlan3dAmacSR:  _plan3dAmacSR.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlan3dImacRS:  _plan3dImacRS.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    case kPlan3dImacSR:  _plan3dImacSR.build_plan(aHashDirSize, aLog2ChunkSize);
                         break;
    default: assert(0 == 1);
  }
}

template<typename T_rel_R, typename T_rel_S, bool Ttest, bool Ttrace>
template<typename Trel>
void
PlanContainer<T_rel_R,T_rel_S,Ttest,Ttrace>
::run_build_part_a(const plan_et aPlanId,
                   const Trel&   aRel,
                   const uint    aNoRep,
                   meas_eval_t&  aMeasBuild) {
  switch(aPlanId) {
    case kPlanChOrigRS:  _planChOrigRS.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlanChOrigSR:  _planChOrigSR.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlanChImprRS:  _planChImprRS.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlanChImprSR:  _planChImprSR.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlanChRpUpkRS: _planChUpkRpRS.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlanChRpUpkSR: _planChUpkRpSR.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlanChRpPkdRS: _planChPkdRpRS.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlanChRpPkdSR: _planChPkdRpSR.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlanChAmacRS:  _planChAmacRS.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlanChAmacSR:  _planChAmacSR.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlanChImacRS:  _planChImacRS.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlanChImacSR:  _planChImacSR.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlan3dOrigRS:  _plan3dOrigRS.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlan3dOrigSR:  _plan3dOrigSR.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlan3dImprRS:  _plan3dImprRS.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlan3dImprSR:  _plan3dImprSR.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlan3dRpUpkRS: _plan3dUpkRpRS.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlan3dRpUpkSR: _plan3dUpkRpSR.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlan3dRpPkdRS: _plan3dPkdRpRS.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlan3dRpPkdSR: _plan3dPkdRpSR.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlan3dAmacRS:  _plan3dAmacRS.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlan3dAmacSR:  _plan3dAmacSR.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlan3dImacRS:  _plan3dImacRS.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    case kPlan3dImacSR:  _plan3dImacSR.run_build_part_a(aRel, aNoRep, aMeasBuild);
                         break;
    default: assert(0 == 1);
  }
}

template<typename T_rel_R, typename T_rel_S, bool Ttest, bool Ttrace>
template<typename Trel>
uint64_t
PlanContainer<T_rel_R,T_rel_S,Ttest,Ttrace>
::run_probe(const plan_et aPlanId,
            const Trel&   aRel,
            const uint    aNoRep,
            meas_eval_t&  aMeasProbe) {
  uint64_t lRes = 0;
  switch(aPlanId) {
    case kPlanChOrigRS:  lRes = _planChOrigRS.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlanChOrigSR:  lRes = _planChOrigSR.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlanChImprRS:  lRes = _planChImprRS.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlanChImprSR:  lRes = _planChImprSR.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlanChRpUpkRS: lRes = _planChUpkRpRS.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlanChRpUpkSR: lRes = _planChUpkRpSR.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlanChRpPkdRS: lRes = _planChPkdRpRS.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlanChRpPkdSR: lRes = _planChPkdRpSR.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlanChAmacRS:  lRes = _planChAmacRS.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlanChAmacSR:  lRes = _planChAmacSR.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlanChImacRS:  lRes = _planChImacRS.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlanChImacSR:  lRes = _planChImacSR.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlan3dOrigRS:  lRes = _plan3dOrigRS.run_probe(aRel, aNoRep, aMeasProbe);
                         _card_sj_sr = _plan3dOrigRS.card_sj_sr();
                         break;
    case kPlan3dOrigSR:  lRes = _plan3dOrigSR.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlan3dImprRS:  lRes = _plan3dImprRS.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlan3dImprSR:  lRes = _plan3dImprSR.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlan3dRpUpkRS: lRes = _plan3dUpkRpRS.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlan3dRpUpkSR: lRes = _plan3dUpkRpSR.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlan3dRpPkdRS: lRes = _plan3dPkdRpRS.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlan3dRpPkdSR: lRes = _plan3dPkdRpSR.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlan3dAmacRS:  lRes = _plan3dAmacRS.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlan3dAmacSR:  lRes = _plan3dAmacSR.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    case kPlan3dImacRS:  lRes = _plan3dImacRS.run_probe(aRel, aNoRep, aMeasProbe);
                         break;    
    case kPlan3dImacSR:  lRes = _plan3dImacSR.run_probe(aRel, aNoRep, aMeasProbe);
                         break;
    default: assert(0 == 1);
  }
  return lRes;
}

template<typename T_rel_R, typename T_rel_S, bool Ttest, bool Ttrace>
void 
PlanContainer<T_rel_R,T_rel_S,Ttest,Ttrace>
::run_build_part_b(const plan_et aPlanId,
                   meas_eval_t&  aMeasBuild) {
  switch(aPlanId) {
    case kPlanChOrigRS:  _planChOrigRS.run_build_part_b(aMeasBuild);
                         break;
    case kPlanChOrigSR:  _planChOrigSR.run_build_part_b(aMeasBuild);
                         break;
    case kPlanChImprRS:  _planChImprRS.run_build_part_b(aMeasBuild);
                         break;
    case kPlanChImprSR:  _planChImprSR.run_build_part_b(aMeasBuild);
                         break;
    case kPlanChRpUpkRS: _planChUpkRpRS.run_build_part_b(aMeasBuild);
                         break;
    case kPlanChRpUpkSR: _planChUpkRpSR.run_build_part_b(aMeasBuild);
                         break;
    case kPlanChRpPkdRS: _planChPkdRpRS.run_build_part_b(aMeasBuild);
                         break;
    case kPlanChRpPkdSR: _planChPkdRpSR.run_build_part_b(aMeasBuild);
                         break;
    case kPlanChAmacRS:  _planChAmacRS.run_build_part_b(aMeasBuild);
                         break;
    case kPlanChAmacSR:  _planChAmacSR.run_build_part_b(aMeasBuild);
                         break;
    case kPlanChImacRS:  _planChImacRS.run_build_part_b(aMeasBuild);
                         break;
    case kPlanChImacSR:  _planChImacSR.run_build_part_b(aMeasBuild);
                         break;
    case kPlan3dOrigRS:  _plan3dOrigRS.run_build_part_b(aMeasBuild);
                         break;
    case kPlan3dOrigSR:  _plan3dOrigSR.run_build_part_b(aMeasBuild);
                         break;
    case kPlan3dImprRS:  _plan3dImprRS.run_build_part_b(aMeasBuild);
                         break;
    case kPlan3dImprSR:  _plan3dImprSR.run_build_part_b(aMeasBuild);
                         break;
    case kPlan3dRpUpkRS: _plan3dUpkRpRS.run_build_part_b(aMeasBuild);
                         break;
    case kPlan3dRpUpkSR: _plan3dUpkRpSR.run_build_part_b(aMeasBuild);
                         break;
    case kPlan3dRpPkdRS: _plan3dPkdRpRS.run_build_part_b(aMeasBuild);
                         break;
    case kPlan3dRpPkdSR: _plan3dPkdRpSR.run_build_part_b(aMeasBuild);
                         break;
    case kPlan3dAmacRS:  _plan3dAmacRS.run_build_part_b(aMeasBuild);
                         break;
    case kPlan3dAmacSR:  _plan3dAmacSR.run_build_part_b(aMeasBuild);
                         break;
    case kPlan3dImacRS:  _plan3dImacRS.run_build_part_b(aMeasBuild);
                         break;
    case kPlan3dImacSR:  _plan3dImacSR.run_build_part_b(aMeasBuild);
                         break;
    default: assert(0 == 1);
  }
}

template<typename T_rel_R, typename T_rel_S, bool Ttest, bool Ttrace>
size_t
PlanContainer<T_rel_R,T_rel_S,Ttest,Ttrace>
::mem_consumption_ht(const plan_et aPlanId) const {
  switch(aPlanId) {
    case kPlanChOrigRS:  return _planChOrigRS.mem_consumption_ht();
    case kPlanChOrigSR:  return _planChOrigSR.mem_consumption_ht();
    case kPlanChImprRS:  return _planChImprRS.mem_consumption_ht();
    case kPlanChImprSR:  return _planChImprSR.mem_consumption_ht();
    case kPlanChRpUpkRS: return _planChUpkRpRS.mem_consumption_ht();
    case kPlanChRpUpkSR: return _planChUpkRpSR.mem_consumption_ht();
    case kPlanChRpPkdRS: return _planChPkdRpRS.mem_consumption_ht();
    case kPlanChRpPkdSR: return _planChPkdRpSR.mem_consumption_ht();
    case kPlanChAmacRS:  return _planChAmacRS.mem_consumption_ht();
    case kPlanChAmacSR:  return _planChAmacSR.mem_consumption_ht();
    case kPlanChImacRS:  return _planChImacRS.mem_consumption_ht();
    case kPlanChImacSR:  return _planChImacSR.mem_consumption_ht();
    case kPlan3dOrigRS:  return _plan3dOrigRS.mem_consumption_ht();
    case kPlan3dOrigSR:  return _plan3dOrigSR.mem_consumption_ht();
    case kPlan3dImprRS:  return _plan3dImprRS.mem_consumption_ht();
    case kPlan3dImprSR:  return _plan3dImprSR.mem_consumption_ht();
    case kPlan3dRpUpkRS: return _plan3dUpkRpRS.mem_consumption_ht();
    case kPlan3dRpUpkSR: return _plan3dUpkRpSR.mem_consumption_ht();
    case kPlan3dRpPkdRS: return _plan3dPkdRpRS.mem_consumption_ht();
    case kPlan3dRpPkdSR: return _plan3dPkdRpSR.mem_consumption_ht();
    case kPlan3dAmacRS:  return _plan3dAmacRS.mem_consumption_ht();
    case kPlan3dAmacSR:  return _plan3dAmacSR.mem_consumption_ht();
    case kPlan3dImacRS:  return _plan3dImacRS.mem_consumption_ht();
    case kPlan3dImacSR:  return _plan3dImacSR.mem_consumption_ht();
    default: assert(0 == 1); return -7;
  }
}

template<typename T_rel_R, typename T_rel_S, bool Ttest, bool Ttrace>
size_t
PlanContainer<T_rel_R,T_rel_S,Ttest,Ttrace>
::ht_num_nodes(const plan_et aPlanId) const {
  switch(aPlanId) {
    case kPlanChOrigRS:  return _planChOrigRS.build()->ht_num_nodes();
    case kPlanChOrigSR:  return _planChOrigSR.build()->ht_num_nodes();
    case kPlanChImprRS:  return _planChImprRS.build()->ht_num_nodes();
    case kPlanChImprSR:  return _planChImprSR.build()->ht_num_nodes();
    case kPlanChRpUpkRS: return _planChUpkRpRS.build()->ht_num_nodes();
    case kPlanChRpUpkSR: return _planChUpkRpSR.build()->ht_num_nodes();
    case kPlanChRpPkdRS: return _planChPkdRpRS.build()->ht_num_nodes();
    case kPlanChRpPkdSR: return _planChPkdRpSR.build()->ht_num_nodes();
    case kPlanChAmacRS:  return _planChAmacRS.build()->ht_num_nodes();
    case kPlanChAmacSR:  return _planChAmacSR.build()->ht_num_nodes();
    case kPlanChImacRS:  return _planChImacRS.build()->ht_num_nodes();
    case kPlanChImacSR:  return _planChImacSR.build()->ht_num_nodes();
    //
    case kPlan3dOrigRS:  return _plan3dOrigRS.build()->ht_num_main_nodes();
    case kPlan3dOrigSR:  return _plan3dOrigSR.build()->ht_num_main_nodes();
    case kPlan3dImprRS:  return _plan3dImprRS.build()->ht_num_main_nodes();
    case kPlan3dImprSR:  return _plan3dImprSR.build()->ht_num_main_nodes();
    case kPlan3dRpUpkRS: return _plan3dUpkRpRS.build()->ht_num_main_nodes();
    case kPlan3dRpUpkSR: return _plan3dUpkRpSR.build()->ht_num_main_nodes();
    case kPlan3dRpPkdRS: return _plan3dPkdRpRS.build()->ht_num_main_nodes();
    case kPlan3dRpPkdSR: return _plan3dPkdRpSR.build()->ht_num_main_nodes();
    case kPlan3dAmacRS:  return _plan3dAmacRS.build()->ht_num_main_nodes();
    case kPlan3dAmacSR:  return _plan3dAmacSR.build()->ht_num_main_nodes();
    case kPlan3dImacRS:  return _plan3dImacRS.build()->ht_num_main_nodes();
    case kPlan3dImacSR:  return _plan3dImacSR.build()->ht_num_main_nodes();
    default: assert(0 == 1); return -7;
  }
}

template<typename T_rel_R, typename T_rel_S, bool Ttest, bool Ttrace>
size_t
PlanContainer<T_rel_R,T_rel_S,Ttest,Ttrace>
::ht_num_nodes_sub(const plan_et aPlanId) const {
  switch(aPlanId) {
    case kPlanChOrigRS:  return 0;
    case kPlanChOrigSR:  return 0;
    case kPlanChImprRS:  return 0;
    case kPlanChImprSR:  return 0;
    case kPlanChRpUpkRS: return 0;
    case kPlanChRpUpkSR: return 0;
    case kPlanChRpPkdRS: return 0;
    case kPlanChRpPkdSR: return 0;
    case kPlanChAmacRS:  return 0;
    case kPlanChAmacSR:  return 0;
    case kPlanChImacRS:  return 0;
    case kPlanChImacSR:  return 0;
    //
    case kPlan3dOrigRS:  return _plan3dOrigRS.build()->ht_num_sub_nodes();
    case kPlan3dOrigSR:  return _plan3dOrigSR.build()->ht_num_sub_nodes();
    case kPlan3dImprRS:  return _plan3dImprRS.build()->ht_num_sub_nodes();
    case kPlan3dImprSR:  return _plan3dImprSR.build()->ht_num_sub_nodes();
    case kPlan3dRpUpkRS: return _plan3dUpkRpRS.build()->ht_num_sub_nodes();
    case kPlan3dRpUpkSR: return _plan3dUpkRpSR.build()->ht_num_sub_nodes();
    case kPlan3dRpPkdRS: return _plan3dPkdRpRS.build()->ht_num_sub_nodes();
    case kPlan3dRpPkdSR: return _plan3dPkdRpSR.build()->ht_num_sub_nodes();
    case kPlan3dAmacRS:  return _plan3dAmacRS.build()->ht_num_sub_nodes();
    case kPlan3dAmacSR:  return _plan3dAmacSR.build()->ht_num_sub_nodes();
    case kPlan3dImacRS:  return _plan3dImacRS.build()->ht_num_sub_nodes();
    case kPlan3dImacSR:  return _plan3dImacSR.build()->ht_num_sub_nodes();
    default: assert(0 == 1); return -7;
  }
}

template<typename T_rel_R, typename T_rel_S, bool Ttest, bool Ttrace>
double
PlanContainer<T_rel_R,T_rel_S,Ttest,Ttrace>
::ht_dir_fill_factor(const plan_et aPlanId) const {
  switch(aPlanId) {
    case kPlanChOrigRS:  return _planChOrigRS.build()->hashtable().dir_fill_factor();
    case kPlanChOrigSR:  return _planChOrigSR.build()->hashtable().dir_fill_factor();
    case kPlanChImprRS:  return _planChImprRS.build()->hashtable().dir_fill_factor();
    case kPlanChImprSR:  return _planChImprSR.build()->hashtable().dir_fill_factor();
    case kPlanChRpUpkRS: return _planChUpkRpRS.build()->hashtable().dir_fill_factor();
    case kPlanChRpUpkSR: return _planChUpkRpSR.build()->hashtable().dir_fill_factor();
    case kPlanChRpPkdRS: return _planChPkdRpRS.build()->hashtable().dir_fill_factor();
    case kPlanChRpPkdSR: return _planChPkdRpSR.build()->hashtable().dir_fill_factor();
    case kPlanChAmacRS:  return _planChAmacRS.build()->hashtable().dir_fill_factor();
    case kPlanChAmacSR:  return _planChAmacSR.build()->hashtable().dir_fill_factor();
    case kPlanChImacRS:  return _planChImacRS.build()->hashtable().dir_fill_factor();
    case kPlanChImacSR:  return _planChImacSR.build()->hashtable().dir_fill_factor();
    //
    case kPlan3dOrigRS:  return _plan3dOrigRS.build()->hashtable().dir_fill_factor();
    case kPlan3dOrigSR:  return _plan3dOrigSR.build()->hashtable().dir_fill_factor();
    case kPlan3dImprRS:  return _plan3dImprRS.build()->hashtable().dir_fill_factor();
    case kPlan3dImprSR:  return _plan3dImprSR.build()->hashtable().dir_fill_factor();
    case kPlan3dRpUpkRS: return _plan3dUpkRpRS.build()->hashtable().dir_fill_factor();
    case kPlan3dRpUpkSR: return _plan3dUpkRpSR.build()->hashtable().dir_fill_factor();
    case kPlan3dRpPkdRS: return _plan3dPkdRpRS.build()->hashtable().dir_fill_factor();
    case kPlan3dRpPkdSR: return _plan3dPkdRpSR.build()->hashtable().dir_fill_factor();
    case kPlan3dAmacRS:  return _plan3dAmacRS.build()->hashtable().dir_fill_factor();
    case kPlan3dAmacSR:  return _plan3dAmacSR.build()->hashtable().dir_fill_factor();
    case kPlan3dImacRS:  return _plan3dImacRS.build()->hashtable().dir_fill_factor();
    case kPlan3dImacSR:  return _plan3dImacSR.build()->hashtable().dir_fill_factor();
    default: assert(0 == 1); return -7;
  }
}


template<typename T_rel_R, typename T_rel_S, bool Ttest, bool Ttrace>
uint64_t
PlanContainer<T_rel_R,T_rel_S,Ttest,Ttrace>
::count_build(const plan_et aPlanId) const {
  switch(aPlanId) {
    case kPlanChOrigRS:  return _planChOrigRS.build()->count();
    case kPlanChOrigSR:  return _planChOrigSR.build()->count();
    case kPlanChImprRS:  return _planChImprRS.build()->count();
    case kPlanChImprSR:  return _planChImprSR.build()->count();
    case kPlanChRpUpkRS: return _planChUpkRpRS.build()->count();
    case kPlanChRpUpkSR: return _planChUpkRpSR.build()->count();
    case kPlanChRpPkdRS: return _planChPkdRpRS.build()->count();
    case kPlanChRpPkdSR: return _planChPkdRpSR.build()->count();
    case kPlanChAmacRS:  return _planChAmacRS.build()->count();
    case kPlanChAmacSR:  return _planChAmacSR.build()->count();
    case kPlanChImacRS:  return _planChImacRS.build()->count();
    case kPlanChImacSR:  return _planChImacSR.build()->count();
    //
    case kPlan3dOrigRS:  return _plan3dOrigRS.build()->count();
    case kPlan3dOrigSR:  return _plan3dOrigSR.build()->count();
    case kPlan3dImprRS:  return _plan3dImprRS.build()->count();
    case kPlan3dImprSR:  return _plan3dImprSR.build()->count();
    case kPlan3dRpUpkRS: return _plan3dUpkRpRS.build()->count();
    case kPlan3dRpUpkSR: return _plan3dUpkRpSR.build()->count();
    case kPlan3dRpPkdRS: return _plan3dPkdRpRS.build()->count();
    case kPlan3dRpPkdSR: return _plan3dPkdRpSR.build()->count();
    case kPlan3dAmacRS:  return _plan3dAmacRS.build()->count();
    case kPlan3dAmacSR:  return _plan3dAmacSR.build()->count();
    case kPlan3dImacRS:  return _plan3dImacRS.build()->count();
    case kPlan3dImacSR:  return _plan3dImacSR.build()->count();
    default: assert(0 == 1); return -7;
  }
}

template<typename T_rel_R, typename T_rel_S, bool Ttest, bool Ttrace>
uint64_t
PlanContainer<T_rel_R,T_rel_S,Ttest,Ttrace>
::count_probe(const plan_et aPlanId) const {
  switch(aPlanId) {
    case kPlanChOrigRS:  return _planChOrigRS.probe()->count();
    case kPlanChOrigSR:  return _planChOrigSR.probe()->count();
    case kPlanChImprRS:  return _planChImprRS.probe()->count();
    case kPlanChImprSR:  return _planChImprSR.probe()->count();
    case kPlanChRpUpkRS: return _planChUpkRpRS.probe()->count();
    case kPlanChRpUpkSR: return _planChUpkRpSR.probe()->count();
    case kPlanChRpPkdRS: return _planChPkdRpRS.probe()->count();
    case kPlanChRpPkdSR: return _planChPkdRpSR.probe()->count();
    case kPlanChAmacRS:  return _planChAmacRS.probe()->count();
    case kPlanChAmacSR:  return _planChAmacSR.probe()->count();
    case kPlanChImacRS:  return _planChImacRS.probe()->count();
    case kPlanChImacSR:  return _planChImacSR.probe()->count();
    //
    case kPlan3dOrigRS:  return _plan3dOrigRS.probe()->count();
    case kPlan3dOrigSR:  return _plan3dOrigSR.probe()->count();
    case kPlan3dImprRS:  return _plan3dImprRS.probe()->count();
    case kPlan3dImprSR:  return _plan3dImprSR.probe()->count();
    case kPlan3dRpUpkRS: return _plan3dUpkRpRS.probe()->count();
    case kPlan3dRpUpkSR: return _plan3dUpkRpSR.probe()->count();
    case kPlan3dRpPkdRS: return _plan3dPkdRpRS.probe()->count();
    case kPlan3dRpPkdSR: return _plan3dPkdRpSR.probe()->count();
    case kPlan3dAmacRS:  return _plan3dAmacRS.probe()->count();
    case kPlan3dAmacSR:  return _plan3dAmacSR.probe()->count();
    case kPlan3dImacRS:  return _plan3dImacRS.probe()->count();
    case kPlan3dImacSR:  return _plan3dImacSR.probe()->count();
    default: assert(0 == 1); return -7;
  }
}

#endif
