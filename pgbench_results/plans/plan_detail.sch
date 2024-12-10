# Schema of plan_detail.rel

table plan_detail {
  string query,                # query graph ID = JOB query
  char   plan_generator,       # plan generator, see below
  int    card_provider,        # cardinality provider ID, see below
  int    cost_function,        # cost function, see below
  double loss_factor,          # plan loss factor
  string plan_generator_name,  # plan generator name, see below
  string card_provider_name,   # cardinality provider name, see below
  uint   no_rel                # number of relations in query
};


## Plan Generators
upper case: run with true cardinalities
lower case: run under independence assumption or with alt cardinality provider
- A/a: DPccp
- B/b: DPccpCout, no commute
- C/c: DPccpCout, applies commutativity during plan extraction
- D/d: TDmcc
- E/e: MVP
- G/g: GOO cost
- H/h: GOO card
- I/i: IKKBZ, original, pure left-deep trees
- K/k: IKKBZ, applies commutativity during plan extraction
- L/l: DPlin, no commute
- M/m: DPlin, applies commutativity
- R/r: SimpliSquared (Datta et al.)
- S/s: GreedyStar000, no commute
- T/t: GreedyStar000, applies commutativity
- U/u: DPccp finding wost plan for BuildPlan only
- V/v: DPccp finding wost plan for CCPs only
- W/w: DPccp finding worst plan


## Cardinality Providers
- 0 = (true): CE_{tru}, cf. paper Sec. IV.
- 1 = ia(M):  CE_{IA-M}, cf. paper Sec. IV-A.
- 2 = ia(S):  CE_{IA-S}, cf. paper Sec. IV-A.
- 3 = ia(L):  CE_{IA-L}, cf. paper Sec. IV-A.
- 5 = CrBse:  CE_{base}, cf. paper Sec. IV-B-1.
- 6 = CrSel:  CE_{sel}, cf. paper Sec. IV-B-1.
- 7 = CJBse:  CE_{j-base}, cf. paper Sec. IV-B-2.
- 8 = CJSel:  CE_{j-sel}, cf. paper Sec. IV-B-2.


## Cost Function
- 1 = BP_{trad} with CF_{tru}
- 2 = BP_{trad} with CF_{est}
