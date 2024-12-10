# Cardinality Estimation

## Files

- `cp_*_error.dat`: One file for each cardinality estimator.
   One line represents one cardinality estimate requested by any plan generator
   during plan generation running all JOB queries.
   Duplicate-free, i.e., each plan class is included only once.
   File contains two columns:
   the plan class size, and the p-error of the respective cardinality estimate.

- `plot_*.pdf`: one boxplot for each cardinality estimator

- `mk_boxplots.R`: R script to make box plots from measurements.


## Cardinality estimator:

- `ia_m`:        CE_{IA-M}
- `ia_l`:        CE_{IA-L}
- `ia_s`:        CE_{IA-S}
- `crude_base`:  CE_{base}
- `crude_sel`:   CE_{sel}
- `join_base`:   CE_{j-base}
- `join_sel`:    CE_{j-sel}
