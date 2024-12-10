# baking-plans-icde2025

Source code for the experiments of our ICDE 2025 paper submission "Baking Plans: Know Thy Ingredients".

Content under construction.

This repository contains source code and results from our experiments
relating to our ICDE 2025 paper submission.


## Directory Structure

`cost_function_measurements`: cost function generation (cf. paper Sec. V)

- `src`:
   C++ code of join implementations (cf. paper Sec. V-A)
   and runtime experiment (cf. paper Sec. V-B-1).
   Compile with `make`, binary is in `./build/0/`
- `run_scripts`: bash scripts to run runtime experiment
- `results`: raw data from runtime measurements


`pgbench_results`: evaluation results

- `cardinality_estimation`:
   results for cardinality estimation accuracy in isolation,
   cf. paper Sec. VI-B.
- `cost_function`:
  cost function errors, table similar to Table III from paper.
- `plans`:
  plan loss for all JOB queries,
  all cardinality estimators,
  all cost functions, and
  all plan generators considered.
