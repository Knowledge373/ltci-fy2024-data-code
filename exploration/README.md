# Fixed exploratory calculations

Run from the repository root:

    python exploration/run_all.py --out /path/to/new/exploration-output

The output directory must not exist. The wrapper runs the existing four-weight,
distribution/influence, institution-size, and fixed top-assessment-exclusion
calculations. Inputs and model specification are fixed; no new model selection,
inference, or threshold tuning is performed.

Install the versions in `requirements.txt`. Calculations are single-threaded.
The input data are insurer-level public administrative aggregates for FY2024.
The 61-coefficient baseline fixture is required by the weight-comparison output
and baseline checks. It is a fixed computational input, not an output directory.

Each analysis writes newly computed tables and figures under the requested output
path. Results are descriptive or point estimates where specified, not causal
estimates. Empty institutional cells remain empty. The proxy weight N*A/T is a
constructed exposure proxy, not an observed ordinary-insured population count.

## Reading coefficient and descriptive outputs

See [生成出力の読み方](../OUTPUT_DICTIONARY_JA.md). In particular, `estimate_pp` is an outcome-scale coefficient per one stored unit of a continuous predictor, not per one percentage-point increase in a share predictor. Monetary outputs use yen and insurer-size outputs use end-of-fiscal-year Type 1 insured counts.
