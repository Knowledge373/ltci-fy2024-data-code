# Fixed-support sensitivity analyses

This execution-only module reconstructs six fixed-support design matrices from the supplied 1,573-row FY2024 input and support memberships. It runs six fractional-logit GLMs, two additional half-SD OLS point fits, and, by default, four sets of 400 freshly computed TOST boundaries (1,600 boundaries / 800 widths). It does not rerun the original raw-source collection or support-selection procedure.

## Run

Use Python 3.12 with `numpy==2.3.5`, and an R executable or launcher providing R/stats 4.6.1, fwildclusterboot 0.14.3 and dqrng 0.4.1. `sha256sum` must be available. The version checks are strict, including in points-only mode.

```
python -m pip install -r sensitivity/requirements.txt
python sensitivity/run_public.py --out results/sensitivity --r-command /path/to/run-r-local.sh
```

`--out` must be a new directory outside this module. Every model and boundary is computed anew. Results, runtime details, logs and completion status are written only there. No historical results are loaded or copied. `COMPLETE_FULL_FIXED_SUPPORT` means all requested computations finished; it does not imply that every TOST width passed or that historical numerical equality was tested. Consult `THRESHOLDS.csv` for PASS/NONPASS/HOLD.

For a partial installation smoke test, add `--points-only`. This runs six GLMs and two half-SD OLS point fits, with independent fixed-coefficient arithmetic checks, but no TOST/bootstrap. A successful points-only run is not evidence of a full TOST run.

## Optional pinned runtime setup

`environment/setup_runtime.py` downloads SHA256-locked official artifacts and builds a directory-local runtime on Debian 13 amd64. It needs an existing GCC/G++/make toolchain, `dpkg-deb`, `sha256sum`, compatible system shared libraries, internet access and sufficient memory/disk. It does not install OS packages. Other platforms must provide the exact pinned R/package versions themselves; no fallback version is selected. The download specification is supplied; runtime binaries and package archives are not included.

```
python sensitivity/environment/setup_runtime.py --out runtime/sensitivity
python sensitivity/run_public.py --out results/sensitivity --r-command runtime/sensitivity/run-r-local.sh
```

## Fixed computation and provenance

The supplied input bytes and each support's X/y/cluster/CSV fingerprints are checked before fitting. The R computational code and diagnostic formulas are preserved from the recovered analysis scripts; only the boundary runners' dependency paths were made relative. The Python entry point removes historical comparisons and derives current width statuses from fresh boundary outputs. The independent reconstruction retains the original rank, logistic and linear-predictor calculations.

The original settings remain: B=99,999; R and dqrng seeds 20261015 reset per boundary; null-imposed WCR11/fnw11; Rademacher weights; R engine; dqrng sampling; one thread; alpha=1/60; grid 0.01–2.00 percentage points; original small-sample corrections; no confidence-interval inversion. P1/P3 enumerate 32,768 signs. Diagnostic HOLD retains priority over the two official one-sided p-value decisions. No support, outcome, model, threshold or scientific rule is newly selected by this module.

The recovered historical archive had 1,567 of 1,600 scalar boundary RDS files. The 33 absent records were P1 j=1 lower and P2 j=1–16 lower/upper. Scalar bit-exact historical equality cannot be established for those absent records. This package supplies fresh execution, not historical comparison evidence, and includes no historical RDS, result tables, manuscript targets or prior run logs. Auxiliary runtime package versions are pinned to the reconstruction environment; historical auxiliary versions were not separately captured. See `licenses/NOTICE.md` and the official sources in the download lock for third-party provenance.
