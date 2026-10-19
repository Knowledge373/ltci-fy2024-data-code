# Fixed support selection

Run from the repository root:

    python support_selection/run_selection.py --out /path/to/new/support-output

The output directory must not exist. This module uses Python's standard library.
It computes the three fixed contrast support sets at the fixed 1-SD and 0.5-SD
calipers. An edge requires the same prefecture and simultaneous satisfaction of
all four covariate calipers; retained observations have at least one edge.

`math_early.py` and `math_p3.py` contain unchanged mathematical functions from the
support-set builders. The adapter uses the final stored covariates and the supplied
monthly institutional classification input. It uses Fraction arithmetic for the
three shares and Decimal arithmetic at precision 80 for the stored log count.
It does not use collection outcomes or monetary values in selection.

This reproduces selection from final covariates, not the unavailable raw-count
source reconstruction. Stored share/log precision can differ from that of original
raw-count calculations; no scales are tuned. Compact retained-ID hash pins verify
the resulting fixed sets. There are no stored membership or pair-decision results
in this module. Those files are computed afresh into the chosen output directory.
