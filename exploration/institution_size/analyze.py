#!/usr/bin/env python3
"""Fixed institutional and insurer-size descriptive calculations."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(os.environ['EXPLORATION_OUTPUT_DIR'])
ANNUAL = list(range(3, 13))
QLEVELS = ["Q1", "Q2", "Q3", "Q4"]
BASE_COLS = ["insref_id", "fiscal_year", "prefecture_code", "annual_installments",
             "provisional_none", "collection_ratio_fraction", "income_low3_share",
             "age_75_84_share", "age_85plus_share", "ln_insured_end_fy_count"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def csv(name, rows):
    tab = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    tab.to_csv(OUT / name, index=False, na_rep="NA", float_format="%.12g")
    return tab


def summary(g, national):
    n = len(g)
    empty = {x: np.nan for x in [
        "insured_sum", "insured_mean", "insured_median", "insured_p10", "insured_p90",
        "insured_min", "insured_max", "insured_share_national_pct",
        "ordinary_assessed_sum", "ordinary_assessed_share_national_pct",
        "special_assessed_sum", "total_assessed_sum", "total_assessed_share_national_pct",
        "ordinary_share_of_total_assessed_pct", "rate_mean_pct", "rate_median_pct",
        "rate_sd_pp", "rate_p10_pct", "rate_p90_pct", "rate_min_pct", "rate_max_pct",
        "rate_exact100_pct", "rate_ge99_5_pct", "rate_ge99_9_pct", "rate_above100_pct"]}
    out = dict(n=n, support="empty" if n == 0 else "sparse_n_lt_10" if n < 10 else "n_ge_10", **empty)
    if not n:
        return out
    rates = g["collection_ratio_fraction"] * 100
    out.update(
        insured_sum=g["N"].sum(), insured_mean=g["N"].mean(), insured_median=g["N"].median(),
        insured_p10=g["N"].quantile(.1), insured_p90=g["N"].quantile(.9),
        insured_min=g["N"].min(), insured_max=g["N"].max(),
        insured_share_national_pct=g["N"].sum()/national["N"]*100,
        ordinary_assessed_sum=g["A"].sum(), ordinary_assessed_share_national_pct=g["A"].sum()/national["A"]*100,
        special_assessed_sum=g["S"].sum(), total_assessed_sum=g["T"].sum(),
        total_assessed_share_national_pct=g["T"].sum()/national["T"]*100,
        ordinary_share_of_total_assessed_pct=g["A"].sum()/g["T"].sum()*100 if g["T"].sum() != 0 else np.nan,
        rate_mean_pct=rates.mean(), rate_median_pct=rates.median(), rate_sd_pp=rates.std(ddof=1),
        rate_p10_pct=rates.quantile(.1), rate_p90_pct=rates.quantile(.9),
        rate_min_pct=rates.min(), rate_max_pct=rates.max(),
        rate_exact100_pct=g["collection_ratio_fraction"].eq(1).mean()*100,
        rate_ge99_5_pct=g["collection_ratio_fraction"].ge(.995).mean()*100,
        rate_ge99_9_pct=g["collection_ratio_fraction"].ge(.999).mean()*100,
        rate_above100_pct=g["collection_ratio_fraction"].gt(1).mean()*100)
    return out


def group_table(df, cols, levels, national):
    rows = []
    for values in itertools.product(*levels):
        m = np.ones(len(df), dtype=bool)
        for col, val in zip(cols, values):
            m &= df[col].eq(val).to_numpy()
        rows.append(dict(zip(cols, values)) | summary(df.loc[m], national))
    return pd.DataFrame(rows)


def contrast(df, label, ma, mb):
    """A minus B. Prelisted two-part symmetric descriptive decomposition."""
    a, b = df.loc[ma], df.loc[mb]
    full_weights = df["size_quartile"].value_counts(normalize=True).reindex(QLEVELS)
    rows = []
    for q in QLEVELS:
        qa, qb = a.loc[a["size_quartile"].eq(q)], b.loc[b["size_quartile"].eq(q)]
        ym_a = qa["collection_ratio_fraction"].mean()*100
        ym_b = qb["collection_ratio_fraction"].mean()*100
        med_a = qa["collection_ratio_fraction"].median()*100
        med_b = qb["collection_ratio_fraction"].median()*100
        wa, wb = len(qa)/len(a), len(qb)/len(b)
        rows.append(dict(contrast=label, size_quartile=q, n_a=len(qa), n_b=len(qb),
                         mean_a_pct=ym_a, mean_b_pct=ym_b, mean_a_minus_b_pp=ym_a-ym_b,
                         median_a_pct=med_a, median_b_pct=med_b, median_a_minus_b_pp=med_a-med_b,
                         weight_a=wa, weight_b=wb, full_population_weight=full_weights[q],
                         within_term_pp=(wa+wb)/2*(ym_a-ym_b),
                         composition_term_pp=(ym_a+ym_b)/2*(wa-wb)))
    tab = pd.DataFrame(rows)
    supported = bool((tab[["n_a", "n_b"]] > 0).all().all())
    raw = (a["collection_ratio_fraction"].mean()-b["collection_ratio_fraction"].mean())*100
    standardized = (tab["full_population_weight"]*tab["mean_a_minus_b_pp"]).sum() if supported else np.nan
    within = tab["within_term_pp"].sum() if supported else np.nan
    comp = tab["composition_term_pp"].sum() if supported else np.nan
    result = dict(contrast=label, direction="A minus B", n_a=len(a), n_b=len(b),
                  status="all_four_size_strata_supported" if supported else "unidentified_missing_group_stratum",
                  raw_mean_gap_pp=raw, raw_median_gap_pp=(a["collection_ratio_fraction"].median()-b["collection_ratio_fraction"].median())*100,
                  full_population_standardized_gap_pp=standardized,
                  symmetric_within_component_pp=within, symmetric_composition_component_pp=comp,
                  reconstruction_error_pp=raw-within-comp if supported else np.nan)
    if supported:
        assert abs(result["reconstruction_error_pp"]) < 1e-9
    return tab, result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    path = Path(cfg["input_path"])
    assert sha(path) == cfg["input_sha256"], "Common input hash mismatch"
    raw = pd.read_csv(path, dtype=str, keep_default_na=False)
    assert raw.shape == (1573, 15)
    assert raw["insref_id"].is_unique
    assert not raw.eq("").any().any()
    base = pd.read_csv(cfg["baseline_path"], dtype=str, keep_default_na=False)
    assert raw[BASE_COLS].equals(base[BASE_COLS]), "Original ten fields differ from baseline"
    df = raw.copy()
    for col in raw.columns.difference(["insref_id", "prefecture_code"]):
        df[col] = pd.to_numeric(raw[col], errors="raise")
    mapping = cfg["columns"]
    for alias, source in mapping.items():
        df[alias] = df[source]
    assert set(mapping) == {"N", "A", "C", "S", "T"}
    assert np.isfinite(df[list(mapping)].to_numpy()).all()
    assert (df["N"] > 0).all()
    assert (df["N"] == df["N"].astype(np.int64)).all()
    assert (df[["A", "S", "T"]] >= 0).all().all()
    assert (df["T"] == df["A"] + df["S"]).all()
    assert (df["A"] > 0).all(), "Nonpositive ordinary assessment requires explicit diagnostics"
    assert np.allclose(df["C"]/df["A"], df["collection_ratio_fraction"], rtol=0, atol=1e-12)
    assert set(df["annual_installments"]) == set(ANNUAL)
    assert set(df["provisional_none"]) == {0, 1}
    cut = df["N"].quantile([.25, .5, .75], interpolation="linear").to_numpy()
    assert len(set(cut)) == 3, "Prelisted four-stratum design unavailable: duplicate thresholds"
    df["size_quartile"] = pd.cut(df["N"], bins=[-np.inf, *cut, np.inf], right=True, labels=QLEVELS).astype(str)
    # Record only the derived assignment, not a competing 15-column input.
    csv("derived_size_quartile_assignments.csv", df[["insref_id", "size_quartile"]])
    national = df[list(mapping)].sum().to_dict()
    tables = {}
    specs = {
        "size_quartiles": (["size_quartile"], [QLEVELS]),
        "provisional_groups": (["provisional_none"], [[0, 1]]),
        "institution_cells20": (["annual_installments", "provisional_none"], [ANNUAL, [0, 1]]),
        "size_by_provisional8": (["size_quartile", "provisional_none"], [QLEVELS, [0, 1]]),
        "size_by_institution_cells80": (["annual_installments", "provisional_none", "size_quartile"], [ANNUAL, [0, 1], QLEVELS]),
    }
    csv("national_summary.csv", [summary(df, national)])
    for name, (cols, levels) in specs.items():
        tables[name] = csv(name+".csv", group_table(df, cols, levels, national))
    assert tables["institution_cells20"]["n"].gt(0).sum() == 16
    assert tables["institution_cells20"]["n"].sum() == 1573
    assert tables["size_by_institution_cells80"]["n"].sum() == 1573
    decomp, contrasts = [], []
    for label, ma, mb in [
        ("provisional_none_minus_yes", df["provisional_none"].eq(1), df["provisional_none"].eq(0)),
        ("annual12_yes_minus_annual8_none", df["annual_installments"].eq(12)&df["provisional_none"].eq(0), df["annual_installments"].eq(8)&df["provisional_none"].eq(1)),
    ]:
        t, r = contrast(df, label, ma, mb)
        decomp.append(t)
        contrasts.append(r)
    dtab = csv("contrasts_by_size_quartile.csv", pd.concat(decomp, ignore_index=True))
    ctab = csv("contrast_decomposition.csv", contrasts)
    qtab = tables["size_quartiles"]
    bounds = []
    for i, q in enumerate(QLEVELS):
        g = df.loc[df["size_quartile"].eq(q)]
        bounds.append(dict(size_quartile=q, n=len(g), observed_n_min=int(g["N"].min()), observed_n_max=int(g["N"].max()),
                           lower_boundary_exclusive=None if i == 0 else float(cut[i-1]),
                           upper_boundary_inclusive=None if i == 3 else float(cut[i])))
    dump(OUT/"quartile_definition.json", {"method":"global original-N quantiles, linear interpolation, right-closed bins, tied N unsplit", "cutpoints":[float(x) for x in cut], "strata":bounds})
    plot(tables, dtab)
    dump(OUT/"software_versions.json", dict(python=sys.version, platform=platform.platform(), numpy=np.__version__, pandas=pd.__version__, matplotlib=matplotlib.__version__, executable=sys.executable))
    dump(OUT/"audit_receipt.json", dict(run_time_utc=datetime.now(timezone.utc).isoformat(), input_path=str(path), input_sha256=sha(path),
        baseline_path=cfg["baseline_path"], baseline_sha256=sha(cfg["baseline_path"]), baseline_original10_exact_equal=True,
        input_dimensions=list(raw.shape), input_ids_unique=True, missing_input_values=0,
        amount_units=cfg["amount_units"], column_mapping=mapping, total_identity_exact=True,
        rate_from_collected_over_assessed_max_abs_error=float(np.max(np.abs(df["C"]/df["A"]-df["collection_ratio_fraction"]))),
        rate_above_one_count=int(df["collection_ratio_fraction"].gt(1).sum()), rate_exact_one_count=int(df["collection_ratio_fraction"].eq(1).sum()),
        possible_institution_cells=20, observed_institution_cells=16, possible_institution_size_cells=80,
        observed_institution_size_cells=int(tables["size_by_institution_cells80"]["n"].gt(0).sum()),
        sparse_nonempty_institution_size_cells=int(tables["size_by_institution_cells80"]["n"].between(1,9).sum()),
        decomposition_reconstruction_max_abs_error_pp=float(ctab["reconstruction_error_pp"].abs().max()),
        analysis_code_sha256=sha(__file__)))
    print("Quartile definitions")
    print(json.dumps(bounds, ensure_ascii=False))
    print("Size summaries")
    print(qtab.to_string(index=False))
    print("Institution summaries")
    print(tables["institution_cells20"].to_string(index=False))
    print("Contrasts")
    print(ctab.to_string(index=False))
    print(dtab.to_string(index=False))


def plot(tables, dtab):
    fig = plt.figure(figsize=(14, 10), layout="constrained")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1])
    ax0 = fig.add_subplot(gs[:, 0])
    tab = tables["size_by_institution_cells80"]
    p = tab.pivot(index=["annual_installments", "provisional_none"], columns="size_quartile", values="n").reindex(columns=QLEVELS)
    counts = p.to_numpy()
    masked = np.ma.masked_where(counts == 0, counts)
    cm = plt.get_cmap("Blues").copy()
    cm.set_bad("#ededed")
    im = ax0.imshow(masked, cmap=cm, aspect="auto", vmin=0, vmax=max(1,counts.max()))
    ax0.set_yticks(np.arange(len(p)), [f"{a} / {'none' if v else 'yes'}" for a,v in p.index])
    ax0.set_xticks(range(4), QLEVELS)
    for i in range(counts.shape[0]):
        for j in range(counts.shape[1]):
            ax0.text(j, i, str(int(counts[i,j])) if counts[i,j] else "empty", ha="center", va="center", fontsize=9,
                     color="white" if counts[i,j] > counts.max()*.65 else "#222222")
    ax0.set_title("Institutional support across insurer-size quartiles", loc="left", fontsize=12)
    ax0.set_xlabel("Quartile of original year-end insured count")
    ax0.set_ylabel("Annual installments / adopted provisional classification")
    fig.colorbar(im, ax=ax0, shrink=.35, label="Insurers (n)")
    ax1 = fig.add_subplot(gs[0, 1])
    q = tables["size_quartiles"].set_index("size_quartile").loc[QLEVELS]
    x = np.arange(4)
    ax1.vlines(x, q["rate_p10_pct"], q["rate_p90_pct"], color="#8da0ae", linewidth=6, label="10th to 90th percentile")
    ax1.scatter(x-.045, q["rate_mean_pct"], color="#b2472d", label="Mean", zorder=3)
    ax1.scatter(x+.045, q["rate_median_pct"], color="#24587e", marker="D", label="Median", zorder=3)
    ax1.set_xticks(x, [f"{qv}\nn={int(q.loc[qv,'n'])}" for qv in QLEVELS])
    ax1.set_ylabel("Ordinary current collection rate (%)")
    ax1.set_title("Collection-rate distribution by insurer size", loc="left", fontsize=12)
    ax1.grid(axis="y", alpha=.2)
    ax1.legend(fontsize=9, loc="upper right")
    ax2 = fig.add_subplot(gs[1, 1])
    for label, color, marker, display in [
        ("provisional_none_minus_yes", "#b2472d", "o", "None minus yes"),
        ("annual12_yes_minus_annual8_none", "#24587e", "s", "12/yes minus 8/none")]:
        t = dtab.loc[dtab["contrast"].eq(label)].set_index("size_quartile").loc[QLEVELS]
        ax2.plot(x, t["mean_a_minus_b_pp"], color=color, marker=marker, label=display)
    ax2.axhline(0, color="#777777", linewidth=.8)
    ax2.set_xticks(x, QLEVELS)
    ax2.set_ylabel("Difference in equal-insurer mean (pp)")
    ax2.set_xlabel("Quartile of original year-end insured count")
    ax2.set_title("Prelisted descriptive contrasts within each quartile", loc="left", fontsize=12)
    ax2.grid(axis="y", alpha=.2)
    ax2.legend(fontsize=9)
    fig.suptitle("FY2024: institutional support and collection distributions across insurer size", fontsize=15)
    fig.supxlabel("Descriptive, not causal adjustment. Annual 8 is the existing model reference. Gray cells have no observations; no rate is imputed.", fontsize=10)
    fig.savefig(OUT/"institution_size_patterns.png", dpi=180)
    fig.savefig(OUT/"institution_size_patterns.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
