"""Unchanged independent QR-projector computations for the fixed 16-cell model."""
from pathlib import Path
import json
import numpy as np
from scipy import linalg, stats
OUT = Path(__file__).resolve().parent
ATOL, RTOL, DF_ATOL, CUTOFF = 1e-10, 1e-8, 1e-8, 1e-12
CHECKS = []

def write_json(name, data):
    (OUT / name).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n")

def check(label, ok, detail=None):
    CHECKS.append({"check": label, "pass": bool(ok), "detail": detail})
    write_json("CHECKS.json", CHECKS)
    if not ok:
        raise RuntimeError(f"HOLD: {label}: {detail}")

def compare(label, a, b, df=False):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    check(label + "_shape", a.shape == b.shape)
    if a.size == 0:
        raise ValueError("Empty numerical comparison")
    diff = np.abs(a - b)
    allowance = (DF_ATOL if df else ATOL) + RTOL * np.maximum(np.abs(a), np.abs(b))
    detail = {"max_absolute_difference": float(np.max(diff)), "max_tolerance_ratio": float(np.max(diff / allowance))}
    check(label, np.all(np.isfinite(a)) and np.all(np.isfinite(b)) and np.all(diff <= allowance), detail)
    return detail

def fit_representation(Z, y, pref, tag):
    singular_values = linalg.svdvals(Z)
    ranks = {str(t): int(np.count_nonzero(singular_values > singular_values[0] * t)) for t in (1e-7, 1e-9, 1e-11)}
    check(tag + "_rank66_all_frozen_thresholds", all(v == 66 for v in ranks.values()), ranks)
    Q, R = linalg.qr(Z, mode="economic")
    beta = linalg.solve_triangular(R, Q.T @ y)
    fitted = Q @ (Q.T @ y)
    e = y - fitted
    clusters = []
    eigenaudit = []
    for p in range(1, 48):
        idx = np.flatnonzero(pref == p)
        Qj = Q[idx]
        Mjj = np.eye(len(idx)) - Qj @ Qj.T
        eigenvalues, vectors = linalg.eigh(Mjj)
        check(f"{tag}_pref{p:02d}_no_material_negative_eigen", eigenvalues.min() >= -1e-10, {"minimum": float(eigenvalues.min())})
        keep = eigenvalues > CUTOFF
        retained = {str(c): int(np.count_nonzero(eigenvalues > c)) for c in (1e-14, 1e-12, 1e-10)}
        check(f"{tag}_pref{p:02d}_retention_invariant", len(set(retained.values())) == 1, retained)
        positive = eigenvalues[keep]
        A = (vectors[:, keep] * (1 / np.sqrt(positive))) @ vectors[:, keep].T
        clusters.append((idx, Qj, A))
        eigenaudit.append({"prefecture": p, "n": len(idx), "minimum_eigenvalue": float(eigenvalues.min()), "minimum_retained_eigenvalue": float(positive.min()), "max_eigenvalue": float(eigenvalues.max()), "suppressed": int(np.count_nonzero(~keep)), "retained_by_cutoff": retained})
    check(tag + "_expected50_suppressed_FE_singleton_directions", sum(r["suppressed"] for r in eigenaudit) == 50, sum(r["suppressed"] for r in eigenaudit))
    write_json(tag + "_EIGEN_RANK_AUDIT.json", {"ranks": ranks, "condition_number": float(singular_values[0] / singular_values[-1]), "cluster_blocks": eigenaudit})
    return {"Q": Q, "R": R, "beta": beta, "fitted": fitted, "e": e, "clusters": clusters}

def infer_from_projector(fit, c, y):
    # A distinct residual-projector route, without the author's GH subtraction.
    Q, R, e = fit["Q"], fit["R"], fit["e"]
    w = linalg.solve_triangular(R.T, c, lower=True)
    q = np.empty((len(y), 47))
    scores = np.empty(47)
    for j, (idx, Qj, A) in enumerate(fit["clusters"]):
        gj = A @ (Qj @ w)
        q[:, j] = -Q @ (Qj.T @ gj)
        q[idx, j] += gj
        scores[j] = gj @ e[idx]
    P = q.T @ q
    variance = float(scores @ scores)
    df = float(np.trace(P) ** 2 / np.sum(P ** 2))
    estimate = float(c @ fit["beta"])
    se = float(np.sqrt(variance))
    critical = float(stats.t.ppf(.975, df))
    lo, hi = estimate - critical * se, estimate + critical * se
    if not np.all(np.isfinite([estimate, variance, se, df, lo, hi])) or not min(variance, se, df) > 0:
        raise RuntimeError("HOLD: invalid finite/positive variance, SE or df")
    # A second identity within this route checks residual-score construction.
    if not np.allclose(q.T @ y, scores, atol=ATOL, rtol=RTOL):
        raise RuntimeError("HOLD: q^T y score identity")
    return {"estimate_fraction": estimate, "variance_fraction2": variance, "SE_fraction": se, "df": df, "critical_t": critical, "CI_lower_fraction": lo, "CI_upper_fraction": hi, "estimate_pp": 100 * estimate, "variance_pp2": 10000 * variance, "SE_pp": 100 * se, "CI_lower_pp": 100 * lo, "CI_upper_pp": 100 * hi, "P_trace": float(np.trace(P)), "P_squared_frobenius": float(np.sum(P ** 2))}, scores
