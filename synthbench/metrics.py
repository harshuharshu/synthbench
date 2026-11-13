# -*- coding: utf-8 -*-
"""
The four quantities synthbench reports, and the weighting they all share.

Every AUC here is computed with per-set sample weights. Reference set sizes span 1,000 to
10,000, so an unweighted pooled AUC would largely report the largest set. Weighting each set
equally within its class makes the benchmark a statement about agreement across sets.
"""
from __future__ import annotations

import numpy as np
from sklearn.metrics import roc_auc_score

MAIN_SETS = ("ZINC_purchasable", "Enamine_on_demand", "ChEMBL_kinase", "THQ_library",
             "BRICS_recombined", "GDB17_enumerated")
EASY_NEGATIVE = "GDB17_enumerated"
HARD_CASE = ("THQ_library", "BRICS_recombined")
CONFIRMED_CASE = ("THQ_made_assayed", "BRICS_recombined")


def set_weights(sets, made):
    """Equal total weight per set within each class, so set sizes cannot drive the result."""
    sets = np.asarray(sets)
    made = np.asarray(made)
    w = np.ones(len(sets), dtype=float)
    for cls in (0, 1):
        in_cls = made == cls
        names = np.unique(sets[in_cls])
        for s in names:
            m = in_cls & (sets == s)
            w[m] = 1.0 / (len(names) * m.sum())
    return w


def weighted_auc(score, made, sets, higher_is_easier=True):
    """Set-weighted AUC for separating made (1) from never-made (0).

    Returns nan when the subset does not contain both classes or no finite scores.
    """
    score = np.asarray(score, dtype=float)
    made = np.asarray(made)
    sets = np.asarray(sets)
    ok = np.isfinite(score)
    if ok.sum() == 0 or len(np.unique(made[ok])) < 2:
        return float("nan")
    s = score[ok] if higher_is_easier else -score[ok]
    return float(roc_auc_score(made[ok], s, sample_weight=set_weights(sets[ok], made[ok])))


def bootstrap_ci(score, made, sets, higher_is_easier=True, n_boot=400, seed=7, alpha=0.05):
    """Percentile CI, resampling within each set.

    This interval describes sampling within these particular reference sets. It does not
    capture uncertainty about which sets were chosen; use leave_one_set_out for that.
    """
    score, made, sets = np.asarray(score, float), np.asarray(made), np.asarray(sets)
    rs = np.random.default_rng(seed)
    idx = {s: np.flatnonzero(sets == s) for s in np.unique(sets)}
    vals = []
    for _ in range(n_boot):
        take = np.concatenate([rs.choice(ix, len(ix), replace=True) for ix in idx.values()])
        a = weighted_auc(score[take], made[take], sets[take], higher_is_easier)
        if np.isfinite(a):
            vals.append(a)
    if not vals:
        return (float("nan"), float("nan"))
    return (float(np.percentile(vals, 100 * alpha / 2)),
            float(np.percentile(vals, 100 * (1 - alpha / 2))))


def _subset(score, made, sets, keep):
    m = np.isin(np.asarray(sets), list(keep))
    return np.asarray(score, float)[m], np.asarray(made)[m], np.asarray(sets)[m]


def leave_one_set_out(score, made, sets, higher_is_easier=True, main_sets=MAIN_SETS):
    """Makeability AUC with each set removed in turn.

    The spread of this mapping is how much a reported AUC depends on the choice of reference
    sets, and in practice is dominated by which negative set is present.
    """
    out = {}
    for s in main_sets:
        keep = [x for x in main_sets if x != s]
        out[s] = weighted_auc(*_subset(score, made, sets, keep), higher_is_easier=higher_is_easier)
    return out


def all_metrics(score, made, sets, higher_is_easier=True, n_boot=400, seed=7):
    """The four reported quantities, plus the two diagnostic sub-comparisons."""
    main = _subset(score, made, sets, MAIN_SETS)
    no_easy = _subset(score, made, sets, [s for s in MAIN_SETS if s != EASY_NEGATIVE])

    auc = weighted_auc(*main, higher_is_easier=higher_is_easier)
    auc_no_easy = weighted_auc(*no_easy, higher_is_easier=higher_is_easier)
    lo, hi = bootstrap_ci(*main, higher_is_easier=higher_is_easier, n_boot=n_boot, seed=seed)
    loso = leave_one_set_out(score, made, sets, higher_is_easier)
    finite = [v for v in loso.values() if np.isfinite(v)]

    return {
        "makeability_auc": auc,
        "ci_low": lo,
        "ci_high": hi,
        "auc_without_easy_negative": auc_no_easy,
        "easy_negative_dependence": auc - auc_no_easy,
        "hard_case_auc": weighted_auc(*_subset(score, made, sets, HARD_CASE),
                                      higher_is_easier=higher_is_easier),
        "confirmed_subset_auc": weighted_auc(*_subset(score, made, sets, CONFIRMED_CASE),
                                             higher_is_easier=higher_is_easier),
        "leave_one_set_out": loso,
        "set_dependence_spread": (max(finite) - min(finite)) if finite else float("nan"),
    }
