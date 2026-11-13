# -*- coding: utf-8 -*-
"""
The entry point: hand synthbench a scoring function, get back what it actually measures.
"""
from __future__ import annotations

import functools

import numpy as np
import pandas as pd

from . import datasets as _data
from . import metrics as _metrics

DIRECTIONS = ("higher_easier", "lower_easier")

# Published scores on the packaged reference molecules, for comparison without installing weights.
_PUBLISHED_DIRECTION = {"sascore": "lower_easier", "scscore": "lower_easier",
                        "rascore": "higher_easier", "brsascore": "lower_easier"}
_BASELINE_DIRECTION = {"catalogue_sim": "higher_easier"}


def _apply(score_fn, smiles_list):
    """Accept either a vectorised callable or a per-molecule one, and a .score method.

    The vectorised path is tried first. A call that returns a sized sequence of the wrong
    length is a bug in the scorer, not a signal to fall back, so it raises rather than
    silently producing an all-nan column that would read as "scored nothing".
    """
    fn = getattr(score_fn, "score", score_fn)
    try:
        out = fn(smiles_list)
    except Exception:
        out = None
    if out is not None and hasattr(out, "__len__") and not isinstance(out, (str, bytes)):
        if len(out) != len(smiles_list):
            raise ValueError(f"scorer returned {len(out)} values for "
                             f"{len(smiles_list)} molecules")
        try:
            return np.asarray(list(out), dtype=float)
        except (TypeError, ValueError):
            pass
    return np.asarray([_safe(fn, s) for s in smiles_list], dtype=float)


def _safe(fn, s):
    try:
        v = fn(s)
        return float("nan") if v is None else float(v)
    except Exception:
        return float("nan")


class Result:
    """What a scorer did on the benchmark. Numbers only; the reading is in summary()."""

    def __init__(self, name, direction, values, frame, metrics):
        self.name, self.direction = name, direction
        self.values, self.frame, self.metrics = values, frame, metrics

    @property
    def n_scored(self):
        return int(np.isfinite(self.values).sum())

    def to_frame(self) -> pd.DataFrame:
        m = dict(self.metrics)
        m.pop("leave_one_set_out", None)
        m.update(scorer=self.name, direction=self.direction, n_scored=self.n_scored)
        cols = ["scorer", "direction", "n_scored", "makeability_auc", "ci_low", "ci_high",
                "auc_without_easy_negative", "easy_negative_dependence", "hard_case_auc",
                "confirmed_subset_auc", "set_dependence_spread"]
        return pd.DataFrame([m])[cols]

    def leave_one_set_out(self) -> pd.DataFrame:
        d = self.metrics["leave_one_set_out"]
        return pd.DataFrame({"set_removed": list(d), "makeability_auc": list(d.values())})

    def compare_to_published(self) -> pd.DataFrame:
        """This scorer beside the four published scores and the nearest-neighbour baseline."""
        out = pd.concat([self.to_frame(), published_table()], ignore_index=True)
        out["kind"] = np.where(out.scorer == self.name, "yours",
                               np.where(out.scorer.isin(_BASELINE_DIRECTION), "baseline", "published"))
        return out.sort_values("makeability_auc", ascending=False).reset_index(drop=True)

    def summary(self) -> str:
        m = self.metrics
        base = _baseline_auc()
        L = [f"synthbench: {self.name}  ({self.direction}, {self.n_scored} of {len(self.values)} scored)",
             "",
             f"  makeability AUC            {m['makeability_auc']:.3f}  "
             f"(95% CI {m['ci_low']:.3f} to {m['ci_high']:.3f})",
             f"  without the easy negative  {m['auc_without_easy_negative']:.3f}  "
             f"(dependence {m['easy_negative_dependence']:+.3f})",
             f"  hard case, THQ vs BRICS    {m['hard_case_auc']:.3f}",
             f"  confirmed subset, n=57     {m['confirmed_subset_auc']:.3f}",
             f"  set-dependence spread      {m['set_dependence_spread']:.3f}",
             "",
             f"  nearest-neighbour baseline {base:.3f}   "
             f"(you are {m['makeability_auc'] - base:+.3f} against it)",
             ""]
        if not np.isfinite(m["makeability_auc"]):
            L.append("  Nothing scored. Check that the scorer returns a float for a SMILES string.")
            return "\n".join(L)
        if m["makeability_auc"] < 0.5:
            L.append("  Below chance: on this material the score ranks never-made molecules")
            L.append("  above made ones. Check the direction argument before reading further.")
        elif m["makeability_auc"] < base:
            L.append("  Below the trivial baseline. Maximum Tanimoto similarity to a purchasable")
            L.append("  catalogue separates made from never-made better than this score does.")
        else:
            L.append("  Above the trivial baseline, which no published score in the panel manages.")
        if m["easy_negative_dependence"] > 0.15:
            L.append(f"  {m['easy_negative_dependence']:.3f} of this comes from one easy negative class")
            L.append("  (GDB-17 strained cages). Without them the AUC is "
                     f"{m['auc_without_easy_negative']:.3f}.")
        if np.isfinite(m["hard_case_auc"]) and m["hard_case_auc"] < 0.5:
            L.append("  Inverted on the hard case: molecules made in one step rank below molecules")
            L.append("  never made. Every score tested so far fails here, the baseline included.")
        return "\n".join(L)

    __repr__ = summary


@functools.lru_cache(maxsize=1)
def published_table():
    """The panel scored on the packaged molecules.

    Served from a packaged snapshot so the first comparison is instant and needs no
    third-party weights. The snapshot is regenerated by tools/refresh_published.py and
    the test suite checks it still matches a live recomputation.
    """
    from importlib import resources
    try:
        with resources.as_file(resources.files("synthbench.data")
                               .joinpath("published_metrics.csv.gz")) as q:
            return pd.read_csv(q)
    except (FileNotFoundError, ModuleNotFoundError):
        pass
    return recompute_published_table()


def recompute_published_table():
    """Recompute the panel from the packaged molecules, bypassing the snapshot.

    Uses all seven sets, like evaluate(), so the 57-molecule confirmed subset is reported
    rather than coming back as nan. The headline metrics still use only the six main sets.
    """
    df = _data.load_reference_sets()
    rows = []
    for col, direction in {**_PUBLISHED_DIRECTION, **_BASELINE_DIRECTION}.items():
        if col not in df:
            continue
        m = _metrics.all_metrics(df[col].values, df["made"].values, df["set"].values,
                                 higher_is_easier=(direction == "higher_easier"))
        rows.append(Result(col, direction, df[col].values, df, m).to_frame())
    return pd.concat(rows, ignore_index=True)


def _baseline_auc():
    t = published_table()
    return float(t.loc[t.scorer == "catalogue_sim", "makeability_auc"].iloc[0])


def evaluate(score_fn, name="scorer", direction="higher_easier", n_boot=400, seed=7) -> Result:
    """Run a scoring function against the benchmark.

    score_fn      callable mapping a SMILES string to a float, or a list of SMILES to a list of
                  floats, or an object with a .score method of either shape. Returning None or
                  raising for a molecule is allowed and recorded as not scored.
    direction     "higher_easier" or "lower_easier". Declaring it wrong will show up as an AUC
                  below chance, which summary() flags.
    """
    if direction not in DIRECTIONS:
        raise ValueError(f"direction must be one of {DIRECTIONS}, got {direction!r}")
    # score all seven sets: the headline metrics use the six main sets, and the
    # 57-molecule confirmed subset is reported separately rather than dropped.
    df = _data.load_reference_sets()
    vals = _apply(score_fn, df["smiles"].tolist())
    if len(vals) != len(df):
        raise ValueError(f"scorer returned {len(vals)} values for {len(df)} molecules")
    m = _metrics.all_metrics(vals, df["made"].values, df["set"].values,
                             higher_is_easier=(direction == "higher_easier"),
                             n_boot=n_boot, seed=seed)
    return Result(name, direction, vals, df, m)
