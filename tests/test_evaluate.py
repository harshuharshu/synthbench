# -*- coding: utf-8 -*-
import numpy as np
import pytest

import synthbench as sb

BOOT = 30


@pytest.fixture(scope="module")
def truth():
    ref = sb.load_reference_sets()
    return dict(zip(ref["smiles"], ref["made"]))


def test_perfect_oracle(truth):
    r = sb.evaluate(lambda s: float(truth.get(s, 0)), name="oracle", n_boot=BOOT)
    assert r.metrics["makeability_auc"] == pytest.approx(1.0)
    assert r.metrics["hard_case_auc"] == pytest.approx(1.0)
    assert r.n_scored == len(sb.load_reference_sets())


def test_inverted_scorer_is_reported_not_hidden(truth):
    r = sb.evaluate(lambda s: 1.0 - float(truth.get(s, 0)), name="inverted", n_boot=BOOT)
    assert r.metrics["makeability_auc"] == pytest.approx(0.0)
    assert "Below chance" in r.summary()


def test_random_scorer_sits_at_chance():
    rs = np.random.default_rng(3)
    r = sb.evaluate(lambda s: float(rs.random()), name="random", n_boot=200)
    assert r.metrics["ci_low"] < 0.5 < r.metrics["ci_high"]


def test_vectorised_and_per_molecule_agree(truth):
    per = sb.evaluate(lambda s: float(len(s)), name="per", direction="lower_easier", n_boot=BOOT)
    vec = sb.evaluate(lambda lst: [float(len(s)) for s in lst], name="vec",
                      direction="lower_easier", n_boot=BOOT)
    assert per.metrics["makeability_auc"] == pytest.approx(vec.metrics["makeability_auc"])


def test_object_with_score_method(truth):
    class Scorer:
        def score(self, s):
            return float(truth.get(s, 0))
    r = sb.evaluate(Scorer(), name="method", n_boot=BOOT)
    assert r.metrics["makeability_auc"] == pytest.approx(1.0)


def test_scorer_that_raises_or_returns_none_is_tolerated():
    def flaky(s):
        if len(s) % 3 == 0:
            raise ValueError("no")
        if len(s) % 3 == 1:
            return None
        return float(len(s))
    r = sb.evaluate(flaky, name="flaky", direction="lower_easier", n_boot=BOOT)
    assert 0 < r.n_scored < len(sb.load_reference_sets())
    assert np.isfinite(r.metrics["makeability_auc"])


def test_confirmed_subset_is_reported(truth):
    r = sb.evaluate(lambda s: float(truth.get(s, 0)), name="oracle", n_boot=BOOT)
    assert np.isfinite(r.metrics["confirmed_subset_auc"]), \
        "the 57-molecule confirmed subset must be scored, not dropped"


def test_bad_direction_rejected():
    with pytest.raises(ValueError):
        sb.evaluate(lambda s: 1.0, direction="sideways", n_boot=BOOT)


def test_wrong_length_rejected():
    with pytest.raises(ValueError):
        sb.evaluate(lambda lst: [1.0, 2.0], name="short", n_boot=BOOT)


def test_compare_to_published_places_the_baseline_above_every_published_score(truth):
    r = sb.evaluate(lambda s: float(truth.get(s, 0)), name="oracle", n_boot=BOOT)
    t = r.compare_to_published().set_index("scorer")
    base = t.loc["catalogue_sim", "makeability_auc"]
    for col in sb.PUBLISHED_SCORES:
        assert t.loc[col, "makeability_auc"] < base, col
    assert t.loc["oracle", "makeability_auc"] > base


def test_summary_mentions_the_baseline_and_is_a_string(truth):
    r = sb.evaluate(lambda s: float(truth.get(s, 0)) * 0.2 + 0.4, name="weak", n_boot=BOOT)
    s = r.summary()
    assert isinstance(s, str) and "nearest-neighbour baseline" in s


def test_to_frame_has_one_row_and_the_documented_columns(truth):
    r = sb.evaluate(lambda s: float(truth.get(s, 0)), name="oracle", n_boot=BOOT)
    f = r.to_frame()
    assert len(f) == 1
    for c in ["scorer", "makeability_auc", "ci_low", "ci_high", "hard_case_auc",
              "easy_negative_dependence", "confirmed_subset_auc", "set_dependence_spread"]:
        assert c in f.columns


def test_scorer_adapters_degrade_with_a_reason():
    avail = sb.scorers.available()
    assert set(avail) == set(sb.scorers.DIRECTION)
    for name, status in avail.items():
        assert status == "ready" or len(status) > 10, (name, status)
