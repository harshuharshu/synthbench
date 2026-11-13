# -*- coding: utf-8 -*-
import numpy as np
import pytest

import synthbench as sb
from synthbench import metrics as M


@pytest.fixture(scope="module")
def bench():
    return sb.main_benchmark()


def test_reference_sets_load_with_documented_counts():
    ref = sb.load_reference_sets()
    for name, meta in sb.PROVENANCE.items():
        assert (ref["set"] == name).sum() == meta["n"], name
    assert ref["made"].isin([0, 1]).all()
    assert ref["smiles"].notna().all()


def test_published_score_columns_have_no_nulls(bench):
    for col in sb.PUBLISHED_SCORES + ("catalogue_sim",):
        assert bench[col].notna().all(), col


def test_no_molecule_carries_contradictory_labels():
    ref = sb.load_reference_sets()
    assert ref.groupby("smiles")["made"].nunique().max() == 1


def test_set_weights_equalise_sets(bench):
    w = M.set_weights(bench["set"].values, bench["made"].values)
    for cls in (0, 1):
        m = bench["made"].values == cls
        totals = {s: w[m & (bench["set"].values == s)].sum()
                  for s in np.unique(bench["set"].values[m])}
        assert np.allclose(list(totals.values()), 1.0 / len(totals))


def test_oracle_random_and_inverted(bench):
    made = bench["made"].values
    sets = bench["set"].values
    assert M.weighted_auc(made.astype(float), made, sets) == pytest.approx(1.0)
    assert M.weighted_auc(1.0 - made, made, sets) == pytest.approx(0.0)
    rs = np.random.default_rng(0)
    rand = rs.random(len(made))
    auc = M.weighted_auc(rand, made, sets)
    lo, hi = M.bootstrap_ci(rand, made, sets, n_boot=120, seed=1)
    assert 0.45 < auc < 0.55
    assert lo < 0.5 < hi, "a random scorer's CI must contain chance"


def test_direction_flag_inverts(bench):
    made, sets = bench["made"].values, bench["set"].values
    up = M.weighted_auc(bench["sascore"].values, made, sets, higher_is_easier=True)
    down = M.weighted_auc(bench["sascore"].values, made, sets, higher_is_easier=False)
    assert up + down == pytest.approx(1.0, abs=1e-9)


def test_bootstrap_is_reproducible_under_a_fixed_seed(bench):
    args = (bench["rascore"].values, bench["made"].values, bench["set"].values)
    assert M.bootstrap_ci(*args, n_boot=60, seed=11) == M.bootstrap_ci(*args, n_boot=60, seed=11)


def test_nan_scores_are_excluded_not_imputed(bench):
    v = bench["sascore"].values.astype(float).copy()
    full = M.weighted_auc(v, bench["made"].values, bench["set"].values, higher_is_easier=False)
    v[::10] = np.nan
    part = M.weighted_auc(v, bench["made"].values, bench["set"].values, higher_is_easier=False)
    assert np.isfinite(part) and abs(part - full) < 0.05


def test_degenerate_inputs_return_nan(bench):
    made, sets = bench["made"].values, bench["set"].values
    assert np.isnan(M.weighted_auc(np.full(len(made), np.nan), made, sets))
    one_class = bench[bench["made"] == 1]
    assert np.isnan(M.weighted_auc(one_class["sascore"].values, one_class["made"].values,
                                   one_class["set"].values))


@pytest.mark.parametrize("col,higher", [("sascore", False), ("scscore", False),
                                        ("rascore", True), ("brsascore", False),
                                        ("catalogue_sim", True)])
def test_published_values_match_the_recorded_snapshot(bench, col, higher):
    from synthbench.evaluate import published_table
    got = M.all_metrics(bench[col].values, bench["made"].values, bench["set"].values,
                        higher_is_easier=higher, n_boot=30)
    snap = published_table().set_index("scorer").loc[col]
    assert got["makeability_auc"] == pytest.approx(snap["makeability_auc"], abs=2e-3)
    assert got["hard_case_auc"] == pytest.approx(snap["hard_case_auc"], abs=2e-3)
    assert got["easy_negative_dependence"] == pytest.approx(snap["easy_negative_dependence"], abs=2e-3)


def test_leave_one_set_out_covers_every_main_set(bench):
    d = M.leave_one_set_out(bench["catalogue_sim"].values, bench["made"].values,
                            bench["set"].values)
    assert set(d) == set(M.MAIN_SETS)
    assert all(np.isfinite(v) for v in d.values())
