# -*- coding: utf-8 -*-
"""
The reference sets, packaged with the library. No network access at import or at load.

Labels come from provenance, never from a score. "made" means there is a record that the
molecule exists. "never made" means no evidence of existence, which is weaker than unmakeable:
no public collection of compounds that chemists attempted and failed to make exists, so this
benchmark has no true negative set and does not claim one.
"""
from __future__ import annotations

import functools
from importlib import resources

import pandas as pd

PUBLISHED_SCORES = ("sascore", "scscore", "rascore", "brsascore")

PROVENANCE = {
    "ZINC_purchasable": dict(
        n=3580, made=1, known="purchasable, therefore made", route_length="unknown",
        source="SYBA easy test set", citation="Vorsilak et al., J. Cheminform. 2020, 12, 35",
        licence="MIT (source repository)"),
    "Enamine_on_demand": dict(
        n=1000, made=1, known="made", route_length="2 to 3 steps by construction",
        source="synformer sample of Enamine make-on-demand",
        citation="Gao et al., arXiv:2410.03494", licence="MIT (source repository)"),
    "ChEMBL_kinase": dict(
        n=2968, made=1, known="made", route_length="unknown",
        source="built from ChEMBL kinase actives",
        citation="Zdrazil et al., Nucleic Acids Res. 2024, 52, D1180", licence="CC BY-SA 3.0"),
    "THQ_library": dict(
        n=10000, made=1, known="made", route_length="one step, demonstrated experimentally",
        source="enumerated aza-Diels-Alder tetrahydroquinoline library",
        citation="Papidocha et al., chemRxiv 2026, doi:10.26434/chemrxiv.15010755/v1",
        licence="CC BY 4.0"),
    "THQ_made_assayed": dict(
        n=57, made=1, known="synthesized and assayed", route_length="one step, experimental",
        source="the subset of the THQ library actually made",
        citation="Papidocha et al., chemRxiv 2026, doi:10.26434/chemrxiv.15010755/v1",
        licence="CC BY 4.0"),
    "BRICS_recombined": dict(
        n=4776, made=0, known="no evidence of existence", route_length="not applicable",
        source="BRICS recombination of ChEMBL JAK inhibitor fragments",
        citation="this work", licence="MIT"),
    "GDB17_enumerated": dict(
        n=3581, made=0, known="no evidence of existence", route_length="not applicable",
        source="SYBA hard test set, drawn from GDB-17",
        citation="Ruddigkeit et al., J. Chem. Inf. Model. 2012, 52, 2864", licence="MIT"),
}

# THQ_made_assayed shares 3 molecules with THQ_library, so the two must not be pooled.
OVERLAPS = {("THQ_library", "THQ_made_assayed"): 3}


@functools.lru_cache(maxsize=1)
def load_reference_sets() -> pd.DataFrame:
    """All seven sets as one frame: set, made, provenance, smiles, catalogue_sim, published scores."""
    with resources.as_file(
            resources.files("synthbench.data").joinpath("reference_sets.csv.gz")) as p:
        df = pd.read_csv(p, compression="gzip")
    return df


def main_benchmark() -> pd.DataFrame:
    """The six sets the headline metrics are computed on, excluding the 57-molecule subset."""
    from .metrics import MAIN_SETS
    return load_reference_sets().query("set in @MAIN_SETS").reset_index(drop=True)


def smiles(main_only: bool = True):
    """The SMILES a scorer will be asked to score, in benchmark order."""
    df = main_benchmark() if main_only else load_reference_sets()
    return df["smiles"].tolist()


def provenance_table() -> pd.DataFrame:
    return pd.DataFrame(PROVENANCE).T.rename_axis("set").reset_index()
