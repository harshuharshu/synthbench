# -*- coding: utf-8 -*-
"""
Optional adapters for the four published scores, for scoring molecules outside the benchmark.

You do not need any of this to use synthbench. The published panel is already computed on the
packaged reference molecules, so `Result.compare_to_published()` works with nothing installed.
These adapters exist for the case where you want to run the same panel on your own molecules.

No third-party weights are redistributed here. `fetch_weights()` downloads them from the
authors' own repositories into a local cache. Each scorer reports itself unavailable with a
reason rather than raising, so a partial install degrades to a smaller panel.

    from synthbench import scorers
    scorers.available()                 # what is ready, and why the rest is not
    scorers.fetch_weights()             # download what is missing
    sa = scorers.load("sascore")
    sa.score("CC(=O)Oc1ccccc1C(=O)O")
"""
from __future__ import annotations

import os
import pathlib

DIRECTION = {"sascore": "lower_easier", "scscore": "lower_easier",
             "rascore": "higher_easier", "brsascore": "lower_easier"}

CITATION = {
    "sascore":   "Ertl & Schuffenhauer, J. Cheminform. 2009, 1, 8",
    "scscore":   "Coley et al., J. Chem. Inf. Model. 2018, 58, 252",
    "rascore":   "Thakkar et al., Chem. Sci. 2021, 12, 3339",
    "brsascore": "Chen et al., J. Cheminform. 2024, 16, 83",
}

# Author repositories the weights come from. Paths are resolved through the GitHub API so that
# a fork or a moved default branch does not silently fetch the wrong file.
WEIGHTS = {
    "scscore":   [("connorcoley", "scscore",
                   "models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz",
                   "scscore/model.ckpt-10654.as_numpy.json.gz")],
    "rascore":   [("reymond-group", "RAscore",
                   "RAscore/models/XGB_chembl_ecfp_counts/model.pkl",
                   "rascore/XGB_chembl_ecfp_counts_model.pkl")],
    "brsascore": [("awadell1", "BR-SAScore", "BRSAScore/pickle/BRScores_uspto_emolecules.pkl.gz",
                   "brsascore/pickle/BRScores_uspto_emolecules.pkl.gz"),
                  ("awadell1", "BR-SAScore", "BRSAScore/pickle/BScores_emolecules.pkl.gz",
                   "brsascore/pickle/BScores_emolecules.pkl.gz"),
                  ("awadell1", "BR-SAScore", "BRSAScore/pickle/RScores_uspto.pkl.gz",
                   "brsascore/pickle/RScores_uspto.pkl.gz")],
}


def cache_dir() -> pathlib.Path:
    """Where fetched weights live. Override with SYNTHBENCH_CACHE."""
    p = os.environ.get("SYNTHBENCH_CACHE")
    d = pathlib.Path(p) if p else pathlib.Path.home() / ".cache" / "synthbench"
    d.mkdir(parents=True, exist_ok=True)
    return d


class _SAscore:
    """RDKit's bundled implementation of Ertl and Schuffenhauer's score. No download needed."""
    name, direction, citation = "sascore", DIRECTION["sascore"], CITATION["sascore"]

    def __init__(self):
        import sys
        from rdkit.Chem import RDConfig
        sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
        import sascorer  # noqa: F401  (provided by the rdkit contrib directory)
        self._f = sascorer.calculateScore

    def score(self, smiles):
        from rdkit import Chem
        m = Chem.MolFromSmiles(str(smiles))
        return None if m is None else float(self._f(m))


class _Missing:
    """Stands in for a scorer that could not be loaded, carrying the reason."""

    def __init__(self, name, reason):
        self.name, self.reason = name, reason
        self.direction, self.citation = DIRECTION.get(name), CITATION.get(name)

    def score(self, smiles):
        raise RuntimeError(f"{self.name} is unavailable: {self.reason}")

    def __repr__(self):
        return f"<unavailable {self.name}: {self.reason}>"


def _missing_files(name):
    return [dst for *_ , dst in WEIGHTS.get(name, []) if not (cache_dir() / dst).exists()]


def load(name):
    """Return a scorer with a .score(smiles) method, or a _Missing placeholder explaining why not."""
    name = name.lower()
    if name not in DIRECTION:
        raise KeyError(f"unknown scorer {name!r}; known: {sorted(DIRECTION)}")
    if name == "sascore":
        try:
            return _SAscore()
        except Exception as exc:
            return _Missing(name, f"rdkit contrib SA_Score not importable ({type(exc).__name__})")
    gone = _missing_files(name)
    if gone:
        return _Missing(name, f"weights not fetched ({len(gone)} file(s)); "
                              f"run synthbench.scorers.fetch_weights({name!r})")
    try:
        return _load_vendored(name)
    except Exception as exc:
        return _Missing(name, f"{type(exc).__name__}: {str(exc)[:160]}")


def _load_vendored(name):
    """Build the scorer from cached weights. Wrapper code is ours; weights are the authors'."""
    from .vendor_impl import build
    return build(name, cache_dir())


def available() -> dict:
    """Map each scorer to 'ready' or the reason it is not, without raising."""
    out = {}
    for n in DIRECTION:
        s = load(n)
        out[n] = s.reason if isinstance(s, _Missing) else "ready"
    return out


def fetch_weights(only=None, token=None, quiet=False):
    """Download the authors' published weights into the local cache.

    Files are resolved through the GitHub API, so a moved default branch fails loudly rather
    than fetching the wrong blob. Set GITHUB_TOKEN to raise the anonymous rate limit.
    """
    import base64
    import json
    import urllib.request

    token = token or os.environ.get("GITHUB_TOKEN")
    hdr = {"Accept": "application/vnd.github+json", "User-Agent": "synthbench"}
    if token:
        hdr["Authorization"] = f"Bearer {token}"

    def api(path):
        req = urllib.request.Request("https://api.github.com" + path, headers=hdr)
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.loads(r.read().decode())

    wanted = [only] if isinstance(only, str) else (only or list(WEIGHTS))
    report = {}
    for name in wanted:
        for owner, repo, src, dst in WEIGHTS.get(name, []):
            target = cache_dir() / dst
            if target.exists():
                report[dst] = "already present"
                continue
            try:
                tree = api(f"/repos/{owner}/{repo}/git/trees/HEAD?recursive=1")
                sha = next(e["sha"] for e in tree["tree"] if e["path"] == src)
                blob = api(f"/repos/{owner}/{repo}/git/blobs/{sha}")
                raw = (base64.b64decode(blob["content"]) if blob.get("encoding") == "base64"
                       else blob["content"].encode())
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(raw)
                report[dst] = f"fetched {len(raw) / 1e6:.1f} MB"
            except Exception as exc:
                report[dst] = f"FAILED {type(exc).__name__}: {str(exc)[:120]}"
            if not quiet:
                print(f"  {report[dst]:28s} {dst}")
    return report
