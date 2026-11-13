# -*- coding: utf-8 -*-
"""
Wrapper code for the three published scores whose weights are fetched rather than bundled.

These are the implementations used to produce the packaged reference scores, so a user
scoring their own molecules gets numbers on the same footing as the benchmark. The
arithmetic follows each paper; only the weight files belong to the original authors.
"""
from __future__ import annotations

import gzip
import json
import math
import pickle

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def _mol(x):
    return x if isinstance(x, Chem.Mol) else Chem.MolFromSmiles(str(x))


class SCScore:
    """Coley et al. 2018. Forward pass as published: alternating weight and bias pairs,
    ReLU on all but the last layer, then 1 + (scale - 1) * sigmoid."""
    name, direction = "scscore", "lower_easier"
    FP_LEN, FP_RAD, SCALE = 1024, 2, 5.0

    def __init__(self, weight_path):
        with gzip.GzipFile(str(weight_path), "r") as fh:
            self.vars = [np.array(v) for v in json.loads(fh.read().decode("utf-8"))]

    def _fp(self, mol):
        from rdkit.Chem import AllChem
        bv = AllChem.GetMorganFingerprintAsBitVect(mol, self.FP_RAD, nBits=self.FP_LEN,
                                                   useChirality=True)
        return np.array(bv, dtype=bool).astype(np.float32)

    def score(self, x):
        mol = _mol(x)
        if mol is None:
            return None
        v = self._fp(mol)
        if v.sum() == 0:
            return None
        for i in range(0, len(self.vars), 2):
            v = np.matmul(v, self.vars[i]) + self.vars[i + 1]
            if i != len(self.vars) - 2:
                v = v * (v > 0)
        s = 1.0 / (1.0 + np.exp(-v))
        return round(float(1.0 + (self.SCALE - 1.0) * s.item()), 4)


class RAScore:
    """Thakkar et al. 2021, probability that a route exists.

    The published model.pkl was serialized by XGBoost 1.x and loads only under xgboost 1.7;
    even there the sklearn wrapper is missing attributes removed after 1.5, so predict_proba
    raises. The raw Booster is used instead. The objective is checked to be binary:logistic,
    for which Booster.predict returns P(class=1), the same quantity predict_proba would give.
    """
    name, direction = "rascore", "higher_easier"

    def __init__(self, weight_path):
        import xgboost
        self._xgb = xgboost
        with open(str(weight_path), "rb") as fh:
            clf = pickle.load(fh)
        self._bst = clf.get_booster()
        obj = json.loads(self._bst.save_config())["learner"]["objective"]["name"]
        if obj != "binary:logistic":
            raise RuntimeError(f"unexpected RAscore objective {obj!r}; "
                               "Booster.predict may not return P(route exists)")

    @staticmethod
    def _ecfp(mol):
        from rdkit.Chem import AllChem
        fp = AllChem.GetMorganFingerprint(mol, 3, useCounts=True, useFeatures=False)
        arr = np.zeros((2048,), np.int32)
        for idx, v in fp.GetNonzeroElements().items():
            arr[idx % 2048] += int(v)
        return arr

    def score(self, x):
        mol = _mol(x)
        if mol is None:
            return None
        d = self._xgb.DMatrix(self._ecfp(mol).reshape(1, -1).astype(np.float32))
        return round(float(self._bst.predict(d)[0]), 4)


class BRSAScore:
    """Chen et al. 2024, building-block and reaction aware synthetic accessibility."""
    name, direction = "brsascore", "lower_easier"

    def __init__(self, weight_path, frag_penalty=-6.0, complexity_buffer=1.0):
        with gzip.open(str(weight_path), "rb") as fh:
            self._fscores = pickle.load(fh)
        self.frag_penalty = frag_penalty
        self.max_score = 0
        self.min_score = frag_penalty - complexity_buffer

    def score(self, x):
        mol = _mol(x)
        if mol is None:
            return None
        bi = {}
        fp = rdMolDescriptors.GetMorganFingerprint(mol, 2, useChirality=True, bitInfo=bi)
        fps = fp.GetNonzeroElements()

        score1, nf = 0.0, 0
        for bit_id, vs in bi.items():
            if vs[0][1] != 2:
                continue
            fscore = self._fscores.get(bit_id, self.frag_penalty)
            if fscore < 0:
                nf += 1
                score1 += fscore
        if nf:
            score1 /= nf

        n_atoms = mol.GetNumAtoms()
        n_chiral = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        n_bridge = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        n_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        ri = mol.GetRingInfo()
        n_macro = sum(1 for r in ri.AtomRings() if len(r) > 6)
        seen = {}
        for r in ri.AtomRings():
            for at in r:
                seen[at] = seen.get(at, 0) + 1
        n_multi = sum(v - 1 for v in seen.values() if v > 1)

        score2 = (-(n_atoms ** 1.005 - n_atoms)
                  - math.log10(n_chiral + 1) - math.log10(n_spiro + 1)
                  - math.log10(n_bridge + 1)
                  - (math.log10(2) if n_macro > 0 else 0)
                  - math.log10(n_multi + 1))
        score3 = math.log(float(n_atoms) / len(fps)) * 0.5 if n_atoms > len(fps) else 0.0

        raw = score1 + score2 + score3
        s = min(max((raw - self.min_score) / (self.max_score - self.min_score), 0.0), 1.0)
        return round(10 - s * 9, 4)


_BUILDERS = {
    "scscore":   (SCScore,  "scscore/model.ckpt-10654.as_numpy.json.gz"),
    "rascore":   (RAScore,  "rascore/XGB_chembl_ecfp_counts_model.pkl"),
    "brsascore": (BRSAScore, "brsascore/pickle/BRScores_uspto_emolecules.pkl.gz"),
}


def build(name, cache):
    cls, rel = _BUILDERS[name]
    return cls(cache / rel)
