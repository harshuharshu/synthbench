# -*- coding: utf-8 -*-
"""Rescore the packaged SMILES with RAscore under xgboost 1.7, so the shipped scores are
reproducible from the shipped structures."""
import sys, json, pickle, numpy as np, pandas as pd, xgboost
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")

W = sys.argv[1]; IN = sys.argv[2]; OUT = sys.argv[3]
with open(W, "rb") as fh:
    bst = pickle.load(fh).get_booster()
assert json.loads(bst.save_config())["learner"]["objective"]["name"] == "binary:logistic"

def ecfp(m):
    fp = AllChem.GetMorganFingerprint(m, 3, useCounts=True, useFeatures=False)
    a = np.zeros((2048,), np.int32)
    for i, v in fp.GetNonzeroElements().items():
        a[i % 2048] += int(v)
    return a

df = pd.read_csv(IN)
vals = []
for s in df.smiles:
    m = Chem.MolFromSmiles(str(s))
    if m is None:
        vals.append(np.nan); continue
    vals.append(round(float(bst.predict(xgboost.DMatrix(ecfp(m).reshape(1, -1).astype(np.float32)))[0]), 4))
df["rascore"] = vals
df.to_csv(OUT, index=False)
print("rescored", len(df), "nulls", int(df.rascore.isna().sum()))
