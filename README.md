# synthbench

Does a synthesizability score predict whether a molecule can be made, or only whether it looks
like something already in a catalogue?

synthbench answers that for any scorer, using six reference sets whose makeability is known from
provenance rather than from a score. On these sets a one-line nearest-neighbour baseline
outperforms all four published scores tested, and all four fall below chance once a single easy
negative class is removed.

```
pip install synthbench
```

## Sixty seconds

```python
import synthbench as sb

def my_score(smiles):
    return some_model.predict(smiles)          # higher means easier to make

res = sb.evaluate(my_score, name="MyScore", direction="higher_easier")
print(res.summary())
res.compare_to_published()
```

`score_fn` may take one SMILES and return a float, or take a list and return a list. An object
with a `.score` method works too. Returning `None` or raising for a molecule is fine and is
recorded as not scored.

## What it reports

| Quantity | Question it answers |
|---|---|
| Makeability AUC | Does the score separate molecules that were made from molecules that were not? |
| Easy-negative dependence | How much of that comes from rejecting one easy class of virtual cages? |
| Hard-case AUC | Can it tell one-step molecules from never-made ones when both are unfamiliar? |
| Baseline margin | Does it beat maximum Tanimoto similarity to a purchasable catalogue? |

Reference set sizes span 1,000 to 10,000, so every AUC is computed with per-set weights that give
each set equal total weight within its class. A result here is a statement about agreement across
sets, not about the largest one.

## The reference sets

| Set | n | Label | What is known | Source |
|---|---|---|---|---|
| ZINC purchasable | 3,580 | made | can be bought | SYBA easy set, Vorsilak 2020 |
| Enamine make-on-demand | 1,000 | made | 2 to 3 steps by construction | synformer sample, Gao 2024 |
| ChEMBL kinase actives | 2,968 | made | route length unknown | built from ChEMBL |
| THQ library | 10,000 | made | one step, demonstrated experimentally | Papidocha 2026, CC BY 4.0 |
| BRICS recombination | 4,776 | never made | no evidence of existence | this work |
| GDB-17 cages | 3,581 | never made | no evidence of existence | SYBA hard set |

A seventh set of 57 THQ members that were **synthesized and assayed** is reported separately. It
is the only set where makeability is confirmed by experiment; at n=57 it is too small to weight
alongside the rest.

Full provenance, citations and licences: `sb.provenance_table()`.

## Where the published scores stand

| Scorer | Makeability AUC | Without GDB-17 | Hard case | Confirmed subset |
|---|---|---|---|---|
| **Nearest-neighbour baseline** | **0.785** | **0.611** | 0.295 | 0.226 |
| BR-SAScore | 0.706 | 0.418 | 0.001 | 0.011 |
| SAscore | 0.691 | 0.401 | 0.017 | 0.044 |
| RAscore | 0.647 | 0.491 | 0.016 | 0.026 |
| SCScore | 0.352 | 0.464 | 0.266 | 0.162 |

Two readings that matter more than the ranking. Removing one negative class moves every published
score below chance, and swapping which negative set is used moves SAscore between 0.401 and 0.981.
A reported AUC for a synthesizability score is therefore a statement about its negative set at
least as much as about the score, and should not be quoted without one.

On the confirmed subset, where one side was actually synthesized and assayed, every method is
inverted: it ranks molecules that were made *below* molecules that have never been made.

## Scoring your own molecules with the published panel

Not needed for the benchmark, which ships the panel precomputed. For other molecules:

```python
from synthbench import scorers
scorers.available()            # what is ready, and why the rest is not
scorers.fetch_weights()        # downloads from the authors' own repositories
scorers.load("sascore").score("CC(=O)Oc1ccccc1C(=O)O")
```

No third-party weights are redistributed here. SAscore needs nothing beyond RDKit. RAscore's
published model was serialized by XGBoost 1.x and loads only under `xgboost<2`; where it does not
load, it reports that rather than raising.

## Limitations, stated plainly

**No set here is proven unmakeable.** "Never made" means no evidence of existence, which is
weaker. A true negative set would need compounds that chemists attempted and failed to make, and
no such collection is public, so this benchmark has none and does not claim one.

**The baseline is a floor, not a solution.** It fails the hard case too (0.295). Nothing tested
here distinguishes molecules made in one step from molecules never made, once both are unfamiliar.
That the baseline wins is a statement about the published scores, not a recommendation to use it.

**The benchmark is diagnostic.** It supplies no better score. A method that wins on it has shown
only that it clears a floor the current generation does not.

**Bootstrap intervals are conditional on these sets.** They describe sampling within the reference
sets, not uncertainty about which sets were chosen. For the latter, read
`Result.leave_one_set_out()` and the spread it reports.

## Reproducing

Every number above is regenerated by `examples/quickstart.ipynb` from the packaged data. The
metric definitions and the full ablation tables are in `benchmark_definition.md`.

## Citation

See `CITATION.cff`. MIT licensed; derived ChEMBL data under CC BY-SA, the THQ reference set reused
under CC BY 4.0 with attribution to its original authors.
