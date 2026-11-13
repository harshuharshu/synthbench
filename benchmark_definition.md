# synthbench: benchmark definition

A synthesizability score claims to predict whether a molecule can be made. This benchmark measures
that claim directly, using reference sets whose makeability is known from provenance rather than from
any score.

## Reference sets

Six sets, 25,905 molecules, deduplicated. The label is binary and comes from the record, not from chemistry.

| Set | n | Label | What is known | Source |
|---|---|---|---|---|
| ZINC purchasable | 3,581 | made | can be bought, therefore exists | SYBA easy set, Vorsilak 2020 |
| Enamine make-on-demand | 1,000 | made | 2 to 3 steps by construction | synformer sample |
| ChEMBL kinase actives | 3,012 | made | made, route length unknown | built here from ChEMBL |
| THQ library | 10,000 | made | one step, demonstrated experimentally | Papidocha 2026, CC-BY 4.0 |
| BRICS recombination | 5,000 | never made | no evidence of existence | generated here |
| GDB-17 cages | 3,581 | never made | no evidence of existence | SYBA hard set |

A seventh set is reported separately rather than included in the main benchmark: 57 THQ members that
were **synthesized and assayed** by their original authors. It is the only set where makeability is
confirmed by experiment rather than inferred, and at n=57 it is too small to weight alongside the rest.

**Wording that binds every use of this benchmark.** No set here is proven unmakeable. "Never made"
means no evidence of existence, nothing stronger. A true negative set would require compounds that
chemists attempted and failed to make, and no such collection is public.

## Why set sizes do not drive the result

Set sizes span 1,000 to 10,000, so an unweighted pooled AUC would largely report the THQ library.
Every quantity below is computed with sample weights that give each set equal total weight within its
class, so the benchmark measures agreement across sets rather than the size of any one of them.

## The four reported quantities

**1. Makeability AUC.** Area under the ROC curve for separating made from never-made, over all six
sets, equally weighted. This is the headline: it asks whether the score does the thing it claims.
0.5 is chance. Below 0.5 means the score is inverted on this material.

**2. Easy-negative dependence.** Makeability AUC minus the same AUC with GDB-17 removed. GDB-17 cages
are strained virtual polycyclics that every method rejects, and they are half the negative class in the
field's conventional benchmark. A large positive dependence means a score's apparent performance is
carried by rejecting one easy class and does not survive on harder material.

**3. Hard-case AUC.** The discriminating comparison: THQ molecules made in one step against BRICS
recombination products never made. Both sets are novel chemotypes absent from catalogues, so this
isolates makeability from familiarity. A score that tracks makeability scores above 0.5 here.

**4. Baseline margin.** Makeability AUC minus that of a nearest-neighbour baseline, the maximum
Tanimoto similarity (Morgan radius 2, 2048 bits) to the purchasable reference set. The baseline is
deliberately trivial. A score that cannot beat it is not earning its complexity.

Bootstrap 95% confidence intervals are computed by resampling within each set, 400 iterations, fixed
seed. These intervals describe sampling within these particular reference sets. They do **not** capture
uncertainty about which sets were chosen; that is what the leave-one-set-out ablation addresses.

## Results for the four published scores

| Scorer | Makeability AUC | 95% CI | Without easy negative | Easy-negative dependence | Hard case |
|---|---|---|---|---|---|
| SAscore | 0.691 | 0.687 to 0.696 | 0.401 | +0.290 | 0.017 |
| SCScore | 0.352 | 0.346 to 0.358 | 0.464 | −0.112 | 0.266 |
| RAscore | 0.647 | 0.642 to 0.651 | 0.491 | +0.156 | 0.016 |
| BR-SAScore | 0.706 | 0.702 to 0.710 | 0.418 | +0.287 | 0.001 |
| **Nearest-neighbour baseline** | **0.785** | **0.781 to 0.789** | **0.611** | +0.174 | 0.295 |

On the experimentally confirmed subset, 57 made-and-assayed molecules against 5,000 never made:
SAscore 0.044, SCScore 0.162, RAscore 0.026, BR-SAScore 0.011, baseline 0.226.

## What these numbers say, and what they do not

The trivial baseline outperforms all four published scores on makeability AUC, and the gap survives
removing the purchasable set that serves as its own reference (0.611 against 0.401 to 0.491).

Three of the four scores lose essentially all discrimination once GDB-17 is removed, landing between
0.401 and 0.491. Their apparent performance on the conventional benchmark is carried by one easy
negative class. SCScore is a separate case: its dependence is negative because it rates GDB-17 cages
as *easier* than real kinase drugs (median 3.29 against 4.29, lower meaning easier), which is why its
overall AUC sits below chance.

On the hard case every method fails, the baseline included. Nothing tested here distinguishes molecules
made in one step from molecules never made, once both are unfamiliar. **The baseline is not a solution
and must not be presented as one.** It is a floor: a measure of how much of a published score's
behaviour is reproduced by asking only whether a molecule looks like something already in a catalogue.

The purpose of this benchmark is diagnostic. It does not supply a better score, and a method that wins
on it has demonstrated only that it clears a floor the current generation does not.

## The baseline family, and why only one of them wins

Three deliberately trivial baselines are reported, so that the floor is not a single lucky choice.

| Baseline | Makeability AUC | 95% CI |
|---|---|---|
| Nearest neighbour: max Tanimoto to the purchasable set | 0.785 | 0.781 to 0.789 |
| Heavy-atom count | 0.310 | 0.305 to 0.314 |
| Three descriptors, logistic, leave-one-set-out | 0.119 | 0.115 to 0.124 |

Only the similarity baseline beats the published scores. Size fails, and a logistic model on heavy-atom
count, fraction sp3 and aromatic ring count fails badly: fitted on five sets it anti-predicts the sixth.
Two consequences follow, and both matter.

The benchmark is not gameable by memorising descriptor ranges. Nothing here can be won by learning that
makeable molecules are small, or flat, or ring-rich, because across these sets they are not consistently
any of those things.

And the published scores are not merely reproducing some generic trivial signal. They are specifically
reproducing *catalogue similarity*, the one trivial measure that works. That is the sharpest available
statement of what these scores have learned.

## Which set decides the answer

Each set is removed in turn and the makeability AUC recomputed. The two negative sets pull in opposite
directions, and the size of that effect is the main finding of this table.

| Scorer | All six | Drop GDB-17 | Drop BRICS | Spread across all leave-one-out runs |
|---|---|---|---|---|
| SAscore | 0.691 | 0.401 | 0.981 | 0.580 |
| SCScore | 0.352 | 0.464 | 0.240 | 0.238 |
| RAscore | 0.647 | 0.491 | 0.803 | 0.318 |
| BR-SAScore | 0.706 | 0.418 | 0.993 | 0.575 |
| Nearest-neighbour baseline | 0.785 | 0.611 | 0.959 | 0.348 |

Remove GDB-17 and the published scores fall to chance. Remove BRICS instead and SAscore reads 0.981 and
BR-SAScore 0.993, close to perfect. The same scores, the same positive sets, and a swing of more than
half a point of AUC depending only on which negative set is used.

This is the methodological point the benchmark exists to make. The field's conventional evaluation uses
computer-enumerated strained cages as its negative class, and against that negative the published scores
are near-perfect. Substitute a negative class built by recombining real fragments, molecules that are
unremarkable to look at and simply have never been made, and the same scores retain nothing. A reported
AUC for a synthesizability score is therefore a statement about its negative set at least as much as
about the score, and should not be quoted without one.

Note that SAscore and BR-SAScore, the two most sensitive to this choice at a spread of 0.580 and 0.575,
are more set-dependent than the trivial similarity baseline at 0.348.

