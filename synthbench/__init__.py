# -*- coding: utf-8 -*-
"""
synthbench: does a synthesizability score predict whether a molecule can be made?

Six reference sets whose makeability is known from provenance rather than from any score,
and four reported quantities that separate predicting makeability from recognising
catalogue chemistry.

    import synthbench as sb
    res = sb.evaluate(my_score, name="MyScore", direction="higher_easier")
    print(res.summary())
    res.compare_to_published()

The benchmark is diagnostic. It supplies no better score, and the trivial baseline it holds
published scores against is a floor rather than a solution: that baseline fails the hard case
too. See benchmark_definition.md for what each number means and what it does not.
"""
from .datasets import (PROVENANCE, PUBLISHED_SCORES, load_reference_sets, main_benchmark,
                       provenance_table, smiles)
from . import scorers
from .evaluate import Result, evaluate
from .metrics import (all_metrics, bootstrap_ci, leave_one_set_out, set_weights, weighted_auc)

__version__ = "1.0.0"

__all__ = ["evaluate", "Result", "scorers", "load_reference_sets", "main_benchmark", "smiles",
           "provenance_table", "PROVENANCE", "PUBLISHED_SCORES", "all_metrics", "weighted_auc",
           "bootstrap_ci", "leave_one_set_out", "set_weights", "__version__"]
