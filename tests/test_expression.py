from __future__ import annotations

import sys
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from primer_build.codon_opt import optimize_codons  # noqa: E402
from primer_build.expression import (  # noqa: E402
    evaluate_expression_metrics,
    rank_expression_candidates,
)
from primer_build.io import load_rules  # noqa: E402


class ExpressionOptimizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.rules = load_rules(PROJECT_ROOT / "config" / "design_rules.json")
        cls.expression_rules = cls.rules["translation_initiation"]
        cls.upstream = "ACGT" * 5 + "AAGAAGGAGATATACAT"

    def test_metrics_include_accessibility_local_folding_and_sd(self) -> None:
        cds = "ATGAGGAGG" + "GCT" * 20
        metrics = evaluate_expression_metrics(
            cds,
            self.upstream,
            self.expression_rules,
        )
        self.assertGreaterEqual(metrics.initiation_unpaired_probability, 0.0)
        self.assertLessEqual(metrics.initiation_unpaired_probability, 1.0)
        self.assertLessEqual(metrics.five_prime_mfe_kcal_mol, 0.0)
        self.assertLessEqual(metrics.local_mfe_min_kcal_mol, 0.0)
        self.assertEqual(len(metrics.sd_like_sites), 1)
        self.assertEqual(metrics.sd_like_sites[0].motif, "AGGAGG")
        self.assertEqual(metrics.sd_like_sites[0].start_nt, 4)

    def test_codon_beam_retains_five_prime_diversity(self) -> None:
        protein = "MVTETATVQLPSPEVQQLLQAGVEAGDNVGLLLR"
        candidates = optimize_codons(
            protein,
            self.rules["codon_optimization"],
            self.rules["input"]["ncaa_codon"],
            self.expression_rules["sd_like_motifs"],
            self.rules["codon_optimization"]["sd_like_search_penalty"],
        )
        self.assertEqual(
            len(candidates),
            self.rules["codon_optimization"]["candidate_count"],
        )
        starts = {candidate.cds[:45] for candidate in candidates}
        self.assertGreater(len(starts), 1)

    def test_expression_ranking_is_deterministic(self) -> None:
        protein = "MVTETATVQLPSPEVQQLLQAGVEAGDNVGLLLR"
        candidates = optimize_codons(
            protein,
            self.rules["codon_optimization"],
            self.rules["input"]["ncaa_codon"],
            self.expression_rules["sd_like_motifs"],
            self.rules["codon_optimization"]["sd_like_search_penalty"],
        )
        first = rank_expression_candidates(
            candidates,
            self.upstream,
            self.expression_rules,
        )
        second = rank_expression_candidates(
            candidates,
            self.upstream,
            self.expression_rules,
        )
        self.assertEqual(
            [(candidate.cds, candidate.combined_score) for candidate in first],
            [(candidate.cds, candidate.combined_score) for candidate in second],
        )


if __name__ == "__main__":
    unittest.main()
