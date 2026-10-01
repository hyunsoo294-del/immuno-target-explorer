"""Alias resolution. Whole tokens only; a prefix is not a gene."""

from __future__ import annotations

import unittest

from immunoscore.src.data.gene_resolver import resolve
from immunoscore.src.genes import gene_identity

_PROTEASOME = {f"PSMA{i}" for i in range(1, 9)}


class AliasResolutionTests(unittest.TestCase):
    def test_psma_resolves_to_folh1(self):
        result = resolve("PSMA")
        self.assertEqual(result.symbol, "FOLH1")
        self.assertEqual(result.route, "curated")
        self.assertNotIn(result.symbol, _PROTEASOME)
        identity = gene_identity("PSMA")
        self.assertEqual(identity["gene_symbol"], "FOLH1")
        self.assertEqual(identity["gene_id"], "ENSG00000086205")
        self.assertEqual(identity["route_label"], "별칭 해석(큐레이션)")
        self.assertFalse(identity["refused"])

    def test_psm_does_not_resolve_to_a_proteasome_gene(self):
        # No curated row and no substring. HGNC lists the whole token PSM on FOLH1.
        result = resolve("PSM")
        self.assertNotIn(result.symbol, _PROTEASOME)
        self.assertEqual(result.symbol, "FOLH1")
        self.assertEqual(result.route, "hgnc")

    def test_prefix_without_an_alias_token_is_no_match(self):
        result = resolve(
            "PSM",
            symbols=("PSMA1", "PSMA2", "PSMA8", "FOLH1"),
            aliases={},
            hgnc_index={},
        )
        self.assertIsNone(result.symbol)
        self.assertFalse(result.refused)
        self.assertIn("사용하지 않았습니다", result.failure_message())
        for name in _PROTEASOME:
            self.assertNotEqual(result.symbol, name)

    def test_gd2_refuses_to_score(self):
        result = resolve("GD2")
        self.assertTrue(result.refused)
        self.assertIsNone(result.symbol)
        self.assertEqual(result.antigen_class, "glycolipid")
        self.assertEqual(result.biosynthesis_genes, ["ST8SIA1", "B4GALNT1"])
        self.assertIn("Ganglioside", result.warning)
        identity = gene_identity("GD2")
        self.assertTrue(identity["refused"])
        self.assertEqual(identity["gene_symbol"], "")
        self.assertIn("Ganglioside", identity["warning"])

    def test_isoform_keeps_the_dot_and_flags_the_parent(self):
        result = resolve("CLDN18.2")
        self.assertEqual(result.symbol, "CLDN18")
        self.assertEqual(result.antigen_class, "isoform_variant")
        self.assertEqual(result.flags, ["ISOFORM_UNRESOLVED"])
        self.assertEqual(result.confidence_override, "LOW")
        self.assertIn("18.1", result.warning)
        stripped = resolve(
            "CLDN182",
            symbols=("CLDN18",),
            aliases={},
            hgnc_index={},
        )
        self.assertIsNone(stripped.symbol)

    def test_exact_column_wins_over_the_alias_table(self):
        result = resolve("STEAP1")
        self.assertEqual(result.symbol, "STEAP1")
        self.assertEqual(result.route, "exact")

    def test_existing_protein_aliases_still_resolve(self):
        self.assertEqual(resolve("HER2").symbol, "ERBB2")
        self.assertEqual(resolve("HER-2").symbol, "ERBB2")
        self.assertEqual(resolve("PD-L1").symbol, "CD274")
        self.assertEqual(resolve("CEA").symbol, "CEACAM5")
        self.assertEqual(resolve("TROP2").symbol, "TACSTD2")


if __name__ == "__main__":
    unittest.main()
