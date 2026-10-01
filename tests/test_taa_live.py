# -*- coding: utf-8 -*-
"""Live public-data checks. Skips when the source cannot be reached."""

from __future__ import annotations

import unittest

from taa_analysis.pipeline import expression_bundle, healthy_blood_reference
from taa_analysis.xena import RetrievalFailed


class LiveBrcaTests(unittest.TestCase):
    def test_erbb2_brca_views_use_real_sources(self):
        try:
            her2 = expression_bundle("HER2", "solid_primary", "BRCA", "PAM50")
            erbb2 = expression_bundle("ERBB2", "solid_primary", "BRCA", "PAM50")
        except (RetrievalFailed, OSError) as exc:
            self.skipTest(f"source unavailable: {exc}")
        self.assertEqual(her2["identity"]["gene_symbol"], "ERBB2")
        self.assertEqual(her2["meta"]["input_checksum"], erbb2["meta"]["input_checksum"])
        self.assertEqual(her2["meta"]["unit_display"], "log2(TPM+1)")
        self.assertEqual(her2["meta"]["unit_source"], "log2(tpm+0.001)")
        brca = her2["summary"][her2["summary"]["cancer_code"] == "BRCA"]
        self.assertEqual(len(brca), 1)
        self.assertGreater(int(brca.iloc[0]["n_patients"]), 100)
        labels = set(her2["subtype_summary"]["subtype"])
        self.assertTrue({"Luminal A", "Luminal B", "HER2-enriched", "Basal-like", "Normal-like"} & labels)
        self.assertNotIn("BRCA", set(her2["subtype_summary"]["subtype"]))
        self.assertEqual(her2["immune"]["status"], "Unavailable")
        self.assertFalse((her2["patients"]["display_value"].dropna() < 0).any())

    def test_hpa_memory_cd8_values_are_healthy_blood(self):
        try:
            table = healthy_blood_reference(["TNFRSF9", "PDCD1"])
        except Exception as exc:
            self.skipTest(f"HPA unavailable: {exc}")
        memory = table[table["cell_type"] == "memory CD8 T-cell"]
        tnfrsf9 = float(memory.loc[memory["gene_symbol"] == "TNFRSF9", "rna_summary_nTPM"].iloc[0])
        pdcd1 = float(memory.loc[memory["gene_symbol"] == "PDCD1", "rna_summary_nTPM"].iloc[0])
        self.assertEqual(tnfrsf9, 9.4)
        self.assertEqual(pdcd1, 12.4)
        self.assertTrue((table["evidence_status"] == "Reference only").all())
        self.assertTrue(table["cancer_type"].isna().all())


if __name__ == "__main__":
    unittest.main()
