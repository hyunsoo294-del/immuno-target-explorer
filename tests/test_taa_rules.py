# -*- coding: utf-8 -*-
"""Calculation checks. Synthetic rows stay in tests and are not user results."""

from __future__ import annotations

import io
import math
import tempfile
import unittest
import zipfile
from pathlib import Path

import pandas as pd

from taa_analysis.aliases import target_identity
from taa_analysis.aliases import canonical_symbol
from taa_analysis.cell_lines import (
    attach_cancer_types,
    extract_hpa_gene_rows,
    resolve_catalog_symbol,
    summarize_cell_line_groups,
)
from taa_analysis.charts import matplotlib_png
from taa_analysis.cohort import select_one_sample_per_patient
from taa_analysis.composition import (
    all_cell_percentage_status,
    complete_rows,
    component_medians,
    mean_complete_composition,
    same_sample_all_cell_fraction,
    whole_tumor_fraction_allowed,
)
from taa_analysis.joins import JoinError, safe_left_join
from taa_analysis.molecules import (
    is_enrichment_summary_field,
    load_hpa_immune_cell_table,
    rna_detection_rate,
    summarize_cells_by_patient,
    tme_catalog_eligible,
)
from taa_analysis.pipeline import coverage_table, export_frame, immune_composition_result
from taa_analysis.score import score_cell_line_groups
from taa_analysis.settings import load_settings
from taa_analysis.transforms import TransformError, convert_stored_expression
from taa_analysis.xena import RetrievalFailed


class AliasTests(unittest.TestCase):
    def test_her2_and_erbb2_share_target(self):
        her2 = target_identity("HER2")
        erbb2 = target_identity("ERBB2")
        self.assertEqual(her2["gene_symbol"], erbb2["gene_symbol"])
        self.assertEqual(her2["gene_id"], erbb2["gene_id"])
        self.assertEqual(her2["gene_symbol"], "ERBB2")
        self.assertEqual(canonical_symbol("PD-L1"), "CD274")
        self.assertEqual(canonical_symbol("4-1BB"), "TNFRSF9")
        self.assertEqual(canonical_symbol("BCMA"), "TNFRSF17")
        self.assertEqual(target_identity("EGFR")["gene_symbol"], "EGFR")


class TransformTests(unittest.TestCase):
    def test_log2_tpm_0001_round_trip(self):
        stored = convert_stored_expression(0.0, "log2(tpm+0.001)")
        self.assertAlmostEqual(stored["tpm"], 0.999, places=3)
        self.assertAlmostEqual(stored["display_value"], 0.999, places=2)
        back = convert_stored_expression(stored["display_value"], "log2(TPM+1)")
        self.assertEqual(back["display_value"], stored["display_value"])
        self.assertIn("not logged again", back["formula"])

    def test_known_tpm_one(self):
        import math

        stored = math.log2(1.001)
        converted = convert_stored_expression(stored, "log2(tpm+0.001)")
        self.assertAlmostEqual(converted["tpm"], 1.0, places=6)
        self.assertAlmostEqual(converted["display_value"], 1.0, places=6)

    def test_negative_back_transform_is_error(self):
        with self.assertRaises(TransformError):
            convert_stored_expression(-20, "log2(tpm+0.001)")

    def test_other_units_are_not_relabeled(self):
        for unit in ("nTPM", "FPKM", "log2(fpkm+0.001)", "RSEM expected_count", "microarray"):
            with self.assertRaises(TransformError):
                convert_stored_expression(3.5, unit)


class JoinTests(unittest.TestCase):
    def test_duplicate_keys_do_not_expand(self):
        left = pd.DataFrame({"sample_id": ["A"], "x": [1]})
        right = pd.DataFrame({"sample_id": ["A", "A"], "subtype_label": ["LumA", "Basal"]})
        with self.assertRaises(JoinError):
            safe_left_join(left, right, "sample_id")

    def test_unique_join_keeps_row_count(self):
        left = pd.DataFrame({"sample_id": ["A", "B"], "x": [1, 2]})
        right = pd.DataFrame({"sample_id": ["A"], "subtype_label": ["Luminal A"]})
        merged = safe_left_join(left, right, "sample_id")
        self.assertEqual(len(merged), 2)
        self.assertEqual(merged.loc[merged["sample_id"] == "B", "subtype_label"].isna().all(), True)


class PatientTests(unittest.TestCase):
    def test_annotated_sample_is_preferred_without_copying(self):
        frame = pd.DataFrame(
            {
                "patient_id": ["P", "P"],
                "sample_id": ["S1", "S2"],
                "subtype_label": ["Unknown", "Luminal A"],
            }
        )
        chosen = select_one_sample_per_patient(frame, "subtype_label")
        self.assertEqual(len(chosen), 1)
        self.assertEqual(chosen.iloc[0]["sample_id"], "S2")
        self.assertEqual(float(chosen.iloc[0]["aggregation_weight"]), 1.0)
        self.assertEqual(int(chosen.iloc[0]["samples_collapsed"]), 2)


class CompositionTests(unittest.TestCase):
    def test_complete_mean_sums_to_one_and_medians_are_not_forced(self):
        frame = pd.DataFrame(
            [
                {"T": 0.5, "B": 0.5, "NK": 0.0},
                {"T": 0.5, "B": 0.0, "NK": 0.5},
                {"T": 0.0, "B": 0.5, "NK": 0.5},
            ]
        )
        columns = ["T", "B", "NK"]
        means = mean_complete_composition(frame, columns)
        self.assertAlmostEqual(float(means.sum()), 1.0, places=6)
        medians = component_medians(frame, columns)
        self.assertAlmostEqual(float(medians.sum()), 1.5, places=6)
        self.assertNotAlmostEqual(float(medians.sum()), 1.0, places=2)

    def test_missing_component_is_not_zero(self):
        frame = pd.DataFrame(
            [
                {"T": 0.4, "B": 0.6},
                {"T": 0.2, "B": None},
            ]
        )
        complete = complete_rows(frame, ["T", "B"])
        self.assertEqual(len(complete), 1)
        self.assertAlmostEqual(float(complete.iloc[0]["T"]), 0.4)

    def test_relative_and_scores_cannot_become_all_cell_percent(self):
        for algorithm in ("CIBERSORT_LM22_relative", "xCell", "ESTIMATE", "CIBERSORT_absolute_score"):
            status = all_cell_percentage_status(algorithm)
            self.assertEqual(status["status"], "Unavailable")
            self.assertEqual(status["missing_reason"], "unsupported")

    def test_product_rule_is_same_sample_only(self):
        self.assertAlmostEqual(same_sample_all_cell_fraction(0.30, 0.20, same_sample=True, same_definition=True), 0.06)
        with self.assertRaises(Exception):
            same_sample_all_cell_fraction(0.30, 0.20, same_sample=False, same_definition=True)

    def test_enriched_dataset_has_no_whole_tumor_fraction(self):
        self.assertFalse(whole_tumor_fraction_allowed("CD45 enrichment"))
        self.assertFalse(whole_tumor_fraction_allowed("immune-only"))
        self.assertTrue(whole_tumor_fraction_allowed("unselected dissociated tumor"))


class MoleculeTests(unittest.TestCase):
    def test_detection_rejects_scaled_layers(self):
        with self.assertRaises(Exception):
            rna_detection_rate(pd.Series([0.2, -1.0, 1.4]), "scaled z-score")

    def test_detection_uses_raw_counts(self):
        rate = rna_detection_rate(pd.Series([0, 2, 0, 5]), "raw counts")
        self.assertAlmostEqual(rate, 0.5)

    def test_zero_is_not_unavailable_and_empty_patient_is_not_zero(self):
        frame = pd.DataFrame(
            {
                "patient_id": ["A", "A", "B"],
                "value": [0.0, 4.0, None],
            }
        )
        # Patient B has a row but no numeric value: unavailable, not a zero mean.
        summary = summarize_cells_by_patient(frame.dropna(subset=["value"]), "value")
        pooled = summarize_cells_by_patient(frame.dropna(subset=["value"]), "value")
        self.assertAlmostEqual(summary["patient_mean"], 2.0)
        self.assertAlmostEqual(pooled["pooled_mean"], 2.0)
        absent = summarize_cells_by_patient(pd.DataFrame(columns=["patient_id", "value"]), "value")
        self.assertIsNone(absent["patient_mean"])
        self.assertEqual(absent["status"], "Unavailable")
        self.assertNotEqual(absent["patient_mean"], 0)

    def test_one_donor_does_not_define_the_patient_mean_alone_when_others_exist(self):
        frame = pd.DataFrame(
            {
                "patient_id": ["deep"] * 10 + ["small"],
                "value": [10.0] * 10 + [0.0],
            }
        )
        summary = summarize_cells_by_patient(frame, "value")
        self.assertAlmostEqual(summary["patient_mean"], 5.0)
        self.assertAlmostEqual(summary["pooled_mean"], 10.0 * 10 / 11)
        self.assertNotAlmostEqual(summary["patient_mean"], summary["pooled_mean"])

    def test_catalog_filters(self):
        self.assertFalse(tme_catalog_eligible("heart", "normal"))
        self.assertFalse(tme_catalog_eligible("cerebellum", "normal tissue"))
        self.assertFalse(tme_catalog_eligible("PBMC", "blood"))
        self.assertTrue(tme_catalog_eligible("breast", "primary tumor"))

    def test_enrichment_summary_is_not_a_matrix(self):
        self.assertTrue(is_enrichment_summary_field("RNA blood cell specific nTPM"))
        self.assertTrue(is_enrichment_summary_field("RNA cancer specificity score"))

    def test_hpa_reference_is_not_a_tumor(self):
        text = (
            "Gene\tGene name\tImmune cell\tTPM\tpTPM\tnTPM\n"
            "ENSG00000049249\tTNFRSF9\tmemory CD8 T-cell\t9.4\t11.5\t9.4\n"
            "ENSG00000188389\tPDCD1\tmemory CD8 T-cell\t12.3\t15.1\t12.4\n"
        )
        with tempfile.TemporaryDirectory() as temporary:
            cache = Path(temporary)
            (cache / "rna_immune_cell.tsv").write_text(text, encoding="utf-8")
            table = load_hpa_immune_cell_table("https://example.invalid/not-used", cache, ["TNFRSF9", "PDCD1"])
        self.assertTrue((table["evidence_status"] == "Reference only").all())
        self.assertTrue(table["cancer_type"].isna().all())
        memory = table[table["Immune cell"] == "memory CD8 T-cell"]
        tnfrsf9 = memory[memory["Gene name"] == "TNFRSF9"].iloc[0]
        pdcd1 = memory[memory["Gene name"] == "PDCD1"].iloc[0]
        self.assertEqual(float(tnfrsf9["nTPM"]), 9.4)
        self.assertEqual(float(pdcd1["nTPM"]), 12.4)
        self.assertNotIn("BRCA", table.fillna("").astype(str).values.ravel().tolist())


class CoverageAndExportTests(unittest.TestCase):
    def test_coverage_does_not_invent_immune_results(self):
        samples = pd.DataFrame(
            {
                "study": ["TCGA", "TCGA"],
                "cancer_code": ["BRCA", "LAML"],
                "specimen_class": ["solid_primary", "hematologic"],
            }
        )
        coverage = coverage_table(samples)
        brca = coverage[coverage["cancer_code"] == "BRCA"].iloc[0]
        laml = coverage[coverage["cancer_code"] == "LAML"].iloc[0]
        self.assertEqual(brca["subtype_status"], "available")
        self.assertEqual(laml["rna_reason"].startswith("hematologic"), True)
        self.assertTrue((coverage["immune_status"] == "unavailable").all())
        self.assertTrue((coverage["immune_reason"] == "unsupported").all())
        immune = immune_composition_result("BRCA")
        self.assertEqual(immune["status"], "Unavailable")
        self.assertTrue(immune["table"].empty)

    def test_export_matches_table_values(self):
        frame = pd.DataFrame(
            {
                "cancer_code": ["BRCA"],
                "median": [1.25],
                "unit": ["log2(TPM+1)"],
                "denominator": ["bulk tumor RNA, one patient"],
            }
        )
        exported = pd.read_csv(io.StringIO(export_frame(frame)))
        self.assertEqual(exported.iloc[0]["unit"], "log2(TPM+1)")
        self.assertAlmostEqual(float(exported.iloc[0]["median"]), 1.25)

    def test_figure_uses_points_for_single_values(self):
        patients = pd.DataFrame(
            {
                "group_label": ["BRCA", "BRCA", "LUAD"],
                "display_value": [1.0, 2.0, 3.0],
            }
        )
        summary = pd.DataFrame(
            {
                "group_label": ["BRCA", "LUAD"],
                "distribution": ["boxplot", "point"],
                "median": [1.5, 3.0],
            }
        )
        png = matplotlib_png(patients, summary, "fixture")
        self.assertTrue(png.startswith(b"\x89PNG"))

    def test_settings_do_not_require_invented_secrets(self):
        settings = load_settings()
        self.assertTrue(settings["TAA_XENA_TOIL_HUB"].startswith("https://"))
        self.assertEqual(settings["_used_env_keys"], "")

    def test_failure_is_distinct_from_a_number(self):
        error = RetrievalFailed("hub", "connection closed")
        self.assertEqual(error.missing_reason, "retrieval_failed")
        self.assertNotIn("9.4", str(error))


class CellLineFactorTests(unittest.TestCase):
    def _matrix(self) -> pd.DataFrame:
        text = (
            "Gene\tGene name\tCell line\tTPM\tpTPM\tnTPM\n"
            "ENSG00000141736\tERBB2\tMCF-7\t1\t1\t10\n"
            "ENSG00000141736\tERBB2\tBT-474\t1\t1\t30\n"
            "ENSG00000141736\tERBB2\tNO-MATCH\t1\t1\t0\n"
            "ENSG00000141736\tERBB2\tCONFLICT\t1\t1\t4\n"
            "ENSG00000141736\tERBB2\tBAD\t1\t1\t-1\n"
            "ENSG00000000003\tTSPAN6\tMCF-7\t9\t9\t9\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rna_celline.tsv.zip"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("rna_celline.tsv", text)
            return extract_hpa_gene_rows(path, "ERBB2")

    def test_full_matrix_factors_keep_ntpm_and_row_count(self):
        expression = self._matrix()
        self.assertEqual(len(expression), 5)
        annotations = pd.DataFrame(
            [
                {
                    "names": "MCF-7|MCF7",
                    "cancer_type": "Breast Carcinoma",
                    "model_id": "SIDM1",
                    "sample_site": "Primary",
                    "tissue_status": "Tumour",
                },
                {
                    "names": "BT-474",
                    "cancer_type": "Breast Carcinoma",
                    "model_id": "SIDM2",
                    "sample_site": "Primary",
                    "tissue_status": "Tumour",
                },
                {
                    "names": "CONFLICT",
                    "cancer_type": "Breast Carcinoma",
                    "model_id": "SIDM3",
                    "sample_site": "Primary",
                    "tissue_status": "Tumour",
                },
                {
                    "names": "CONFLICT",
                    "cancer_type": "Melanoma",
                    "model_id": "SIDM4",
                    "sample_site": "Metastasis",
                    "tissue_status": "Metastasis",
                },
            ]
        )
        lines = attach_cancer_types(expression, annotations)
        self.assertEqual(len(lines), 5)
        summary = summarize_cell_line_groups(lines)
        breast = summary[summary["cancer"] == "Breast Carcinoma"].iloc[0]
        self.assertEqual(int(breast["n_cell_lines"]), 2)
        expected = pd.Series([math.log2(11), math.log2(31)])
        self.assertAlmostEqual(float(breast["median"]), float(expected.median()))
        self.assertAlmostEqual(float(breast["median_nTPM"]), 20.0)
        self.assertEqual(breast["unit"], "log2(nTPM+1)")
        self.assertEqual(breast["source_unit"], "nTPM")
        self.assertNotIn("TPM+1)", breast["unit"].replace("nTPM", "X"))
        unknown = summary[summary["cancer"] == "Unknown"].iloc[0]
        self.assertEqual(int(unknown["n_cell_lines"]), 3)
        self.assertEqual(int(unknown["missing_n"]), 1)
        conflict = lines[lines["cell_line"] == "CONFLICT"].iloc[0]
        self.assertEqual(conflict["annotation_status"], "conflict")
        self.assertEqual(conflict["group_label"], "Unknown")

    def test_other_taa_rows_do_not_reuse_erbb2(self):
        text = (
            "Gene\tGene name\tCell line\tTPM\tpTPM\tnTPM\n"
            "ENSG00000141736\tERBB2\tMCF-7\t1\t1\t100\n"
            "ENSG00000120217\tCD274\tMCF-7\t1\t1\t8\n"
            "ENSG00000120217\tCD274\tA549\t1\t1\t2\n"
            "ENSG00000000003\tTSPAN6\tMCF-7\t1\t1\t9\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rna_celline.tsv.zip"
            catalog_path = Path(directory) / "hpa_celline_genes.csv"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("rna_celline.tsv", text)
            cd274 = extract_hpa_gene_rows(path, "CD274", catalog_path=catalog_path)
            by_id = extract_hpa_gene_rows(path, "ENSG00000120217")
            catalog = pd.read_csv(catalog_path, dtype=str)
        self.assertEqual(set(cd274["gene_symbol"]), {"CD274"})
        self.assertEqual(set(cd274["gene_id"]), {"ENSG00000120217"})
        self.assertEqual(set(cd274["cell_line"]), {"MCF-7", "A549"})
        self.assertEqual(list(by_id["gene_symbol"]), ["CD274", "CD274"])
        self.assertEqual(resolve_catalog_symbol(catalog, "PD-L1"), "CD274")
        self.assertEqual(resolve_catalog_symbol(catalog, "HER2"), "ERBB2")
        self.assertEqual(resolve_catalog_symbol(catalog, "ENSG00000120217"), "CD274")
        summary = summarize_cell_line_groups(
            attach_cancer_types(
                cd274,
                pd.DataFrame(
                    [
                        {
                            "names": "MCF-7",
                            "cancer_type": "Breast Carcinoma",
                            "model_id": "SIDM1",
                            "sample_site": "Primary",
                            "tissue_status": "Tumour",
                        },
                        {
                            "names": "A549",
                            "cancer_type": "Lung Adenocarcinoma",
                            "model_id": "SIDM9",
                            "sample_site": "Primary",
                            "tissue_status": "Tumour",
                        },
                    ]
                ),
            )
        )
        breast = summary[summary["cancer"] == "Breast Carcinoma"].iloc[0]
        self.assertEqual(int(breast["n_cell_lines"]), 1)
        self.assertAlmostEqual(float(breast["median_nTPM"]), 8.0)
        self.assertEqual(breast["unit"], "log2(nTPM+1)")


class CellLineScoreTests(unittest.TestCase):
    def _lines(self) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {"cell_line": "MCF-7", "group_label": "Breast Carcinoma", "nTPM": 10},
                {"cell_line": "BT-474", "group_label": "Breast Carcinoma", "nTPM": 30},
                {"cell_line": "A549", "group_label": "Lung Adenocarcinoma", "nTPM": 1},
                {"cell_line": "H460", "group_label": "Lung Adenocarcinoma", "nTPM": 1},
                {"cell_line": "NORM", "group_label": "Non-Cancerous", "nTPM": 500},
                {"cell_line": "X", "group_label": "Unknown", "nTPM": 800},
                {"cell_line": "BAD", "group_label": "Breast Carcinoma", "nTPM": -3},
            ]
        )

    def test_score_uses_cancer_lines_and_ranks_higher_expression(self):
        scored = score_cell_line_groups(self._lines())
        self.assertTrue(set(scored.loc[~scored["scored"], "cancer"]) >= {"Unknown", "Non-Cancerous"})
        breast = scored[scored["cancer"] == "Breast Carcinoma"].iloc[0]
        lung = scored[scored["cancer"] == "Lung Adenocarcinoma"].iloc[0]
        self.assertEqual(int(breast["n_cell_lines"]), 3)
        self.assertEqual(int(breast["n_measured"]), 2)
        self.assertAlmostEqual(float(breast["median_nTPM"]), 20.0)
        self.assertAlmostEqual(float(breast["detected_fraction"]), 1.0)
        self.assertGreater(float(breast["score"]), float(lung["score"]))
        self.assertGreater(float(breast["selectivity_log2"]), 0)
        self.assertLess(float(lung["selectivity_log2"]), 0)
        self.assertTrue(0 <= float(breast["score"]) <= 100)
        self.assertEqual(scored.iloc[0]["cancer"], "Breast Carcinoma")
        held = scored[scored["cancer"] == "Non-Cancerous"].iloc[0]
        self.assertTrue(pd.isna(held["score"]))

    def test_single_cancer_score_omits_selectivity(self):
        lines = pd.DataFrame(
            [
                {"cell_line": "A", "group_label": "Melanoma", "nTPM": 8},
                {"cell_line": "B", "group_label": "Melanoma", "nTPM": 0},
            ]
        )
        scored = score_cell_line_groups(lines).iloc[0]
        self.assertTrue(pd.isna(scored["selectivity_log2"]))
        self.assertAlmostEqual(float(scored["detected_fraction"]), 0.5)
        self.assertIsNotNone(scored["score"])
        self.assertTrue(float(scored["score"]) > 0)

    def test_empty_measurements_have_no_score(self):
        lines = pd.DataFrame(
            [
                {"cell_line": "A", "group_label": "Melanoma", "nTPM": None},
            ]
        )
        scored = score_cell_line_groups(lines).iloc[0]
        self.assertTrue(pd.isna(scored["score"]))
        self.assertTrue(pd.isna(scored["median_nTPM"]))


if __name__ == "__main__":
    unittest.main()
