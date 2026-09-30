"""Weight sums, bounds, missing values, and effector gating."""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from immunoscore.src.config_loader import load_configs
from immunoscore.src.scoring import score_models


def _configs():
    load_configs.cache_clear()
    return load_configs()


def _fixture():
    configs = _configs()
    ids = [f"ACH-{i:06d}" for i in range(1, 5)]
    # Rising ERBB2, varied adhesion and mucin.
    expr = pd.DataFrame(
        {
            "ERBB2": [2.0, 6.0, 10.0, np.nan],
            "ICAM1": [1.0, 3.0, 8.0, 4.0],
            "CD58": [1.0, 2.0, 7.0, 3.0],
            "MUC1": [1.0, 2.0, 3.0, 8.0],
            "MUC16": [0.2, 0.4, 1.0, 9.0],
            "CD274": [1.0, 4.0, 2.0, 6.0],
            "HLA-A": [0.2, 5.0, 6.0, 4.0],
            "HLA-B": [0.2, 5.0, 6.0, 4.0],
            "HLA-C": [0.2, 5.0, 6.0, 4.0],
            "B2M": [0.1, 5.0, 6.0, 4.0],
            "CLTC": [3.0, 3.0, 3.0, 3.0],
        },
        index=ids,
    )
    expr.index.name = "ModelID"
    models = pd.DataFrame(
        {
            "ModelID": ids,
            "CellLineName": ["LOW", "MID", "HIGH", "MISSING"],
            "StrippedCellLineName": ["LOW", "MID", "HIGH", "MISSING"],
            "cancer": ["Breast", "Breast", "Gastric", "Breast"],
            "OncotreeLineage": ["Breast", "Breast", "Esophagus/Stomach", "Breast"],
            "OncotreePrimaryDisease": ["Breast", "Breast", "Gastric", "Breast"],
            "GrowthPattern": ["Adherent", "Suspension", "Adherent", "Adherent"],
        }
    )
    cfg = {
        "weights": configs["weights"],
        "effectors": configs["effectors"],
        "params": configs["params"],
        "genes": configs["genes"],
        "mucins": configs["mucins"],
        "expression": expr,
        "models": models,
        "assay_duration_h": 24,
        "abc50_base": configs["params"]["f1_expression"]["abc50_base"],
        "morphology": configs["morphology"],
        "hill_h_override": None,
        "jurkat_pd1": False,
        "taa_entry": configs["taa"]["TAAs"]["ERBB2"],
        "surfaceome": None,
    }
    return models, cfg


class WeightTests(unittest.TestCase):
    def test_every_combination_sums_to_100(self):
        combinations = _configs()["weights"]["combinations"]
        self.assertGreaterEqual(len(combinations), 6)
        for name, spec in combinations.items():
            total = sum(spec["weights"].values())
            self.assertAlmostEqual(total, 100.0, places=2, msg=name)
            self.assertNotIn("F1_expression", spec["weights"], name)
            self.assertEqual(spec.get("f1_mode"), "gate", name)

    def test_gene_set_weights_sum_to_one(self):
        genes = _configs()["genes"]
        for name in ("checkpoint", "adhesion"):
            total = sum(genes[name].values())
            self.assertAlmostEqual(total, 1.0, places=6, msg=name)


class ScoreTests(unittest.TestCase):
    def test_scores_stay_in_range_and_rank_expression(self):
        models, cfg = _fixture()
        result = score_models(models, "ERBB2", "CD3", "PBMC", cfg)
        self.assertEqual(len(result), 4)
        finite = result["score"].dropna()
        self.assertTrue(((finite >= 0) & (finite <= 100)).all())
        high = result.loc[result["cell_line"] == "HIGH", "F1_expression"].iloc[0]
        low = result.loc[result["cell_line"] == "LOW", "F1_expression"].iloc[0]
        self.assertGreater(high, low)

    def test_missing_expression_is_not_zero(self):
        models, cfg = _fixture()
        result = score_models(models, "ERBB2", "CD3", "PBMC", cfg)
        missing = result.loc[result["cell_line"] == "MISSING"].iloc[0]
        self.assertTrue(pd.isna(missing["F1_expression"]))
        self.assertNotEqual(missing["F1_expression"], 0)

    def test_zero_expression_is_a_real_zero(self):
        models, cfg = _fixture()
        cfg["expression"].loc["ACH-000001", "ERBB2"] = 0.0
        result = score_models(models.iloc[:1], "ERBB2", "CD3", "PBMC", cfg)
        self.assertEqual(result.iloc[0]["F1_expression"], 0.0)

    def test_hek_checkpoint_is_inert_and_adhesion_is_scaled(self):
        models, cfg = _fixture()
        pbmc = score_models(models, "ERBB2", "CD3", "PBMC", cfg)
        hek = score_models(models, "ERBB2", "CD40", "HEK", cfg)
        self.assertTrue((hek["F4_checkpoint_tier"] == "not_applicable").all())
        self.assertTrue(hek["F4_checkpoint"].isna().all())
        raw = pbmc.set_index("cell_line")["F5_adhesion"]
        scaled = hek.set_index("cell_line")["F5_adhesion"]
        competence = cfg["effectors"]["effectors"]["HEK"]["synapse_competence"]
        for line in ("LOW", "MID", "HIGH"):
            self.assertAlmostEqual(scaled.loc[line], raw.loc[line] * competence, places=4)
        self.assertIn("EFFECTOR_FACTOR_INERT", hek.iloc[0]["flags"])

    def test_pbmc_keeps_checkpoint_weight(self):
        models, cfg = _fixture()
        result = score_models(models, "ERBB2", "CD3", "PBMC", cfg)
        self.assertTrue((result["F4_checkpoint_tier"] == "rna_only").all())
        self.assertNotIn("EFFECTOR_FACTOR_INERT", result.iloc[0]["flags"])

    def test_suspension_pair_and_low_confidence(self):
        models, cfg = _fixture()
        result = score_models(models, "ERBB2", "CD3", "PBMC", cfg)
        mid = result.loc[result["cell_line"] == "MID"].iloc[0]
        self.assertIn("SUSPENSION_PAIR", mid["flags"])
        self.assertIn("LOW_CONFIDENCE", mid["flags"])
        self.assertEqual(mid["confidence"], "LOW")

    def test_constant_factor_is_removed_from_the_sum(self):
        import json

        models, cfg = _fixture()
        result = score_models(models, "ERBB2", "CD3", "PBMC", cfg)
        dropped = result.iloc[0]["zero_variance_factors"].split(",")
        self.assertIn("F7_epitope_proximity", dropped)
        effective = json.loads(result.iloc[0]["effective_weights"])
        self.assertNotIn("F7_epitope_proximity", effective)
        self.assertAlmostEqual(sum(effective.values()), 100.0, places=6)
        text = result.iloc[0]["flag_text"]
        self.assertIn("무변별인자(", text)
        self.assertIn("균일성", text)
        self.assertIn("에피톱", text)

    def test_a_single_line_keeps_constant_factors(self):
        models, cfg = _fixture()
        result = score_models(models.iloc[:1], "ERBB2", "CD3", "PBMC", cfg)
        self.assertEqual(result.iloc[0]["zero_variance_factors"], "")
        self.assertNotIn("ZERO_VARIANCE_FACTOR", result.iloc[0]["flags"])

    def test_cell_size_stays_out_until_diameter_and_flag_exist(self):
        from immunoscore.src.data.morphology import effective_accessibility

        configs = _configs()
        params = configs["params"]
        self.assertFalse(params["cell_size"]["use_cell_size"])
        self.assertIsNone(params["cell_size"]["d_ref_um"])
        plain, source = effective_accessibility("Calu-3", "CALU3", "Adherent", configs["morphology"], params)
        self.assertEqual(plain, 1.0)
        self.assertEqual(source, "curated")
        enabled = {"cell_size": {"use_cell_size": True, "d_ref_um": 15}}
        still, _source = effective_accessibility("Calu-3", "CALU3", "Adherent", configs["morphology"], enabled)
        self.assertEqual(still, 1.0)
        morphology = {"cell_lines": {"LOW": {"accessibility": 1.0, "mean_diameter_um": 30}}, "default": {}}
        scaled, _source = effective_accessibility("LOW", "LOW", "Adherent", morphology, enabled)
        self.assertAlmostEqual(scaled, 4.0, places=6)

    def test_accessibility_spacing_is_not_steepened(self):
        lines = _configs()["morphology"]["cell_lines"]
        self.assertAlmostEqual(lines["Calu-3"]["accessibility"], 1.00)
        self.assertAlmostEqual(lines["HCC1954"]["accessibility"], 0.80)
        self.assertAlmostEqual(lines["NCI-N87"]["accessibility"], 0.60)
        self.assertAlmostEqual(lines["SK-BR-3"]["accessibility"], 0.55)
        self.assertAlmostEqual(lines["MDA-MB-231"]["accessibility"], 0.90)

    def test_fitted_jurkat_weights_are_not_copied(self):
        combinations = _configs()["weights"]["combinations"]
        fitted = combinations["4-1BB__Jurkat_NFkB"]["weights"]
        for name, spec in combinations.items():
            if name == "4-1BB__Jurkat_NFkB":
                self.assertNotIn("not fitted", spec["fitted_on"])
                continue
            self.assertNotEqual(spec["weights"], fitted, name)
            self.assertIn("not fitted", spec["fitted_on"])


class GateTests(unittest.TestCase):
    def test_hill_gate_at_jurkat_nfkb_abc50(self):
        abc50 = 60000.0
        hill = 2.0

        def s1(abc: float) -> float:
            return 100.0 * abc**hill / (abc50**hill + abc**hill)

        self.assertAlmostEqual(s1(150000), 86.2, places=1)
        self.assertAlmostEqual(s1(500000), 98.6, places=1)
        self.assertAlmostEqual(s1(900000), 99.6, places=1)
        self.assertAlmostEqual(s1(15000), 5.9, places=1)

    def test_her2_jurkat_refit_matches_measured_order(self):
        from pathlib import Path

        parquet = Path(__file__).resolve().parents[1] / "data" / "processed" / "expression_log2tpm.parquet"
        if not parquet.exists():
            self.skipTest("DepMap parquet is not in this checkout")
        from immunoscore.src.data.abc import measured_abc_table
        from immunoscore.src.genes import gene_identity
        from immunoscore.src.scoring import build_cfg, score_models

        measured_abc_table.cache_clear()
        load_configs.cache_clear()
        identity = gene_identity("HER2")
        cfg = build_cfg({"taa_entry": identity["curated_entry"]})
        names = ["Calu-3", "HCC1954", "NCI-N87", "SK-BR-3", "MDA-MB-231"]
        chosen = cfg["models"][cfg["models"]["CellLineName"].isin(names)]
        result = score_models(chosen, identity["gene_symbol"], "4-1BB", "Jurkat_NFkB", cfg)
        expected = {
            "Calu-3": 39.2,
            "HCC1954": 35.0,
            "NCI-N87": 32.0,
            "SK-BR-3": 28.8,
            "MDA-MB-231": 2.7,
        }
        got = {}
        for name, score in expected.items():
            value = round(float(result.loc[result["cell_line"] == name, "score"].iloc[0]), 1)
            got[name] = value
            self.assertEqual(value, score, name)
        order = [name for name, _score in sorted(got.items(), key=lambda item: -item[1])]
        self.assertEqual(order, names)
        import json

        effective = json.loads(result.iloc[0]["effective_weights"])
        labels = load_configs()["weights"]["factor_labels"]
        shown = {labels[factor]: round(weight, 1) for factor, weight in effective.items()}
        self.assertEqual(shown, {"접촉": 47.6, "내재화": 28.6, "당질층": 23.8})
        self.assertIn("균일성", result.iloc[0]["flag_text"])
        self.assertIn("에피톱", result.iloc[0]["flag_text"])


if __name__ == "__main__":
    unittest.main()
