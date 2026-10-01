"""Emax calibration. The line assigns magnitude and does not re-rank."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from immunoscore.src.calibration import (
    _fit_from_cache,
    attach_calibration,
    fit_log_emax,
    predict_emax,
)
from immunoscore.src.config_loader import load_configs
from immunoscore.src.scoring import score_models


SPEC_SCORES = [39.2, 35.1, 32.0, 28.8, 2.7]
SPEC_EMAX = [25000, 15000, 8000, 7000, 900]
SPEC_LINES = ["Calu-3", "HCC1954", "NCI-N87", "SK-BR-3", "MDA-MB-231"]
LIVE_PRED = [18820, 12143, 10117, 9922, 824]
LIVE_RATIO = [0.75, 0.81, 1.26, 1.42, 0.92]
LIVE_MEASURED = [25000, 15000, 8000, 7000, 900]


def _fixture():
    load_configs.cache_clear()
    configs = load_configs()
    ids = [f"ACH-{i:06d}" for i in range(1, 5)]
    expr = pd.DataFrame({"ERBB2": [2.0, 6.0, 10.0, 4.0]}, index=ids)
    expr.index.name = "ModelID"
    models = pd.DataFrame(
        {
            "ModelID": ids,
            "CellLineName": ["LOW", "MID", "HIGH", "MISSING"],
            "StrippedCellLineName": ["LOW", "MID", "HIGH", "MISSING"],
            "cancer": ["Breast"] * 4,
            "OncotreeLineage": ["Breast"] * 4,
            "OncotreePrimaryDisease": ["Breast"] * 4,
            "GrowthPattern": ["Adherent"] * 4,
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


class LogisticTests(unittest.TestCase):
    def test_median_is_half_scale_and_p90_is_near_88(self):
        from immunoscore.src.factors.common import gene_logistic

        values = pd.Series([5.0, 7.0, 3.0])
        score = gene_logistic(values, median=5.0, k=1.0, scale=100.0)
        self.assertAlmostEqual(float(score.iloc[0]), 50.0, places=6)
        self.assertAlmostEqual(float(score.iloc[1]), 100.0 / (1.0 + np.exp(-2.0)), places=6)
        self.assertAlmostEqual(float(score.iloc[2]), 100.0 / (1.0 + np.exp(2.0)), places=6)


class CalibrationMathTests(unittest.TestCase):
    def test_spec_line_on_change1_scores(self):
        fit = fit_log_emax(SPEC_SCORES, SPEC_EMAX)
        # The prompt prints a=2.815, b=0.0377. OLS on those same 1-decimal scores is 2.8168 and 0.03768.
        self.assertAlmostEqual(fit["a"], 2.8168, places=4)
        self.assertAlmostEqual(fit["b"], 0.03768, places=5)
        self.assertAlmostEqual(fit["r2"], 0.974, places=3)
        for score in SPEC_SCORES:
            pred = predict_emax(fit, score)
            self.assertLess(pred["emax_lo"], pred["emax_pred"])
            self.assertGreater(pred["emax_hi"], pred["emax_pred"])
            self.assertFalse(pred["extrapolated"])
        outside = predict_emax(fit, 80)
        self.assertTrue(outside["extrapolated"])

    def test_refuses_below_four_points(self):
        with self.assertRaises(ValueError):
            fit_log_emax([1, 2, 3], [10, 20, 30])

    def test_unfitted_combination_leaves_emax_blank(self):
        models, cfg = _fixture()
        result = score_models(models, "ERBB2", "CD3", "PBMC", cfg)
        out, meta = attach_calibration(result, "ERBB2", "CD3", "PBMC", cfg)
        self.assertFalse(meta["ok"])
        self.assertEqual(meta["caption"], "캘리브레이션 없음 — 순위만 유효")
        self.assertTrue((out["emax_text"] == "").all())
        self.assertTrue(out["emax_pred"].isna().all())

    def test_cache_invalid_when_yaml_or_scores_change(self):
        payload = {
            "a": 1.0,
            "b": 0.1,
            "r2": 0.9,
            "n": 4,
            "dof": 2,
            "s2": 0.01,
            "x_mean": 10.0,
            "sxx": 4.0,
            "score_min": 1.0,
            "score_max": 4.0,
            "fit_date": "2026-09-30",
            "scores": {"A": 1.0, "B": 2.0},
            "yaml_hash": "abc",
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "calibration.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            self.assertIsNotNone(_fit_from_cache(path, "abc", {"A": 1.0, "B": 2.0}))
            self.assertIsNone(_fit_from_cache(path, "other", {"A": 1.0, "B": 2.0}))
            self.assertIsNone(_fit_from_cache(path, "abc", {"A": 1.0, "B": 2.2}))


class Her2CalibrationTests(unittest.TestCase):
    def test_live_panel_matches_spec_predictions(self):
        parquet = Path(__file__).resolve().parents[1] / "data" / "processed" / "expression_log2tpm.parquet"
        if not parquet.exists():
            self.skipTest("DepMap parquet is not in this checkout")
        from immunoscore.src.data.abc import measured_abc_table
        from immunoscore.src.genes import gene_identity
        from immunoscore.src.scoring import build_cfg

        measured_abc_table.cache_clear()
        load_configs.cache_clear()
        identity = gene_identity("HER2")
        cfg = build_cfg({"taa_entry": identity["curated_entry"]})
        chosen = cfg["models"][cfg["models"]["CellLineName"].isin(SPEC_LINES)]
        result = score_models(chosen, identity["gene_symbol"], "4-1BB", "Jurkat_NFkB", cfg)
        out, meta = attach_calibration(result, identity["gene_symbol"], "4-1BB", "Jurkat_NFkB", cfg)
        self.assertTrue(meta["ok"], meta.get("reason"))
        self.assertAlmostEqual(meta["a"], 2.8213, places=4)
        self.assertAlmostEqual(meta["b"], 0.03936, places=5)
        self.assertAlmostEqual(meta["r2"], 0.9518, places=4)
        self.assertEqual(meta["n"], 5)
        score_order = out.sort_values("score", ascending=False)["cell_line"].tolist()
        emax_order = out.sort_values("emax_pred", ascending=False)["cell_line"].tolist()
        self.assertEqual(score_order, emax_order)
        preds = []
        for name, expected, measured, ratio in zip(SPEC_LINES, LIVE_PRED, LIVE_MEASURED, LIVE_RATIO):
            pred = float(out.loc[out["cell_line"] == name, "emax_pred"].iloc[0])
            preds.append(pred)
            self.assertEqual(round(pred), expected, name)
            self.assertAlmostEqual(pred / measured, ratio, places=2, msg=name)
            self.assertFalse(bool(out.loc[out["cell_line"] == name, "emax_extrapolated"].iloc[0]))
        self.assertAlmostEqual(preds[0] / preds[1], 1.55, places=2)
        self.assertIn("2.8213", meta["caption"])
        self.assertIn("R²=0.952", meta["caption"])
        self.assertNotIn("EXTRAPOLATED", out.iloc[0]["flags"])


if __name__ == "__main__":
    unittest.main()
