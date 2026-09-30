# Validation 2026-09-30

Predicted score versus dose-response Emax. This is a check, not a claim of validity.

## HER2 × 4-1BB × Jurkat_NFkB

**n = 5**

n is below 8, so no p-value is reported.

- Spearman vs log10(Emax): 1.000
- Pearson vs log10(Emax): 0.987
- Rank agreement: 5/5
- Predicted order: Calu-3, HCC1954, NCI-N87, SK-BR-3, MDA-MB-231
- Measured order: Calu-3, HCC1954, NCI-N87, SK-BR-3, MDA-MB-231
- Predicted dynamic range: 14.34x
- Measured dynamic range: 27.78x

| cell line | score | predicted Emax | measured Emax | ratio |
|---|---:|---:|---:|---:|
| Calu-3 | 39.2 | 19,703 | 25000 | 0.79 |
| HCC1954 | 35.0 | 13,734 | 15000 | 0.92 |
| NCI-N87 | 32.0 | 10,561 | 8000 | 1.32 |
| SK-BR-3 | 28.8 | 7,975 | 7000 | 1.14 |
| MDA-MB-231 | 2.7 | 829 | 900 | 0.92 |

log10(Emax) = 2.8156 + 0.03776 × 점수 · n=5 · R²=0.974 · 2026-09-30

Calibration is a monotone transform of the score. It does not re-rank.

Per-factor Spearman with log10(Emax). The antigen-positive column keeps lines with F1 gate ≥ 30, which is the set where adhesion was expected to carry the rank. Undefined when the factor does not vary. No p-values at this n.

| factor | all lines | antigen-positive |
|---|---:|---:|
| F1_expression | 0.000 | -1.000 |
| F2_heterogeneity | undefined | undefined |
| F3_internalization | 0.300 | 0.200 |
| F4_checkpoint | undefined | undefined |
| F5_adhesion | 0.000 | 1.000 |
| F6_glycocalyx | -0.300 | 0.400 |
| F7_epitope_proximity | undefined | undefined |
| accessibility | 0.400 | 1.000 |

