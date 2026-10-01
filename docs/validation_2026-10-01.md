# Validation 2026-10-01

Predicted score versus dose-response Emax. This is a check, not a claim of validity.

## HER2 × 4-1BB × Jurkat_NFkB

**n = 5**

n is below 8, so no p-value is reported.

- Spearman vs log10(Emax): 1.000
- Pearson vs log10(Emax): 0.976
- Rank agreement: 5/5
- Predicted order: Calu-3, HCC1954, NCI-N87, SK-BR-3, MDA-MB-231
- Measured order: Calu-3, HCC1954, NCI-N87, SK-BR-3, MDA-MB-231
- Predicted dynamic range: 15.37x
- Measured dynamic range: 27.78x

| cell line | score | predicted Emax | measured Emax | ratio |
|---|---:|---:|---:|---:|
| Calu-3 | 36.9 | 18,820 | 25000 | 0.75 |
| HCC1954 | 32.1 | 12,143 | 15000 | 0.81 |
| NCI-N87 | 30.1 | 10,117 | 8000 | 1.26 |
| SK-BR-3 | 29.9 | 9,922 | 7000 | 1.42 |
| MDA-MB-231 | 2.4 | 824 | 900 | 0.92 |

log10(Emax) = 2.8213 + 0.03936 × 점수 · n=5 · R²=0.952 · 2026-10-01

Calibration is a monotone transform of the score. It does not re-rank.

Per-factor Spearman with log10(Emax). The antigen-positive column keeps lines with F1 gate ≥ 30, which is the set where adhesion was expected to carry the rank. Undefined when the factor does not vary. No p-values at this n.

| factor | all lines | antigen-positive |
|---|---:|---:|
| F1_expression | 0.000 | -1.000 |
| F2_heterogeneity | undefined | undefined |
| F3_internalization | 0.300 | 0.200 |
| F4_checkpoint | undefined | undefined |
| F5_adhesion | 0.000 | 1.000 |
| F6_glycocalyx | -0.400 | 0.200 |
| F7_epitope_proximity | undefined | undefined |
| accessibility | 0.400 | 1.000 |

