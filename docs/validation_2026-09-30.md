# Validation 2026-09-30

Predicted score versus dose-response Emax. This is a check, not a claim of validity.

## HER2 × 4-1BB × Jurkat_NFkB

**n = 5**

n is below 8, so no p-value is reported.

- Spearman vs log10(Emax): 1.000
- Pearson vs log10(Emax): 0.967
- Rank agreement: 5/5
- Predicted order: Calu-3, HCC1954, NCI-N87, SK-BR-3, MDA-MB-231
- Measured order: Calu-3, HCC1954, NCI-N87, SK-BR-3, MDA-MB-231
- Predicted dynamic range: 14.41x
- Measured dynamic range: 27.78x

| cell line | score | Emax |
|---|---:|---:|
| Calu-3 | 42.4 | 25000 |
| HCC1954 | 40.3 | 15000 |
| NCI-N87 | 37.9 | 8000 |
| SK-BR-3 | 35.2 | 7000 |
| MDA-MB-231 | 2.9 | 900 |

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

