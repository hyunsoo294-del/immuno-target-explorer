"""Factor modules. Each one is a heuristic, not a validated predictor."""

from immunoscore.src.factors.f1_expression import score as score_f1
from immunoscore.src.factors.f2_heterogeneity import score as score_f2
from immunoscore.src.factors.f3_internalization import score as score_f3
from immunoscore.src.factors.f4_checkpoint import score as score_f4
from immunoscore.src.factors.f5_adhesion import score as score_f5
from immunoscore.src.factors.f6_glycocalyx import score as score_f6
from immunoscore.src.factors.f7_proximity import score as score_f7

FACTORS = [
    ("F1_expression", score_f1),
    ("F2_heterogeneity", score_f2),
    ("F3_internalization", score_f3),
    ("F4_checkpoint", score_f4),
    ("F5_adhesion", score_f5),
    ("F6_glycocalyx", score_f6),
    ("F7_epitope_proximity", score_f7),
]
