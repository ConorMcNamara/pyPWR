"""Method for calculating the power of multiple statistical tests

This library contains functionality for calculating the power,
sample size, effect size or significance level of a wide array of statistical
tests.
"""

from effect_size import (
    es_h,
    es_w1,
    es_w2,
    cohen_es,
)
from pwr_tests import (
    pwr_2p_test,
    pwr_2p2n_test,
    pwr_anova_test,
    pwr_chisq_test,
    pwr_f2_test,
    pwr_norm_test,
    pwr_p_test,
    pwr_r_test,
    pwr_t_test,
    pwr_t2n_test,
)

__all__ = [
    "es_h",
    "es_w1",
    "es_w2",
    "cohen_es",
    "pwr_2p_test",
    "pwr_2p2n_test",
    "pwr_anova_test",
    "pwr_chisq_test",
    "pwr_f2_test",
    "pwr_norm_test",
    "pwr_p_test",
    "pwr_r_test",
    "pwr_t_test",
    "pwr_t2n_test",
]


def __dir__():
    """List of methods."""
    return __all__
