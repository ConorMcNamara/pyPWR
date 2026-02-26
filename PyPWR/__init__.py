"""Methods for calculating the power of multiple statistical tests.

This library contains functionality for calculating the power, sample size,
effect size or significance level of a wide array of statistical tests.

The package provides Python equivalents to R's pwr library functions.
"""

from PyPWR.effect_size import (
    cohen_es,
    es_h,
    es_w1,
    es_w2,
)
from PyPWR.pwr_tests import (
    pwr_2p2n_test,
    pwr_2p_test,
    pwr_anova_test,
    pwr_chisq_test,
    pwr_f2_test,
    pwr_norm_test,
    pwr_p_test,
    pwr_r_test,
    pwr_t2n_test,
    pwr_t_test,
)

__version__ = "1.0.0"

__all__ = [
    # Effect size functions
    "cohen_es",
    "es_h",
    "es_w1",
    "es_w2",
    "pwr_2p2n_test",
    # Power test functions
    "pwr_2p_test",
    "pwr_anova_test",
    "pwr_chisq_test",
    "pwr_f2_test",
    "pwr_norm_test",
    "pwr_p_test",
    "pwr_r_test",
    "pwr_t2n_test",
    "pwr_t_test",
]


def __dir__() -> list[str]:
    """Return list of public methods and attributes."""
    return __all__
