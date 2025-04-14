from math import asin, sqrt
from typing import Sequence, Union, Dict

import numpy as np
import numpy.typing as npt


def es_h(p1: float, p2: float) -> float:
    """Compute effect size h for two proportions.

    Parameters
    ----------
    p1: float
        First proportion
    p2: float
        Second proportion

    Returns
    -------
    The corresponding effect size
    """
    return 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))


def es_w1(
    p0: Union[npt.ArrayLike, Sequence], p1: Union[npt.ArrayLike, Sequence]
) -> float:
    """Compute effect size w for two sets of k probabilities P0 (null hypothesis) and P1 (alternative hypothesis)

    Parameters
    ----------
    p0: array-like
        First set of k probabilities
    p1: array-like
        Second set of k probabilities

    Returns
    -------
    The corresponding effect size w
    """
    p0, p1 = np.array(p0), np.array(p1)
    return sqrt(np.sum(np.power(p1 - p0, 2) / p0))


def es_w2(p: Union[Sequence, npt.ArrayLike]) -> float:
    """Compute effect size w for a two-way probability table corresponding to the alternative hypothesis in the chi-
    squared test of association in two-way contingency tables.

    Parameters
    ----------
    p: 2x2 array
        A two-way probability table

    Returns
    -------
    The corresponding effect size w
    """
    p_i = np.sum(p, axis=1)[np.newaxis]
    p_j = np.sum(p, axis=0)[np.newaxis]
    p0 = np.dot(np.transpose(p_i), p_j)
    return sqrt(np.sum(np.power(p - p0, 2) / p0))


def cohen_es(test: str = "p", size: str = "small") -> Dict:
    """Give the conventional effect size (small, medium, large) for the tests available in this package.

    Parameters
    ----------
    test: {'p', 't', 'r', 'anov', 'chisq', 'f2'}
        The statistical test of interest
    size: {'small', 'medium', 'large'}
        The effect size

    Returns
    -------
    The corresponding effect size
    """
    test, size = test.casefold(), size.casefold()
    if test in ["p", "t"]:
        sizes = {"small": 0.2, "medium": 0.5, "large": 0.8}
    elif test == "r":
        sizes = {"small": 0.1, "medium": 0.3, "large": 0.5}
    elif test in ["anov", "anova"]:
        sizes = {"small": 0.1, "medium": 0.25, "large": 0.4}
    elif test in ["chisq", "chisquare"]:
        sizes = {"small": 0.1, "medium": 0.3, "large": 0.5}
    elif test == "f2":
        sizes = {"small": 0.02, "medium": 0.15, "large": 0.35}
    else:
        raise ValueError("Cannot identify statistical test")
    es = sizes[size]
    return {
        "test": test,
        "size": size,
        "effect_size": es,
        "method": "Conventional effect size from Cohen (1982)",
    }
