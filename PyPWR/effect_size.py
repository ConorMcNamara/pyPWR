"""Effect size calculation functions for statistical power analysis."""

from collections.abc import Sequence
from math import asin, sqrt

import numpy as np
import numpy.typing as npt


def es_h(p1: float, p2: float) -> float:
    """Compute effect size h for two proportions.

    Parameters
    ----------
    p1 : float
        First proportion (0 <= p1 <= 1)
    p2 : float
        Second proportion (0 <= p2 <= 1)

    Returns
    -------
    float
        The corresponding effect size h
    """
    return 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))


def es_w1(p0: npt.ArrayLike | Sequence[float], p1: npt.ArrayLike | Sequence[float]) -> float:
    """Compute effect size w for two sets of k probabilities.

    Calculates the effect size for goodness of fit chi-squared tests comparing
    the null hypothesis probabilities P0 with alternative hypothesis probabilities P1.

    Parameters
    ----------
    p0 : array-like
        First set of k probabilities (null hypothesis)
    p1 : array-like
        Second set of k probabilities (alternative hypothesis)

    Returns
    -------
    float
        The corresponding effect size w
    """
    p0_arr, p1_arr = np.array(p0), np.array(p1)
    return float(sqrt(np.sum(np.power(p1_arr - p0_arr, 2) / p0_arr)))


def es_w2(p: Sequence[Sequence[float]] | npt.ArrayLike) -> float:
    """Compute effect size w for a two-way probability table.

    Calculates effect size for chi-squared test of association in two-way
    contingency tables using the alternative hypothesis probability table.

    Parameters
    ----------
    p : 2D array-like
        A two-way probability table

    Returns
    -------
    float
        The corresponding effect size w
    """
    p_arr = np.array(p)
    p_i = np.sum(p_arr, axis=1)[np.newaxis]
    p_j = np.sum(p_arr, axis=0)[np.newaxis]
    p0 = np.dot(np.transpose(p_i), p_j)
    return float(sqrt(np.sum(np.power(p_arr - p0, 2) / p0)))


def cohen_es(test: str = "p", size: str = "small") -> dict[str, str | float]:
    """Get conventional effect size values for statistical tests.

    Returns the conventional effect size (small, medium, large) for various
    statistical tests as defined by Cohen (1988).

    Parameters
    ----------
    test : {'p', 't', 'r', 'anov', 'anova', 'chisq', 'chisquare', 'f2'}
        The statistical test of interest
    size : {'small', 'medium', 'large'}
        The effect size category

    Returns
    -------
    dict[str, str | float]
        Dictionary containing test name, size, effect_size value, and method

    Raises
    ------
    ValueError
        If the test type is not recognized

    References
    ----------
    Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences
    (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum Associates.
    """
    test_lower, size_lower = test.casefold(), size.casefold()

    size_mappings: dict[str, dict[str, float]] = {
        "p_t": {"small": 0.2, "medium": 0.5, "large": 0.8},
        "r": {"small": 0.1, "medium": 0.3, "large": 0.5},
        "anova": {"small": 0.1, "medium": 0.25, "large": 0.4},
        "chisq": {"small": 0.1, "medium": 0.3, "large": 0.5},
        "f2": {"small": 0.02, "medium": 0.15, "large": 0.35},
    }

    if test_lower in ["p", "t"]:
        sizes = size_mappings["p_t"]
    elif test_lower == "r":
        sizes = size_mappings["r"]
    elif test_lower in ["anov", "anova"]:
        sizes = size_mappings["anova"]
    elif test_lower in ["chisq", "chisquare"]:
        sizes = size_mappings["chisq"]
    elif test_lower == "f2":
        sizes = size_mappings["f2"]
    else:
        raise ValueError(f"Cannot identify statistical test: {test}")

    if size_lower not in sizes:
        raise ValueError(f"Invalid size: {size}. Must be 'small', 'medium', or 'large'")

    return {
        "test": test_lower,
        "size": size_lower,
        "effect_size": sizes[size_lower],
        "method": "Conventional effect size from Cohen (1988)",
    }
