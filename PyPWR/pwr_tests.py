"""User-facing functions for statistical power calculations.

This module provides convenient wrapper functions for calculating power,
sample size, effect size, or significance level for various statistical tests.
"""

from PyPWR.power_classes import (
    pwr_2p,
    pwr_2p2n,
    pwr_anova,
    pwr_chisq,
    pwr_f2,
    pwr_norm,
    pwr_p,
    pwr_r,
    pwr_t,
    pwr_t2n,
)


def pwr_2p_test(
    h: float | None = None,
    n: int | None = None,
    sig_level: float | None = None,
    power: float | None = None,
    alternative: str = "two-sided",
    print_pretty: bool = True,
) -> dict[str, float | int | str]:
    """Compute power for two proportions test with equal sample sizes.

    Computes power of test, or determines parameters to obtain target power.
    Similar to R's pwr.2p.test function.

    Parameters
    ----------
    h : float | None, default=None
        Effect size (arcsine transformation)
    n : int | None, default=None
        Number of observations per sample
    sig_level : float | None, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power : float | None, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
        The alternative hypothesis
    print_pretty : bool, default=True
        Whether to print the results in a formatted output

    Returns
    -------
    dict[str, float | int | str]
        Dictionary containing h, n, sig_level, power, alternative, method, and note

    Raises
    ------
    ValueError
        If not exactly one of h, n, sig_level, or power is None, or if
        parameter values are out of valid range
    """
    if not any(v is None for v in [h, n, sig_level, power]):
        raise ValueError("One of h, n, sig_level or power must be None")
    if sum([v is None for v in [h, n, sig_level, power]]) > 1:
        raise ValueError("Only one of h, n, sig_level or power may be None")
    if sig_level is not None and (sig_level < 0 or sig_level > 1):
        raise ValueError("sig_level must be between 0 and 1")
    if power is not None and (power < 0 or power > 1):
        raise ValueError("power must be between 0 and 1")
    alternative = alternative.casefold()
    if alternative == "two-sided" and h is not None:
        h = abs(h)
    pwr = pwr_2p(h, n, sig_level, power, alternative).pwr_test()
    if print_pretty:
        str_print = (
            "\t"
            + pwr["method"]
            + "\n" * 2
            + "\t" * 2
            + " " * 2
            + f"h = {pwr['effect_size']}"
            + "\n"
            + "\t" * 2
            + " " * 2
            + f"n = {pwr['n']}"
            + "\n"
            + "\t"
            + " " * 2
            + f"sig_level = {pwr['sig_level']}"
            + "\n"
            + "\t"
            + " " * 6
            + f"power = {pwr['power']}"
            + "\n"
            + "\t"
            + f"alternative = {pwr['alternative']}"
            + "\n" * 2
            + f"NOTE: {pwr['note']}"
        )
        print(str_print)
    return pwr


def pwr_2p2n_test(
    h: float | None = None,
    n1: int | None = None,
    n2: int | None = None,
    sig_level: float | None = None,
    power: float | None = None,
    alternative: str = "two-sided",
    print_pretty: bool = True,
) -> dict[str, float | int | str]:
    """Compute power for two proportions test with unequal sample sizes.

    Computes power of test, or determines parameters to obtain target power.
    Similar to R's pwr.2p2n.test function.

    Parameters
    ----------
    h : float | None, default=None
        Effect size (arcsine transformation)
    n1 : int | None, default=None
        Number of observations in the first sample
    n2 : int | None, default=None
        Number of observations in the second sample
    sig_level : float | None, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power : float | None, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
        The alternative hypothesis
    print_pretty : bool, default=True
        Whether to print the results in a formatted output

    Returns
    -------
    dict[str, float | int | str]
        Dictionary containing h, n1, n2, sig_level, power, alternative, method, and note

    Raises
    ------
    ValueError
        If not exactly one of h, n1, n2, sig_level, or power is None, or if
        parameter values are out of valid range
    """
    if not any(v is None for v in [h, n1, n2, sig_level, power]):
        raise ValueError("One of h, n1, n2, sig_level or power must be None")
    if sum([v is None for v in [h, n1, n2, sig_level, power]]) > 1:
        raise ValueError("Only one of h, n1, n2, sig_level or power may be None")
    if n1 is not None and n1 < 2:
        raise ValueError("Number of observations in the first group must be at least 2")
    if n2 is not None and n2 < 2:
        raise ValueError("Number of observations in the second group must be at least 2")
    if sig_level is not None and (sig_level < 0 or sig_level > 1):
        raise ValueError("sig_level must be between 0 and 1")
    if power is not None and (power < 0 or power > 1):
        raise ValueError("power must be between 0 and 1")
    alternative = alternative.casefold()
    if alternative == "two-sided" and h is not None:
        h = abs(h)
    pwr = pwr_2p2n(h, n1, n2, sig_level, power, alternative).pwr_test()
    if print_pretty:
        str_print = (
            "\t"
            + pwr["method"]
            + "\n" * 2
            + "\t" * 2
            + " " * 2
            + f"h = {pwr['effect_size']}"
            + "\n"
            + "\t" * 2
            + " "
            + f"n1 = {pwr['n1']}"
            + "\n"
            + "\t" * 2
            + " "
            + f"n2 = {pwr['n2']}"
            + "\n"
            + "\t"
            + " " * 2
            + f"sig_level = {pwr['sig_level']}"
            + "\n"
            + "\t"
            + " " * 6
            + f"power = {pwr['power']}"
            + "\n"
            + "\t"
            + f"alternative = {pwr['alternative']}"
            + "\n" * 2
            + f"NOTE: {pwr['note']}"
        )
        print(str_print)
    return pwr


def pwr_anova_test(
    k: int | None = None,
    n: int | None = None,
    f: float | None = None,
    sig_level: float | None = None,
    power: float | None = None,
    print_pretty: bool = True,
) -> dict[str, int | float | str]:
    """Compute power for balanced one-way ANOVA.

    Computes power of test or determines parameters to obtain target power.
    Similar to R's pwr.anova.test function.

    Parameters
    ----------
    k : int | None, default=None
        Number of groups
    n : int | None, default=None
        Number of observations per group
    f : float | None, default=None
        Effect size
    sig_level : float | None, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power : float | None, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    print_pretty : bool, default=True
        Whether to print the results in a formatted output

    Returns
    -------
    dict[str, int | float | str]
        Dictionary containing k, n, effect_size (f), sig_level, power, method, and note

    Raises
    ------
    ValueError
        If not exactly one of k, n, f, sig_level, or power is None, or if
        parameter values are out of valid range
    """
    if not any(v is None for v in [k, n, f, sig_level, power]):
        raise ValueError("One of k, n, f, sig_level or power must be None")
    if sum([v is None for v in [k, n, f, sig_level, power]]) > 1:
        raise ValueError("Only one of k, n, f, sig_level or power may be None")
    if f is not None and f < 0:
        raise ValueError("f must be positive")
    if k is not None and k < 2:
        raise ValueError("Number of groups must be at least 2")
    if n is not None and n < 2:
        raise ValueError("Number of observations must be at least 2")
    if sig_level is not None and (sig_level < 0 or sig_level > 1):
        raise ValueError("sig_level must be between 0 and 1")
    if power is not None and (power < 0 or power > 1):
        raise ValueError("power must be between 0 and 1")
    pwr = pwr_anova(k, n, f, sig_level, power).pwr_test()
    if print_pretty:
        str_print = (
            "\t"
            + pwr["method"]
            + "\n" * 2
            + "\t" * 2
            + " " * 2
            + f"k = {pwr['k']}"
            + "\n"
            + "\t" * 2
            + " " * 2
            + f"n = {pwr['n']}"
            + "\n"
            + "\t" * 2
            + " " * 2
            + f"f = {pwr['effect_size']}"
            + "\n"
            + "\t"
            + " " * 2
            + f"sig_level = {pwr['sig_level']}"
            + "\n"
            + "\t"
            + " " * 6
            + f"power = {pwr['power']}"
            + "\n" * 2
            + f"NOTE: {pwr['note']}"
        )
        print(str_print)
    return pwr


def pwr_chisq_test(
    w: float | None = None,
    n: int | None = None,
    df: int = 1,
    sig_level: float | None = None,
    power: float | None = None,
    print_pretty: bool = True,
) -> dict[str, int | float | str]:
    """Compute power for chi-squared test.

    Computes power of test or determines parameters to obtain target power.
    Similar to R's pwr.chisq.test function.

    Parameters
    ----------
    w : float | None, default=None
        Effect size
    n : int | None, default=None
        Total number of observations
    df : int, default=1
        Degrees of freedom
    sig_level : float | None, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power : float | None, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    print_pretty : bool, default=True
        Whether to print the results in a formatted output

    Returns
    -------
    dict[str, int | float | str]
        Dictionary containing effect_size (w), n, df, sig_level, power, method, and note

    Raises
    ------
    ValueError
        If not exactly one of w, n, sig_level, or power is None, or if
        parameter values are out of valid range
    """
    if not any(v is None for v in [w, n, sig_level, power]):
        raise ValueError("One of w, n, sig_level or power must be None")
    if sum([v is None for v in [w, n, sig_level, power]]) > 1:
        raise ValueError("Only one of w, n, sig_level or power may be None")
    if w is not None and w < 0:
        raise ValueError("w must be positive")
    if n is not None and n < 1:
        raise ValueError("Number of observations must be at least 1")
    if sig_level is not None and (sig_level < 0 or sig_level > 1):
        raise ValueError("sig_level must be between 0 and 1")
    if power is not None and (power < 0 or power > 1):
        raise ValueError("power must be between 0 and 1")
    pwr = pwr_chisq(w, n, df, sig_level, power).pwr_test()
    if print_pretty:
        str_print = (
            "\t"
            + pwr["method"]
            + "\n" * 2
            + "\t" * 2
            + " " * 2
            + f"w = {pwr['effect_size']}"
            + "\n"
            + "\t" * 2
            + " " * 2
            + f"n = {pwr['n']}"
            + "\n"
            + "\t" * 2
            + " "
            + f"df = {pwr['df']}"
            + "\n"
            + "\t"
            + " " * 2
            + f"sig_level = {pwr['sig_level']}"
            + "\n"
            + "\t"
            + " " * 6
            + f"power = {pwr['power']}"
            + "\n" * 2
            + f"NOTE: {pwr['note']}"
        )
        print(str_print)
    return pwr


def pwr_f2_test(
    u: int | None = None,
    v: int | None = None,
    f2: float | None = None,
    sig_level: float | None = None,
    power: float | None = None,
    print_pretty: bool = True,
) -> dict[str, int | float | str]:
    """Compute power for general linear model F-test.

    Computes power of test or determines parameters to obtain target power.
    Similar to R's pwr.f2.test function.

    Parameters
    ----------
    u : int | None, default=None
        Degrees of freedom for the numerator
    v : int | None, default=None
        Degrees of freedom for the denominator
    f2 : float | None, default=None
        Effect size (f-squared)
    sig_level : float | None, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power : float | None, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    print_pretty : bool, default=True
        Whether to print the results in a formatted output

    Returns
    -------
    dict[str, int | float | str]
        Dictionary containing u, v, effect_size (f2), sig_level, power, and method

    Raises
    ------
    ValueError
        If not exactly one of u, v, f2, sig_level, or power is None, or if
        parameter values are out of valid range
    """
    if not any(x is None for x in [u, v, f2, sig_level, power]):
        raise ValueError("One of u, v, f2, sig_level or power must be None")
    if sum([x is None for x in [u, v, f2, sig_level, power]]) > 1:
        raise ValueError("Only one of u, v, f2, sig_level or power may be None")
    if f2 is not None and f2 < 0:
        raise ValueError("f2 must be positive")
    if u is not None and u < 1:
        raise ValueError("Degrees of freedom u for numerator must be at least 1")
    if v is not None and v < 1:
        raise ValueError("Degrees of freedom v for denominator must be at least 1")
    if sig_level is not None and (sig_level < 0 or sig_level > 1):
        raise ValueError("sig_level must be between 0 and 1")
    if power is not None and (power < 0 or power > 1):
        raise ValueError("power must be between 0 and 1")
    pwr = pwr_f2(u, v, f2, sig_level, power).pwr_test()
    if print_pretty:
        str_print = (
            "\t"
            + pwr["method"]
            + "\n" * 2
            + "\t" * 2
            + " " * 2
            + f"u = {pwr['u']}"
            + "\n"
            + "\t" * 2
            + " " * 2
            + f"v = {pwr['v']}"
            + "\n"
            + "\t" * 2
            + " "
            + f"f2 = {pwr['effect_size']}"
            + "\n"
            + "\t"
            + " " * 2
            + f"sig_level = {pwr['sig_level']}"
            + "\n"
            + "\t"
            + " " * 6
            + f"power = {pwr['power']}"
        )
        print(str_print)
    return pwr


def pwr_norm_test(
    d: float | None = None,
    n: int | None = None,
    sig_level: float | None = None,
    power: float | None = None,
    alternative: str = "two-sided",
    print_pretty: bool = True,
) -> dict[str, float | int | str]:
    """Compute power for normal distribution with known variance.

    Computes power of test or determines parameters to obtain target power.
    Similar to R's pwr.norm.test function.

    Parameters
    ----------
    d : float | None, default=None
        Effect size (standardized mean difference)
    n : int | None, default=None
        Number of observations
    sig_level : float | None, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power : float | None, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
        The alternative hypothesis
    print_pretty : bool, default=True
        Whether to print the results in a formatted output

    Returns
    -------
    dict[str, float | int | str]
        Dictionary containing effect_size (d), n, sig_level, power, alternative, and method

    Raises
    ------
    ValueError
        If not exactly one of d, n, sig_level, or power is None, or if
        parameter values are out of valid range
    """
    if not any(x is None for x in [d, n, sig_level, power]):
        raise ValueError("One of d, n, sig_level or power must be None")
    if sum([x is None for x in [d, n, sig_level, power]]) > 1:
        raise ValueError("Only one of d, n, sig_level or power may be None")
    if n is not None and n < 1:
        raise ValueError("Number of observations in each group must be at least 1")
    if sig_level is not None and (sig_level < 0 or sig_level > 1):
        raise ValueError("sig_level must be between 0 and 1")
    if power is not None and (power < 0 or power > 1):
        raise ValueError("power must be between 0 and 1")
    alternative = alternative.casefold()
    if alternative == "two-sided" and d is not None:
        d = abs(d)
    pwr = pwr_norm(d, n, sig_level, power, alternative).pwr_test()
    if print_pretty:
        str_print = (
            "\t"
            + pwr["method"]
            + "\n" * 2
            + "\t" * 2
            + " " * 2
            + f"d = {pwr['effect_size']}"
            + "\n"
            + "\t" * 2
            + " " * 2
            + f"n = {pwr['n']}"
            + "\n"
            + "\t"
            + " " * 2
            + f"sig_level = {pwr['sig_level']}"
            + "\n"
            + "\t"
            + " " * 6
            + f"power = {pwr['power']}"
            + "\n"
            + "\t"
            + f"alternative = {pwr['alternative']}"
        )
        print(str_print)
    return pwr


def pwr_p_test(
    h: float | None = None,
    n: int | None = None,
    sig_level: float | None = None,
    power: float | None = None,
    alternative: str = "two-sided",
    print_pretty: bool = True,
) -> dict[str, float | int | str]:
    """Compute power for one-sample proportion test.

    Computes power of test or determines parameters to obtain target power.
    Similar to R's pwr.p.test function.

    Parameters
    ----------
    h : float | None, default=None
        Effect size (arcsine transformation)
    n : int | None, default=None
        Number of observations
    sig_level : float | None, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power : float | None, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
        The alternative hypothesis
    print_pretty : bool, default=True
        Whether to print the results in a formatted output

    Returns
    -------
    dict[str, float | int | str]
        Dictionary containing effect_size (h), n, sig_level, power, alternative, and method

    Raises
    ------
    ValueError
        If not exactly one of h, n, sig_level, or power is None, or if
        parameter values are out of valid range
    """
    if not any(x is None for x in [h, n, sig_level, power]):
        raise ValueError("One of h, n, sig_level or power must be None")
    if sum([x is None for x in [h, n, sig_level, power]]) > 1:
        raise ValueError("Only one of h, n, sig_level or power may be None")
    if n is not None and n < 1:
        raise ValueError("Number of observations in each group must be at least 1")
    if sig_level is not None and (sig_level < 0 or sig_level > 1):
        raise ValueError("sig_level must be between 0 and 1")
    if power is not None and (power < 0 or power > 1):
        raise ValueError("power must be between 0 and 1")
    alternative = alternative.casefold()
    if h is not None and alternative == "two-sided":
        h = abs(h)
    pwr = pwr_p(h, n, sig_level, power, alternative).pwr_test()
    if print_pretty:
        str_print = (
            "\t"
            + pwr["method"]
            + "\n" * 2
            + "\t" * 2
            + " " * 2
            + f"h = {pwr['effect_size']}"
            + "\n"
            + "\t" * 2
            + " " * 2
            + f"n = {pwr['n']}"
            + "\n"
            + "\t"
            + " " * 2
            + f"sig_level = {pwr['sig_level']}"
            + "\n"
            + "\t"
            + " " * 6
            + f"power = {pwr['power']}"
            + "\n"
            + "\t"
            + f"alternative = {pwr['alternative']}"
        )
        print(str_print)
    return pwr


def pwr_r_test(
    n: int | None = None,
    r: float | None = None,
    sig_level: float | None = None,
    power: float | None = None,
    alternative: str = "two-sided",
    print_pretty: bool = True,
) -> dict[str, float | int | str]:
    """Compute power for correlation test.

    Computes approximate power of test or determines parameters to obtain target
    power using the arctanh transformation. Similar to R's pwr.r.test function.

    Parameters
    ----------
    n : int | None, default=None
        Number of observations
    r : float | None, default=None
        Linear correlation coefficient
    sig_level : float | None, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power : float | None, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
        The alternative hypothesis
    print_pretty : bool, default=True
        Whether to print the results in a formatted output

    Returns
    -------
    dict[str, float | int | str]
        Dictionary containing effect_size (r), n, sig_level, power, alternative, and method

    Raises
    ------
    ValueError
        If not exactly one of r, n, sig_level, or power is None, or if
        parameter values are out of valid range
    """
    if not any(x is None for x in [r, n, sig_level, power]):
        raise ValueError("One of r, n, sig_level or power must be None")
    if sum([x is None for x in [r, n, sig_level, power]]) > 1:
        raise ValueError("Only one of r, n, sig_level or power may be None")
    if n is not None and n < 4:
        raise ValueError("Number of observations must be at least 4")
    if sig_level is not None and (sig_level < 0 or sig_level > 1):
        raise ValueError("sig_level must be between 0 and 1")
    if power is not None and (power < 0 or power > 1):
        raise ValueError("power must be between 0 and 1")
    alternative = alternative.casefold()
    if alternative == "two-sided" and r is not None:
        r = abs(r)
    pwr = pwr_r(r, n, sig_level, power, alternative).pwr_test()
    if print_pretty:
        str_print = (
            "\t"
            + pwr["method"]
            + "\n" * 2
            + "\t" * 2
            + " " * 2
            + f"r = {pwr['effect_size']}"
            + "\n"
            + "\t" * 2
            + " " * 2
            + f"n = {pwr['n']}"
            + "\n"
            + "\t"
            + " " * 2
            + f"sig_level = {pwr['sig_level']}"
            + "\n"
            + "\t"
            + " " * 6
            + f"power = {pwr['power']}"
            + "\n"
            + "\t"
            + f"alternative = {pwr['alternative']}"
        )
        print(str_print)
    return pwr


def pwr_t_test(
    n: int | None = None,
    d: float | None = None,
    sig_level: float | None = None,
    power: float | None = None,
    test_type: str = "paired",
    alternative: str = "two-sided",
    print_pretty: bool = True,
) -> dict[str, int | float | str | None]:
    """Compute power for t-test.

    Computes power of test or determines parameters to obtain target power.
    Similar to R's pwr.t.test function.

    Parameters
    ----------
    n : int | None, default=None
        Number of observations (per sample for two-sample, pairs for paired)
    d : float | None, default=None
        Effect size (Cohen's d) - standardized mean difference
    sig_level : float | None, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power : float | None, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    test_type : {'two-sample', 'one-sample', 'paired'}, default='paired'
        Type of t-test
    alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
        The alternative hypothesis
    print_pretty : bool, default=True
        Whether to print the results in a formatted output

    Returns
    -------
    dict[str, int | float | str | None]
        Dictionary containing n, effect_size (d), sig_level, power, alternative,
        method, and optionally note

    Raises
    ------
    ValueError
        If not exactly one of n, d, sig_level, or power is None, or if
        parameter values are out of valid range
    """
    if not any(x is None for x in [n, d, sig_level, power]):
        raise ValueError("One of n, d, sig_level or power must be None")
    if sum([x is None for x in [n, d, sig_level, power]]) > 1:
        raise ValueError("Only one of n, d, sig_level or power may be None")
    if n is not None and n < 2:
        raise ValueError("Number of observations must be at least 2")
    if sig_level is not None and (sig_level < 0 or sig_level > 1):
        raise ValueError("sig_level must be between 0 and 1")
    if power is not None and (power < 0 or power > 1):
        raise ValueError("power must be between 0 and 1")
    pwr = pwr_t(n, d, sig_level, power, test_type, alternative).pwr_test()
    if print_pretty:
        if "note" in pwr.keys():
            str_print = (
                "\t"
                + pwr["method"]
                + "\n" * 2
                + "\t" * 2
                + " " * 2
                + f"d = {pwr['effect_size']}"
                + "\n"
                + "\t" * 2
                + " " * 2
                + f"n = {pwr['n']}"
                + "\n"
                + "\t"
                + " " * 2
                + f"sig_level = {pwr['sig_level']}"
                + "\n"
                + "\t"
                + " " * 6
                + f"power = {pwr['power']}"
                + "\n"
                + "\t"
                + f"alternative = {pwr['alternative']}"
                + "\n" * 2
                + f"NOTE: {pwr['note']}"
            )
        else:
            str_print = (
                "\t"
                + pwr["method"]
                + "\n" * 2
                + "\t" * 2
                + " " * 2
                + f"d = {pwr['effect_size']}"
                + "\n"
                + "\t" * 2
                + " " * 2
                + f"n = {pwr['n']}"
                + "\n"
                + "\t"
                + " " * 2
                + f"sig_level = {pwr['sig_level']}"
                + "\n"
                + "\t"
                + " " * 6
                + f"power = {pwr['power']}"
                + "\n"
                + "\t"
                + f"alternative = {pwr['alternative']}"
            )
        print(str_print)
    return pwr


def pwr_t2n_test(
    n1: int | None = None,
    n2: int | None = None,
    d: float | None = None,
    sig_level: float | None = None,
    power: float | None = None,
    alternative: str = "two-sided",
    print_pretty: bool = True,
) -> dict[str, int | float | str]:
    """Compute power for two-sample t-test with unequal sample sizes.

    Computes power of test or determines parameters to obtain target power.
    Similar to R's pwr.t2n.test function.

    Parameters
    ----------
    n1 : int | None, default=None
        Number of observations in the first sample
    n2 : int | None, default=None
        Number of observations in the second sample
    d : float | None, default=None
        Effect size (Cohen's d)
    sig_level : float | None, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power : float | None, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
        The alternative hypothesis
    print_pretty : bool, default=True
        Whether to print the results in a formatted output

    Returns
    -------
    dict[str, int | float | str]
        Dictionary containing effect_size (d), n1, n2, sig_level, power,
        alternative, method, and note

    Raises
    ------
    ValueError
        If not exactly one of n1, n2, d, sig_level, or power is None, or if
        parameter values are out of valid range
    """
    if not any(x is None for x in [n1, n2, d, sig_level, power]):
        raise ValueError("One of n1, n2, d sig_level or power must be None")
    if sum([x is None for x in [n1, n2, d, sig_level, power]]) > 1:
        raise ValueError("Only one of n1, n2, d, sig_level or power may be None")
    if n1 is not None and n1 < 2:
        raise ValueError("Number of observations in the first group must be at least 2")
    if n2 is not None and n2 < 2:
        raise ValueError("Number of observations in the second group must be at least 2")
    if sig_level is not None and (sig_level < 0 or sig_level > 1):
        raise ValueError("sig_level must be between 0 and 1")
    if power is not None and (power < 0 or power > 1):
        raise ValueError("power must be between 0 and 1")
    alternative = alternative.casefold()
    if alternative == "two-sided" and d is not None:
        d = abs(d)
    pwr = pwr_t2n(d, n1, n2, sig_level, power, alternative).pwr_test()
    if print_pretty:
        str_print = (
            "\t"
            + pwr["method"]
            + "\n" * 2
            + "\t" * 2
            + " " * 2
            + f"d = {pwr['effect_size']}"
            + "\n"
            + "\t" * 2
            + " "
            + f"n1 = {pwr['n1']}"
            + "\n"
            + "\t" * 2
            + " "
            + f"n2 = {pwr['n2']}"
            + "\n"
            + "\t"
            + " " * 2
            + f"sig_level = {pwr['sig_level']}"
            + "\n"
            + "\t"
            + " " * 6
            + f"power = {pwr['power']}"
            + "\n"
            + "\t"
            + f"alternative = {pwr['alternative']}"
            + "\n" * 2
            + f"NOTE: {pwr['note']}"
        )
        print(str_print)
    return pwr
