from src.power_classes import *


def pwr_2p_test(h=None, n=None, sig_level=None, power=None, alternative="two-sided"):
    """Compute power of test, or determine parameters to obtain target power (similar to power.prop.test).

    Parameters
    ----------
    h: float, default=None
        The effect size
    n: int, default=None
        Number of observations (per sample)
    sig_level: float, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power: float, default=None
        Power of test (1 minus Type II error probability). Must be betwen 0 and 1
    alternative: {'two-sided', 'greater', 'less'}
        A character string specifying the alternative hypothesis

    Returns
    -------
    A dict containing our h, n, sig_level, power and alternative hypothesis
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
    pwr = pwr_2p(h, n, sig_level, power, alternative)
    return pwr.pwr_test()


def pwr_2p2n_test(
    h=None, n1=None, n2=None, sig_level=None, power=None, alternative="two-sided"
):
    """Compute power of test, or determine parameters to obtain target power.

    Parameters
    ----------
    h: float, default=None
        The effect size
    n1: int, default=None
        Number of observations in the first sample
    n2: int, default=None
        Number of observations in the second sample
    sig_level: float, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power: float, default=None
        Power of the test (1 minus Type II error probability). Must be between 0 and 1
    alternative: {'two-sided', 'greater', 'less'}
        A character string specifying the alternative hypothesis

    Returns
    -------
    A dict containing our h, n1, n2, sig_level, power and alternative hypothesis
    """
    if not any(v is None for v in [h, n1, n2, sig_level, power]):
        raise ValueError("One of h, n1, n2, sig_level or power must be None")
    if sum([v is None for v in [h, n1, n2, sig_level, power]]) > 1:
        raise ValueError("Only one of h, n1, n2, sig_level or power may be None")
    if n1 is not None and n1 < 2:
        raise ValueError("Number of observations in the first group must be at least 2")
    if n2 is not None and n2 < 2:
        raise ValueError(
            "Number of observations in the second group must be at least 2"
        )
    if sig_level is not None and (sig_level < 0 or sig_level > 1):
        raise ValueError("sig_level must be between 0 and 1")
    if power is not None and (power < 0 or power > 1):
        raise ValueError("power must be between 0 and 1")
    alternative = alternative.casefold()
    if alternative == "two-sided" and h is not None:
        h = abs(h)
    pwr = pwr_2p2n(h, n1, n2, sig_level, power, alternative)
    return pwr.pwr_test()


def pwr_anova_test(k=None, n=None, f=None, sig_level=None, power=None):
    """Compute power of test or determine parameters to obtain target power (same as power.anova.test).

    Parameters
    ----------
    k: int, default=None
        Number of groups
    n: int, default=None
        Number of observations (per group)
    f: float, default=None
        Effect size
    sig_level: float, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power: float, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1

    Returns
    -------
    A dict containing our k, n, f, sig_level, and power
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
    pwr = pwr_anova(k, n, f, sig_level, power)
    return pwr.pwr_test()


def pwr_chisq_test(w=None, n=None, df=1, sig_level=None, power=None):
    """Compute power of test or determine parameters to obtain target power (same as power.anova.test).

    Parameters
    ----------
    w: float, default=None
        Effect size
    n: int, default=None
        Total number of observations
    df: int, default=1
        Degrees of freedom
    sig_level: float, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power: float, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1

    Returns
    -------
    A dict containing our w, n, df, sig_level and power
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
    pwr = pwr_chisq(w, n, df, sig_level, power)
    return pwr.pwr_test()


def pwr_f2_test(u=None, v=None, f2=None, sig_level=None, power=None):
    """Compute power of test or determine parameters to obtain target power (same as power.anova.test).

    Parameters
    ----------
    u: int, default=None
        Degrees of freedom for the numerator
    v: int, default=None
        Degrees of freedom for the denominator
    f2: float, default=None
        Effect size
    sig_level: float, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power: float, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1

    Returns
    -------
    A dict containing our u, v, f2, sig_level and power
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
    pwr = pwr_f2(u, v, f2, sig_level, power)
    return pwr.pwr_test()


def pwr_norm_test(d=None, n=None, sig_level=None, power=None, alternative="two-sided"):
    """Compute power of test or determine parameters to obtain target power (same as power.anova.test).

    Parameters
    ----------
    d: float, default=None
        Effect size
    n: int, default=None
        Number of observations
    sig_level: float, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power: float, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    alternative: {'two-sided', 'greater', 'less'}
        A character string specifying the alternative hypothesis

    Returns
    -------
    A dict containing our d, n, sig_level, power and alternative hypothesis
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
    pwr = pwr_norm(d, n, sig_level, power, alternative)
    return pwr.pwr_test()


def pwr_p_test(h=None, n=None, sig_level=None, power=None, alternative="two-sided"):
    """Compute power of test or determine parameters to obtain target power (same as power.anova.test).

    Parameters
    ----------
    h: float, default=None
        Effect size
    n: int, default=None
        Number of observations
    sig_level: float, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power: float, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    alternative: {'two-sided', 'greater', 'less'}
        A character string specifying the alternative hypothesis

    Returns
    -------
    A dict containing our h, n, sig_level, power and alternative hypothesis
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
    pwr = pwr_p(h, n, sig_level, power, alternative)
    return pwr.pwr_test()


def pwr_r_test(n=None, r=None, sig_level=None, power=None, alternative="two-sided"):
    """Compute power of test or determine parameters to obtain target power (same as power.anova.test).

    Parameters
    ----------
    n: int, default=None
        Number of observations
    r: float, default=none
        Linear correlation coefficient
    sig_level: float, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power: float, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    alternative: {'two-sided', 'greater', 'less'}
        A character string specifying the alternative hypothesis

    Returns
    -------
    A dict containing our h, n, sig_level, power and alternative hypothesis
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
    pwr = pwr_r(r, n, sig_level, power, alternative)
    return pwr.pwr_test()


def pwr_t_test(n=None, d=None, sig_level=None, power=None, type="paired", alternative="two-sided"):
    """Compute power of tests or determine parameters to obtain target power (similar to as power.t.test)

    Parameters
    ----------
    n: int, default=None
        Number of observations (per sample)
    d: float
        Effect size (Cohen’s d) - difference between the means divided by the pooled standard deviation
    sig_level: float, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power: float, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    type: {'two-sample', 'one-sample', 'paired'}
        Type of t-test: One sample, two sample or paired sample
    alternative: {'two-sided', 'greater', 'less'}
        A character string specifying the alternative hypothesis

    Returns
    -------
    A dict containing our n, d, sig_level, power and alternative hypothesis
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
    pwr = pwr_t(n, d, sig_level, power, type, alternative)
    return pwr.pwr_test()


def pwr_t2n_test(n1=None, n2=None, d=None, sig_level=None, power=None, alternative="two-sided"):
    """Compute power of tests or determine parameters to obtain target power (similar to as power.t.test)

    Parameters
    ----------
    n1: int, default=None
        Number of observations in the first sample
    n2: int, default=None
        Number of observations in the second sample
    d: float, default=None
        Effect size
    sig_level: float, default=None
        Significance level (Type I error probability). Must be between 0 and 1
    power: float, default=None
        Power of test (1 minus Type II error probability). Must be between 0 and 1
    alternative: {'two-sided', 'greater', 'less'}
        A character string specifying the alternative hypothesis

    Returns
    -------
    A dict containing our n1, n2, d, sig_level, power and alternative hypothesis
    """
    if not any(x is None for x in [n1, n2, d, sig_level, power]):
        raise ValueError("One of n1, n2, d sig_level or power must be None")
    if sum([x is None for x in [n1, n2, d, sig_level, power]]) > 1:
        raise ValueError("Only one of n1, n2, d, sig_level or power may be None")
    if n1 is not None and n1 < 2:
        raise ValueError("Number of observations in the first group must be at least 2")
    if n2 is not None and n2 < 2:
        raise ValueError(
            "Number of observations in the second group must be at least 2"
        )
    if sig_level is not None and (sig_level < 0 or sig_level > 1):
        raise ValueError("sig_level must be between 0 and 1")
    if power is not None and (power < 0 or power > 1):
        raise ValueError("power must be between 0 and 1")
    alternative = alternative.casefold()
    if alternative == "two-sided" and d is not None:
        d = abs(d)
    pwr = pwr_t2n(d, n1, n2, sig_level, power, alternative)
    return pwr.pwr_test()