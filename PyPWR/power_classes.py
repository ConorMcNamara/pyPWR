"""Power calculation classes for various statistical tests.

This module provides abstract base classes and concrete implementations for
calculating statistical power, sample size, effect size, or significance level
for various hypothesis tests.
"""

import abc
from math import atanh, sqrt

import numpy as np
from scipy.optimize import brentq, bisect
from scipy.stats import chi2, f as f_dist, ncf, nct, ncx2, norm, t as t_dist


class pwr_1n(abc.ABC):
    """Abstract base class for power tests with one sample size parameter.

    This class provides the common interface and logic for power calculations
    that involve a single sample size (or equal sample sizes per group).
    """

    def __init__(
        self,
        effect_size: float | None,
        n: int | None,
        sig_level: float | None,
        power: float | None,
        alternative: str = "two-sided",
    ) -> None:
        """Initialize the power test with one sample size.

        Parameters
        ----------
        effect_size : float | None
            The effect size
        n : int | None
            Number of observations (per sample)
        sig_level : float | None
            Significance level (Type I error probability). Must be between 0 and 1
        power : float | None
            Power of test (1 minus Type II error probability). Must be between 0 and 1
        alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
            A character string specifying the alternative hypothesis
        """
        self.effect_size = effect_size
        self.n = n
        self.sig_level = sig_level
        self.power = power
        self.alternative = alternative.casefold()
        self.method = "Difference of proportion power calculation for binomial distribution (arcsine transformation)"
        self.note: str | None = "Same sample sizes"

    @abc.abstractmethod
    def _get_power(self) -> float:
        """Calculate power given effect_size, n, and sig_level."""
        pass

    @abc.abstractmethod
    def _get_effect_size(self, effect_size: float) -> float:
        """Calculate difference from target power for given effect_size."""
        pass

    @abc.abstractmethod
    def _get_n(self, n: int) -> float:
        """Calculate difference from target power for given n."""
        pass

    @abc.abstractmethod
    def _get_sig_level(self, sig_level: float) -> float:
        """Calculate difference from target power for given sig_level."""
        pass

    def pwr_test(self) -> dict[str, float | int | str | None]:
        if self.power is None:
            self.power = self._get_power()
        elif self.effect_size is None:
            if self.alternative == "two-sided":
                self.effect_size = brentq(self._get_effect_size, 1e-10, 10)
            elif self.alternative == "greater":
                self.effect_size = brentq(self._get_effect_size, -5, 10)
            else:
                self.effect_size = brentq(self._get_effect_size, -10, 5)
        elif self.n is None:
            self.n = np.ceil(brentq(self._get_n, 2 + 1e-10, 1e09))
        else:
            self.sig_level = brentq(self._get_sig_level, 1e-10, 1 - 1e-10)
        if self.note is not None:
            return {
                "effect_size": self.effect_size,
                "n": self.n,
                "sig_level": self.sig_level,
                "power": self.power,
                "alternative": self.alternative,
                "method": self.method,
                "note": self.note,
            }
        else:
            return {
                "effect_size": self.effect_size,
                "n": self.n,
                "sig_level": self.sig_level,
                "power": self.power,
                "alternative": self.alternative,
                "method": self.method,
            }


class pwr_2n(abc.ABC):
    """Abstract base class for power tests with two different sample sizes.

    This class provides the common interface and logic for power calculations
    that involve two independent samples with potentially different sizes.
    """

    def __init__(
        self,
        effect_size: float | None,
        n1: int | None,
        n2: int | None,
        sig_level: float | None,
        power: float | None,
        alternative: str = "two-sided",
    ) -> None:
        """Initialize the power test with two sample sizes.

        Parameters
        ----------
        effect_size : float | None
            The effect size
        n1 : int | None
            Number of observations in the first sample
        n2 : int | None
            Number of observations in the second sample
        sig_level : float | None
            Significance level (Type I error probability). Must be between 0 and 1
        power : float | None
            Power of the test (1 minus Type II error probability). Must be between 0 and 1
        alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
            A character string specifying the alternative hypothesis
        """
        self.effect_size = effect_size
        self.n1 = n1
        self.n2 = n2
        self.sig_level = sig_level
        self.power = power
        self.alternative = alternative.casefold()
        self.method = "Difference of proportion power calculation for binomial distribution (arcsine transformation)"
        self.note: str | None = "Different sample sizes"

    @abc.abstractmethod
    def _get_power(self) -> float:
        """Calculate power given effect_size, n1, n2, and sig_level."""
        pass

    @abc.abstractmethod
    def _get_effect_size(self, effect_size: float) -> float:
        """Calculate difference from target power for given effect_size."""
        pass

    @abc.abstractmethod
    def _get_n1(self, n1: int) -> float:
        """Calculate difference from target power for given n1."""
        pass

    @abc.abstractmethod
    def _get_n2(self, n2: int) -> float:
        """Calculate difference from target power for given n2."""
        pass

    @abc.abstractmethod
    def _get_sig_level(self, sig_level: float) -> float:
        """Calculate difference from target power for given sig_level."""
        pass

    def pwr_test(self) -> dict[str, float | int | str | None]:
        if self.power is None:
            self.power = self._get_power()
        elif self.effect_size is None:
            if self.alternative == "two-sided":
                self.effect_size = brentq(self._get_effect_size, 1e-10, 10)
            elif self.alternative == "greater":
                self.effect_size = brentq(self._get_effect_size, -5, 10)
            else:
                self.effect_size = brentq(self._get_effect_size, -10, 5)
        elif self.n1 is None:
            self.n1 = np.ceil(brentq(self._get_n1, 2 + 1e-10, 500))
        elif self.n2 is None:
            self.n2 = np.ceil(brentq(self._get_n2, 2 + 1e-10, 500))
        else:
            self.sig_level = brentq(self._get_sig_level, 1e-10, 1 - 1e-10)
        if self.note is not None:
            return {
                "effect_size": self.effect_size,
                "n1": self.n1,
                "n2": self.n2,
                "sig_level": self.sig_level,
                "power": self.power,
                "alternative": self.alternative,
                "method": self.method,
                "note": "Different sample sizes",
            }
        else:
            return {
                "effect_size": self.effect_size,
                "n1": self.n1,
                "n2": self.n2,
                "sig_level": self.sig_level,
                "power": self.power,
                "alternative": self.alternative,
                "method": self.method,
            }


class pwr_2p(pwr_1n):
    """Power calculation for two proportions test (equal sample sizes).

    Computes power calculations for testing the difference between two proportions
    using the arcsine transformation, with equal sample sizes per group.
    """

    def __init__(
        self,
        h: float | None,
        n: int | None,
        sig_level: float | None,
        power: float | None,
        alternative: str = "two-sided",
    ) -> None:
        """Initialize two proportions test with equal sample sizes.

        Parameters
        ----------
        h : float | None
            Effect size (arcsine transformation)
        n : int | None
            Sample size per group
        sig_level : float | None
            Significance level (Type I error probability). Must be between 0 and 1
        power : float | None
            Power of test. Must be between 0 and 1
        alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
            The alternative hypothesis
        """
        super().__init__(h, n, sig_level, power, alternative)

    def _get_power(self) -> float:
        if self.alternative == "two-sided":
            power = norm.sf(norm.isf(self.sig_level / 2) - self.effect_size * sqrt(self.n / 2)) + norm.cdf(
                norm.ppf(self.sig_level / 2) - self.effect_size * sqrt(self.n / 2)
            )
        elif self.alternative == "greater":
            power = norm.sf(norm.isf(self.sig_level) - self.effect_size * sqrt(self.n / 2))
        else:
            power = norm.cdf(norm.ppf(self.sig_level) - self.effect_size * sqrt(self.n / 2))
        return power

    def _get_effect_size(self, h: float) -> float:
        if self.alternative == "two-sided":
            h = (
                norm.sf(norm.isf(self.sig_level / 2) - h * sqrt(self.n / 2))
                + norm.cdf(norm.ppf(self.sig_level / 2) - h * sqrt(self.n / 2))
                - self.power
            )
        elif self.alternative == "greater":
            h = norm.sf(norm.isf(self.sig_level) - h * sqrt(self.n / 2)) - self.power
        else:
            h = norm.cdf(norm.ppf(self.sig_level) - h * sqrt(self.n / 2)) - self.power
        return h

    def _get_n(self, n: int) -> float:
        if self.alternative == "two-sided":
            n = (
                norm.sf(norm.isf(self.sig_level / 2) - self.effect_size * sqrt(n / 2))
                + norm.cdf(norm.ppf(self.sig_level / 2) - self.effect_size * sqrt(n / 2))
                - self.power
            )
        elif self.alternative == "greater":
            n = norm.sf(norm.isf(self.sig_level) - self.effect_size * sqrt(n / 2)) - self.power
        else:
            n = norm.cdf(norm.ppf(self.sig_level) - self.effect_size * sqrt(n / 2)) - self.power
        return n

    def _get_sig_level(self, sig_level: float) -> float:
        if self.alternative == "two-sided":
            sig_level = (
                norm.sf(norm.isf(sig_level / 2) - self.effect_size * sqrt(self.n / 2))
                + norm.cdf(norm.ppf(sig_level / 2) - self.effect_size * sqrt(self.n / 2))
                - self.power
            )
        elif self.alternative == "greater":
            sig_level = norm.sf(norm.isf(sig_level) - self.effect_size * sqrt(self.n / 2)) - self.power
        else:
            sig_level = norm.cdf(norm.ppf(sig_level) - self.h * sqrt(self.n / 2)) - self.power
        return sig_level


class pwr_2p2n(pwr_2n):
    """Power calculation for two proportions test (unequal sample sizes).

    Computes power calculations for testing the difference between two proportions
    using the arcsine transformation, with different sample sizes per group.
    """

    def __init__(
        self,
        h: float | None,
        n1: int | None,
        n2: int | None,
        sig_level: float | None,
        power: float | None,
        alternative: str = "two-sided",
    ) -> None:
        """Initialize two proportions test with unequal sample sizes.

        Parameters
        ----------
        h : float | None
            Effect size (arcsine transformation)
        n1 : int | None
            Sample size for the first group
        n2 : int | None
            Sample size for the second group
        sig_level : float | None
            Significance level (Type I error probability). Must be between 0 and 1
        power : float | None
            Power of test. Must be between 0 and 1
        alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
            The alternative hypothesis
        """
        super().__init__(h, n1, n2, sig_level, power, alternative)

    def _get_power(self) -> float:
        if self.alternative == "two-sided":
            power = norm.sf(
                norm.isf(self.sig_level / 2) - self.effect_size * sqrt(self.n1 * self.n2 / (self.n1 + self.n2))
            ) + norm.cdf(
                norm.ppf(self.sig_level / 2) - self.effect_size * sqrt(self.n1 * self.n2 / (self.n1 + self.n2))
            )
        elif self.alternative == "greater":
            power = norm.sf(norm.isf(self.sig_level) - self.effect_size * sqrt(self.n1 * self.n2 / (self.n1 + self.n2)))
        else:
            power = norm.cdf(
                norm.ppf(self.sig_level) - self.effect_size * sqrt(self.n1 * self.n2 / (self.n1 + self.n2))
            )
        return power

    def _get_effect_size(self, h: float) -> float:
        if self.alternative == "two-sided":
            h = (
                norm.sf(norm.isf(self.sig_level / 2) - h * sqrt(self.n1 * self.n2 / (self.n1 + self.n2)))
                + norm.cdf(norm.ppf(self.sig_level / 2) - h * sqrt(self.n1 * self.n2 / (self.n1 + self.n2)))
                - self.power
            )
        elif self.alternative == "greater":
            h = norm.sf(norm.isf(self.sig_level) - h * sqrt(self.n1 * self.n2 / (self.n1 + self.n2))) - self.power
        else:
            h = norm.cdf(norm.ppf(self.sig_level) - h * sqrt(self.n1 * self.n2 / (self.n1 + self.n2))) - self.power
        return h

    def _get_n1(self, n1: int) -> float:
        if self.alternative == "two-sided":
            n1 = (
                norm.sf(norm.isf(self.sig_level / 2) - self.effect_size * sqrt(n1 * self.n2 / (n1 + self.n2)))
                + norm.cdf(norm.ppf(self.sig_level / 2) - self.effect_size * sqrt(n1 * self.n2 / (n1 + self.n2)))
                - self.power
            )
        elif self.alternative == "greater":
            n1 = norm.sf(norm.isf(self.sig_level) - self.effect_size * sqrt(n1 * self.n2 / (n1 + self.n2))) - self.power
        else:
            n1 = (
                norm.cdf(norm.ppf(self.sig_level) - self.effect_size * sqrt(n1 * self.n2 / (n1 + self.n2))) - self.power
            )
        return n1

    def _get_n2(self, n2: int) -> float:
        if self.alternative == "two-sided":
            n2 = (
                norm.sf(norm.isf(self.sig_level / 2) - self.effect_size * sqrt(self.n1 * n2 / (self.n1 + n2)))
                + norm.cdf(norm.ppf(self.sig_level / 2) - self.effect_size * sqrt(self.n1 * n2 / (self.n1 + n2)))
                - self.power
            )
        elif self.alternative == "greater":
            n2 = norm.sf(norm.isf(self.sig_level) - self.effect_size * sqrt(self.n1 * n2 / (self.n1 + n2))) - self.power
        else:
            n2 = (
                norm.cdf(norm.ppf(self.sig_level) - self.effect_size * sqrt(self.n1 * n2 / (self.n1 + n2))) - self.power
            )
        return n2

    def _get_sig_level(self, sig_level: float) -> float:
        if self.alternative == "two-sided":
            sig_level = (
                norm.sf(norm.isf(sig_level / 2) - self.effect_size * sqrt(self.n1 * self.n2 / (self.n1 + self.n2)))
                + norm.cdf(norm.ppf(sig_level / 2) - self.effect_size * sqrt(self.n1 * self.n2 / (self.n1 + self.n2)))
                - self.power
            )
        elif self.alternative == "greater":
            sig_level = (
                norm.sf(norm.isf(sig_level) - self.effect_size * sqrt(self.n1 * self.n2 / (self.n1 + self.n2)))
                - self.power
            )
        else:
            sig_level = (
                norm.cdf(norm.ppf(sig_level) - self.effect_size * sqrt(self.n1 * self.n2 / (self.n1 + self.n2)))
                - self.power
            )
        return sig_level


class pwr_anova:
    """Power calculation for balanced one-way ANOVA.

    Computes power calculations for balanced one-way analysis of variance tests.
    """

    def __init__(
        self,
        k: int | None,
        n: int | None,
        f: float | None,
        sig_level: float | None,
        power: float | None,
    ) -> None:
        """Initialize balanced one-way ANOVA power test.

        Parameters
        ----------
        k : int | None
            Number of groups
        n : int | None
            Number of samples per group
        f : float | None
            Effect size
        sig_level : float | None
            Significance level (Type I error probability). Must be between 0 and 1
        power : float | None
            Power of test. Must be between 0 and 1
        """
        self.k = k
        self.n = n
        self.f = f
        self.sig_level = sig_level
        self.power = power
        self.method = "Balanced one-way analysis of variance power calculation"
        self.note = "n is number in each group"

    def _get_power(self) -> float:
        l_var = self.k * self.n * pow(self.f, 2)
        power = ncf.sf(
            f_dist.isf(self.sig_level, self.k - 1, (self.n - 1) * self.k),
            self.k - 1,
            (self.n - 1) * self.k,
            l_var,
        )
        return power

    def _get_k(self, k: int) -> float:
        l_var = k * self.n * pow(self.f, 2)
        k = (
            ncf.sf(
                f_dist.isf(self.sig_level, k - 1, (self.n - 1) * k),
                k - 1,
                (self.n - 1) * k,
                l_var,
            )
            - self.power
        )
        return k

    def _get_n(self, n: int) -> float:
        l_var = self.k * n * pow(self.f, 2)
        n = (
            ncf.sf(
                f_dist.isf(self.sig_level, self.k - 1, (n - 1) * self.k),
                self.k - 1,
                (n - 1) * self.k,
                l_var,
            )
            - self.power
        )
        return n

    def _get_effect_size(self, f: float) -> float:
        l_var = self.k * self.n * pow(f, 2)
        f = (
            ncf.sf(
                f_dist.isf(self.sig_level, self.k - 1, (self.n - 1) * self.k),
                self.k - 1,
                (self.n - 1) * self.k,
                l_var,
            )
            - self.power
        )
        return f

    def _get_sig_level(self, sig_level: float) -> float:
        l_var = self.k * self.n * pow(self.f, 2)
        sig_level = (
            ncf.sf(
                f_dist.isf(sig_level, self.k - 1, (self.n - 1) * self.k),
                self.k - 1,
                (self.n - 1) * self.k,
                l_var,
            )
            - self.power
        )
        return sig_level

    def pwr_test(self) -> dict[str, int | float | str]:
        """Perform power calculation for ANOVA test."""
        if self.power is None:
            self.power = self._get_power()
        elif self.k is None:
            self.k = np.ceil(brentq(self._get_k, 2 + 1e-10, 100))
        elif self.n is None:
            self.n = np.ceil(brentq(self._get_n, 2 + 1e-10, 1e09))
        elif self.f is None:
            self.f = bisect(self._get_effect_size, 1e-07, 1e07)
        else:
            self.sig_level = brentq(self._get_sig_level, 1e-10, 1 - 1e-10)
        return {
            "k": self.k,
            "n": self.n,
            "effect_size": self.f,
            "sig_level": self.sig_level,
            "power": self.power,
            "method": self.method,
            "note": self.note,
        }


class pwr_chisq:
    """Power calculation for chi-squared test.

    Computes power calculations for chi-squared tests of association
    or goodness of fit.
    """

    def __init__(
        self,
        w: float | None,
        n: int | None,
        df: int | None,
        sig_level: float | None,
        power: float | None,
    ) -> None:
        """Initialize chi-squared test power calculation.

        Parameters
        ----------
        w : float | None
            Effect size
        n : int | None
            Total sample size
        df : int | None
            Degrees of freedom
        sig_level : float | None
            Significance level (Type I error probability). Must be between 0 and 1
        power : float | None
            Power of test. Must be between 0 and 1
        """
        self.w = w
        self.n = n
        self.df = df
        self.sig_level = sig_level
        self.power = power
        self.method = "Chi squared power calculation"
        self.note = "N is the number of observations"

    def _get_power(self) -> float:
        k = chi2.isf(self.sig_level, self.df)
        power = ncx2.sf(k, self.df, self.n * pow(self.w, 2))
        return power

    def _get_effect_size(self, w: float) -> float:
        k = chi2.isf(self.sig_level, self.df)
        w = ncx2.sf(k, self.df, self.n * pow(w, 2)) - self.power
        return w

    def _get_n(self, n: int) -> float:
        k = chi2.isf(self.sig_level, self.df)
        n = ncx2.sf(k, self.df, n * pow(self.w, 2)) - self.power
        return n

    def _get_sig_level(self, sig_level: float) -> float:
        k = chi2.isf(sig_level, self.df)
        sig_level = ncx2.sf(k, self.df, self.n * pow(self.w, 2)) - self.power
        return sig_level

    def pwr_test(self) -> dict[str, int | float | str]:
        """Perform power calculation for chi-squared test."""
        if self.power is None:
            self.power = self._get_power()
        elif self.w is None:
            self.w = brentq(self._get_effect_size, 1e-10, 1e08)
        elif self.n is None:
            self.n = np.ceil(brentq(self._get_n, 1 + 1e-10, 1e09))
        else:
            self.sig_level = brentq(self._get_sig_level, 1e-10, 1 - 1e-10)
        return {
            "effect_size": self.w,
            "n": self.n,
            "df": self.df,
            "sig_level": self.sig_level,
            "power": self.power,
            "method": self.method,
            "note": self.note,
        }


class pwr_f2:
    """Power calculation for general linear model (F-test).

    Computes power calculations for the general linear model, including
    multiple regression and other F-tests.
    """

    def __init__(
        self,
        u: int | None,
        v: int | None,
        f2: float | None,
        sig_level: float | None,
        power: float | None,
    ) -> None:
        """Initialize general linear model power calculation.

        Parameters
        ----------
        u : int | None
            Degrees of freedom for the numerator
        v : int | None
            Degrees of freedom for the denominator
        f2 : float | None
            Effect size (f-squared)
        sig_level : float | None
            Significance level (Type I error probability). Must be between 0 and 1
        power : float | None
            Power of test. Must be between 0 and 1
        """
        self.u = u
        self.v = v
        self.f2 = f2
        self.sig_level = sig_level
        self.power = power
        self.method = "Multiple regression power calculator"

    def _get_power(self) -> float:
        l_var = self.f2 * (self.u + self.v + 1)
        power = ncf.sf(f_dist.isf(self.sig_level, self.u, self.v), self.u, self.v, l_var)
        return power

    def _get_u(self, u: int) -> float:
        l_var = self.f2 * (u + self.v + 1)
        u = ncf.sf(f_dist.isf(self.sig_level, u, self.v), u, self.v, l_var) - self.power
        return u

    def _get_v(self, v: int) -> float:
        l_var = self.f2 * (self.u + v + 1)
        v = ncf.sf(f_dist.isf(self.sig_level, self.u, v), self.u, v, l_var) - self.power
        return v

    def _get_effect_size(self, f2: float) -> float:
        l_var = f2 * (self.u + self.v + 1)
        f2 = ncf.sf(f_dist.isf(self.sig_level, self.u, self.v), self.u, self.v, l_var) - self.power
        return f2

    def _get_sig_level(self, sig_level: float) -> float:
        l_var = self.f2 * (self.u + self.v + 1)
        sig_level = ncf.sf(f_dist.isf(sig_level, self.u, self.v), self.u, self.v, l_var) - self.power
        return sig_level

    def pwr_test(self) -> dict[str, int | float | str]:
        """Perform power calculation for general linear model."""
        if self.power is None:
            self.power = self._get_power()
        elif self.u is None:
            self.u = np.ceil(brentq(self._get_u, 1 + 1e-10, 100))
        elif self.v is None:
            self.v = np.ceil(brentq(self._get_v, 1 + 1e-10, 1e09))
        elif self.f2 is None:
            self.f2 = bisect(self._get_effect_size, 1e-07, 1e07)
        else:
            self.sig_level = brentq(self._get_sig_level, 1e-10, 1 - 1e-10)
        return {
            "u": self.u,
            "v": self.v,
            "effect_size": self.f2,
            "sig_level": self.sig_level,
            "power": self.power,
            "method": self.method,
        }


class pwr_norm(pwr_1n):
    """Power calculation for normal distribution with known variance.

    Computes power calculations for testing the mean of a normal distribution
    when the variance is known.
    """

    def __init__(
        self,
        d: float | None,
        n: int | None,
        sig_level: float | None,
        power: float | None,
        alternative: str = "two-sided",
    ) -> None:
        """Initialize normal distribution power calculation.

        Parameters
        ----------
        d : float | None
            Effect size (standardized mean difference)
        n : int | None
            Sample size
        sig_level : float | None
            Significance level (Type I error probability). Must be between 0 and 1
        power : float | None
            Power of test. Must be between 0 and 1
        alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
            The alternative hypothesis
        """
        super().__init__(d, n, sig_level, power, alternative)
        self.method = "Mean power calculation for normal distribution with known variance"
        self.note = None

    def _get_power(self) -> float:
        if self.alternative == "two-sided":
            power = norm.sf(norm.isf(self.sig_level / 2) - self.effect_size * sqrt(self.n)) + norm.cdf(
                norm.ppf(self.sig_level / 2) - self.effect_size * sqrt(self.n)
            )
        elif self.alternative == "greater":
            power = norm.sf(norm.isf(self.sig_level) - self.effect_size * sqrt(self.n))
        else:
            power = norm.cdf(norm.ppf(self.sig_level) - self.effect_size * sqrt(self.n))
        return power

    def _get_effect_size(self, effect_size: float) -> float:
        if self.alternative == "two-sided":
            effect_size = (
                norm.sf(norm.isf(self.sig_level / 2) - effect_size * sqrt(self.n))
                + norm.cdf(norm.ppf(self.sig_level / 2) - effect_size * sqrt(self.n))
                - self.power
            )
        elif self.alternative == "greater":
            effect_size = norm.sf(norm.isf(self.sig_level) - effect_size * sqrt(self.n)) - self.power
        else:
            effect_size = norm.cdf(norm.ppf(self.sig_level) - effect_size * sqrt(self.n)) - self.power
        return effect_size

    def _get_n(self, n: int) -> float:
        if self.alternative == "two-sided":
            n = (
                norm.sf(norm.isf(self.sig_level / 2) - self.effect_size * sqrt(n))
                + norm.cdf(norm.ppf(self.sig_level / 2) - self.effect_size * sqrt(n))
                - self.power
            )
        elif self.alternative == "greater":
            n = norm.sf(norm.isf(self.sig_level) - self.effect_size * sqrt(n)) - self.power
        else:
            n = norm.cdf(norm.ppf(self.sig_level) - self.effect_size * sqrt(n)) - self.power
        return n

    def _get_sig_level(self, sig_level: float) -> float:
        if self.alternative == "two-sided":
            sig_level = (
                norm.sf(norm.isf(sig_level / 2) - self.effect_size * sqrt(self.n))
                + norm.cdf(norm.ppf(sig_level / 2) - self.effect_size * sqrt(self.n))
                - self.power
            )
        elif self.alternative == "greater":
            sig_level = norm.sf(norm.isf(sig_level) - self.effect_size * sqrt(self.n)) - self.power
        else:
            sig_level = norm.cdf(norm.ppf(sig_level) - self.effect_size * sqrt(self.n)) - self.power
        return sig_level


class pwr_p(pwr_norm):
    """Power calculation for one-sample proportion test.

    Computes power calculations for testing a single proportion using
    the arcsine transformation.
    """

    def __init__(
        self,
        h: float | None,
        n: int | None,
        sig_level: float | None,
        power: float | None,
        alternative: str = "two-sided",
    ) -> None:
        """Initialize one-sample proportion test power calculation.

        Parameters
        ----------
        h : float | None
            Effect size (arcsine transformation)
        n : int | None
            Sample size
        sig_level : float | None
            Significance level (Type I error probability). Must be between 0 and 1
        power : float | None
            Power of test. Must be between 0 and 1
        alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
            The alternative hypothesis
        """
        super().__init__(h, n, sig_level, power, alternative)
        self.method = "Proportion power calculation for binomial distribution (arcsine transformation)"
        self.note = None


class pwr_r(pwr_1n):
    """Power calculation for correlation test.

    Computes approximate power calculations for testing correlation coefficients
    using the arctanh (Fisher's Z) transformation.
    """

    def __init__(
        self,
        r: float | None,
        n: int | None,
        sig_level: float | None,
        power: float | None,
        alternative: str = "two-sided",
    ) -> None:
        """Initialize correlation test power calculation.

        Parameters
        ----------
        r : float | None
            Correlation coefficient
        n : int | None
            Sample size
        sig_level : float | None
            Significance level (Type I error probability). Must be between 0 and 1
        power : float | None
            Power of test. Must be between 0 and 1
        alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
            The alternative hypothesis
        """
        super().__init__(r, n, sig_level, power, alternative)
        self.method = "Approximate correlation power calculation (arctanh transformation)"
        self.note = None

    def _get_power(self) -> float:
        if self.alternative == "less":
            self.effect_size *= -1
        if self.alternative == "two-sided":
            sig_level = self.sig_level / 2
        else:
            sig_level = self.sig_level
        ttt = t_dist.isf(sig_level, df=self.n - 2)
        rc = sqrt(pow(ttt, 2) / (pow(ttt, 2) + self.n - 2))
        zr = atanh(self.effect_size) + self.effect_size / (2 * (self.n - 1))
        zrc = atanh(rc)
        if self.alternative == "two-sided":
            power = norm.cdf((zr - zrc) * sqrt(self.n - 3)) + norm.cdf((-zr - zrc) * sqrt(self.n - 3))
        else:
            power = norm.cdf((zr - zrc) * sqrt(self.n - 3))
        return power

    def _get_effect_size(self, effect_size: float) -> float:
        if self.alternative == "less":
            effect_size *= -1
        if self.alternative == "two-sided":
            sig_level = self.sig_level / 2
        else:
            sig_level = self.sig_level
        ttt = t_dist.isf(sig_level, df=self.n - 2)
        zr = atanh(effect_size) + effect_size / (2 * (self.n - 1))
        rc = sqrt(pow(ttt, 2) / (pow(ttt, 2) + self.n - 2))
        zrc = atanh(rc)
        if self.alternative == "two-sided":
            effect_size = (
                norm.cdf((zr - zrc) * sqrt(self.n - 3)) + norm.cdf((-zr - zrc) * sqrt(self.n - 3)) - self.power
            )
        else:
            effect_size = norm.cdf((zr - zrc) * sqrt(self.n - 3)) - self.power
        return effect_size

    def _get_n(self, n: int) -> float:
        if self.alternative == "less":
            self.effect_size *= -1
        if self.alternative == "two-sided":
            sig_level = self.sig_level / 2
        else:
            sig_level = self.sig_level
        ttt = t_dist.isf(sig_level, df=n - 2)
        zr = atanh(self.effect_size) + self.effect_size / (2 * (n - 1))
        rc = sqrt(pow(ttt, 2) / (pow(ttt, 2) + n - 2))
        zrc = atanh(rc)
        if self.alternative == "two-sided":
            n = norm.cdf((zr - zrc) * sqrt(n - 3)) + norm.cdf((-zr - zrc) * sqrt(n - 3)) - self.power
        else:
            n = norm.cdf((zr - zrc) * sqrt(n - 3)) - self.power
        return n

    def _get_sig_level(self, sig_level: float) -> float:
        if self.alternative == "less":
            self.effect_size *= -1
        if self.alternative == "two-sided":
            sig_level /= 2
        ttt = t_dist.isf(sig_level, df=self.n - 2)
        zr = atanh(self.effect_size) + self.effect_size / (2 * (self.n - 1))
        rc = sqrt(pow(ttt, 2) / (pow(ttt, 2) + self.n - 2))
        zrc = atanh(rc)
        if self.alternative == "two-sided":
            sig_level = norm.cdf((zr - zrc) * sqrt(self.n - 3)) + norm.cdf((-zr - zrc) * sqrt(self.n - 3)) - self.power
        else:
            sig_level = norm.cdf((zr - zrc) * sqrt(self.n - 3)) - self.power
        return sig_level

    def pwr_test(self) -> dict[str, float | int | str]:
        """Perform power calculation for correlation test."""
        if self.power is None:
            self.power = self._get_power()
        elif self.effect_size is None:
            if self.alternative == "two-sided":
                self.effect_size = brentq(self._get_effect_size, 1e-10, 1 - 1e-10)
            else:
                self.effect_size = brentq(self._get_effect_size, -1 + 1e-10, 1 - 1e-10)
        elif self.n is None:
            self.n = np.ceil(brentq(self._get_n, 4 + 1e-10, 1e09))
        else:
            self.sig_level = brentq(self._get_sig_level, 1e-10, 1 - 1e-10)
        return {
            "effect_size": self.effect_size,
            "n": self.n,
            "sig_level": self.sig_level,
            "power": self.power,
            "alternative": self.alternative,
            "method": self.method,
        }


class pwr_t:
    """Power calculation for t-test.

    Computes power calculations for one-sample, two-sample, or paired t-tests.
    """

    def __init__(
        self,
        n: int | None,
        d: float | None,
        sig_level: float | None,
        power: float | None,
        type: str = "two-sample",
        alternative: str = "two-sided",
    ) -> None:
        """Initialize t-test power calculation.

        Parameters
        ----------
        n : int | None
            Sample size (per group for two-sample, pairs for paired)
        d : float | None
            Effect size (Cohen's d)
        sig_level : float | None
            Significance level (Type I error probability). Must be between 0 and 1
        power : float | None
            Power of test. Must be between 0 and 1
        type : {'two-sample', 'one-sample', 'paired'}, default='two-sample'
            Type of t-test
        alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
            The alternative hypothesis
        """
        self.n = n
        self.d = d
        self.sig_level = sig_level
        self.power = power
        self.type = type.casefold()
        self.alternative = alternative.casefold()
        if self.type == "one-sample":
            self.method = "One Sample"
            self.note: str | None = None
            self.t_sample = 1
        elif self.type == "paired":
            self.method = "Paired Sample"
            self.note = "n is number of pairs"
            self.t_sample = 1
        else:
            self.method = "Two Sample"
            self.note = "n is the number in each group"
            self.t_sample = 2

    def _get_power(self) -> float:
        nu = (self.n - 1) * self.t_sample
        if self.alternative == "two-sided":
            qu = t_dist.isf(self.sig_level / 2, nu)
            power = nct.sf(qu, nu, sqrt(self.n / self.t_sample) * self.d) + nct.cdf(
                -qu, nu, sqrt(self.n / self.t_sample) * self.d
            )
        elif self.alternative == "greater":
            power = nct.sf(
                t_dist.isf(self.sig_level, nu),
                nu,
                sqrt(self.n / self.t_sample) * self.d,
            )
        else:
            power = nct.cdf(
                t_dist.ppf(self.sig_level, nu),
                nu,
                sqrt(self.n / self.t_sample) * self.d,
            )
        return power

    def _get_effect_size(self, effect_size: float) -> float:
        nu = (self.n - 1) * self.t_sample
        if self.alternative == "two-sided":
            qu = t_dist.isf(self.sig_level / 2, nu)
            effect_size = (
                nct.sf(qu, nu, sqrt(self.n / self.t_sample) * effect_size)
                + nct.cdf(-qu, nu, sqrt(self.n / self.t_sample) * effect_size)
                - self.power
            )
        elif self.alternative == "greater":
            effect_size = (
                nct.sf(
                    t_dist.isf(self.sig_level, nu),
                    nu,
                    sqrt(self.n / self.t_sample) * effect_size,
                )
                - self.power
            )
        else:
            effect_size = (
                nct.cdf(
                    t_dist.ppf(self.sig_level, nu),
                    nu,
                    sqrt(self.n / self.t_sample) * effect_size,
                )
                - self.power
            )
        return effect_size

    def _get_n(self, n: int) -> float:
        nu = (n - 1) * self.t_sample
        if self.alternative == "two-sided":
            qu = t_dist.isf(self.sig_level / 2, nu)
            n = (
                nct.sf(qu, nu, sqrt(n / self.t_sample) * self.d)
                + nct.cdf(-qu, nu, sqrt(n / self.t_sample) * self.d)
                - self.power
            )
        elif self.alternative == "greater":
            n = nct.sf(t_dist.isf(self.sig_level, nu), nu, sqrt(n / self.t_sample) * self.d) - self.power
        else:
            n = nct.cdf(t_dist.ppf(self.sig_level, nu), nu, sqrt(n / self.t_sample) * self.d) - self.power
        return n

    def _get_sig_level(self, sig_level: float) -> float:
        nu = (self.n - 1) * self.t_sample
        if self.alternative == "two-sided":
            qu = t_dist.isf(sig_level / 2, nu)
            sig_level = (
                nct.sf(qu, nu, sqrt(self.n / self.t_sample) * self.d)
                + nct.cdf(-qu, nu, sqrt(self.n / self.t_sample) * self.d)
                - self.power
            )
        elif self.alternative == "greater":
            sig_level = nct.sf(t_dist.isf(sig_level, nu), nu, sqrt(self.n / self.t_sample) * self.d) - self.power
        else:
            sig_level = nct.cdf(t_dist.ppf(sig_level, nu), nu, sqrt(self.n / self.t_sample) * self.d) - self.power
        return sig_level

    def pwr_test(self) -> dict[str, int | float | str | None]:
        """Perform power calculation for t-test."""
        if self.power is None:
            self.power = self._get_power()
        elif self.d is None:
            if self.alternative == "two-sided":
                self.d = brentq(self._get_effect_size, 1e-07, 10)
            elif self.alternative == "greater":
                self.d = brentq(self._get_effect_size, -5, 10)
            else:
                self.d = brentq(self._get_effect_size, -10, 5)
        elif self.n is None:
            self.n = np.ceil(brentq(self._get_n, 2 + 1e-10, 1e09))
        else:
            self.sig_level = brentq(self._get_sig_level, 1e-10, 1 - 1e-10)
        if self.note is not None:
            return {
                "n": self.n,
                "effect_size": self.d,
                "sig_level": self.sig_level,
                "power": self.power,
                "alternative": self.alternative,
                "method": "{} t test power calculation".format(self.method),
                "note": self.note,
            }
        else:
            return {
                "n": self.n,
                "effect_size": self.d,
                "sig_level": self.sig_level,
                "power": self.power,
                "alternative": self.alternative,
                "method": "{} t test power calculation".format(self.method),
            }


class pwr_t2n(pwr_2n):
    """Power calculation for two-sample t-test with unequal sample sizes.

    Computes power calculations for two-sample t-tests when the two groups
    have different sample sizes.
    """

    def __init__(
        self,
        d: float | None,
        n1: int | None,
        n2: int | None,
        sig_level: float | None,
        power: float | None,
        alternative: str = "two-sided",
    ) -> None:
        """Initialize two-sample t-test with unequal sizes.

        Parameters
        ----------
        d : float | None
            Effect size (Cohen's d)
        n1 : int | None
            Sample size of the first group
        n2 : int | None
            Sample size of the second group
        sig_level : float | None
            Significance level (Type I error probability). Must be between 0 and 1
        power : float | None
            Power of test. Must be between 0 and 1
        alternative : {'two-sided', 'greater', 'less'}, default='two-sided'
            The alternative hypothesis
        """
        super().__init__(d, n1, n2, sig_level, power, alternative)
        self.method = "T test power calculation"

    def _get_power(self) -> float:
        nu = self.n1 + self.n2 - 2
        if self.alternative == "two-sided":
            qu = t_dist.isf(self.sig_level / 2, nu)
            power = nct.sf(qu, nu, self.effect_size * (1 / sqrt(1 / self.n1 + 1 / self.n2))) + nct.cdf(
                -qu, nu, self.effect_size * (1 / sqrt(1 / self.n1 + 1 / self.n2))
            )
        elif self.alternative == "greater":
            power = nct.sf(
                t_dist.isf(self.sig_level, nu),
                nu,
                self.effect_size * (1 / sqrt(1 / self.n1 + 1 / self.n2)),
            )
        else:
            power = nct.cdf(
                t_dist.ppf(self.sig_level, nu),
                nu,
                self.effect_size * (1 / sqrt(1 / self.n1 + 1 / self.n2)),
            )
        return power

    def _get_effect_size(self, effect_size: float) -> float:
        nu = self.n1 + self.n2 - 2
        if self.alternative == "two-sided":
            qu = t_dist.isf(self.sig_level / 2, nu)
            effect_size = (
                nct.sf(qu, nu, effect_size * (1 / sqrt(1 / self.n1 + 1 / self.n2)))
                + nct.cdf(-qu, nu, effect_size * (1 / sqrt(1 / self.n1 + 1 / self.n2)))
                - self.power
            )
        elif self.alternative == "greater":
            effect_size = (
                nct.sf(
                    t_dist.isf(self.sig_level, nu),
                    nu,
                    effect_size * (1 / sqrt(1 / self.n1 + 1 / self.n2)),
                )
                - self.power
            )
        else:
            effect_size = (
                nct.cdf(
                    t_dist.ppf(self.sig_level, nu),
                    nu,
                    effect_size * (1 / sqrt(1 / self.n1 + 1 / self.n2)),
                )
                - self.power
            )
        return effect_size

    def _get_n1(self, n1: int) -> float:
        nu = n1 + self.n2 - 2
        if self.alternative == "two-sided":
            qu = t_dist.isf(self.sig_level / 2, nu)
            n1 = (
                nct.sf(qu, nu, self.effect_size * (1 / sqrt(1 / n1 + 1 / self.n2)))
                + nct.cdf(-qu, nu, self.effect_size * (1 / sqrt(1 / n1 + 1 / self.n2)))
                - self.power
            )
        elif self.alternative == "greater":
            n1 = (
                nct.sf(
                    t_dist.isf(self.sig_level, nu),
                    nu,
                    self.effect_size * (1 / sqrt(1 / n1 + 1 / self.n2)),
                )
                - self.power
            )
        else:
            n1 = (
                nct.cdf(
                    t_dist.ppf(self.sig_level, nu),
                    nu,
                    self.effect_size * (1 / sqrt(1 / n1 + 1 / self.n2)),
                )
                - self.power
            )
        return n1

    def _get_n2(self, n2: int) -> float:
        nu = self.n1 + n2 - 2
        if self.alternative == "two-sided":
            qu = t_dist.isf(self.sig_level / 2, nu)
            n2 = (
                nct.sf(qu, nu, self.effect_size * (1 / sqrt(1 / self.n1 + 1 / n2)))
                + nct.cdf(-qu, nu, self.effect_size * (1 / sqrt(1 / self.n1 + 1 / n2)))
                - self.power
            )
        elif self.alternative == "greater":
            n2 = (
                nct.sf(
                    t_dist.isf(self.sig_level, nu),
                    nu,
                    self.effect_size * (1 / sqrt(1 / self.n1 + 1 / n2)),
                )
                - self.power
            )
        else:
            n2 = (
                nct.cdf(
                    t_dist.ppf(self.sig_level, nu),
                    nu,
                    self.effect_size * (1 / sqrt(1 / self.n1 + 1 / n2)),
                )
                - self.power
            )
        return n2

    def _get_sig_level(self, sig_level: float) -> float:
        nu = self.n1 + self.n2 - 2
        if self.alternative == "two-sided":
            qu = t_dist.isf(sig_level / 2, nu)
            sig_level = (
                nct.sf(qu, nu, self.effect_size * (1 / sqrt(1 / self.n1 + 1 / self.n2)))
                + nct.cdf(-qu, nu, self.effect_size * (1 / sqrt(1 / self.n1 + 1 / self.n2)))
                - self.power
            )
        elif self.alternative == "greater":
            sig_level = (
                nct.sf(
                    t_dist.isf(sig_level, nu),
                    nu,
                    self.effect_size * (1 / sqrt(1 / self.n1 + 1 / self.n2)),
                )
                - self.power
            )
        else:
            sig_level = (
                nct.cdf(
                    t_dist.ppf(sig_level, nu),
                    nu,
                    self.effect_size * (1 / sqrt(1 / self.n1 + 1 / self.n2)),
                )
                - self.power
            )
        return sig_level

    def pwr_test(self) -> dict[str, int | float | str | None]:
        """Perform power calculation for two-sample t-test with unequal sizes."""
        if self.power is None:
            self.power = self._get_power()
        elif self.effect_size is None:
            if self.alternative == "two-sided":
                self.effect_size = brentq(self._get_effect_size, 1e-10, 10)
            elif self.alternative == "greater":
                self.effect_size = brentq(self._get_effect_size, -1, 10)
            else:
                self.effect_size = brentq(self._get_effect_size, -10, 1)
        elif self.n1 is None:
            self.n1 = np.ceil(brentq(self._get_n1, 2 + 1e-10, 500))
        elif self.n2 is None:
            self.n2 = np.ceil(brentq(self._get_n2, 2 + 1e-10, 500))
        else:
            self.sig_level = brentq(self._get_sig_level, 1e-10, 1 - 1e-10)
        if self.note is not None:
            return {
                "effect_size": self.effect_size,
                "n1": self.n1,
                "n2": self.n2,
                "sig_level": self.sig_level,
                "power": self.power,
                "alternative": self.alternative,
                "method": self.method,
                "note": "Different sample sizes",
            }
        else:
            return {
                "effect_size": self.effect_size,
                "n1": self.n1,
                "n2": self.n2,
                "sig_level": self.sig_level,
                "power": self.power,
                "alternative": self.alternative,
                "method": self.method,
            }
