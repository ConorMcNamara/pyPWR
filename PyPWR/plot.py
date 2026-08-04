"""Power-versus-sample-size plots for :mod:`PyPWR` results.

This module provides a Plotly equivalent of R ``pwr``'s ``plot.power.htest``.
Given a result dictionary returned by any of the ``pwr_*_test`` functions, it
sweeps the sample size, recomputes power at each point, and renders an
interactive power curve with the optimal sample size highlighted.

The heavy lifting (the power computations) is delegated back to the public
``pwr_*_test`` wrappers, so the curve is always consistent with the value the
user originally solved for.
"""

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from PyPWR.pwr_tests import (
    pwr_2p2n_test,
    pwr_2p_test,
    pwr_anova_test,
    pwr_chisq_test,
    pwr_norm_test,
    pwr_p_test,
    pwr_r_test,
    pwr_t2n_test,
    pwr_t_test,
)

if TYPE_CHECKING:
    import plotly.graph_objects as go

# Number of intervals along the sample-size axis (matches R's ``breaks``).
_BREAKS = 20


@dataclass
class PowerCurve:
    """Data backing a power-versus-sample-size plot.

    Attributes
    ----------
    sample_sizes : list[float]
        Sample sizes swept along the x-axis.
    powers : list[float]
        Power at each sample size; ``nan`` where the combination is invalid
        (e.g. a two-sample split leaving fewer than two observations in a group).
    optimal_n : int
        The originally solved sample size, marked with a vertical reference line.
    step : float
        Spacing between successive sample sizes.
    title : str
        Plot title (the test's ``method`` string).
    legend_lines : list[str]
        Parameter summary lines (effect size, alpha, ...).
    optimal_label : str
        Annotation text describing the optimal sample size.
    power_at_optimal : float
        Power reported by the original result, used for annotation placement.
    """

    sample_sizes: list[float]
    powers: list[float]
    optimal_n: int
    step: float
    title: str
    legend_lines: list[str]
    optimal_label: str
    power_at_optimal: float


@dataclass(frozen=True)
class _Spec:
    """Dispatch configuration for a single plottable test."""

    fn: Callable[..., dict[str, Any]]
    effect_kwarg: str
    two_sample: bool
    uses_alternative: bool
    effect_label: str


# Effect-size label as shown in the R legend, keyed by internal test kind.
_SPECS: dict[str, _Spec] = {
    "t": _Spec(pwr_t_test, "d", False, True, "effect size d"),
    "t2n": _Spec(pwr_t2n_test, "d", True, True, "effect size d"),
    "2p": _Spec(pwr_2p_test, "h", False, True, "effect size h"),
    "2p2n": _Spec(pwr_2p2n_test, "h", True, True, "effect size h"),
    "p": _Spec(pwr_p_test, "h", False, True, "effect size h"),
    "r": _Spec(pwr_r_test, "r", False, True, "r"),
    "norm": _Spec(pwr_norm_test, "d", False, True, "effect size d"),
    "anova": _Spec(pwr_anova_test, "f", False, False, "effect size f"),
    "chisq": _Spec(pwr_chisq_test, "w", False, False, "effect size w"),
}


def _identify(result: dict[str, Any]) -> str:
    """Map a result dictionary to an internal test kind.

    Parameters
    ----------
    result : dict[str, Any]
        A dictionary returned by one of the ``pwr_*_test`` functions.

    Returns
    -------
    str
        The internal kind key used to index :data:`_SPECS`.

    Raises
    ------
    ValueError
        If the result is missing a ``method`` field or corresponds to a test
        that has no single sample-size axis to plot (e.g. ``pwr_f2_test``).
    """
    method_raw = result.get("method")
    if not isinstance(method_raw, str):
        raise ValueError("result is missing a 'method' field; pass a pwr_*_test result dictionary")
    method = method_raw.casefold()
    two_sample = "n1" in result and "n2" in result

    if "multiple regression" in method:
        raise ValueError(
            "plot_power does not support pwr_f2_test results: the F2 test has no single "
            "sample-size axis (it is parameterized by numerator/denominator degrees of freedom)"
        )
    if two_sample:
        return "2p2n" if "proportion" in method else "t2n"
    if "difference of proportion" in method:
        return "2p"
    if "proportion" in method:
        return "p"
    if "analysis of variance" in method:
        return "anova"
    if "chi squared" in method:
        return "chisq"
    if "normal distribution" in method:
        return "norm"
    if "correlation" in method:
        return "r"
    if "t test" in method:
        return "t"
    raise ValueError(f"plot_power does not support results with method {method_raw!r}")


def _t_type(method: str) -> str:
    """Recover the ``test_type`` argument from a t-test ``method`` string."""
    lowered = method.casefold()
    if "one sample" in lowered:
        return "one-sample"
    if "two sample" in lowered:
        return "two-sample"
    return "paired"


def _axis(n_total: float) -> tuple[list[float], float]:
    """Build the sample-size sweep from 10 to ``max(n*1.5, n+30)``."""
    n_upper = max(n_total * 1.5, n_total + 30)
    step = (n_upper - 10) / _BREAKS
    sizes = [10 + i * step for i in range(_BREAKS + 1)]
    return sizes, step


def power_curve(result: dict[str, Any]) -> PowerCurve:
    """Compute the power-versus-sample-size curve for a test result.

    This is the pure, plotting-library-independent core of :func:`plot_power`.
    It recomputes power across a range of sample sizes by calling back into the
    matching ``pwr_*_test`` function, holding the effect size, significance
    level, and alternative fixed.

    Parameters
    ----------
    result : dict[str, Any]
        A dictionary returned by one of the ``pwr_*_test`` functions.

    Returns
    -------
    PowerCurve
        The swept sample sizes, corresponding powers, and annotation metadata.

    Raises
    ------
    ValueError
        If the result is not a supported, plottable test (see :func:`_identify`).
    """
    kind = _identify(result)
    spec = _SPECS[kind]
    effect = result["effect_size"]
    sig_level = result["sig_level"]

    base: dict[str, Any] = {spec.effect_kwarg: effect, "sig_level": sig_level, "print_pretty": False}
    if spec.uses_alternative:
        base["alternative"] = result["alternative"]
    if kind == "anova":
        base["k"] = result["k"]
    elif kind == "chisq":
        base["df"] = result["df"]
    elif kind == "t":
        base["test_type"] = _t_type(result["method"])

    legend_lines: list[str] = []
    if spec.uses_alternative:
        legend_lines.append(f"tails = {result['alternative']}")
    if kind == "anova":
        legend_lines.append(f"groups k = {result['k']}")
    legend_lines.append(f"{spec.effect_label} = {effect}")
    if kind == "chisq":
        legend_lines.append(f"df = {result['df']}")

    powers: list[float] = []
    if spec.two_sample:
        n1, n2 = int(result["n1"]), int(result["n2"])
        total = n1 + n2
        n_rel = n1 / total
        sizes, step = _axis(total)
        for size in sizes:
            group1 = math.ceil(size * n_rel)
            group2 = size - group1
            if group1 < 2 or group2 < 2:
                powers.append(math.nan)
                continue
            powers.append(float(spec.fn(**base, n1=group1, n2=group2)["power"]))
        optimal_n = total
        optimal_label = f"optimal sample size<br>n = {n1} + {n2} = {total}"
        legend_lines.append(f"n1/n2 = {round(n_rel, 2)}")
    else:
        n = result["n"]
        sizes, step = _axis(n)
        for size in sizes:
            powers.append(float(spec.fn(**base, n=size)["power"]))
        optimal_n = math.ceil(n)
        optimal_label = f"optimal sample size<br>n = {optimal_n}"
        note = result.get("note")
        if note:
            optimal_label += f"<br>{note}"

    legend_lines.append(f"alpha = {sig_level}")

    return PowerCurve(
        sample_sizes=sizes,
        powers=powers,
        optimal_n=optimal_n,
        step=step,
        title=str(result["method"]),
        legend_lines=legend_lines,
        optimal_label=optimal_label,
        power_at_optimal=float(result["power"]),
    )


def plot_power(
    result: dict[str, Any],
    *,
    title: str | None = None,
    xlab: str = "sample size",
    ylab: str = "test power = 1 - β",
) -> "go.Figure":  # type: ignore[name-defined]  # plotly ships no type stubs
    """Plot statistical power as a function of sample size.

    A Plotly equivalent of R ``pwr``'s ``plot.power.htest``. The power curve is
    drawn as a red line with markers, the originally solved sample size is marked
    with a dashed vertical line, and the test parameters and optimal sample size
    are shown as annotations. Hovering a marker reveals the exact sample size and
    power.

    Parameters
    ----------
    result : dict[str, Any]
        A dictionary returned by one of the ``pwr_*_test`` functions. The
        ``pwr_f2_test`` result is not supported (it has no sample-size axis).
    title : str | None, default=None
        Plot title. Defaults to the test's ``method`` string.
    xlab : str, default="sample size"
        Label for the x-axis.
    ylab : str, default="test power = 1 - β"
        Label for the y-axis.

    Returns
    -------
    plotly.graph_objects.Figure
        The interactive power-versus-sample-size figure.

    Raises
    ------
    ImportError
        If Plotly is not installed. Install it with ``pip install "pwr-py[plot]"``.
    ValueError
        If ``result`` is not a supported, plottable test.

    Examples
    --------
    >>> from PyPWR import pwr_t_test, plot_power
    >>> result = pwr_t_test(d=0.5, power=0.8, sig_level=0.05, test_type="two-sample")
    >>> fig = plot_power(result)  # doctest: +SKIP
    >>> fig.show()  # doctest: +SKIP
    """
    try:
        import plotly.graph_objects as go
    except ImportError as exc:  # pragma: no cover - exercised only without plotly
        raise ImportError('plot_power requires plotly. Install it with: pip install "pwr-py[plot]"') from exc

    curve = power_curve(result)

    valid_powers = [p for p in curve.powers if not math.isnan(p)]
    min_power = min(valid_powers) if valid_powers else 0.0

    fig = go.Figure()  # type: ignore[attr-defined]  # plotly ships no type stubs
    fig.add_trace(
        go.Scatter(  # type: ignore[attr-defined]  # plotly ships no type stubs
            x=curve.sample_sizes,
            y=curve.powers,
            mode="lines+markers",
            line={"color": "red", "width": 1},
            marker={"color": "black", "size": 6},
            connectgaps=False,
            name="power",
            hovertemplate="n = %{x:.0f}<br>power = %{y:.1%}<extra></extra>",
        )
    )
    fig.add_vline(x=curve.optimal_n, line={"color": "darkblue", "dash": "dot", "width": 1.5})

    # Legend block: top-left when the curve dips low, bottom-left otherwise, so it
    # stays clear of the line (mirrors the R placement heuristic).
    legend_top = min_power < 0.6
    fig.add_annotation(
        x=curve.sample_sizes[0],
        y=1.0 if legend_top else 0.0,
        xref="x",
        yref="y",
        text="<br>".join(curve.legend_lines),
        showarrow=False,
        align="left",
        xanchor="left",
        yanchor="top" if legend_top else "bottom",
        font={"size": 12},
    )
    # Optimal-n callout: opposite vertical anchor from the legend based on the
    # solved power, again to avoid overlapping the curve.
    optimal_top = curve.power_at_optimal < 0.5
    fig.add_annotation(
        x=curve.optimal_n + curve.step,
        y=1.0 if optimal_top else 0.0,
        xref="x",
        yref="y",
        text=curve.optimal_label,
        showarrow=False,
        align="left",
        xanchor="left",
        yanchor="top" if optimal_top else "bottom",
        font={"size": 12, "color": "darkblue"},
    )

    fig.update_xaxes(title_text=xlab)
    fig.update_yaxes(title_text=ylab, range=[0, 1], tickformat=".0%")
    fig.update_layout(title=title if title is not None else curve.title, template="plotly_white", showlegend=False)
    return fig
