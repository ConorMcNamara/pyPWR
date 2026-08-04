"""Tests for the power-versus-sample-size plotting module."""

import math
from itertools import pairwise

import pytest

from PyPWR import (
    plot_power,
    power_curve,
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

# One representative result per plottable test kind, each solved for sample size.
_SINGLE_N_RESULTS = [
    pwr_t_test(d=0.5, power=0.8, sig_level=0.05, test_type="two-sample", print_pretty=False),
    pwr_t_test(d=0.5, power=0.8, sig_level=0.05, test_type="one-sample", print_pretty=False),
    pwr_2p_test(h=0.3, power=0.8, sig_level=0.05, print_pretty=False),
    pwr_p_test(h=0.3, power=0.8, sig_level=0.05, print_pretty=False),
    pwr_r_test(r=0.3, power=0.8, sig_level=0.05, print_pretty=False),
    pwr_norm_test(d=0.3, power=0.8, sig_level=0.05, print_pretty=False),
    pwr_anova_test(k=4, f=0.25, power=0.8, sig_level=0.05, print_pretty=False),
    pwr_chisq_test(w=0.3, df=4, power=0.8, sig_level=0.05, print_pretty=False),
]

_TWO_SAMPLE_RESULTS = [
    pwr_t2n_test(n1=40, n2=60, d=0.5, sig_level=0.05, print_pretty=False),
    pwr_2p2n_test(h=0.3, n1=80, n2=120, sig_level=0.05, print_pretty=False),
]


class TestPowerCurve:
    """Tests for the pure ``power_curve`` data generator."""

    @pytest.mark.parametrize("result", _SINGLE_N_RESULTS + _TWO_SAMPLE_RESULTS)
    def test_curve_shape(self, result):
        curve = power_curve(result)
        # 20 breaks -> 21 points.
        assert len(curve.sample_sizes) == 21
        assert len(curve.powers) == 21
        # Every non-nan power is a valid probability.
        for p in curve.powers:
            assert math.isnan(p) or 0.0 <= p <= 1.0
        # Sample sizes are strictly increasing starting at 10.
        assert curve.sample_sizes[0] == 10
        assert all(b > a for a, b in pairwise(curve.sample_sizes))

    @pytest.mark.parametrize("result", _SINGLE_N_RESULTS)
    def test_optimal_n_matches_single(self, result):
        curve = power_curve(result)
        assert curve.optimal_n == math.ceil(result["n"])

    @pytest.mark.parametrize("result", _TWO_SAMPLE_RESULTS)
    def test_optimal_n_matches_two_sample(self, result):
        curve = power_curve(result)
        assert curve.optimal_n == result["n1"] + result["n2"]

    def test_power_increases_with_sample_size(self):
        # Power is monotonically non-decreasing in n for a fixed effect/alpha.
        curve = power_curve(pwr_t_test(d=0.5, power=0.8, sig_level=0.05, test_type="two-sample", print_pretty=False))
        powers = [p for p in curve.powers if not math.isnan(p)]
        assert all(b >= a - 1e-9 for a, b in pairwise(powers))

    def test_legend_contains_effect_and_alpha(self):
        curve = power_curve(pwr_anova_test(k=4, f=0.25, power=0.8, sig_level=0.05, print_pretty=False))
        joined = "\n".join(curve.legend_lines)
        assert "groups k = 4" in joined
        assert "effect size f = 0.25" in joined
        assert "alpha = 0.05" in joined

    def test_f2_result_is_rejected(self):
        result = pwr_f2_test(u=5, v=89, f2=0.1, sig_level=0.05, print_pretty=False)
        with pytest.raises(ValueError, match="pwr_f2_test"):
            power_curve(result)

    def test_missing_method_is_rejected(self):
        with pytest.raises(ValueError, match="method"):
            power_curve({"effect_size": 0.5, "n": 30, "sig_level": 0.05, "power": 0.8})


class TestPlotPower:
    """Tests for the Plotly rendering wrapper."""

    def test_returns_figure(self):
        go = pytest.importorskip("plotly.graph_objects")
        result = pwr_t_test(d=0.5, power=0.8, sig_level=0.05, test_type="two-sample", print_pretty=False)
        fig = plot_power(result)
        assert isinstance(fig, go.Figure)

    def test_figure_has_curve_and_optimal_line(self):
        pytest.importorskip("plotly")
        result = pwr_2p2n_test(h=0.3, n1=80, n2=120, sig_level=0.05, print_pretty=False)
        fig = plot_power(result)
        # One scatter trace for the power curve.
        assert len(fig.data) == 1
        assert fig.data[0].mode == "lines+markers"
        # A dashed vertical line (a shape) marks the optimal sample size.
        assert any(shape.type == "line" for shape in fig.layout.shapes)
        # Two annotations: parameter legend and optimal-n callout.
        assert len(fig.layout.annotations) == 2

    def test_title_override(self):
        pytest.importorskip("plotly")
        result = pwr_r_test(r=0.3, power=0.8, sig_level=0.05, print_pretty=False)
        fig = plot_power(result, title="Custom Title")
        assert fig.layout.title.text == "Custom Title"
