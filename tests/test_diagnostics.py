"""Known functional laws and independent dimer data constrain the diagnostics."""
import json

import numpy as np
import pytest

from rvb.diagnostics import diagnose_decay
from rvb.orthogonal import square_plane_benchmark


def model(result, name, window_index=0):
    return result["windows"][window_index]["models"][name]


def test_power_and_exponential_have_distinct_effective_signatures():
    r = np.geomspace(2., 80., 24)
    power = diagnose_decay(r, -.7 * r**-1.3, windows=[(2, 80)], floor=0)
    assert model(power, "power_law")["p"] == pytest.approx(1.3, abs=1e-12)
    assert model(power, "power_law")["log_rms"] < 1e-12
    assert [v["p_eff"] for v in power["local"]] == pytest.approx([1.3] * 23)
    assert power["local"][-1]["inverse_xi_eff"] < power["local"][0]["inverse_xi_eff"]
    exponential = diagnose_decay(r, .2*np.exp(-r/7), windows=[(2, 80)], floor=0)
    assert model(exponential, "exponential")["xi"] == pytest.approx(7, abs=1e-12)
    assert [v["inverse_xi_eff"] for v in exponential["local"]] == pytest.approx([1/7] * 23)
    assert exponential["local"][-1]["p_eff"] > exponential["local"][0]["p_eff"]
    for value in exponential["local"]:
        assert value["p_eff"] == pytest.approx(value["inverse_xi_eff"] * value["r_logarithmic"])
    assert power["status"] == exponential["status"] == "descriptive_only"


def test_mixed_cutoff_recovers_parameters_but_does_not_classify_bulk():
    r = np.geomspace(2., 200., 40)
    result = diagnose_decay(r, 3*r**-1.2*np.exp(-r/60), windows=[(2, 200)], floor=0)
    fit = model(result, "power_exponential")
    assert fit["p"] == pytest.approx(1.2, abs=1e-12)
    assert fit["xi"] == pytest.approx(60, abs=1e-10)
    assert fit["physical_decay"]
    assert fit["log_rms"] < 1e-12
    assert model(result, "power_law")["log_rms"] > .1
    assert model(result, "exponential")["log_rms"] > .1
    assert result["status"] == "descriptive_only"


def test_amplitude_and_distance_units_do_not_change_physics_or_condition():
    r = np.geomspace(2., 70., 32)
    c = r**-.8*np.exp(-r/20)
    original = diagnose_decay(r, c, windows=[(2, 70)], floor=0)
    changed = diagnose_decay(r*1e8, c*1e-20, windows=[(2e8, 70e8)], floor=0)
    first, second = model(original, "power_exponential"), model(changed, "power_exponential")
    assert second["p"] == pytest.approx(first["p"], rel=1e-10)
    assert second["xi"] / 1e8 == pytest.approx(first["xi"], rel=1e-10)
    assert second["condition"] == pytest.approx(first["condition"], rel=1e-10)
    assert [v["p_eff"] for v in changed["local"]] == pytest.approx([v["p_eff"] for v in original["local"]])


def test_zeros_floor_sign_and_error_split_without_bridging():
    r = np.arange(1., 13.)
    c = np.array([1., .5, 0., .2, .1, -1., -.5, -1e-15, -.2, -.1, -.05, -.025])
    error = np.zeros(12)
    error[10] = .02
    result = diagnose_decay(r, c, windows=[(1, 12)], absolute_error=error)
    edges = {(v["r_left"], v["r_right"]) for v in result["local"]}
    assert edges == {(1., 2.), (4., 5.), (6., 7.), (9., 10.)}
    assert result["points"][2]["reason"] == "zero_unresolved"
    assert result["points"][7]["reason"] == "below_absolute_floor"
    assert result["points"][10]["reason"] == "below_signal_to_error"
    assert result["status"] == "insufficient_resolved_data"
    assert all(entry["npoints"] <= 2 for entry in result["windows"])


def test_explicit_parity_handles_staggering_but_preserves_breaks():
    r = np.arange(1., 15.)
    c = (-1.)**r / r**2
    unbranched = diagnose_decay(r, c, windows=[(1, 14)], floor=0)
    assert not unbranched["local"]
    branched = diagnose_decay(r, c, windows=[(1, 14)], branch=r.astype(int)%2, floor=0)
    assert len(branched["segments"]) == 2
    assert all(v["p_eff"] == pytest.approx(2) for v in branched["local"])
    c[5] = 0  # r=6 breaks only the even branch.
    split = diagnose_decay(r, c, windows=[(1, 14)], branch=r.astype(int)%2, floor=0)
    assert (4., 8.) not in {(v["r_left"], v["r_right"]) for v in split["local"]}
    assert (5., 7.) in {(v["r_left"], v["r_right"]) for v in split["local"]}


def test_short_span_flags_unidentifiable_mixed_parameters():
    r = 1 + np.arange(10)*1e-8
    c = r**-1.2*np.exp(-r/2)
    result = diagnose_decay(r, c, windows=[(r[0], r[-1])], floor=0)
    fit = model(result, "power_exponential")
    assert fit["status"] in {"ill_conditioned", "rank_deficient"}
    assert fit["condition"] > 1e6
    assert model(result, "power_law")["status"] == "ok"


def test_exact_zero_and_plateau_have_no_invented_correlation_length():
    r = np.arange(1., 9.)
    zeros = diagnose_decay(r, np.zeros(8), windows=[(1, 8)], exact_zero=True)
    assert zeros["status"] == "insufficient_resolved_data"
    assert zeros["local"] == zeros["windows"] == []
    assert all(point["reason"] == "exact_zero_declared" for point in zeros["points"])
    plateau = diagnose_decay(r, np.full(8, .02), windows=[(1, 8)])
    fit = model(plateau, "exponential")
    assert fit["inverse_xi"] == pytest.approx(0, abs=1e-14)
    assert not fit["physical_decay"]
    assert fit["xi"] is None
    assert plateau["status"] == "descriptive_only"
    json.dumps(zeros, allow_nan=False)
    json.dumps(plateau, allow_nan=False)


def test_independent_square_dimer_integral_has_known_power_tail():
    reference = square_plane_benchmark(128)
    r = np.arange(1., 129.)
    c = np.asarray(reference["connected_correlation"])[1:]
    result = diagnose_decay(r, c, windows=[(40, 128)], branch=r.astype(int)%2)
    fits = [entry["models"]["power_law"] for entry in result["windows"]]
    assert len(fits) == 2
    assert all(fit["p"] == pytest.approx(2, abs=.01) for fit in fits)
    assert result["status"] == "descriptive_only"


def test_independent_transverse_bonds_have_different_parity_exponents():
    reference = square_plane_benchmark(128, orientation="vertical")
    r = np.arange(1., 129.)
    c = np.asarray(reference["connected_correlation"])[1:]
    result = diagnose_decay(r, c, windows=[(40, 128)], branch=r.astype(int)%2)
    powers = {entry["branch"]: entry["models"]["power_law"]["p"]
              for entry in result["windows"]}
    # A single sign-agnostic fit would erase this physical cancellation.
    assert powers["0"] == pytest.approx(4, abs=.02)
    assert powers["1"] == pytest.approx(2, abs=.01)


@pytest.mark.parametrize("kwargs", [
    {"distance": [0, 1, 2]}, {"distance": [1, 1, 2]}, {"distance": [3, 2, 1]},
    {"distance": [1, 2, np.inf]}, {"distance": [1, 2, 3j]},
    {"correlation": [1, np.nan, .1]}, {"correlation": [1, .5, .1j]},
    {"correlation": [.1]}, {"floor": -1}, {"floor": np.nan},
    {"absolute_error": -.1}, {"absolute_error": [0, np.inf, 0]},
    {"absolute_error": [0, 0]}, {"min_signal_to_error": 0},
    {"branch": [0, 1]}, {"windows": [(3, 1)]}, {"windows": []},
])
def test_invalid_input_is_rejected(kwargs):
    args = {"distance": [1, 2, 3], "correlation": [1, .5, .1], "windows": [(1, 3)]}
    args.update(kwargs)
    with pytest.raises(ValueError):
        diagnose_decay(**args)
