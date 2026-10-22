"""Descriptive decay diagnostics for signed, finite-accuracy correlations.

These calculations do not identify a thermodynamic phase. Establishing bulk
behavior also requires convergence in contraction dimension, circumference,
iteration tolerance, and the subtraction used for connected correlations.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Iterable

import numpy as np


_CAVEATS = [
    "All fits describe the supplied finite-distance data; no bulk classification or confidence interval is inferred.",
    "Compare windows and independently converged chi, circumference, and tolerance values before interpreting a tail.",
    "The absolute floor and supplied errors must include relevant subtraction error; chi differences are not rigorous error bounds.",
    "Residuals are unweighted in log absolute correlation. A smaller residual for the three-parameter model alone does not establish that model.",
    "A plateau can reflect order, a sector choice, finite geometry, or a subtraction offset. No fitted constant offset is introduced.",
    "Zeros below numerical resolution are not evidence of exactly vanishing correlations. Exact zeros require independent justification.",
    "Oscillatory branches must be specified explicitly. Nodes and sign changes terminate a segment; this module does not infer an oscillation period.",
]


def _real_array(values, name: str) -> np.ndarray:
    raw = np.asarray(values)
    if np.iscomplexobj(raw):
        raise ValueError(f"{name} must be real")
    try:
        result = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain real numbers") from exc
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result


def _scalar(value, name: str, *, positive=False) -> float:
    array = _real_array(value, name)
    if array.ndim != 0:
        raise ValueError(f"{name} must be a scalar")
    result = float(array)
    if result < 0 or (positive and result == 0):
        raise ValueError(f"{name} must be {'positive' if positive else 'nonnegative'}")
    return result


def _fit_model(r: np.ndarray, y: np.ndarray, model: str) -> dict:
    mixed = model == "power_exponential"
    minimum = 6 if mixed else 4
    result = {"status": "insufficient_points", "minimum_points": minimum,
              "p": None, "inverse_xi": None, "xi": None, "log_amplitude": None,
              "log_rms": None, "condition": None, "rank": None,
              "physical_decay": None}
    if len(r) < minimum:
        return {**result, "reason": f"requires at least {minimum} resolved points in one constant-sign segment"}
    features = ([np.log(r), r] if mixed else
                [np.log(r)] if model == "power_law" else [r])
    means = np.array([np.sum(feature / len(feature)) for feature in features])
    # Center and scale before solving: raw condition numbers depend on the
    # distance unit and obscure collinearity between log r and r.
    centered = np.column_stack([feature - mean for feature, mean in zip(features, means)])
    maxima = np.max(np.abs(centered), axis=0)
    if np.any(maxima == 0):
        return {**result, "status": "rank_deficient", "reason": "fit coordinates have no resolved span"}
    # RMS scaling avoids squaring very large distances before normalizing.
    scales = maxima * np.sqrt(np.mean((centered / maxima) ** 2, axis=0))
    if np.any(scales == 0):
        return {**result, "status": "rank_deficient", "reason": "fit coordinates have no resolved span"}
    design = np.column_stack([np.ones(len(r)), centered / scales])
    y_mean = float(y.mean())
    coefficients, _, rank, singular = np.linalg.lstsq(design, y - y_mean, rcond=None)
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else None
    with np.errstate(over="ignore", invalid="ignore"):
        slopes = coefficients[1:] / scales
        intercept = float(y_mean + coefficients[0] - slopes @ means)
    if not np.all(np.isfinite(slopes)) or not np.isfinite(intercept):
        return {**result, "status": "numerical_range", "rank": int(rank),
                "condition": condition, "reason": "fit parameters exceed floating-point range"}
    if mixed:
        p, inverse_xi = -float(slopes[0]), -float(slopes[1])
    elif model == "power_law":
        p, inverse_xi = -float(slopes[0]), 0.0
    else:
        p, inverse_xi = 0.0, -float(slopes[0])
    residual = y - y_mean - design @ coefficients
    status = ("rank_deficient" if rank < design.shape[1] else
              "ill_conditioned" if condition is None or condition > 1e6 else "ok")
    # r*d[-log C]/dr is dimensionless. The threshold avoids treating numerical
    # roundoff in a flat curve as physical decay; it is not an uncertainty.
    physical_decay = bool(np.all(p + inverse_xi * r > 1e-10))
    xi = float(1 / inverse_xi) if inverse_xi > 0 else None
    if xi is not None and not np.isfinite(xi):
        xi = None
    result.update(status=status, p=p, inverse_xi=inverse_xi, xi=xi,
                  log_amplitude=intercept, log_rms=float(np.sqrt(np.mean(residual**2))),
                  condition=condition, rank=int(rank), physical_decay=physical_decay)
    if status != "ok":
        result["reason"] = "parameters are not numerically identifiable in this window"
    return result


def diagnose_decay(distance, correlation, *, windows: Iterable,
                   floor=1e-13, branch=None, absolute_error=None,
                   min_signal_to_error=5., exact_zero=False) -> dict:
    """Return JSON-friendly local slopes and descriptive log-space fits.

    Distances must be positive and strictly increasing. ``branch`` is an
    optional label per point, e.g. ``distance % 2`` for a staggered correlator.
    Branch labels are returned as strings. A masked/zero point or a sign change
    splits its branch, so neither local slopes nor fits bridge that break.

    ``absolute_error`` accepts a scalar or one nonnegative value per point.
    A point is resolved only when ``abs(C)>floor`` and, if errors are supplied,
    ``abs(C)>=min_signal_to_error*absolute_error``. This is a resolution gate,
    not a statistical weighting or an inferred error model. ``exact_zero=True``
    records the caller's independent justification for zero-valued entries; it
    does not establish that fact or assign a correlation length to those zeros.

    Windows are inclusive ``(minimum_distance, maximum_distance)`` pairs.
    Each constant-sign segment is fitted separately within a window. Pure
    power/exponential fits need four points; power times exponential needs six.
    The latter uses ``log|C|=a-p*log(r)-inverse_xi*r``. ``condition`` refers to
    the centered, scaled design, and values above 1e6 flag unstable parameters.
    ``physical_decay`` only says that the fitted curve decreases over the used
    interval. It does not validate the model or establish asymptotic behavior.
    """
    r = _real_array(distance, "distance")
    c = _real_array(correlation, "correlation")
    if r.ndim != 1 or c.ndim != 1 or len(r) != len(c) or not len(r):
        raise ValueError("distance and correlation must be nonempty one-dimensional arrays of equal length")
    if np.any(r <= 0) or np.any(np.diff(r) <= 0):
        raise ValueError("distance must be positive and strictly increasing")
    floor = _scalar(floor, "floor")
    min_signal_to_error = _scalar(min_signal_to_error, "min_signal_to_error", positive=True)
    if not isinstance(exact_zero, (bool, np.bool_)):
        raise ValueError("exact_zero must be a boolean")
    errors = None
    if absolute_error is not None:
        errors = _real_array(absolute_error, "absolute_error")
        if errors.ndim == 0:
            errors = np.full(len(r), float(errors))
        if errors.ndim != 1 or len(errors) != len(r) or np.any(errors < 0):
            raise ValueError("absolute_error must be nonnegative and scalar or match the data length")
    if branch is None:
        labels = ["all"] * len(r)
    else:
        raw_labels = np.asarray(branch)
        if raw_labels.ndim != 1 or len(raw_labels) != len(r):
            raise ValueError("branch must provide one label per point")
        if np.iscomplexobj(raw_labels):
            raise ValueError("branch labels must be real or strings")
        if np.issubdtype(raw_labels.dtype, np.number) and not np.all(np.isfinite(raw_labels)):
            raise ValueError("numeric branch labels must be finite")
        labels = [str(label) for label in raw_labels]
    checked_windows = []
    try:
        requested_windows = list(windows)
    except TypeError as exc:
        raise ValueError("windows must be an iterable of distance pairs") from exc
    for window in requested_windows:
        bounds = _real_array(window, "window")
        if bounds.shape != (2,) or bounds[0] <= 0 or bounds[1] < bounds[0]:
            raise ValueError("each window must be a positive, ordered pair of distances")
        checked_windows.append([float(bounds[0]), float(bounds[1])])
    if not checked_windows:
        raise ValueError("provide at least one fit window")

    points, by_branch = [], defaultdict(list)
    for index, (position, value, label) in enumerate(zip(r, c, labels)):
        reason = ("exact_zero_declared" if exact_zero else "zero_unresolved") if value == 0 else None
        if reason is None and abs(value) <= floor:
            reason = "below_absolute_floor"
        if reason is None and errors is not None and abs(value) / min_signal_to_error < errors[index]:
            reason = "below_signal_to_error"
        points.append({"index": index, "distance": float(position), "correlation": float(value),
                       "branch": label, "resolved": reason is None, "reason": reason,
                       "absolute_error": None if errors is None else float(errors[index])})
        by_branch[label].append(index)

    segments = []
    for label, indices in by_branch.items():
        current = []
        def finish():
            if current:
                segments.append({"id": len(segments), "branch": label,
                                 "sign": 1 if c[current[0]] > 0 else -1,
                                 "indices": current.copy(), "distances": [float(r[i]) for i in current]})
                current.clear()
        for index in indices:
            if not points[index]["resolved"]:
                finish()
                continue
            if current and np.sign(c[current[-1]]) != np.sign(c[index]):
                finish()
            current.append(index)
        finish()

    local, fitted_windows = [], []
    for segment in segments:
        indices = segment["indices"]
        for left, right in zip(indices, indices[1:]):
            delta_r = float(r[right] - r[left])
            delta_log_r = float(np.log1p(delta_r / r[left]))
            delta_log_c = float(np.log(abs(c[right])) - np.log(abs(c[left])))
            p_eff = -delta_log_c / delta_log_r
            inverse_xi = -delta_log_c / delta_r
            if not np.isfinite(p_eff) or not np.isfinite(inverse_xi):
                # The supplied points remain present, but no unrepresentable
                # local coefficient is emitted into the JSON diagnostics.
                continue
            xi_eff = float(1 / inverse_xi) if inverse_xi > 0 else None
            if xi_eff is not None and not np.isfinite(xi_eff):
                xi_eff = None
            local.append({"segment": segment["id"], "branch": segment["branch"],
                          "r_left": float(r[left]), "r_right": float(r[right]),
                          "r_geometric": float(np.exp((np.log(r[left]) + np.log(r[right])) / 2)),
                          "r_logarithmic": delta_r / delta_log_r,
                          "p_eff": p_eff, "inverse_xi_eff": inverse_xi, "xi_eff": xi_eff})
        for window in checked_windows:
            used = [i for i in indices if window[0] <= r[i] <= window[1]]
            if not used:
                continue
            models = {model: _fit_model(r[used], np.log(np.abs(c[used])), model)
                      for model in ["power_law", "exponential", "power_exponential"]}
            fitted_windows.append({"segment": segment["id"], "branch": segment["branch"],
                                   "requested_window": window.copy(),
                                   "used_distances": [float(r[i]) for i in used],
                                   "npoints": len(used), "models": models})
    have_fit = any(model["status"] == "ok" for window in fitted_windows for model in window["models"].values())
    return {"status": "descriptive_only" if have_fit else "insufficient_resolved_data",
            "settings": {"floor": floor, "min_signal_to_error": min_signal_to_error,
                         "exact_zero_declared": bool(exact_zero), "requested_windows": checked_windows,
                         "fit_weighting": "unweighted log absolute correlation", "condition_limit": 1e6},
            "points": points, "segments": segments, "local": local,
            "windows": fitted_windows, "caveats": _CAVEATS.copy()}
