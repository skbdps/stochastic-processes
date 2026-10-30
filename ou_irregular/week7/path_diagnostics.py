"""Exploratory, immutable diagnostics of a fixed selection of headline paths.

The original Week 3/4 estimators, source addresses, covariance gates and archived
records are composed unchanged.  A diagnostic never supplies a replacement fit
or interval.  In particular an invalid point's finite boundary candidate is
preserved in the original record, but is not treated as an interior estimate.

The independent profile and exact full NLL use ``longdouble`` arithmetic.  This
is the platform's available precision, not arbitrary precision.  Full Hessian
diagnostics lift the *registered float64 normalized sample*, whereas the profile
normalizes the original observations in longdouble using the recorded observed
scales.  Those representations and their differences are retained explicitly.
"""
from collections.abc import Mapping
from dataclasses import dataclass
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from ..experiments.design import ESTIMATORS
from ..experiments.execution import ENGINE_FILES, registered_wald, source_records
from ..experiments.provenance import (ROOT, canonical, capture, digest_bytes,
                                      file_digest, utc_now, verify)
from ..experiments.storage import CellStore, atomic_bytes, atomic_json, run_lock
from ..experiments.workflow import load_saved_design
from ..inference.adapters import OuLikelihood
from ..inference.coordinates import ObservationCoordinates
from ..inference.hessian import FiniteDifferenceHessian
from ..week4.simulation import OuSampleSource, array_digest
from .selection import build_selection


FILES = (*ENGINE_FILES, "ou_irregular/experiments/audit.py",
         "ou_irregular/experiments/workflow.py", "ou_irregular/experiments/analysis.py",
         "ou_irregular/week5/config.py", "ou_irregular/week6/config.py",
         "configs/week5_e1.yaml", "configs/week6_e2.yaml",
         "ou_irregular/week7/__init__.py", "ou_irregular/week7/selection.py",
         "ou_irregular/week7/path_diagnostics.py", "docs/week7_exploratory_protocol.md",
         "docs/week7_mathematics.md")
SETTINGS = {
    "replay_rtol": 2e-10, "replay_atol": 2e-12,
    "relative_steps": [.001, .0005, .00025, .000125],
    "log_theta_lower_cap": -40., "log_theta_upper_cap": 40.,
    "grid_max_spacing": .05,
    "bounded_options": {"xatol": 1e-12, "maxiter": 500},
    "start_offsets": [-2., -1., 0., 1., 2.], "simplex_width": .1,
    "nelder_mead_options": {"xatol": 1e-10, "fatol": 1e-10, "maxiter": 2000,
                            "maxfev": 4000, "adaptive": False},
    "resolution_absolute": 1e-7, "resolution_relative": 1e-9,
    "replacement_fits_or_intervals": False,
}
_LD = np.longdouble
_LOG_2PI = np.log(_LD(2) * np.arccos(_LD(-1)))
_OPEN_STATES = {"running_incomplete", "incomplete"}
_CLOSED_STATES = {"completed_diagnostic_evidence", "diagnostic_review_required"}


def precision_record():
    """Record significand precision separately from storage width."""
    info = np.finfo(np.longdouble)
    return {"dtype": np.dtype(np.longdouble).name, "storage_bytes": info.dtype.itemsize,
            "fraction_bits": int(info.nmant), "significand_bits": int(info.nmant + 1),
            "decimal_precision": int(info.precision), "epsilon": str(info.eps),
            "extended_beyond_float64": bool(info.nmant > np.finfo(float).nmant),
            "arbitrary_precision": False}


def _safe(value):
    """JSON values with no NaN/Infinity; extended scalars are handled explicitly."""
    if isinstance(value, Mapping):
        return {str(k): _safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list, np.ndarray)):
        return [_safe(v) for v in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (float, np.floating)):
        converted = float(value)
        return converted if math.isfinite(converted) else None
    return value


def _extended(value):
    value = _LD(value)
    return str(value) if np.isfinite(value) else None


def _scalar_record(values):
    """Retain a plotting float and the unrounded longdouble decimal together."""
    result = {}
    for name, value in values.items():
        result[name] = _safe(value)
        result[name + "_extended"] = _extended(value)
    return result


def _data(values, gaps):
    y, d = np.asarray(values, dtype=np.longdouble), np.asarray(gaps, dtype=np.longdouble)
    if (y.ndim != 1 or d.ndim != 1 or len(y) != len(d) + 1 or len(d) < 2
            or not np.isfinite(y).all() or not np.isfinite(d).all() or np.any(d <= 0)):
        raise ValueError("finite observations and matching strictly positive gaps required")
    return y, d


def _components(log_theta, gaps):
    theta = np.exp(_LD(log_theta))
    if not np.isfinite(theta) or theta <= 0:
        raise ValueError("theta must be positive and representable")
    # expm1 prevents small-theta cancellation.  a is innovation variance per
    # sigma^2, and b is the conditional-mean coefficient of the level mu.
    b = -np.expm1(-theta * gaps)
    a = -np.expm1(-2 * theta * gaps) / (2 * theta)
    if not np.isfinite(a).all() or np.any(a <= 0):
        raise ValueError("innovation variances are not representable")
    return theta, b, a


def profile_at(log_theta, values, gaps):
    """Independent fixed-theta GLS profile, with variance divisor K=n-1.

    For z=dy+b*y_previous and A=diag(a), solve the one-dimensional normal
    equation mu=(b' A^-1 z)/(b' A^-1 b); then sigma^2=RSS/K.  The return values
    remain longdouble, so numerical comparisons need not round them to float64.
    """
    y, d = _data(values, gaps)
    theta, b, a = _components(log_theta, d)
    z = np.diff(y) + b * y[:-1]
    denominator = np.sum(b * b / a, dtype=np.longdouble)
    if not np.isfinite(denominator) or denominator <= 0:
        raise ValueError("fixed-theta level normal equation is degenerate")
    mu = np.sum(b * z / a, dtype=np.longdouble) / denominator
    residual = z - b * mu
    sigma2 = np.mean(residual * residual / a, dtype=np.longdouble)
    if not np.isfinite(sigma2) or sigma2 <= 0:
        raise ValueError("profile residual variance is nonpositive or nonfinite")
    nll = (len(d) * (_LOG_2PI + 1 + np.log(sigma2)) + np.sum(np.log(a))) / 2
    if not np.isfinite(nll) or not np.isfinite(mu):
        raise ValueError("profile parameters or objective are not finite")
    return {"log_theta": _LD(log_theta), "theta": theta, "mu": mu,
            "sigma2": sigma2, "nll": nll}


def profile_limits(values, gaps):
    """The profiled Brownian-with-free-drift and iid closure likelihoods.

    In the lower limit alpha=theta*mu can stay finite while mu diverges.  Setting
    mu=0 before taking the limit would incorrectly impose zero drift.  In the
    upper limit sigma^2/(2*theta), rather than sigma itself, remains finite.
    """
    y, d = _data(values, gaps)
    dy = np.diff(y)
    alpha = np.sum(dy) / np.sum(d)
    sigma2 = np.mean((dy - alpha * d) ** 2 / d)
    mu = np.mean(y[1:])
    variance = np.mean((y[1:] - mu) ** 2)
    if sigma2 <= 0 or variance <= 0 or not np.isfinite([sigma2, variance]).all():
        raise ValueError("a limiting residual variance is nonpositive or nonfinite")
    brownian_nll = (len(d) * (_LOG_2PI + 1 + np.log(sigma2)) + np.sum(np.log(d))) / 2
    iid_nll = len(d) * (_LOG_2PI + 1 + np.log(variance)) / 2
    return {"brownian": {"alpha": alpha, "sigma2": sigma2, "nll": brownian_nll},
            "iid": {"mu": mu, "variance": variance, "nll": iid_nll}}


def full_exact_nll(working_parameters, values, gaps):
    """Independent full three-parameter conditional NLL, without an X0 term."""
    y, d = _data(values, gaps)
    eta = np.asarray(working_parameters, dtype=np.longdouble)
    if eta.shape != (3,) or not np.isfinite(eta).all():
        raise ValueError("three finite working parameters required")
    _, b, a = _components(eta[0], d)
    log_variance = 2 * eta[2] + np.log(a)
    residual = np.diff(y) + b * (y[:-1] - eta[1])
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        nll = np.sum(_LOG_2PI + log_variance + residual * residual * np.exp(-log_variance)) / 2
    if not np.isfinite(nll):
        raise ValueError("full conditional NLL is nonfinite")
    return nll


def central_hessian(objective, point, steps, extended=False):
    """Full central stencil; registered float64 or independent longdouble.

    Only the float64 route delegates to the historical provider's exact stencil.
    In the extended route all points, function values, subtractions and quotients
    remain longdouble until a *completed* matrix is used for linear algebra.
    """
    dtype = np.longdouble if extended else float
    point, steps = np.asarray(point, dtype=dtype), np.asarray(steps, dtype=dtype)
    if (point.shape != (3,) or steps.shape != (3,) or not np.isfinite(point).all()
            or not np.isfinite(steps).all() or np.any(steps <= 0)):
        raise ValueError("finite three-dimensional point and positive steps required")
    center = dtype(objective(point))
    if not extended:
        return FiniteDifferenceHessian._central(objective, point, steps, center)
    hessian, gradient = np.empty((3, 3), dtype=dtype), np.empty(3, dtype=dtype)
    for i in range(3):
        ei = np.zeros(3, dtype=dtype)
        ei[i] = steps[i]
        plus, minus = dtype(objective(point + ei)), dtype(objective(point - ei))
        gradient[i] = (plus - minus) / (2 * steps[i])
        hessian[i, i] = ((plus - center) + (minus - center)) / steps[i] ** 2
        for j in range(i):
            ej = np.zeros(3, dtype=dtype)
            ej[j] = steps[j]
            pp = dtype(objective(point + ei + ej))
            pm = dtype(objective(point + ei - ej))
            mp = dtype(objective(point - ei + ej))
            mm = dtype(objective(point - ei - ej))
            hessian[i, j] = hessian[j, i] = ((pp - pm) - (mp - mm)) / (4 * steps[i] * steps[j])
    return hessian, gradient


def _matrix_diagnostics(hessian, gradient):
    hessian, gradient = np.asarray(hessian, dtype=float), np.asarray(gradient, dtype=float)
    result = {"positive_definite": False, "eigenvalues": None,
              "condition_number": None, "standardized_score": None,
              "theta_hessian_diagonal": None, "theta_schur_complement_at_stored_point": None,
              "schur_status": "unavailable",
              "schur_profile_interpretation": "only_at_nuisance_stationarity_on_same_representation"}
    if not np.isfinite(hessian).all():
        return dict(result, status="nonfinite_completed_derivatives")
    try:
        result["theta_hessian_diagonal"] = float(hessian[0, 0])
        try:
            schur = hessian[0, 0] - hessian[0, 1:] @ np.linalg.solve(hessian[1:, 1:], hessian[1:, 0])
            result["theta_schur_complement_at_stored_point"] = _safe(schur)
            result["schur_status"] = "evaluated_at_stored_point" if math.isfinite(schur) else "nonfinite"
        except np.linalg.LinAlgError:
            result["schur_status"] = "singular_nuisance_block"
        eigenvalues = np.linalg.eigvalsh(hessian)
        result["eigenvalues"] = eigenvalues.tolist()
        result["condition_number"] = _safe(np.linalg.cond(hessian))
        if eigenvalues[0] <= 0:
            return dict(result, status="nonpositive_hessian")
        np.linalg.cholesky(hessian)
        result["positive_definite"] = True
        if not np.isfinite(gradient).all():
            return dict(result, status="nonfinite_completed_derivatives")
        score2 = float(gradient @ np.linalg.solve(hessian, gradient))
        if score2 < 0 or not math.isfinite(score2):
            return dict(result, status="invalid_score_quadratic")
        result["standardized_score"] = math.sqrt(score2)
        return dict(result, status="evaluated")
    except np.linalg.LinAlgError as exc:
        return dict(result, status="linear_algebra_failure", error=str(exc))


def _adjacent(coarse, fine):
    result = {"coarse_relative_step": coarse["relative_step"],
              "fine_relative_step": fine["relative_step"],
              "relative_frobenius": None, "whitened_spectral": None}
    if coarse.get("hessian") is None or fine.get("hessian") is None:
        return dict(result, status="derivatives_unavailable")
    a, b = np.asarray(coarse["hessian"], float), np.asarray(fine["hessian"], float)
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        return dict(result, status="nonfinite_completed_derivatives")
    denominator = np.linalg.norm(b, "fro")
    if denominator > 0:
        result["relative_frobenius"] = _safe(np.linalg.norm(a - b, "fro") / denominator)
    if not fine["positive_definite"]:
        return dict(result, status="fine_hessian_not_positive_definite")
    try:
        lower = np.linalg.cholesky(b)
        left = np.linalg.solve(lower, a - b)
        whitened = np.linalg.solve(lower, left.T).T
        whitened = (whitened + whitened.T) / 2
        result["whitened_spectral"] = _safe(np.max(np.abs(np.linalg.eigvalsh(whitened))))
        return dict(result, status="evaluated")
    except np.linalg.LinAlgError as exc:
        return dict(result, status="linear_algebra_failure", error=str(exc))


def hessian_ladder(objective, point, extended=False):
    """All four fixed steps, including failures; no replacement covariance."""
    point = np.asarray(point, dtype=float)
    if point.shape != (3,) or not np.isfinite(point).all():
        raise ValueError("a finite original working point is required")
    entries = []
    for relative in SETTINGS["relative_steps"]:
        # The registered float64 step vector is lifted unchanged for LD.  This
        # separates precision of the stencil from changing its target spacing.
        steps = relative * np.maximum(1., np.abs(point))
        entry = {"relative_step": relative, "steps": steps.tolist(),
                 "hessian": None, "gradient": None, "positive_definite": False}
        try:
            hessian, gradient = central_hessian(objective, point, steps, extended=extended)
            if extended:
                entry.update(hessian_extended=[[_extended(v) for v in row] for row in hessian],
                             gradient_extended=[_extended(v) for v in gradient])
            completed_hessian = np.asarray(hessian, dtype=float)
            completed_gradient = np.asarray(gradient, dtype=float)
            entry.update(hessian=_safe(completed_hessian), gradient=_safe(completed_gradient),
                         **_matrix_diagnostics(completed_hessian, completed_gradient))
        except Exception as exc:
            entry.update(status="derivative_evaluation_failure",
                         error=f"{type(exc).__name__}: {exc}", eigenvalues=None,
                         condition_number=None, standardized_score=None,
                         theta_hessian_diagonal=None, theta_schur_complement_at_stored_point=None,
                         schur_status="unavailable",
                         schur_profile_interpretation="only_at_nuisance_stationarity_on_same_representation")
        entries.append(entry)
    return {"status": "evaluated", "precision": "longdouble" if extended else "registered_float64",
            "working_point": point.tolist(), "entries": entries,
            "adjacent": [_adjacent(a, b) for a, b in zip(entries[:-1], entries[1:])],
            "replacement_interval_constructed": False}


def _profile_bounds(gaps):
    d = np.asarray(gaps, dtype=np.longdouble)
    raw_lower = min(_LD(-30), np.log(_LD("1e-10")) - np.log(np.max(d)))
    raw_upper = max(_LD(12), np.log(_LD(50)) - np.log(np.min(d)))
    lower_ld, upper_ld = max(_LD(-40), raw_lower), min(_LD(40), raw_upper)
    lower, upper = float(lower_ld), float(upper_ld)
    # Round outwards only inside the caps.  The sampled interval must not miss
    # a specified tail merely because conversion to binary64 rounded inwards.
    if _LD(lower) > lower_ld and lower > -40:
        lower = float(np.nextafter(lower, -np.inf))
    if _LD(upper) < upper_ld and upper < 40:
        upper = float(np.nextafter(upper, np.inf))
    if lower > -40 and np.exp(_LD(lower)) * np.max(d) > _LD("1e-10"):
        lower = max(-40., float(np.nextafter(lower, -np.inf)))
    if upper < 40 and np.exp(_LD(upper)) * np.min(d) < 50:
        upper = min(40., float(np.nextafter(upper, np.inf)))
    return {"lower": lower, "upper": upper,
            "lower_cap_limited": bool(raw_lower < -40), "upper_cap_limited": bool(raw_upper > 40),
            "lower_tail_reached": bool(np.exp(_LD(lower)) * np.max(d) <= _LD("1e-10")),
            "upper_tail_reached": bool(np.exp(_LD(upper)) * np.min(d) >= 50),
            "uncapped_lower_extended": str(raw_lower), "uncapped_upper_extended": str(raw_upper)}


def _profile_search(values, gaps, eta):
    bounds = _profile_bounds(gaps)
    lower, upper = bounds["lower"], bounds["upper"]
    intervals = int(math.ceil((upper - lower) / SETTINGS["grid_max_spacing"]))
    while True:
        grid = np.linspace(lower, upper, intervals + 1)
        if np.max(np.diff(grid)) <= SETTINGS["grid_max_spacing"]:
            break
        intervals += 1

    def evaluate(log_theta):
        try:
            return profile_at(_LD(log_theta), values, gaps)
        except (ValueError, FloatingPointError, OverflowError):
            return None

    def objective(log_theta):
        result = evaluate(float(np.asarray(log_theta).reshape(-1)[0]))
        return float(result["nll"]) if result is not None else float("inf")

    grid_values = [evaluate(value) for value in grid]
    grid_nll = [value["nll"] if value is not None else _LD(np.inf) for value in grid_values]
    candidates, all_candidates, brackets, starts = [], [], [], []

    def candidate(kind, at, **extra):
        value = evaluate(at)
        record = {"kind": kind, "status": "finite" if value is not None else "nonfinite",
                  "log_theta": float(at), **extra}
        if value is not None:
            record.update(_scalar_record(value))
            candidates.append((value["nll"], record))
        all_candidates.append(record)
        return record

    candidate("stored_theta_profiled_nuisances", eta[0])
    candidate("grid_lower_endpoint", lower)
    candidate("grid_upper_endpoint", upper)
    finite_indices = [i for i, value in enumerate(grid_nll) if np.isfinite(value)]
    if not finite_indices:
        raise ValueError("the complete fixed profile grid is nonfinite")
    best_grid = min(finite_indices, key=lambda index: grid_nll[index])
    candidate("grid_minimum", grid[best_grid], grid_index=best_grid)
    for index in range(1, len(grid) - 1):
        # <= is intentional: every finite sampled plateau/tie is a bracket.
        if (np.isfinite(grid_nll[index-1:index+2]).all()
                and grid_nll[index] <= grid_nll[index-1]
                and grid_nll[index] <= grid_nll[index+1]):
            bracket = [float(grid[index-1]), float(grid[index+1])]
            search = {"grid_index": index, "bounds": bracket, "method": "bounded",
                      "options": dict(SETTINGS["bounded_options"])}
            try:
                result = minimize_scalar(objective, bounds=bracket, method="bounded",
                                         options=dict(SETTINGS["bounded_options"]))
                record = candidate("polished_grid_bracket", result.x, grid_index=index)
                search.update(success=bool(result.success), message=str(result.message),
                              nfev=int(result.nfev), nit=int(result.nit), candidate=record)
            except Exception as exc:
                # A solver exception does not delete the fixed attempted bracket
                # or prevent the remaining starts.  The retained grid point is
                # explicitly not claimed to be a returned optimizer solution.
                search.update(success=False, message=f"{type(exc).__name__}: {exc}",
                              nfev=None, nit=None, exception=True,
                              candidate=candidate("failed_bracket_grid_center", grid[index], grid_index=index))
            brackets.append(search)
    for offset in SETTINGS["start_offsets"]:
        x0 = float(np.clip(float(eta[0]) + offset, lower, upper))
        second = x0 + SETTINGS["simplex_width"]
        if second > upper:
            second = x0 - SETTINGS["simplex_width"]
        simplex = np.array([[x0], [second]], dtype=float)
        options = {**SETTINGS["nelder_mead_options"], "initial_simplex": simplex}
        search = {"offset": offset, "x0": [x0], "initial_simplex": simplex.tolist(),
                  "bounds": [[lower, upper]], "method": "Nelder-Mead",
                  "options": dict(SETTINGS["nelder_mead_options"])}
        try:
            result = minimize(objective, x0=np.array([x0]), method="Nelder-Mead",
                              bounds=[(lower, upper)], options=options)
            record = candidate("nelder_mead", result.x[0], offset=offset)
            search.update(success=bool(result.success), status=int(result.status),
                          message=str(result.message), nfev=int(result.nfev), nit=int(result.nit), candidate=record)
        except Exception as exc:
            search.update(success=False, status=None, message=f"{type(exc).__name__}: {exc}",
                          nfev=None, nit=None, exception=True,
                          candidate=candidate("failed_nelder_mead_initial_point", x0, offset=offset))
        starts.append(search)
    best_nll, best = min(candidates, key=lambda item: item[0])
    stored = profile_at(eta[0], values, gaps)
    full_nll = full_exact_nll(eta, values, gaps)
    limits = profile_limits(values, gaps)
    resolution = _LD("1e-7") + _LD("1e-9") * abs(best_nll)
    differences = {"stored_profile_minus_best": stored["nll"] - best_nll,
                   "stored_full_minus_best": full_nll - best_nll,
                   "stored_full_minus_stored_profile": full_nll - stored["nll"],
                   "brownian_limit_minus_best": limits["brownian"]["nll"] - best_nll,
                   "iid_limit_minus_best": limits["iid"]["nll"] - best_nll}
    tail_unresolved = (not bounds["lower_tail_reached"] or not bounds["upper_tail_reached"])
    limiting_contender = any(value["nll"] <= best_nll + resolution for value in limits.values())
    if tail_unresolved:
        disposition = "cap_limited_tail_unresolved"
    elif limiting_contender:
        disposition = "limiting_likelihood_competitive_or_numerically_unresolved"
    elif differences["stored_full_minus_best"] > resolution:
        disposition = "finite_search_improves_stored_full_nll_above_diagnostic_resolution"
    else:
        disposition = "stored_full_fit_matches_finite_search_within_diagnostic_resolution"
    return {"status": "evaluated", "bounds": bounds, "grid_log_theta": grid.tolist(),
            "grid_nll": _safe(grid_nll), "grid_nll_extended": [_extended(v) for v in grid_nll],
            "maximum_grid_spacing": float(np.max(np.diff(grid))),
            "polished_brackets": brackets, "local_searches": starts,
            "candidates": all_candidates, "best_finite": best,
            "stored_profile": _scalar_record(stored),
            "stored_full": _scalar_record({"nll": full_nll}),
            "limits": {name: _scalar_record(value) for name, value in limits.items()},
            "differences": _scalar_record(differences), "resolution": _safe(resolution),
            "resolution_extended": str(resolution), "disposition": disposition,
            "global_optimum_proved": False, "verified_software_defect": False,
            "replacement_fit_or_interval": False}


def _ar1_diagnostic(values, gaps):
    y, d = _data(values, gaps)
    x, z = y[:-1], y[1:]
    xc, zc = x - np.mean(x), z - np.mean(z)
    denominator = np.sum(xc * xc)
    # The float64 observed near-regular criterion mirrors the historical solver
    # but is reported as a numerical criterion, not an exact mathematical proof
    # for unequal gaps.  The proof itself requires exactly constant gaps.
    d64 = np.asarray(gaps, float)
    deviations = np.max(np.abs(d64 / float(np.mean(d64.astype(np.longdouble))) - 1))
    tolerance = 128 * np.finfo(float).eps * len(d64)
    exact_regular = bool(np.all(d == d[0]))
    result = {"exactly_equal_gaps": exact_regular,
              "historical_near_regular_criterion": bool(deviations <= tolerance),
              "maximum_relative_gap_deviation": float(deviations),
              "historical_near_regular_tolerance": float(tolerance),
              "proof_scope": "PFML at realized mean gap; exact only for equal gaps",
              "exact_boundary_proof_applicable": exact_regular}
    if denominator <= 0:
        return dict(result, status="degenerate_predictor", phi=None, intercept=None)
    phi = np.sum(xc * zc) / denominator
    intercept = np.mean(z) - phi * np.mean(x)
    variance = np.mean((z - intercept - phi * x) ** 2)
    if variance <= 0:
        conclusion = "zero_residual_variance_degeneracy"
    elif phi <= 0:
        conclusion = "PFML_supremum_at_phi_zero_no_finite_interior_estimate"
    elif phi >= 1:
        conclusion = "PFML_supremum_at_phi_one_lower_theta_boundary"
    else:
        conclusion = "PFML_interior_positive_correlation"
    return dict(result, status="evaluated", **_scalar_record({"phi": phi, "intercept": intercept,
                                                             "innovation_variance": variance}),
                conclusion=conclusion)


def analyze_sample(sample, stored_fits):
    """Analyze one original source after its complete three-fit replay passed."""
    rows = (dict(stored_fits) if isinstance(stored_fits, Mapping)
            else {row["estimator"]: row for row in stored_fits})
    if set(rows) != set(ESTIMATORS) or (not isinstance(stored_fits, Mapping)
                                     and len(stored_fits) != len(ESTIMATORS)):
        raise ValueError("exactly one original record per registered estimator required")
    coordinates = ObservationCoordinates.from_sample(sample)
    normalized = coordinates.normalized_sample
    registered_y = np.asarray(normalized.values, dtype=np.longdouble)
    registered_d = np.asarray(np.diff(normalized.times), dtype=np.longdouble)
    profile_y = ((np.asarray(sample.values, np.longdouble) - _LD(coordinates.offset))
                 / _LD(coordinates.state_scale))
    profile_d = np.diff(np.asarray(sample.times, np.longdouble)) / _LD(coordinates.time_scale)
    representations = {"registered_ladder": "lift_exact_float64_normalized_values_and_diff_times",
                       "profile": "normalize_original_float64_values_and_times_in_longdouble",
                       **_scalar_record({"maximum_absolute_value_difference": np.max(np.abs(profile_y - registered_y)),
                                         "maximum_absolute_gap_difference": np.max(np.abs(profile_d - registered_d))})}
    fits = {}
    for method in ESTIMATORS:
        row = rows[method]
        record = {"original_point_valid": row["point_valid"],
                  "original_point_status": row["point_status"],
                  "original_interval_valid": row["interval_valid"],
                  "original_interval_status": row["interval_status"]}
        if row["point_valid"] is not True:
            fits[method] = dict(record, status="skipped_original_point_invalid",
                                reason="retained candidate is not an interior fitted point")
            continue
        eta = coordinates.to_working(np.array([row[p] for p in ("theta", "mu", "sigma")]))
        likelihood = OuLikelihood(method)
        record.update(status="evaluated", working_point=eta.tolist(),
                      float64_ladder=hessian_ladder(lambda p: likelihood.nll(p, normalized), eta))
        if method == "exact":
            record["independent_longdouble_ladder"] = hessian_ladder(
                lambda p: full_exact_nll(p, registered_y, registered_d), eta, extended=True)
            representations.update(_scalar_record({
                "independent_full_nll_registered_representation": full_exact_nll(eta, registered_y, registered_d),
                "independent_full_nll_profile_representation": full_exact_nll(eta, profile_y, profile_d),
                "registered_float64_full_nll": _LD(likelihood.nll(eta, normalized)),
                "stored_natural_nll_minus_normalization_constant": _LD(row["nll"]) - len(profile_d) * np.log(_LD(coordinates.state_scale)),
            }))
        fits[method] = record
    exact = rows["exact"]
    if exact["point_valid"] is True:
        eta = np.asarray(fits["exact"]["working_point"])
        profile = _profile_search(profile_y, profile_d, eta)
    else:
        profile = {"status": "skipped_original_point_invalid",
                   "reason": "profile search is restricted to original valid exact points"}
    return {"status": "evaluated", "precision": precision_record(),
            "coordinates": {"state_offset": coordinates.offset,
                            "state_scale": coordinates.state_scale, "time_scale": coordinates.time_scale,
                            "source": "observed_data_only"},
            "representation_differences": representations,
            "ar1": _ar1_diagnostic(profile_y, profile_d), "fits": fits,
            "exact_profile": profile, "headline_results_replaced": False}


@dataclass(frozen=True)
class _DiagnosticDesign:
    serialized: str

    def to_dict(self):
        return json.loads(self.serialized)

    @property
    def digest(self):
        return digest_bytes(self.serialized.encode())


def _read_selection(selection):
    """Preserve literal selection bytes in addition to their semantic identity."""
    if isinstance(selection, Mapping):
        selected = dict(selection)
        return selected, (canonical(selected) + "\n").encode(), None
    path = Path(selection).resolve()
    payload = path.read_bytes()
    selected = json.loads(payload)
    if not isinstance(selected, dict):
        raise ValueError("selection must be a JSON object")
    try:
        relative = str(path.relative_to(ROOT.resolve()))
    except ValueError:
        relative = None
    return selected, payload, relative


def _diagnostic_design(selection, selection_bytes):
    return _DiagnosticDesign(canonical({
        "name": "week7_selected_path_diagnostics", "scope": "exploratory",
        "selection": selection, "selection_bytes_sha256": digest_bytes(selection_bytes),
        "settings": SETTINGS, "precision": precision_record(),
        "source_policy": "unchanged_original_headline_address_and_three_fits"}))


def _compare_replay(expected, actual, trail="record"):
    """Every field and missing value, with the fixed tolerance only for floats.

    Nested serialized diagnostics are parsed and checked recursively.  Matching
    just the displayed estimates or the set of successful records is inadequate:
    numerical candidates, rejected intervals and their reasons are evidence too.
    """
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(expected) != set(actual):
            raise ValueError(f"replay object keys changed at {trail}")
        for key in expected:
            if key.endswith("_json") and isinstance(expected[key], str):
                if not isinstance(actual[key], str):
                    raise ValueError(f"replay serialized diagnostic type changed at {trail}.{key}")
                try:
                    left, right = json.loads(expected[key]), json.loads(actual[key])
                except (ValueError, TypeError) as exc:
                    raise ValueError(f"invalid serialized diagnostic at {trail}.{key}") from exc
                _compare_replay(left, right, trail + "." + key)
            else:
                _compare_replay(expected[key], actual[key], trail + "." + key)
        return
    if isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise ValueError(f"replay sequence changed at {trail}")
        for i, (left, right) in enumerate(zip(expected, actual)):
            _compare_replay(left, right, f"{trail}[{i}]")
        return
    if expected is None or isinstance(expected, (bool, str, int)):
        if type(actual) is not type(expected) or actual != expected:
            raise ValueError(f"replay discrete field changed at {trail}")
        return
    if isinstance(expected, float):
        if (isinstance(actual, bool) or not isinstance(actual, (int, float))
                or not math.isfinite(expected) or not math.isfinite(actual)
                or not np.isclose(actual, expected, rtol=SETTINGS["replay_rtol"],
                                  atol=SETTINGS["replay_atol"])):
            raise ValueError(f"replay numerical field changed at {trail}")
        return
    raise ValueError(f"unsupported replay value at {trail}")


def _source_key(item):
    if (item.get("experiment") not in {"E1", "E2"}
            or not isinstance(item.get("cell_id"), str) or len(item["cell_id"]) != 64
            or any(c not in "0123456789abcdef" for c in item["cell_id"])
            or type(item.get("rep")) is not int or item["rep"] < 0
            or not isinstance(item.get("reasons"), list)
            or not item["reasons"] or not all(isinstance(v, str) for v in item["reasons"])
            or item["reasons"] != sorted(set(item["reasons"]))):
        raise ValueError("invalid selected source identity or reasons")
    return item["experiment"], item["cell_id"], item["rep"]


def _filename(item):
    experiment, cid, rep = _source_key(item)
    return f"{experiment}-{cid}-rep-{rep:06d}.json"


def _selected_originals(selection, e1, e2):
    """Load complete original three-method records from validated partitions."""
    wanted = {_source_key(item) for item in selection["sources"]}
    if len(wanted) != len(selection["sources"]):
        raise ValueError("duplicate selected source")
    originals = {key: {} for key in wanted}
    for label, directory in (("E1", e1), ("E2", e2)):
        design = load_saved_design(directory)
        metadata = json.loads((Path(directory) / "run_metadata.json").read_text())
        cells = {cell.identity: cell for cell in design.cells}
        wanted_cells = {cid for experiment, cid, _ in wanted if experiment == label}
        if not wanted_cells <= set(cells):
            raise ValueError("selected cell is absent from its original design")
        for cid in sorted(wanted_cells):
            store = CellStore(Path(directory) / "cells" / cid, cells[cid], design,
                              metadata["provenance"]["identity_sha256"], resume=True)
            if store.next_replication != design.replications or store.orphan_stop is not None:
                raise ValueError("selected input cell is incomplete")
            for row in store.records():
                key = (label, cid, row["rep"])
                if key not in wanted:
                    continue
                method = row["estimator"]
                if method not in ESTIMATORS or method in originals[key]:
                    raise ValueError("duplicate or unexpected original selected fit")
                originals[key][method] = row
    if any(set(rows) != set(ESTIMATORS) for rows in originals.values()):
        raise ValueError("the selected originals lack the full three-method Cartesian set")
    return originals


def _check_output_layout(output, selection, *, require_complete=False):
    output = Path(output)
    if not output.exists():
        return
    if not output.is_dir() or output.is_symlink():
        raise ValueError("output must be an ordinary directory")
    allowed = {".execution.lock", "run_metadata.json", "diagnostic_design.json",
               "selection.json", "sources", "audit.json"}
    if {p.name for p in output.iterdir()} - allowed:
        raise ValueError("unknown artifact in selected diagnostic output")
    for path in output.iterdir():
        if path.is_symlink() or (path.name != "sources" and not path.is_file()):
            raise ValueError("unexpected output object or symbolic link")
    sources = output / "sources"
    expected = {_filename(item) for item in selection["sources"]}
    if len(expected) != len(selection["sources"]):
        raise ValueError("duplicate selected source filename")
    actual = set()
    if sources.exists():
        if not sources.is_dir():
            raise ValueError("sources must be a directory")
        for path in sources.iterdir():
            if path.is_symlink() or not path.is_file() or path.name not in expected:
                raise ValueError("unknown or nonregular selected-source artifact")
            actual.add(path.name)
    if require_complete and actual != expected:
        raise ValueError("selected-source artifact set is incomplete")


def _check_frozen(provenance, design, fresh, require_frozen):
    verify(provenance, expected_design=design)
    if provenance["source_files"] != fresh["source_files"]:
        raise ValueError("selected diagnostics dependency manifest membership or hashes changed")
    if require_frozen:
        if provenance.get("uncommitted_dependency_files") != []:
            raise ValueError("selected diagnostics were not frozen before execution")
        commit = provenance.get("code_commit", "")
        if not isinstance(commit, str) or len(commit) != 40 or any(c not in "0123456789abcdef" for c in commit):
            raise ValueError("invalid recorded freeze commit")
        for path, expected in provenance["source_files"].items():
            recorded = subprocess.run(["git", "show", f"{commit}:{path}"], cwd=ROOT,
                                      capture_output=True, check=False)
            if recorded.returncode or digest_bytes(recorded.stdout) != expected:
                raise ValueError(f"recorded freeze commit does not contain dependency: {path}")


def _payload_identity(payload, item, originals, identity):
    if (type(payload.get("schema_version")) is not int or payload["schema_version"] != 1
            or payload.get("run_identity") != identity or payload.get("headline_results_replaced") is not False
            or canonical(payload.get("original_records")) != canonical(originals)):
        raise ValueError("selected-source payload identity or retained originals changed")
    _compare_replay(item, payload.get("selection"), "selection")
    if payload.get("status") not in {"evaluated", "original_replay_failed", "diagnostic_failure"}:
        raise ValueError("unknown selected-source diagnostic status")
    replay = payload.get("original_replay", {})
    if payload["status"] == "original_replay_failed":
        if replay.get("status") != "failed" or not isinstance(replay.get("error"), str) or not replay["error"]:
            raise ValueError("replay-failed payload has inconsistent replay state")
        if "diagnostic" in payload:
            raise ValueError("a failed original replay must not supply diagnostic analysis")
    else:
        if replay.get("status") != "passed":
            raise ValueError("diagnostics must follow a complete successful original replay")
        _compare_replay(originals, replay.get("records"), "retained_full_fit_replay")
        if replay.get("rtol") != SETTINGS["replay_rtol"] or replay.get("atol") != SETTINGS["replay_atol"]:
            raise ValueError("original replay tolerance changed")
        diagnostics, hashes = replay.get("source_diagnostics"), replay.get("array_hashes")
        if not isinstance(diagnostics, dict) or not isinstance(hashes, dict) or set(hashes) != {"times_sha256", "observations_sha256"}:
            raise ValueError("complete source replay diagnostics/hashes are required")
        if not set(hashes) <= set(diagnostics):
            raise ValueError("source replay diagnostics omit observation hashes")
        for method in ESTIMATORS:
            if not set(diagnostics) <= set(originals[method]):
                raise ValueError("unexpected source diagnostic field")
            _compare_replay({key: originals[method][key] for key in diagnostics}, diagnostics, "source_replay")
            _compare_replay({key: originals[method][key] for key in hashes}, hashes, "source_hashes")
    if payload["status"] == "evaluated":
        _validate_diagnostic(payload.get("diagnostic"), originals)
    if payload["status"] == "diagnostic_failure":
        diagnostic = payload.get("diagnostic", {})
        if (diagnostic.get("status") != "failed" or not isinstance(diagnostic.get("error"), str)
                or not diagnostic["error"]):
            raise ValueError("diagnostic-failure payload must retain its error")


def _finite_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _numeric_shape(value, shape):
    if not isinstance(value, list) or len(value) != shape[0]:
        return False
    if len(shape) == 1:
        return all(v is None or _finite_number(v) for v in value)
    return all(_numeric_shape(v, shape[1:]) for v in value)


def _require_scalars(record, names):
    if not isinstance(record, dict):
        raise ValueError("complete numeric and extended scalar records are required")
    for name in names:
        if not _finite_number(record.get(name)) or not isinstance(record.get(name + "_extended"), str):
            raise ValueError("numeric scalar or extended precision record is missing: " + name)
        try:
            value = _LD(record[name + "_extended"])
        except (ValueError, TypeError) as exc:
            raise ValueError("invalid extended numeric record: " + name) from exc
        if not np.isfinite(value) or float(value) != record[name]:
            raise ValueError("numeric/extended scalar representation mismatch: " + name)


def _validate_ladder(ladder, point, extended):
    if (not isinstance(ladder, dict) or ladder.get("status") != "evaluated"
            or ladder.get("precision") != ("longdouble" if extended else "registered_float64")
            or ladder.get("working_point") != point
            or ladder.get("replacement_interval_constructed") is not False
            or not isinstance(ladder.get("entries"), list) or len(ladder["entries"]) != 4
            or not isinstance(ladder.get("adjacent"), list) or len(ladder["adjacent"]) != 3):
        raise ValueError("complete four-step Hessian ladder is required")
    for relative, entry in zip(SETTINGS["relative_steps"], ladder["entries"]):
        expected_steps = (relative * np.maximum(1., np.abs(point))).tolist()
        required = {"relative_step", "steps", "status", "hessian", "gradient", "positive_definite",
                    "eigenvalues", "condition_number", "standardized_score", "theta_hessian_diagonal",
                    "theta_schur_complement_at_stored_point", "schur_status", "schur_profile_interpretation"}
        if (not isinstance(entry, dict) or not required <= set(entry) or entry["relative_step"] != relative
                or entry["steps"] != expected_steps or type(entry["positive_definite"]) is not bool):
            raise ValueError("Hessian step record is incomplete or changed")
        if entry["status"] == "derivative_evaluation_failure":
            if (entry["hessian"] is not None or entry["gradient"] is not None
                    or not isinstance(entry.get("error"), str) or not entry["error"]):
                raise ValueError("failed Hessian step must retain its evaluation error")
        else:
            if (entry["status"] not in {"evaluated", "nonpositive_hessian", "invalid_score_quadratic",
                                        "linear_algebra_failure", "nonfinite_completed_derivatives"}
                    or not _numeric_shape(entry["hessian"], (3, 3))
                    or not _numeric_shape(entry["gradient"], (3,))):
                raise ValueError("completed full derivative arrays are missing")
            if extended:
                hs, gs = entry.get("hessian_extended"), entry.get("gradient_extended")
                if (not isinstance(hs, list) or len(hs) != 3 or any(not isinstance(row, list) or len(row) != 3 for row in hs)
                        or not isinstance(gs, list) or len(gs) != 3
                        or any(v is not None and not isinstance(v, str) for row in hs for v in row)
                        or any(v is not None and not isinstance(v, str) for v in gs)):
                    raise ValueError("extended derivative decimal records are missing")
        if entry["eigenvalues"] is not None and not _numeric_shape(entry["eigenvalues"], (3,)):
            raise ValueError("Hessian eigenvalue record has wrong shape")
        if not entry["positive_definite"] and entry["standardized_score"] is not None:
            raise ValueError("nonpositive Hessian must not supply an inverse-based score")
        if entry["status"] == "evaluated" and (not entry["positive_definite"] or not _finite_number(entry["standardized_score"])):
            raise ValueError("evaluated positive Hessian requires its standardized score")
    for i, comparison in enumerate(ladder["adjacent"]):
        required = {"coarse_relative_step", "fine_relative_step", "relative_frobenius", "whitened_spectral", "status"}
        if (not isinstance(comparison, dict) or not required <= set(comparison)
                or comparison["coarse_relative_step"] != SETTINGS["relative_steps"][i]
                or comparison["fine_relative_step"] != SETTINGS["relative_steps"][i + 1]
                or comparison["status"] not in {"evaluated", "derivatives_unavailable", "nonfinite_completed_derivatives",
                                                "fine_hessian_not_positive_definite", "linear_algebra_failure"}):
            raise ValueError("adjacent Hessian comparisons are incomplete or changed")
        if not ladder["entries"][i + 1]["positive_definite"] and comparison["whitened_spectral"] is not None:
            raise ValueError("nonpositive fine Hessian cannot supply a whitened comparison")


def _validate_profile(profile, point):
    required = {"status", "bounds", "grid_log_theta", "grid_nll", "grid_nll_extended", "maximum_grid_spacing",
                "polished_brackets", "local_searches", "candidates", "best_finite", "stored_profile", "stored_full",
                "limits", "differences", "resolution", "resolution_extended", "disposition", "global_optimum_proved",
                "verified_software_defect", "replacement_fit_or_interval"}
    if (not isinstance(profile, dict) or not required <= set(profile) or profile["status"] != "evaluated"
            or profile["global_optimum_proved"] is not False or profile["verified_software_defect"] is not False
            or profile["replacement_fit_or_interval"] is not False):
        raise ValueError("complete nonreplacement exact profile evidence is required")
    grid, nll, strings = profile["grid_log_theta"], profile["grid_nll"], profile["grid_nll_extended"]
    bounds = profile["bounds"]
    if (not isinstance(grid, list) or len(grid) < 3 or not all(_finite_number(v) for v in grid)
            or not isinstance(nll, list) or len(nll) != len(grid) or not isinstance(strings, list) or len(strings) != len(grid)
            or any(v is not None and not _finite_number(v) for v in nll)
            or any(v is not None and not isinstance(v, str) for v in strings)
            or not isinstance(bounds, dict) or grid[0] != bounds.get("lower") or grid[-1] != bounds.get("upper")
            or min(np.diff(grid)) <= 0 or max(np.diff(grid)) > SETTINGS["grid_max_spacing"]):
        raise ValueError("complete fixed profile grid is required")
    for key in ("lower_cap_limited", "upper_cap_limited", "lower_tail_reached", "upper_tail_reached"):
        if type(bounds.get(key)) is not bool:
            raise ValueError("profile tail/cap disposition is missing")
    searches = profile["local_searches"]
    if not isinstance(searches, list) or len(searches) != 5:
        raise ValueError("all five genuine local-search starts must be retained")
    for offset, search in zip(SETTINGS["start_offsets"], searches):
        x0 = float(np.clip(point[0] + offset, bounds["lower"], bounds["upper"]))
        second = x0 + SETTINGS["simplex_width"]
        if second > bounds["upper"]:
            second = x0 - SETTINGS["simplex_width"]
        if (search.get("offset") != offset or search.get("x0") != [x0]
                or search.get("initial_simplex") != [[x0], [second]]
                or search.get("bounds") != [[bounds["lower"], bounds["upper"]]]
                or search.get("method") != "Nelder-Mead" or search.get("options") != SETTINGS["nelder_mead_options"]
                or type(search.get("success")) is not bool or not isinstance(search.get("candidate"), dict)):
            raise ValueError("local-search settings/start/outcome changed")
    # Count brackets from saved extended objectives, so float64 conversion cannot
    # erase distinctions or invent plateaus in this structural readback check.
    extended = [_LD(v) if v is not None else _LD(np.inf) for v in strings]
    indices = [i for i in range(1, len(grid)-1) if np.isfinite(extended[i-1:i+2]).all()
               and extended[i] <= extended[i-1] and extended[i] <= extended[i+1]]
    if not isinstance(profile["polished_brackets"], list) or [row.get("grid_index") for row in profile["polished_brackets"]] != indices:
        raise ValueError("every finite sampled local-minimum bracket must be retained")
    for row in profile["polished_brackets"]:
        i = row["grid_index"]
        if (row.get("bounds") != [grid[i-1], grid[i+1]] or row.get("method") != "bounded"
                or row.get("options") != SETTINGS["bounded_options"] or type(row.get("success")) is not bool
                or not isinstance(row.get("candidate"), dict)):
            raise ValueError("profile bracket settings/outcome changed")
    if not isinstance(profile["limits"], dict) or set(profile["limits"]) != {"brownian", "iid"}:
        raise ValueError("both profile closure limits are required")
    _require_scalars(profile["limits"]["brownian"], ("alpha", "sigma2", "nll"))
    _require_scalars(profile["limits"]["iid"], ("mu", "variance", "nll"))
    _require_scalars(profile["stored_profile"], ("log_theta", "theta", "mu", "sigma2", "nll"))
    _require_scalars(profile["stored_full"], ("nll",))
    _require_scalars(profile["differences"], ("stored_profile_minus_best", "stored_full_minus_best",
                     "stored_full_minus_stored_profile", "brownian_limit_minus_best", "iid_limit_minus_best"))
    if (not isinstance(profile["candidates"], list) or not profile["candidates"]
            or not isinstance(profile["best_finite"], dict)
            or not _finite_number(profile["best_finite"].get("nll"))):
        raise ValueError("profile finite candidate evidence is missing")


def _validate_ar1(ar1):
    required = {"status", "exactly_equal_gaps", "historical_near_regular_criterion",
                "maximum_relative_gap_deviation", "historical_near_regular_tolerance",
                "proof_scope", "exact_boundary_proof_applicable", "phi", "intercept"}
    if not isinstance(ar1, dict) or not required <= set(ar1):
        raise ValueError("complete AR1 evidence and proof scope are required")
    if (not all(type(ar1[key]) is bool for key in ("exactly_equal_gaps", "historical_near_regular_criterion",
                                                  "exact_boundary_proof_applicable"))
            or ar1["exact_boundary_proof_applicable"] is not ar1["exactly_equal_gaps"]
            or ar1["proof_scope"] != "PFML at realized mean gap; exact only for equal gaps"
            or not _finite_number(ar1["maximum_relative_gap_deviation"])
            or not _finite_number(ar1["historical_near_regular_tolerance"])):
        raise ValueError("AR1 regularity evidence or boundary proof scope changed")
    if ar1["status"] == "evaluated":
        _require_scalars(ar1, ("phi", "intercept", "innovation_variance"))
        if ar1.get("conclusion") not in {"zero_residual_variance_degeneracy",
                "PFML_supremum_at_phi_zero_no_finite_interior_estimate",
                "PFML_supremum_at_phi_one_lower_theta_boundary", "PFML_interior_positive_correlation"}:
            raise ValueError("AR1 boundary conclusion is missing")
    elif ar1["status"] != "degenerate_predictor" or ar1["phi"] is not None or ar1["intercept"] is not None:
        raise ValueError("AR1 degeneracy state is inconsistent")


def _validate_diagnostic(diagnostic, originals):
    required = {"status", "precision", "coordinates", "representation_differences", "ar1", "fits",
                "exact_profile", "headline_results_replaced"}
    if (not isinstance(diagnostic, dict) or not required <= set(diagnostic)
            or diagnostic["status"] != "evaluated" or diagnostic["headline_results_replaced"] is not False
            or diagnostic["precision"] != precision_record() or set(diagnostic.get("fits", {})) != set(ESTIMATORS)):
        raise ValueError("evaluated payload lacks the complete three-method diagnostics")
    coordinates = diagnostic["coordinates"]
    if (not isinstance(coordinates, dict) or coordinates.get("source") != "observed_data_only"
            or not all(_finite_number(coordinates.get(k)) for k in ("state_offset", "state_scale", "time_scale"))
            or coordinates["state_scale"] <= 0 or coordinates["time_scale"] <= 0):
        raise ValueError("observed normalization record is incomplete")
    representations = diagnostic["representation_differences"]
    representation_keys = {"registered_ladder", "profile", "maximum_absolute_value_difference",
                           "maximum_absolute_value_difference_extended", "maximum_absolute_gap_difference",
                           "maximum_absolute_gap_difference_extended"}
    if not isinstance(representations, dict) or not representation_keys <= set(representations):
        raise ValueError("the two observed data representations must be recorded")
    _validate_ar1(diagnostic["ar1"])
    for method in ESTIMATORS:
        row, fit = originals[method], diagnostic["fits"][method]
        for field in ("point_valid", "point_status", "interval_valid", "interval_status"):
            _compare_replay(row[field], fit.get("original_" + field), method + "." + field)
        if row["point_valid"] is not True:
            if fit.get("status") != "skipped_original_point_invalid" or "float64_ladder" in fit or "independent_longdouble_ladder" in fit:
                raise ValueError("an invalid original point must not be promoted to a derivative target")
            continue
        point = fit.get("working_point")
        if fit.get("status") != "evaluated" or not _numeric_shape(point, (3,)) or any(v is None for v in point):
            raise ValueError("valid original point lacks its working coordinates")
        _validate_ladder(fit.get("float64_ladder"), point, False)
        if method == "exact":
            _validate_ladder(fit.get("independent_longdouble_ladder"), point, True)
            _validate_profile(diagnostic["exact_profile"], point)
            for key in ("independent_full_nll_registered_representation", "independent_full_nll_profile_representation",
                        "registered_float64_full_nll", "stored_natural_nll_minus_normalization_constant"):
                if not _finite_number(representations.get(key)) or not isinstance(representations.get(key + "_extended"), str):
                    raise ValueError("full-NLL representation comparison is missing")
    if originals["exact"]["point_valid"] is not True:
        if diagnostic["exact_profile"].get("status") != "skipped_original_point_invalid":
            raise ValueError("invalid original exact point must not enter the profile search")


class ManifestDriftError(RuntimeError):
    """Unexpected root state is preserved, never overwritten by cached state."""


class _ManifestGuard:
    def __init__(self, path, expected):
        self.path, self.expected, self.drifted = Path(path), expected, False

    def _read(self):
        return self.path.read_bytes() if self.path.is_file() else None

    def check(self, stage):
        observed = self._read()
        if observed != self.expected:
            self._incident(stage, self.expected, observed)

    def _incident(self, stage, expected, observed):
        self.drifted = True
        incident = Path(tempfile.mkdtemp(prefix=self.path.parent.name + "-manifest-incident-",
                                         dir=self.path.parent.parent))
        if expected is not None:
            atomic_bytes(incident / "expected_run_metadata.json", expected, immutable=True)
        if observed is not None:
            atomic_bytes(incident / "observed_run_metadata.json", observed, immutable=True)
        atomic_json(incident / "incident.json", {
            "status": "unresolved_external_manifest_drift", "stage": stage, "recorded_utc": utc_now(),
            "manifest_path": str(self.path), "expected_exists": expected is not None,
            "observed_exists": observed is not None,
            "expected_sha256": digest_bytes(expected) if expected is not None else None,
            "observed_sha256": digest_bytes(observed) if observed is not None else None,
            "observed_manifest_overwritten_by_handler": False}, immutable=True)
        raise ManifestDriftError(f"manifest drift at {stage}; preserved evidence in {incident}")

    def write(self, metadata, stage):
        self.check(stage + ":before_write")
        payload = (canonical(metadata) + "\n").encode()
        atomic_bytes(self.path, payload)
        observed = self._read()
        if observed != payload:
            self._incident(stage + ":after_write", payload, observed)
        self.expected = payload


def _check_incident_block(output):
    """An unresolved external-state incident cannot become an automatic resume."""
    if not output.parent.exists():
        return
    prefix = output.name + "-manifest-incident-"
    for incident in output.parent.iterdir():
        if not incident.name.startswith(prefix):
            continue
        record_path = incident / "incident.json"
        if incident.is_symlink() or not incident.is_dir() or not record_path.is_file():
            raise ValueError("unidentified manifest-incident evidence blocks selected diagnostics")
        try:
            record = json.loads(record_path.read_bytes())
        except (ValueError, OSError) as exc:
            raise ValueError("unreadable manifest-incident evidence blocks selected diagnostics") from exc
        if not isinstance(record, dict) or not isinstance(record.get("manifest_path"), str):
            raise ValueError("unbound manifest-incident evidence blocks selected diagnostics")
        if Path(record["manifest_path"]).resolve() == output / "run_metadata.json":
            raise ValueError(f"preserved manifest incident blocks automatic mutation: {incident}")


def _reject_archive_ancestors(output):
    """Do not place new evidence inside any previously saved experiment."""
    markers = {"run_metadata.json", "cell_manifest.json", "resolved_design.json", "diagnostic_design.json",
               "table_manifest.json", "figure_manifest.json", "artifact_manifest.json", "run_manifest.json"}
    for ancestor in output.parents:
        if any((ancestor / name).is_file() for name in markers):
            raise ValueError(f"diagnostic output is nested in a saved experiment: {ancestor}")


def _resume_preflight(output, selection, selection_bytes, design, originals, fresh, require_frozen):
    _check_output_layout(output, selection)
    path = output / "run_metadata.json"
    if not path.is_file():
        raise ValueError("cannot resume without a selected diagnostics manifest")
    raw_metadata = path.read_bytes()
    metadata = json.loads(raw_metadata)
    if (raw_metadata != (canonical(metadata) + "\n").encode()
            or type(metadata.get("schema_version")) is not int or metadata["schema_version"] != 1
            or metadata.get("project") != "ou-week7-path-diagnostics"):
        raise ValueError("unknown selected diagnostics manifest")
    if metadata.get("status") in _CLOSED_STATES:
        raise ValueError("closed selected diagnostics are immutable; resume is not permitted")
    if metadata.get("status") not in _OPEN_STATES:
        raise ValueError("unknown selected diagnostics resume state")
    if metadata.get("require_frozen") is not require_frozen:
        raise ValueError("resume freeze policy changed")
    _check_frozen(metadata["provenance"], design, fresh, require_frozen)
    if (output / "selection.json").read_bytes() != selection_bytes:
        raise ValueError("literal selected input bytes changed")
    expected_design = (design.serialized + "\n").encode()
    if (output / "diagnostic_design.json").read_bytes() != expected_design:
        raise ValueError("persisted diagnostic descriptor changed")
    files = metadata.get("source_files")
    if not isinstance(files, dict):
        raise ValueError("missing selected-source digest map")
    ordered = [_filename(item) for item in selection["sources"]]
    if set(files) != set(ordered[:len(files)]):
        raise ValueError("committed selected sources must be the deterministic prefix")
    existing = {p.name for p in (output / "sources").iterdir()} if (output / "sources").exists() else set()
    if not set(files) <= existing:
        raise ValueError("committed selected-source evidence is missing")
    orphans = existing - set(files)
    if orphans and (len(orphans) != 1 or len(files) == len(ordered) or orphans != {ordered[len(files)]}):
        raise ValueError("only the next uncommitted selected-source orphan can be considered")
    for item in selection["sources"]:
        filename = _filename(item)
        if filename in existing:
            source_path = output / "sources" / filename
            if filename in files and file_digest(source_path) != files[filename]:
                raise ValueError("committed selected-source digest mismatch")
            _payload_identity(json.loads(source_path.read_text()), item,
                              originals[_source_key(item)], metadata["provenance"]["identity_sha256"])
    audit_path = output / "audit.json"
    if metadata.get("audit_sha256") is not None:
        if not audit_path.is_file() or file_digest(audit_path) != metadata["audit_sha256"]:
            raise ValueError("bound selected-source audit changed")
    if audit_path.exists():
        audit = json.loads(audit_path.read_text())
        if (audit.get("run_identity") != metadata["provenance"]["identity_sha256"]
                or audit.get("source_files") != files):
            raise ValueError("existing selected-source audit binding changed")
    return metadata


def _one_selected(item, original, design, run_identity):
    """Replay the identical original source and fits before any diagnostic math."""
    payload = {"schema_version": 1, "run_identity": run_identity, "selection": item,
               "original_records": original, "headline_results_replaced": False}
    attempted = {}
    try:
        cell = next(cell for cell in design.cells if cell.identity == item["cell_id"])
        source = OuSampleSource(design.seed)
        sample, source_diagnostics = source.sample(cell.to_dict(), item["rep"])
        source_diagnostics = _safe(source_diagnostics)
        actual_hashes = {"times_sha256": array_digest(sample.times),
                         "observations_sha256": array_digest(sample.values)}
        attempted.update(source_diagnostics=source_diagnostics, array_hashes=actual_hashes)
        for method in ESTIMATORS:
            row = original[method]
            _compare_replay({name: row[name] for name in source_diagnostics}, source_diagnostics,
                            method + ".source_diagnostics")
            _compare_replay({name: row[name] for name in actual_hashes}, actual_hashes,
                            method + ".array_hashes")
        replay_rows = source_records(cell, item["rep"], source,
                                     {name: registered_wald(name) for name in ESTIMATORS})
        if len(replay_rows) != len(ESTIMATORS):
            raise ValueError("original replay lost planned estimator records")
        replay = {row["estimator"]: _safe(row) for row in replay_rows}
        attempted["records"] = replay
        _compare_replay(original, replay)
        payload["original_replay"] = {"status": "passed", "records": replay,
                                      "source_diagnostics": source_diagnostics,
                                      "array_hashes": actual_hashes,
                                      "rtol": SETTINGS["replay_rtol"], "atol": SETTINGS["replay_atol"]}
    except Exception as exc:
        payload.update(status="original_replay_failed", original_replay={
            "status": "failed", "error": f"{type(exc).__name__}: {exc}", "attempted": attempted})
        return _safe(payload)
    try:
        payload["diagnostic"] = analyze_sample(sample, original)
        payload["status"] = "evaluated"
    except Exception as exc:
        payload.update(status="diagnostic_failure", diagnostic={
            "status": "failed", "error": f"{type(exc).__name__}: {exc}"})
    return _safe(payload)


def _audit_selected_files(output, selection, originals, metadata):
    """Readback identity/replay audit; does not pretend to rederive all math."""
    output = Path(output)
    _check_output_layout(output, selection, require_complete=True)
    descriptor = metadata["provenance"]["design"]
    descriptor_bytes = (canonical(descriptor) + "\n").encode()
    if (output / "diagnostic_design.json").read_bytes() != descriptor_bytes:
        raise ValueError("audited diagnostic descriptor changed")
    if file_digest(output / "selection.json") != descriptor["selection_bytes_sha256"]:
        raise ValueError("audited literal selection bytes changed")
    if json.loads((output / "selection.json").read_bytes()) != selection:
        raise ValueError("audited semantic selection changed")
    expected = {_filename(item) for item in selection["sources"]}
    if set(metadata["source_files"]) != expected:
        raise ValueError("source manifest is not the complete selected Cartesian set")
    failures, replayed, evaluated = [], 0, 0
    for item in selection["sources"]:
        filename = _filename(item)
        path = output / "sources" / filename
        if file_digest(path) != metadata["source_files"][filename]:
            raise ValueError("audited selected-source digest mismatch")
        payload = json.loads(path.read_text())
        _payload_identity(payload, item, originals[_source_key(item)],
                          metadata["provenance"]["identity_sha256"])
        if payload["original_replay"]["status"] == "passed":
            replayed += len(ESTIMATORS)
        if payload["status"] == "evaluated":
            evaluated += 1
        else:
            failures.append({"selection": item, "status": payload["status"]})
    return {"status": "review_required" if failures else "passed",
            "scope": "selected_source_and_full_fit_replay_then_immutable_evidence_readback",
            "independent_recalculation_of_all_diagnostic_math": False,
            "run_identity": metadata["provenance"]["identity_sha256"],
            "diagnostic_design_sha256": file_digest(output / "diagnostic_design.json"),
            "selection_bytes_sha256": file_digest(output / "selection.json"),
            "planned_sources": len(selection["sources"]), "retained_sources": len(expected),
            "original_fit_records": len(expected) * len(ESTIMATORS),
            "verified_fit_replays": replayed, "evaluated_sources": evaluated,
            "source_files": dict(metadata["source_files"]), "failures": failures,
            "headline_results_replaced": False}


def run_selected(selection, e1, e2, output, require_frozen=True, resume=False, progress=True):
    """Durably inspect the fixed selection; full runs require a repository freeze.

    ``require_frozen=False`` exists solely for distinct synthetic/development
    fixtures.  It rejects any headline design, and is not an escape hatch for
    exploratory replay of headline observations before publication.  Closed
    states are immutable; failures in committed source files are never retried.
    """
    if type(require_frozen) is not bool or type(resume) is not bool:
        raise ValueError("freeze and resume flags must be booleans")
    selected, selection_bytes, selection_dependency = _read_selection(selection)
    if require_frozen and selection_dependency != "artifacts/week7_evidence/selection.json":
        raise ValueError("frozen selected replay requires the fixed repository selection file")
    if not isinstance(selected.get("sources"), list) or not selected["sources"]:
        raise ValueError("a nonempty fixed selected source list is required")
    if len({_source_key(item) for item in selected["sources"]}) != len(selected["sources"]):
        raise ValueError("duplicate selected source")
    designs = {label: load_saved_design(directory) for label, directory in (("E1", e1), ("E2", e2))}
    if not require_frozen and any(design.scope == "headline" for design in designs.values()):
        raise ValueError("headline source replay cannot run with unfrozen dependencies")
    output = Path(output).resolve()
    _check_incident_block(output)
    _reject_archive_ancestors(output)
    for directory in (Path(e1).resolve(), Path(e2).resolve()):
        if output == directory or output.is_relative_to(directory) or directory.is_relative_to(output):
            raise ValueError("diagnostic output must be disjoint from original headline directories")
    design = _diagnostic_design(selected, selection_bytes)
    extra_files = (*FILES, *((selection_dependency,) if selection_dependency else ()))
    # Capture and verify the full frozen file set BEFORE output creation or any
    # regeneration of observations.  Reading the sealed inputs below is inert.
    fresh = capture(design, extra_files, require_frozen=require_frozen)
    _check_frozen(fresh, design, fresh, require_frozen)
    current = build_selection(e1, e2)
    if canonical(current) != canonical(selected):
        raise ValueError("the supplied fixed selection differs from sealed inputs")
    originals = _selected_originals(selected, e1, e2)
    _check_output_layout(output, selected)
    if output.exists() and any(output.iterdir()):
        if not resume:
            raise FileExistsError("existing diagnostic evidence requires explicit resume")
        _resume_preflight(output, selected, selection_bytes, design, originals, fresh, require_frozen)
    elif resume:
        raise ValueError("cannot resume an absent selected diagnostics manifest")
    output.mkdir(parents=True, exist_ok=True)
    with run_lock(output):
        _check_incident_block(output)
        path = output / "run_metadata.json"
        if path.exists():
            if not resume:
                raise FileExistsError("diagnostic manifest appeared before execution lock")
            metadata = _resume_preflight(output, selected, selection_bytes, design,
                                         originals, fresh, require_frozen)
            guard = _ManifestGuard(path, (canonical(metadata) + "\n").encode())
            metadata["resume_events"].append({"utc": utc_now()})
        else:
            if any(p.name != ".execution.lock" for p in output.iterdir()):
                raise ValueError("unidentified selected diagnostic output")
            metadata = {"schema_version": 1, "project": "ou-week7-path-diagnostics",
                        "status": "running_incomplete", "started_utc": utc_now(),
                        "provenance": fresh, "require_frozen": require_frozen,
                        "source_files": {}, "resume_events": [], "recovered_orphans": []}
            guard = _ManifestGuard(path, None)
            atomic_bytes(output / "selection.json", selection_bytes, immutable=True)
            atomic_bytes(output / "diagnostic_design.json", (design.serialized + "\n").encode(), immutable=True)
        metadata["status"] = "running_incomplete"
        guard.write(metadata, "attempt_start")
        try:
            for index, item in enumerate(selected["sources"]):
                filename = _filename(item)
                guard.check("before_selected_source:" + filename)
                if filename in metadata["source_files"]:
                    continue
                payload = _one_selected(item, originals[_source_key(item)], designs[item["experiment"]],
                                        metadata["provenance"]["identity_sha256"])
                source_path = output / "sources" / filename
                guard.check("before_source_evidence:" + filename)
                created = atomic_json(source_path, payload, immutable=True)
                if not created:
                    metadata["recovered_orphans"].append({"filename": filename, "utc": utc_now(),
                        "method": "exact_byte_equal_deterministic_replay"})
                metadata["source_files"][filename] = file_digest(source_path)
                guard.write(metadata, "source_commit:" + filename)
                if progress:
                    print(f"Week 7 selected source {index + 1}/{len(selected['sources'])}: {payload['status']}", flush=True)
            _check_frozen(metadata["provenance"], design, fresh, require_frozen)
            if canonical(build_selection(e1, e2)) != canonical(selected):
                raise ValueError("sealed headline inputs changed during selected diagnostics")
            report = _audit_selected_files(output, selected, originals, metadata)
            guard.check("before_audit_commit")
            atomic_json(output / "audit.json", report, immutable=True)
            metadata.update(status="completed_diagnostic_evidence" if report["status"] == "passed"
                            else "diagnostic_review_required", audit_sha256=file_digest(output / "audit.json"),
                            finished_utc=utc_now())
        except BaseException as exc:
            if not guard.drifted:
                metadata.update(status="incomplete", termination_type=type(exc).__name__, termination_reason=str(exc))
                guard.write(metadata, "attempt_failure")
            raise
        else:
            guard.write(metadata, "attempt_finish")
    return metadata
