"""Derivative checks use objectives independent of the production likelihood."""
import numpy as np
import pytest

from ou_irregular.inference.hessian import FiniteDifferenceHessian


def test_quadratic_hessian_and_gradient_match_closed_form():
    matrix = np.array([[5.0, 0.4, -0.2], [0.4, 2.0, 0.1], [-0.2, 0.1, 3.0]])
    center = np.array([0.2, -1.0, 2.5])
    point = np.array([-0.5, 0.4, 3.0])
    objective = lambda p: 0.5 * (p - center) @ matrix @ (p - center) + 1.7
    result = FiniteDifferenceHessian().evaluate(objective, point)
    assert result.valid, result.reason
    np.testing.assert_allclose(result.hessian, matrix, atol=2e-8)
    np.testing.assert_allclose(result.gradient, matrix @ (point - center), atol=1e-10)
    assert result.diagnostics["whitened_relative_difference"] < 1e-7


def test_nonquadratic_derivatives_match_analytic_expression():
    point = np.array([0.3, -0.2, 0.5])
    objective = lambda p: np.exp(p[0]) + p[1] ** 4 + p[1] ** 2 + 2 * p[2] ** 2
    result = FiniteDifferenceHessian().evaluate(objective, point)
    assert result.valid, result.reason
    np.testing.assert_allclose(result.hessian, np.diag([np.exp(point[0]), 12 * point[1] ** 2 + 2, 4]), rtol=1e-6, atol=1e-7)
    np.testing.assert_allclose(result.gradient, [np.exp(point[0]), 4 * point[1] ** 3 + 2 * point[1], 4 * point[2]], atol=3e-7)


@pytest.mark.parametrize("point", [[0, 1], [0, np.nan, 0], [0, np.inf, 1]])
def test_invalid_coordinates_are_explicit(point):
    result = FiniteDifferenceHessian().evaluate(lambda p: np.dot(p, p), point)
    assert not result.valid
    assert result.status == "invalid_coordinates"


@pytest.mark.parametrize("objective", [lambda p: np.inf, lambda p: np.nan])
def test_nonfinite_objective_produces_invalid_derivatives(objective):
    result = FiniteDifferenceHessian().evaluate(objective, np.zeros(3))
    assert not result.valid
    assert result.status == "nonfinite_derivatives"


def test_unstable_curvature_is_not_silently_used():
    # At the origin f=x^4 gives an apparent positive second difference
    # proportional to h^2; halving h changes it fourfold.
    result = FiniteDifferenceHessian().evaluate(lambda p: np.sum(p ** 4), np.zeros(3))
    assert not result.valid
    assert result.status == "unstable_hessian"


def test_whitened_gate_detects_weak_direction_hidden_by_global_norm():
    # Strong y,z curvature makes the global discrepancy ~1e-6; the x
    # curvature is dominated by discretization and changes ~threefold.
    objective = lambda p: 0.5 * (1e-8 * p[0] ** 2 + p[1] ** 2 + p[2] ** 2) + p[0] ** 4
    result = FiniteDifferenceHessian().evaluate(objective, np.zeros(3))
    assert result.relative_difference < 0.005
    assert result.diagnostics["whitened_relative_difference"] > 0.01
    assert not result.valid
    assert result.status == "unstable_hessian"


@pytest.mark.parametrize("diagonal", [[0, 1, 1], [-1, 1, 1]])
def test_nonpositive_hessian_rejected(diagonal):
    matrix = np.diag(diagonal)
    result = FiniteDifferenceHessian().evaluate(lambda p: p @ matrix @ p / 2, np.zeros(3))
    assert not result.valid
    assert result.status == "nonpositive_hessian"


@pytest.mark.parametrize("field", ["relative_step", "max_relative_difference", "max_whitened_difference"])
@pytest.mark.parametrize("value", [0, -1, np.nan, np.inf, True])
def test_hessian_settings_reject_invalid_values(field, value):
    with pytest.raises(ValueError):
        FiniteDifferenceHessian(**{field: value})


def test_derivative_result_arrays_and_diagnostics_are_read_only():
    result = FiniteDifferenceHessian().evaluate(lambda p: np.dot(p, p), np.zeros(3))
    with pytest.raises(ValueError):
        result.hessian[0, 0] = 1
    with pytest.raises(ValueError):
        result.gradient[0] = 1
    with pytest.raises(TypeError):
        result.diagnostics["relative_step"] = 1
