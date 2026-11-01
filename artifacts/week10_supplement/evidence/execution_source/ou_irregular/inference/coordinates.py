"""Observed-data coordinates for numerically stable, unit-equivariant inference.

Write x=a+b*y and t=m*u, using the observed state mean/RMS and mean
gap. Then theta_x=theta_y/m, mu_x=a+b*mu_y and
sigma_x=b*sigma_y/sqrt(m). No simulation truth enters this transformation.
The likelihood changes only by an additive data Jacobian, so differentiating
the normalized conditional likelihood yields the same statistical information.
"""
from dataclasses import dataclass

import numpy as np

from .contracts import Sample


@dataclass(frozen=True)
class ObservationCoordinates:
    offset: float
    state_scale: float
    time_scale: float
    normalized_sample: Sample

    @classmethod
    def from_sample(cls, sample: Sample):
        values = sample.values.astype(np.longdouble)
        offset = float(values.mean())
        centered = values - np.longdouble(offset)
        amplitude = np.max(np.abs(centered))
        if not np.isfinite(amplitude) or amplitude <= 0:
            raise ValueError("constant observed path has no positive state scale")
        scale = amplitude * np.sqrt(np.mean((centered / amplitude) ** 2))
        if not np.isfinite(scale) or not 0 < scale <= np.finfo(float).max:
            raise ValueError("observed state scale is not representable")
        scale = float(scale)
        if scale == 0:
            raise ValueError("observed state scale underflowed")
        # Normalization cannot recover information already lost by rounding a
        # very small fluctuation onto a much larger measurement offset.
        resolution = np.max(np.abs(np.spacing(np.abs(sample.values))))
        if resolution / scale > 1e-6:
            raise ValueError("insufficient floating-point state resolution")
        gaps = np.diff(sample.times)
        mean_gap = float(gaps.astype(np.longdouble).mean())
        if not np.isfinite(mean_gap) or mean_gap <= 0:
            raise ValueError("observed mean gap is not positive and finite")
        # Subtract the time origin in extended precision before division. This
        # retains the given gaps, including any declared spacing floor.
        times = np.asarray((sample.times.astype(np.longdouble) - sample.times[0])
                           / np.longdouble(mean_gap), dtype=float)
        normalized = Sample(np.asarray(centered / scale, dtype=float), times)
        return cls(offset, scale, mean_gap, normalized)

    def to_working(self, parameters):
        """Natural input -> (log theta_y, mu_y, log sigma_y), without products."""
        parameters = np.asarray(parameters, dtype=float)
        if (parameters.shape != (3,) or not np.isfinite(parameters).all()
                or parameters[0] <= 0 or parameters[2] <= 0):
            raise ValueError("finite theta, mu, sigma with theta,sigma>0 required")
        theta, mu, sigma = parameters
        working = np.array([
            np.log(theta) + np.log(self.time_scale),
            float((np.longdouble(mu) - self.offset) / self.state_scale),
            np.log(sigma) + 0.5 * np.log(self.time_scale) - np.log(self.state_scale),
        ])
        if not np.isfinite(working).all():
            raise ValueError("working parameters are not representable")
        return working

    def to_natural(self, working):
        working = np.asarray(working, dtype=float)
        if working.shape != (3,) or not np.isfinite(working).all():
            raise ValueError("three finite working parameters required")
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            parameters = np.array([
                np.exp(working[0] - np.log(self.time_scale)),
                float(np.longdouble(self.offset)
                      + np.longdouble(self.state_scale) * working[1]),
                np.exp(working[2] + np.log(self.state_scale)
                       - 0.5 * np.log(self.time_scale)),
            ])
        if (not np.isfinite(parameters).all()
                or parameters[0] <= 0 or parameters[2] <= 0):
            raise ValueError("natural parameters are not representable")
        return parameters

    def natural_jacobian(self, working):
        """Delta-method derivative d(theta_x,mu_x,sigma_x)/d eta_y."""
        theta, _, sigma = self.to_natural(working)
        return np.diag([theta, self.state_scale, sigma])
