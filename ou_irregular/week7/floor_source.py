"""Paired epsilon sensitivity through the unchanged exact OU source.

The physical cell and its stored identity retain their actual floor.  Only the
cell ID used to address random streams is replaced by the same-CV registered-
floor anchor.  The three variants therefore share raw Gamma draws, the normal
draw for stationary X0 and the subsequent standard normal innovations.  These
are coupled exact OU realizations; variants coincide when their clipped gap
vectors coincide.  They are not one continuous path observed on three grids.

This adapter implements only the twelve physical cells fixed in the dated
Week 7 protocol.  A run design, rather than this source, controls its seed,
replication count, validation namespace and completed-key accounting.
"""
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
import math

from ..experiments.design import Cell
from ..week4.config import canonical_json
from ..week4.simulation import OuSampleSource


ANCHOR_FLOOR_FRACTION = 1e-6
FLOOR_FRACTIONS = (1e-8, 1e-6, 1e-4)
COEFFICIENTS_OF_VARIATION = (0., 1.5, 2., 3.)
COUPLING_POLICY = "week7_paired_floor_v1_shared_gamma_stationary_x0_and_standard_normal_innovations"
_CANONICAL_CELLS = tuple(
    Cell(n=1000, q=.5, cv=cv, theta=1., mu=0., sigma=.5,
         floor_fraction=epsilon, floor_cv_threshold=1.5)
    for cv in COEFFICIENTS_OF_VARIATION for epsilon in FLOOR_FRACTIONS)
_CELL_IDENTITIES = frozenset(cell.identity for cell in _CANONICAL_CELLS)
_FIELDS = frozenset(_CANONICAL_CELLS[0].to_dict())
_FLOAT_FIELDS = _FIELDS - {"cell_id", "n"}


def validate_cell(cell):
    """Require an exact canonical actual Cell mapping in the frozen 12-cell grid.

    Numeric coercion is deliberately absent.  ``Cell.to_dict()`` emits a Python
    int for n and Python floats for all scientific scalars; accepting booleans,
    strings or ambiguous aliases here could make a recorded ID describe a
    different physical source.  The derived mean gap and supplied scientific ID
    must match the reconstructed Cell, including its *actual* floor fraction.
    """
    if not isinstance(cell, Mapping):
        raise ValueError("paired-floor source requires a canonical Cell mapping")
    supplied = dict(cell)
    if set(supplied) != _FIELDS:
        raise ValueError("paired-floor Cell fields must exactly match the canonical schema")
    if type(supplied["n"]) is not int or type(supplied["cell_id"]) is not str:
        raise ValueError("canonical n must be an integer and cell_id must be a string")
    for name in sorted(_FLOAT_FIELDS):
        if type(supplied[name]) is not float or not math.isfinite(supplied[name]):
            raise ValueError(f"canonical scientific scalar {name} must be a finite Python float")
    actual = Cell(n=supplied["n"], q=supplied["theta_mean_gap"], cv=supplied["cv"],
                  theta=supplied["true_theta"], mu=supplied["true_mu"], sigma=supplied["true_sigma"],
                  floor_fraction=supplied["floor_fraction"],
                  floor_cv_threshold=supplied["floor_cv_threshold"])
    # Canonical JSON compares the exact identity representation (including
    # signed zero); Python equality alone equates e.g. +0.0 and -0.0.
    if canonical_json(supplied) != canonical_json(actual.to_dict()):
        raise ValueError("actual Cell ID, mean gap or canonical representation is inconsistent")
    if actual.identity not in _CELL_IDENTITIES:
        raise ValueError("actual Cell is outside the twelve fixed Week 7 floor variants")
    return actual


def anchor_cell(cell):
    """Canonical epsilon=1e-6 seed anchor, holding every other cell field fixed."""
    actual = validate_cell(cell.to_dict() if isinstance(cell, Cell) else cell)
    return replace(actual, floor_fraction=ANCHOR_FLOOR_FRACTION)


@dataclass(frozen=True)
class PairedFloorSource:
    seed: int
    _source: OuSampleSource = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        # Reuse the historical seed contract and source rather than importing
        # simulator internals or copying their numerical transition recursion.
        object.__setattr__(self, "_source", OuSampleSource(self.seed))

    def sample(self, cell, rep):
        actual = validate_cell(cell)
        anchor = anchor_cell(actual)
        addressed = actual.to_dict()
        # OuSampleSource reads cell_id solely for SeedSequence addressing.  It
        # still receives the actual floor, n, mean gap, CV and truth parameters.
        # This mixed internal mapping is NOT a new Cell identity and is never
        # exported as the scientific cell stored in a result row.
        addressed["cell_id"] = anchor.identity
        sample, original_diagnostics = self._source.sample(addressed, rep)
        # A diagnostic cell_id would overwrite the actual row's identity when
        # experiments.execution.source_records composes its dictionaries.  Fail
        # closed on that unexpected upstream schema rather than losing evidence.
        reserved = _FIELDS | {"rep", "actual_cell_id", "seed_anchor_cell_id",
                              "seed_anchor_floor_fraction", "source_coupling_policy"}
        if not isinstance(original_diagnostics, Mapping) or set(original_diagnostics) & reserved:
            raise ValueError("underlying source diagnostics collide with paired-floor identity fields")
        diagnostics = dict(original_diagnostics)
        diagnostics.update(actual_cell_id=actual.identity, seed_anchor_cell_id=anchor.identity,
                           seed_anchor_floor_fraction=ANCHOR_FLOOR_FRACTION,
                           source_coupling_policy=COUPLING_POLICY)
        return sample, diagnostics
