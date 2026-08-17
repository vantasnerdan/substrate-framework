"""Exact stabilizer Lie algebras for standard 2+1D and 3+1D momenta.

Authority status: implements accepted claim C-LOR-002 in release v0.160.0.
The ledgers establish only matrix membership, fixed-vector conditions, and
commutators in the displayed mostly-plus conventions. They do not classify
global unitary representations or derive spin, anyons, continuous-spin states,
helicity, parity pairing, photons, gravitons, or interacting dynamics.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


MatrixTuple = tuple[sp.ImmutableMatrix, ...]


@dataclass(frozen=True)
class LorentzLittleGroupLedger:
    """Named exact generators stabilizing standard massive and null vectors."""

    spacetime_dimension: int
    minkowski_metric: sp.ImmutableMatrix
    standard_massive_momentum: sp.ImmutableMatrix
    standard_null_momentum: sp.ImmutableMatrix
    massive_generator_names: tuple[str, ...]
    massive_generators: MatrixTuple
    nonstabilizing_boost_names: tuple[str, ...]
    nonstabilizing_boosts: MatrixTuple
    massless_generator_names: tuple[str, ...]
    massless_generators: MatrixTuple


def _immutable(matrix: sp.Matrix) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(matrix)


def _boost(dimension: int, axis: int) -> sp.ImmutableMatrix:
    matrix = sp.zeros(dimension)
    matrix[0, axis] = 1
    matrix[axis, 0] = 1
    return _immutable(matrix)


def _validate_ledger(ledger: LorentzLittleGroupLedger) -> None:
    eta = ledger.minkowski_metric
    for generator in ledger.massive_generators + ledger.massless_generators:
        if sp.simplify(eta * generator + (eta * generator).T) != sp.zeros(
            ledger.spacetime_dimension
        ):
            raise AssertionError("stabilizer generator is not Lorentz-algebra valued")
    for generator in ledger.massive_generators:
        if generator * ledger.standard_massive_momentum != sp.zeros(
            ledger.spacetime_dimension, 1
        ):
            raise AssertionError("massive generator does not fix rest momentum")
    for boost in ledger.nonstabilizing_boosts:
        if boost * ledger.standard_massive_momentum == sp.zeros(
            ledger.spacetime_dimension, 1
        ):
            raise AssertionError("listed nonstabilizing boost fixes rest momentum")
    for generator in ledger.massless_generators:
        if generator * ledger.standard_null_momentum != sp.zeros(
            ledger.spacetime_dimension, 1
        ):
            raise AssertionError("massless generator does not fix null momentum")


def little_group_algebra_2plus1() -> LorentzLittleGroupLedger:
    """Return the displayed ``so(2)`` and one-generator null stabilizers."""

    eta = _immutable(sp.diag(-1, 1, 1))
    rotation = _immutable(sp.Matrix([[0, 0, 0], [0, 0, 1], [0, -1, 0]]))
    boost1 = _boost(3, 1)
    boost2 = _boost(3, 2)
    null_rotation = _immutable(boost2 + rotation)
    ledger = LorentzLittleGroupLedger(
        spacetime_dimension=3,
        minkowski_metric=eta,
        standard_massive_momentum=_immutable(sp.Matrix([1, 0, 0])),
        standard_null_momentum=_immutable(sp.Matrix([1, 1, 0])),
        massive_generator_names=("J12",),
        massive_generators=(rotation,),
        nonstabilizing_boost_names=("K1", "K2"),
        nonstabilizing_boosts=(boost1, boost2),
        massless_generator_names=("M=K2+J12",),
        massless_generators=(null_rotation,),
    )
    _validate_ledger(ledger)
    return ledger


def little_group_algebra_3plus1() -> LorentzLittleGroupLedger:
    """Return the displayed ``so(3)`` and ``iso(2)`` stabilizer generators."""

    eta = _immutable(sp.diag(-1, 1, 1, 1))
    j1 = _immutable(
        sp.Matrix(
            [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]]
        )
    )
    j2 = _immutable(
        sp.Matrix(
            [[0, 0, 0, 0], [0, 0, 0, -1], [0, 0, 0, 0], [0, 1, 0, 0]]
        )
    )
    j3 = _immutable(
        sp.Matrix(
            [[0, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0], [0, 0, 0, 0]]
        )
    )
    boosts = tuple(_boost(4, axis) for axis in (1, 2, 3))
    t1 = _immutable(boosts[0] + j2)
    t2 = _immutable(boosts[1] - j1)
    ledger = LorentzLittleGroupLedger(
        spacetime_dimension=4,
        minkowski_metric=eta,
        standard_massive_momentum=_immutable(sp.Matrix([1, 0, 0, 0])),
        standard_null_momentum=_immutable(sp.Matrix([1, 0, 0, 1])),
        massive_generator_names=("J1", "J2", "J3"),
        massive_generators=(j1, j2, j3),
        nonstabilizing_boost_names=("K1", "K2", "K3"),
        nonstabilizing_boosts=boosts,
        massless_generator_names=("T1=K1+J2", "T2=K2-J1", "J3"),
        massless_generators=(t1, t2, j3),
    )
    _validate_ledger(ledger)
    if t1 * t2 - t2 * t1 != sp.zeros(4):
        raise AssertionError("null translations do not commute")
    if j3 * t1 - t1 * j3 != -t2:
        raise AssertionError("first iso(2) commutator has the wrong sign")
    if j3 * t2 - t2 * j3 != t1:
        raise AssertionError("second iso(2) commutator has the wrong sign")
    return ledger
