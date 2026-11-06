"""Frozen dimensionless fixtures and exact-comparison/report helpers."""
import json
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
REGISTRATION = json.loads((HERE / 'registration.json').read_text())
A = sp.Matrix([[1, 1+sp.I, 2-sp.I], [1-sp.I, 2, sp.I], [2+sp.I, -sp.I, -1]])
B = sp.Matrix([[2, 2-sp.I, 1+2*sp.I], [2+sp.I, -1, 1-sp.I], [1-2*sp.I, 1+sp.I, 3]])
OPERATORS = (A, B)
ENERGIES = ((sp.S.Zero, sp.log(2), sp.log(5)),
            (sp.S.Zero, sp.S.Zero, sp.log(3)))
BETA = sp.S.One
HBAR = sp.Integer(2)
Z = (2+sp.I)/3


def zero(value):
    entries = list(value) if isinstance(value, sp.MatrixBase) else [value]
    return all(sp.simplify(sp.expand_log(x, force=True)) == 0 for x in entries)


def serial(value):
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(x) for x in value]
    if isinstance(value, sp.MatrixBase):
        return [[str(sp.simplify(value[i, j])) for j in range(value.cols)]
                for i in range(value.rows)]
    if isinstance(value, sp.Basic):
        return str(value)
    return value
