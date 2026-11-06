"""Exact residual predicates and evidence ledger; no floating tolerances."""
import sympy as s


def reduce(value):
    if isinstance(value, s.MatrixBase):
        return value.applyfunc(s.simplify)
    return s.simplify(value)


def zero(value):
    value = reduce(value)
    return value == s.zeros(*value.shape) if isinstance(value, s.MatrixBase) else value == 0


class Evidence:
    def __init__(self):
        self.rows = []

    def identity(self, name, residual):
        result = reduce(residual)
        self.rows.append({"name": name, "kind": "exact_identity", "residual": str(result)})
        assert zero(result), (name, result)

    def exposed(self, name, residual):
        result = reduce(residual)
        self.rows.append({"name": name, "kind": "nonzero_counterexample", "witness": str(result)})
        assert not zero(result), (name, result)


def native(hp, sp, energies):
    """Inputs are matrices already expressed in the reference eigenbasis."""
    return hp - sp * energies
