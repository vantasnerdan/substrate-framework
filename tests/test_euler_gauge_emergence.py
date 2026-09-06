import pytest
import sympy as sp

from substrate_framework.euler_gauge_emergence import (
    clebsch_shift_kernel,
    transverse_scalar_source,
    uniform_symbol_ledger,
)


def test_uniform_euler_and_maxwell_characteristics_are_distinct():
    lam, k, adv, c = sp.symbols("lambda k adv c", nonzero=True)
    ledger = uniform_symbol_ledger(lam, k, adv, c)
    assert ledger.euler_characteristic == sp.expand((lam + sp.I * adv) ** 2)
    assert ledger.maxwell_characteristic == sp.expand((lam**2 + c**2 * k**2) ** 2)
    assert ledger.euler_tag_characteristic == sp.expand((lam + sp.I * adv) ** 3)
    assert ledger.maxwell_tag_characteristic == sp.expand(
        (lam + sp.I * adv) * (lam**2 + c**2 * k**2) ** 2
    )
    assert ledger.maxwell_temporal_frequencies == (-c * k, c * k)


def test_local_isotropic_scalar_source_is_longitudinal():
    projector, projected = transverse_scalar_source([1, 2, 3], sp.Symbol("q"))
    assert projector * sp.Matrix([1, 2, 3]) == sp.zeros(3, 1)
    assert projected == sp.zeros(3, 1)


def test_clebsch_shift_is_a_physical_kernel_direction():
    assert clebsch_shift_kernel([1, 2, 3], sp.Symbol("fp")) == sp.zeros(3, 1)
    with pytest.raises(ValueError):
        clebsch_shift_kernel([1, 2], 1)
    with pytest.raises(ValueError):
        transverse_scalar_source([0, 0, 0], 1)
