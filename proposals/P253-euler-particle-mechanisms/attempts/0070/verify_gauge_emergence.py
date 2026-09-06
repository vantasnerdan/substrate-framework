import sympy as sp

from substrate_framework.euler_gauge_emergence import (
    clebsch_shift_kernel,
    transverse_scalar_source,
    uniform_symbol_ledger,
)


lam, k, adv, c, q, fp = sp.symbols("lambda k adv c q fp", nonzero=True)
ledger = uniform_symbol_ledger(lam, k, adv, c)

assert ledger.euler_characteristic == sp.expand((lam + sp.I * adv) ** 2)
assert ledger.maxwell_characteristic == sp.expand((lam**2 + c**2 * k**2) ** 2)
assert ledger.euler_temporal_frequencies == (adv, adv)
assert ledger.maxwell_temporal_frequencies == (-c * k, c * k)

assert sp.simplify(
    ledger.maxwell_characteristic.subs(lam, -sp.I * adv)
) != 0

projector, source = transverse_scalar_source([1, 2, 3], q)
assert projector * sp.Matrix([1, 2, 3]) == sp.zeros(3, 1)
assert source == sp.zeros(3, 1)

assert clebsch_shift_kernel([1, 2, 3], fp) == sp.zeros(3, 1)

print("P253/0070 exact checks: 8 passed")
