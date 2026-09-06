import numpy as np
import pytest
import sympy as sp

from substrate_framework.euler_schwinger_hopf import (
    apply_mode_mixing,
    canonical_doublet,
    hopf_identity_residual,
    reduced_kks_area,
    stokes_vector,
    total_action,
)


def test_symbolic_stokes_poisson_algebra_and_casimir():
    q1, p1, q2, p2, B = sp.symbols("q1 p1 q2 p2 B", positive=True, real=True)
    root = sp.sqrt(B / 2)
    z1 = root * (q1 + sp.I * p1)
    z2 = root * (q2 + sp.I * p2)
    sx = sp.re(sp.conjugate(z1) * z2).expand(complex=True)
    sy = sp.im(sp.conjugate(z1) * z2).expand(complex=True)
    sz = sp.expand((z1 * sp.conjugate(z1) - z2 * sp.conjugate(z2)) / 2)
    action = sp.expand((z1 * sp.conjugate(z1) + z2 * sp.conjugate(z2)) / 2)
    variables = ((q1, p1), (q2, p2))

    def bracket(left, right):
        return sp.simplify(
            sum(
                sp.diff(left, q) * sp.diff(right, p)
                - sp.diff(left, p) * sp.diff(right, q)
                for q, p in variables
            )
            / B
        )

    assert sp.simplify(bracket(sx, sy) - sz) == 0
    assert sp.simplify(bracket(sy, sz) - sx) == 0
    assert sp.simplify(bracket(sz, sx) - sy) == 0
    assert sp.simplify(sx**2 + sy**2 + sz**2 - action**2) == 0
    assert all(sp.simplify(bracket(action, component)) == 0 for component in (sx, sy, sz))


def test_hopf_map_and_unitary_mixing_preserve_total_action():
    z = canonical_doublet((1.0, -2.0, 0.5, 3.0), kks_scale=2.5)
    hadamard = np.array([[1, 1], [-1, 1]], dtype=complex) / np.sqrt(2)
    mixed = apply_mode_mixing(z, hadamard)
    assert total_action(mixed) == pytest.approx(total_action(z))
    assert hopf_identity_residual(z) == pytest.approx(0.0, abs=1e-12)
    assert hopf_identity_residual(mixed) == pytest.approx(0.0, abs=1e-12)
    assert np.linalg.norm(stokes_vector(mixed)) == pytest.approx(total_action(z))


def test_common_phase_is_exact_hopf_fibre():
    z = np.array([1 + 2j, -3 + 0.5j])
    rotated = np.exp(0.731j) * z
    assert stokes_vector(rotated) == pytest.approx(stokes_vector(z))
    assert total_action(rotated) == pytest.approx(total_action(z))


def test_domain_and_degenerate_sphere_rejections():
    with pytest.raises(ValueError):
        canonical_doublet((1, 2, 3, 4), kks_scale=0)
    with pytest.raises(ValueError):
        reduced_kks_area(0)
    with pytest.raises(ValueError):
        apply_mode_mixing((1, 0), [[1, 1], [0, 1]])
    with pytest.raises(ValueError):
        total_action((1, complex(float("nan"), 0)))
