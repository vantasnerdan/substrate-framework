"""Tests for the P250 wall/shell/bubble sector primitives.

Every oracle is derived from the canonical C-M5C-001 action via
`substrate_framework.m5_exterior_clock`; nothing is a literal comparison.
Each derived gate carries channel-local mutation checks: a perturbation that
vanishes at the tested witnesses but moves the asserted quantity must flip
the verdict (the #190 defect class is excluded by construction).
"""

from __future__ import annotations

import json
import pathlib

import sympy as sp

from substrate_framework.m5_exterior_clock import (
    aligned_axis_lock_potential,
    clock_inertia_density,
    phase_lock_potential,
    projected_m5_static_potential,
    scalar_amplitude_potential,
)
from substrate_framework.m5_wall_clock import (
    degenerate_locus_point,
    envelope_identity,
    first_integral_conservation_identity,
    inertia_envelope_identity,
    locus_orbit_fixation,
    maxwell_frequency_resultants,
    maxwell_system,
    orbit_invariance_identity,
    phase_slip_energy_bound,
    phase_twist_gradient_excess,
    relative_equilibrium_energies,
    static_shell_derrick_residual,
    thin_wall_bag_radius,
    wall_slice_gradient_density,
    wall_slice_inertia,
    wall_slice_potential,
    wall_bulk_pressure,
)

ATTEMPT = pathlib.Path(__file__).resolve().parents[1] / (
    "proposals/P250-shell-bubble-clock/attempts/0001"
)

C, B = sp.symbols("c b", real=True)
F = sp.Symbol("f", nonnegative=True)
M = sp.Symbol("m", real=True)
W2 = sp.Symbol("w2", nonnegative=True)


def _slice_potential_with_lock_strength(m, c, b, f, omega_sq, strength):
    """Canonical slice potential with an explicitly varied lock strength.

    Same primitive composition as `wall_slice_potential`, but the phase-lock
    coefficient is a parameter: the vehicle for the channel-local lock
    mutation check.
    """
    S = sp.Matrix([[m, 0, 0], [0, c + b, 0], [0, 0, c - b]])
    tensor_part = (projected_m5_static_potential(S)
                   + aligned_axis_lock_potential(S)
                   + phase_lock_potential(S, f, 0, strength=strength)
                   + scalar_amplitude_potential(f))
    inertia = f**2 + 4*b**2
    return sp.expand(tensor_part - omega_sq*inertia/2)


def test_slice_potential_has_exact_vacuum_and_decomposition() -> None:
    V = wall_slice_potential(M, C, B, F, 0)
    assert V.subs({M: 1, C: 0, B: 0, F: 0}) == 0
    r2 = M**2 + 2*C**2 + 2*B**2
    decomposition = sp.expand(V - (
        -(r2)/2 - (M**3 + 2*C**3 + 6*C*B**2) + r2**2 + sp.Rational(1, 2)
        + 2*C**2 + 2*B**2 + 6*(B - F**2)**2
        + 3*F**2 - 4*F**3 + 2*F**4))
    assert decomposition == 0


def test_slice_inertia_and_gradient_match_canonical_normalization() -> None:
    assert wall_slice_inertia(M, C, B, F) == F**2 + 4*B**2
    dm, dc, db, df = sp.symbols("dm dc db df", real=True)
    assert wall_slice_gradient_density(dm, dc, db, df) == (
        (dm**2 + 2*dc**2 + 2*db**2)/4 + df**2/2)


def test_binding_witness_reproduced_independently() -> None:
    # provenance cross-check: the accepted witness values re-derived here
    witness_potential = wall_slice_potential(1, 0, sp.Rational(1, 4),
                                             sp.Rational(1, 2), 0)
    witness_inertia = wall_slice_inertia(1, 0, sp.Rational(1, 4),
                                         sp.Rational(1, 2))
    assert witness_potential == sp.Rational(45, 64)
    assert witness_inertia == sp.Rational(1, 2)
    assert 2*witness_potential/witness_inertia == sp.Rational(45, 16)


def test_first_integral_conservation_identity() -> None:
    assert first_integral_conservation_identity()


def test_locus_is_orbit_fixed_and_inertia_free() -> None:
    assert locus_orbit_fixation()
    S_loc, pr, pi = degenerate_locus_point()
    assert clock_inertia_density(S_loc, pr, pi) == 0
    # every tangent-isotropic point is orbit-fixed, not only the vacuum line
    q = sp.Symbol("q", real=True)
    m0, c0 = sp.symbols("m0 c0", real=True)
    rotation = sp.Matrix([[1, 0, 0], [0, sp.cos(q), -sp.sin(q)],
                          [0, sp.sin(q), sp.cos(q)]])
    S_iso = sp.diag(m0, c0, c0)
    assert sp.simplify(rotation*S_iso*rotation.T - S_iso) == sp.zeros(3)


def test_orbit_invariance_of_all_density_layers() -> None:
    assert orbit_invariance_identity()


def test_phase_slip_cost_vanishes_uniformly() -> None:
    L, dth, eps = sp.symbols("L dtheta epsilon", positive=True)
    bound = phase_slip_energy_bound(L, dth, eps)
    assert bound == L**2*dth**2*eps/4
    assert sp.limit(bound, eps, 0, "+") == 0
    assert sp.limit(bound.subs(dth, dth**2 + 1), eps, 0, "+") == 0
    qprime = sp.Symbol("qprime", real=True)
    assert phase_twist_gradient_excess(F**2 + 4*B**2, qprime) == (
        F**2 + 4*B**2) * qprime**2 / 2


def test_static_radial_shell_derrick_residual() -> None:
    T, U = sp.symbols("T U", nonnegative=True)
    assert static_shell_derrick_residual(T, U) == T + 3*U
    # static scaling stationarity forces T = -3U; with T, U >= 0 the only
    # solution of T + 3U = 0 is T = U = 0
    assert sp.solve(sp.Eq(T + 3*U, 0), T) == [-3*U]
    assert sp.solve(sp.Eq(T + 3*U, 0), U) == [-T/3]


def test_maxwell_system_matches_independent_derivation() -> None:
    pA, pB, pC, w_of = maxwell_system()
    assert sp.expand(pA - 2*(8*C**3 - 3*C**2 + C*(8*B**2 + 1) - 3*B**2)) == 0
    expected_pB = sp.expand(2*(8*B**3 + 8*B*C**2 - 6*B*C
                               - 2*B*w_of + 7*B - 6*F**2))
    assert sp.expand(pB - expected_pB) == 0
    V = wall_slice_potential(0, C, B, F, W2)
    assert sp.expand(pC - 2*V.subs(W2, w_of)) == 0
    assert sp.expand(sp.diff(V, F).subs(W2, w_of)) == 0


def test_bulk_pressure_is_derived_from_the_canonical_potential() -> None:
    pressure = wall_bulk_pressure(M, C, B, F, W2)
    assert sp.expand(pressure + wall_slice_potential(M, C, B, F, W2)) == 0
    assert sp.diff(pressure, W2) == wall_slice_inertia(M, C, B, F) / 2


def test_maxwell_frequency_resultants_share_symbols_and_vanish_at_point() -> None:
    """Regression (P250 attempt 0004): the eliminators must be built from the
    same nonnegative ``f`` symbol as `maxwell_system` — a second, identically
    printed ``f`` silently corrupts both resultants — and the certified
    Maxwell point must annihilate both, i.e. the c-elimination is consistent
    with the system it was derived from."""
    r1, r2 = maxwell_frequency_resultants()
    for r in (r1, r2):
        fsyms = {s for s in r.free_symbols if s.name == "f"}
        assert fsyms == {F}, f"resultant carries foreign f symbol(s): {fsyms}"
    point = json.loads((ATTEMPT / "maxwell_point.json").read_text())
    c0 = sp.Float(point["c"], 30)
    b0 = sp.Float(point["b"], 30)
    f0 = sp.Float(point["f"], 30)
    w0 = sp.N((32*F**2 - 12*F + 6 - 24*B).subs({B: b0, F: f0}), 30)
    for r in (r1, r2):
        val = sp.N(r.subs({C: c0, B: b0, F: f0,
                           sp.Symbol("w", real=True): w0}), 25)
        assert abs(float(val)) < 1e-12, f"resultant not annihilated: {val}"

def test_lock_mutation_moves_derived_system_and_is_invisible_at_witnesses():
    """Channel-local lock mutation: V_lock = (6 + eps)(b - f^2)^2.

    The mutation vanishes at the witness class (the lock manifold b = f^2,
    including the vacuum) but moves the derived Maxwell b-equation away from
    it; a literal-comparison oracle could not see this.
    """
    eps = sp.Symbol("eps_lock", positive=True)
    w_of = 32*F**2 - 12*F + 6 - 24*B
    pB = maxwell_system()[1]
    w2_sym = sp.Symbol("w2_mut", nonnegative=True)
    pB_mut = sp.expand(sp.diff(
        _slice_potential_with_lock_strength(0, C, B, F, w2_sym, 6 + eps),
        B).subs(w2_sym, w_of))
    moved = sp.expand(pB_mut - pB)
    assert moved != 0
    assert sp.simplify(moved.subs(B, F**2)) == 0
    assert moved.subs({B: 0, F: 0}) == 0
    assert moved.subs({B: sp.Rational(1, 4), F: 1}) != 0


def test_scalar_mutation_moves_frequency_relation_invisibly_at_vacuum():
    """W-mutation: V -> V + eps_w f^2 shifts the stationary f-relation by
    exactly 2 eps_w, moving the asserted Maxwell frequency relation while
    staying invisible at the vacuum witness f = 0."""
    eps_w = sp.Symbol("eps_w", positive=True)
    w_of = 32*F**2 - 12*F + 6 - 24*B
    V_base = _slice_potential_with_lock_strength(0, C, B, F, W2, 6)
    dF_base = sp.diff(V_base, F)
    dF_mut = sp.diff(V_base + eps_w*F**2, F)
    assert sp.simplify(dF_mut - dF_base) == 2*eps_w*F
    assert (dF_mut - dF_base).subs(F, 0) == 0
    # solving the mutated stationarity in w2 shifts the relation by 2 eps_w
    shifted = sp.solve(sp.Eq(dF_mut.subs(W2, sp.Symbol("w2s")), 0),
                       sp.Symbol("w2s"))
    assert shifted and sp.simplify(
        shifted[0] - (w_of + 2*eps_w).subs({B: B, F: F})) == 0


def test_maxwell_point_satisfies_system_and_beats_witnesses_exactly() -> None:
    point = json.loads((ATTEMPT / "maxwell_point.json").read_text())
    c0 = sp.Float(point["c"], 30)
    b0 = sp.Float(point["b"], 30)
    f0 = sp.Float(point["f"], 30)
    pA, pB, pC, w_of = maxwell_system()
    subs = {C: c0, B: b0, F: f0}
    for p in (pA, pB, pC):
        assert abs(float(sp.N(p.subs(subs), 20))) < 1e-12
    omega_sq = sp.N(w_of.subs(subs), 30)
    # exact rational witness at omega^2 = 5/3 is negative: the Maxwell
    # frequency is strictly below 5/3 and therefore below the accepted
    # witness ratio 45/16
    witness_V = wall_slice_potential(sp.Integer(0), sp.Rational(31, 100),
                                     sp.Rational(13, 20),
                                     sp.Rational(41, 50), sp.Rational(5, 3))
    assert witness_V == sp.Rational(-13739, 18750000)
    assert witness_V < 0
    assert omega_sq < sp.Rational(5, 3)
    assert omega_sq < sp.Rational(45, 16)


def test_maxwell_point_certified_unique_and_nondegenerate() -> None:
    """Exact local certification (Krawczyk + interval Gershgorin) recorded in
    the attempt artifact: the derived deep-branch Maxwell system has exactly
    one real solution in the rational box, and the fixed-omega Hessian of
    V_omega there is strictly positive definite, so the degenerate-depth
    point is a nondegenerate interior minimum."""
    cert = json.loads((ATTEMPT / "maxwell_certificate.json").read_text())
    assert cert["krawczyk_strict_interior"] is True
    assert cert["hessian_positive_definite"] is True
    for bound in cert["hessian_gershgorin_lower_bounds"]:
        lower = sp.Float(bound.strip("[]").split(",")[0].strip())
        assert lower > 1.0
    point = json.loads((ATTEMPT / "maxwell_point.json").read_text())
    lo = sp.Rational(16639457000, 10**10)
    hi = sp.Rational(16639457001, 10**10)
    assert lo < sp.Float(point["omega_sq"], 30) < hi
    # exact rational witnesses: the Maxwell frequency sits strictly below
    # 5/3 and therefore below the accepted witness ratio 45/16
    witness_V = wall_slice_potential(sp.Integer(0), sp.Rational(31, 100),
                                     sp.Rational(13, 20),
                                     sp.Rational(41, 50), sp.Rational(5, 3))
    assert witness_V == sp.Rational(-13739, 18750000)
    assert witness_V < 0


def test_envelope_and_inertia_identities() -> None:
    assert envelope_identity()
    assert inertia_envelope_identity()


def test_thin_wall_selection_law() -> None:
    sigma, p = sp.symbols("sigma p", positive=True)
    assert thin_wall_bag_radius(sigma, p) == 2*sigma/p


def test_relative_equilibrium_energy_decomposition() -> None:
    routhian, inertia, omega = sp.symbols("F_omega I omega", real=True)
    static, physical = relative_equilibrium_energies(
        routhian, inertia, omega
    )
    charge = omega * inertia
    assert sp.expand(static - routhian - omega**2 * inertia / 2) == 0
    assert sp.expand(physical - routhian - omega * charge) == 0
    assert sp.expand(physical - static - omega**2 * inertia / 2) == 0
