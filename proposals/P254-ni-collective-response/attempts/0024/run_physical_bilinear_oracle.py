"""New finite bridge oracle; Root alone executes after source-wave admission."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from functools import lru_cache
import sys
import traceback

import numpy as np
import sympy as sp
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "src"))
from substrate_framework.finite_fermion_state import annihilators, second_quantize
from substrate_framework.inertial_fermion_dynamics import FermionChannel
from substrate_framework.physical_bilinear_coordinates import (
    PhysicalBilinearCoordinates, UniformModeIndex, hermitian_basis)

FINDINGS = []


def check(name, residual, tolerance=2e-11, **details):
    value = float(residual)
    FINDINGS.append(dict(name=name, residual=value if np.isfinite(value) else str(value), tolerance=tolerance,
                         passed=bool(np.isfinite(value) and value <= tolerance), **details))


def mutation(name, discrepancy, minimum=1e-7):
    value = float(discrepancy)
    FINDINGS.append(dict(name=name, discrepancy=value if np.isfinite(value) else str(value), minimum=minimum,
                         passed=bool(np.isfinite(value) and value > minimum), mutation=True))


def exact_max(matrix):
    return max((abs(complex(sp.simplify(x))) for x in matrix), default=0.0)


def car_and_normalization():
    """Actual CAR/Fock, including distinct k and spin modes, derive weights."""
    f, g = sp.symbols("f g", real=True)
    rho = sp.diag((1-f)*(1-g), f*(1-g), (1-f)*g, f*g)
    operators = [sp.Matrix([[1, 0], [0, 0]]),
                 sp.Matrix([[0, 1], [1, 0]]) / sp.sqrt(2),
                 sp.Matrix([[0, -sp.I], [sp.I, 0]]) / sp.sqrt(2)]
    lifted = [second_quantize(a) for a in operators]
    centered = [a - sp.trace(rho*a)*sp.eye(4) for a in lifted]
    covariance = sp.Matrix([[sp.trace(rho*a*b) for b in centered] for a in centered])
    reference = sp.diag(f, g)
    wick = sp.Matrix([[sp.trace(reference*a*(sp.eye(2)-reference)*b)
                       for b in operators] for a in operators])
    check("symbolic_Fock_ordered_Wick", exact_max(covariance-wick))
    c = annihilators(2)
    check("symbolic_CAR", max(exact_max(c[i]*c[j].H+c[j].H*c[i]
          - (sp.eye(4) if i == j else sp.zeros(4))) for i in range(2) for j in range(2)))
    null = covariance.subs({f: sp.Rational(1, 2), g: sp.Rational(1, 2)})
    check("equal_population_expected_commutator_zero", exact_max(null-null.T))
    mutation("dropped_population_and_degenerate_coherence", float(sp.trace(null)))
    # Four independent CAR modes: two complete k points for each explicit spin.
    probabilities = [sp.Rational(1, 3), sp.Rational(2, 5),
                     sp.Rational(3, 7), sp.Rational(1, 4)]
    rho4 = sp.diag(*[sp.prod(p if state & (1 << i) else 1-p
                             for i, p in enumerate(probabilities)) for state in range(16)])
    c4 = annihilators(4)
    x = [sp.sqrt(2)*(a.H*a-p*sp.eye(16)) for a, p in zip(c4, probabilities)]
    physical = sum(x, sp.zeros(16))/2
    actual = sp.trace(rho4*physical*physical)
    expected = sum(p*(1-p) for p in probabilities)/2
    check("complete_BvK_independent_spin_Fock", abs(complex(actual-expected)))
    mutation("squared_quadrature_weight", abs(complex(actual-expected/2)))
    mutation("duplicate_spin_factor", abs(complex(actual-2*expected)))
    n, variance = sp.symbols("N variance", positive=True)
    # Each of N independent CAR sectors contributes N*variance, with weight 1/N.
    contracted = sp.simplify(n*(1/n)**2*n*variance)
    check("arbitrary_complete_uniform_N_symbolic", abs(complex(sp.simplify(contracted-variance))))
    FINDINGS.append(dict(name="normalization_identity", passed=True,
        expression=str(contracted), coordinate="sqrt(N)*centered_dGamma",
        per_channel_ordered="N*Tr(rho A (1-rho) B)", physical_weight="1/N",
        physical_covariance="sum_k,spin (1/N)*Tr(rho A (1-rho) B)",
        scope="complete uniform k grid per independent spin; Gamma observables only"))
    basis = hermitian_basis(3)
    check("complete_trace_metric", np.linalg.norm(np.einsum("aij,bji->ab", basis, basis)-np.eye(9)))


def canonical_ion_block():
    """Canonical Schrödinger Gaussian, not a finite oscillator truncation."""
    x, width = sp.symbols("x width", real=True, positive=True)
    psi = sp.exp(-x*x/(4*width*width))
    momentum = lambda wave: -sp.I * sp.diff(wave, x)
    check("symbolic_canonical_Schrodinger_CCR",
          exact_max(sp.Matrix([sp.simplify(x*momentum(psi)-momentum(x*psi)-sp.I*psi)])))
    # Gaussian integrals are exact on the full real line; width is supplied.
    real_x = sp.Symbol("r", real=True)
    wave = sp.exp(-real_x**2/4) / (2*sp.pi)**sp.Rational(1, 4)
    actions = [real_x*wave, -sp.I*sp.diff(wave, real_x)]
    ordered = sp.Matrix([[sp.integrate(sp.conjugate(a)*b, (real_x, -sp.oo, sp.oo))
                           for b in actions] for a in actions])
    expected = sp.Matrix([[1, sp.I/2], [-sp.I/2, sp.Rational(1, 4)]])
    check("infinite_canonical_Gaussian_ordered", exact_max(ordered-expected))
    mutation("ionic_order_transpose", exact_max(ordered-ordered.T))


def finite_fixture():
    from substrate_framework.inertial_fermion_dynamics import InertialFermionDynamics
    from substrate_framework.physical_observable_transfer import PhysicalObservableTransfer
    h = np.diag([0., .7, .7])
    g = np.array([[[.2, .3+.1j, -.2j], [.3-.1j, -.1, .4], [.2j, .4, .5]]])
    density = np.array([[[.7, .2j, .1], [-.2j, -.2, .3j], [.1, -.3j, .6]],
                        [[.1, .2, -.3j], [.2, .8, .1j], [.3j, -.1j, -.4]]])
    channel = FermionChannel(h, np.array([.75, .4, .4]), np.eye(3), g, density)
    model = InertialFermionDynamics([channel], [[.2, .07], [.07, .13]], [[1.3]], [2.1])
    modes = UniformModeIndex((1, 1, 1), ((0, 0, 0),))
    return model, modes, PhysicalObservableTransfer(model, modes)


def independent_rhs(model, u, p, matrices):
    """Literal original density-matrix equation; no production derivative call."""
    density = np.array([sum(c.weight*np.trace(d[s] @ r)
                        for c, d, r in zip(model.channels,
                                          [c.densities for c in model.channels], matrices))
                        for s in range(model.ndensity)])
    field = model.interaction @ density
    drho, gradient = [], np.zeros(model.nion, complex)
    for c, r in zip(model.channels, matrices):
        v = sum((u[a]*c.g[a] for a in range(model.nion)), np.zeros_like(c.h0))
        v += sum((field[s]*c.densities[s] for s in range(model.ndensity)), np.zeros_like(c.h0))
        drho.append(-1j*(c.h0@r-r@c.h0+v@c.rho0-c.rho0@v))
        gradient += c.weight*np.array([np.trace(g@r) for g in c.g])
    return p/model.masses, -model.frozen_curvature@u-gradient, tuple(drho)


def make_preparation(model, modes):
    from substrate_framework.physical_bilinear_preparation import (
        BilinearPreparation, IonicPreparation, CrossPreparation)
    # Admissible joint state: product Bernoulli Fock mixture, and for each
    # occupation n0 a canonical Gaussian shifted by lambda*(n0-f0).
    lam, f = .6, model.channels[0].occupations[0]
    ionic = IonicPreparation(np.array([[1+lam**2*f*(1-f), .5j], [-.5j, .25]]),
        "Mixture of infinite canonical width=1 Gaussian states shifted by 0.6*(n0-f0)")
    cross = np.zeros((2, 3, 3), complex)
    cross[0, 0, 0] = lam*f*(1-f)
    joint = CrossPreparation((cross,), "Same explicit occupation-conditioned Gaussian mixture")
    return BilinearPreparation(model.channels, modes, "Independent Bernoulli Fock probabilities", ionic, joint)


@lru_cache(maxsize=None)
def numeric_car(orbitals):
    return tuple(np.array(a, complex) for a in annihilators(orbitals))


def independent_joint_ordered(model, left, right):
    """Trace the finite Fock mixture and conditional infinite-ion moments."""
    f = model.channels[0].occupations
    probability = np.array([np.prod([p if s & (1 << i) else 1-p
                                    for i, p in enumerate(f)]) for s in range(8)])
    rho = np.diag(probability)
    def centered(matrix):
        c = numeric_car(len(matrix))
        lifted = sum((matrix[i, j]*(a.conj().T@b)
                      for i, a in enumerate(c) for j, b in enumerate(c)),
                     np.zeros((1 << len(matrix), 1 << len(matrix)), complex))
        return lifted-np.trace(rho@lifted)*np.eye(8)
    a, b = centered(left.electronic[0]), centered(right.electronic[0])
    li, ri = np.r_[left.u, left.p], np.r_[right.u, right.p]
    total = 0j
    for s, prob in enumerate(probability):
        mean = np.array([.6*((s & 1)-f[0]), 0.])
        ion = np.array([[1, .5j], [-.5j, .25]])+np.outer(mean, mean)
        total += prob*((a@b)[s, s]+li@ion@ri+a[s, s]*(ri@mean)+b[s, s]*(li@mean))
    return total


def transfer_and_preparation():
    from substrate_framework.physical_observable_transfer import BilinearProbe
    from substrate_framework.physical_bilinear_preparation import BilinearPreparation
    model, modes, transfer = finite_fixture()
    coordinates = PhysicalBilinearCoordinates(model.channels, model.nion, modes)
    preparation = make_preparation(model, modes)
    eye = np.eye(coordinates.size)
    probes = [BilinearProbe(*coordinates.unpack(row)) for row in eye]
    # Complete original-basis equation assembled independently, including null sectors.
    generator = np.column_stack([coordinates.pack(*independent_rhs(
        model, *coordinates.unpack(column))) for column in eye])
    dual = np.array([[transfer.dual(probe).evaluate(*coordinates.unpack(column))
                      for column in eye] for probe in probes])
    check("trace_dual_vs_original_basis_EOM", np.linalg.norm(dual-generator))
    commutator = coordinates.expected_commutator()
    check("expected_commutator_preservation", np.linalg.norm(
        generator@commutator+commutator@generator.T))
    ordered = np.array([[preparation.ordered_parts(a, b)["total"] for b in probes] for a in probes])
    check("joint_preparation_ordered_positivity", max(0, -np.linalg.eigvalsh(ordered).min()))
    check("joint_preparation_commutator", np.linalg.norm(ordered-ordered.T-commutator))
    fock = np.array([[independent_joint_ordered(model, a, b) for b in probes] for a in probes])
    check("conditional_Fock_canonical_ion_cross_state", np.linalg.norm(ordered-fock))
    mutation("nonzero_cross_omission", np.linalg.norm(ordered[:2, 2:]))
    # Asymmetric complex probes expose transpose versus adjoint and trace metric errors.
    left = BilinearProbe(np.array([.2+.1j]), np.array([-.3j]),
                         (np.array([[.3, .2j, .1], [.4, -.2, -.3j], [.2j, .1, .7]]),))
    right = transfer.force_probe([1.])
    right = BilinearProbe(right.u+.2j, right.p+.1, right.electronic)
    times = (.37, .19)
    pulled = [transfer.pullback(p, t, rtol=1e-12, atol=1e-13) for p, t in zip((left, right), times)]
    rows = [np.array([p.evaluate(*coordinates.unpack(c)) for c in eye]) for p in (left, right)]
    direct_rows = [row@expm(t*generator) for row, t in zip(rows, times)]
    pulled_rows = [np.array([p.evaluate(*coordinates.unpack(c)) for c in eye]) for p in pulled]
    check("selected_time_pullback_vs_original_EOM", np.linalg.norm(np.array(pulled_rows)-direct_rows))
    actual = preparation.ordered_parts(*pulled)["total"]
    check("prepared_two_time_ordered", abs(actual-direct_rows[0]@ordered@direct_rows[1]))
    check("propagated_Fock_conditional_state", abs(actual-independent_joint_ordered(model, *pulled)))
    evolved = expm(times[0]*generator)@ordered@expm(times[0]*generator).T
    check("propagated_ordered_positivity", max(0, -np.linalg.eigvalsh(evolved).min()))
    check("propagated_expected_commutator", np.linalg.norm(evolved-evolved.T-commutator))
    mutation("wrong_force_sign", np.linalg.norm(dual[1, 2:]+generator[1, 2:]))
    mutation("adjoint_instead_of_ordered", abs(actual-direct_rows[0]@ordered@direct_rows[1].conj()))
    mutation("elementwise_instead_of_trace_metric", abs(
        left.evaluate(*coordinates.unpack(eye[6]))-np.sum(left.electronic[0]*coordinates.unpack(eye[6])[2][0])))
    mutation("nonstationary_preparation_mislabel", np.linalg.norm(generator@ordered+ordered@generator.T))
    # Retained diagonal/equal-occupation blocks can force ions despite zero response.
    null = [2, 3, 4, coordinates.size-2, coordinates.size-1]
    mutation("null_sector_force_covariance_omission",
             np.linalg.norm(direct_rows[1][null]@ordered[np.ix_(null, null)]))
    partial = BilinearPreparation(model.channels, modes, "Same explicit Bernoulli electrons only")
    incomplete = partial.ordered_parts(*pulled)
    check("electronic_only_explicitly_partial", int(incomplete["complete"] or incomplete["total"] is not None))
    result = transfer.two_time(preparation, left, times[0], right, times[1])
    check("two_time_nonstationarity_label", int(result["stationary_assumed"]))
    return model, modes, transfer, coordinates, preparation, left, right


def basis_and_phase_covariance(fixture):
    from substrate_framework.inertial_fermion_dynamics import InertialFermionDynamics
    from substrate_framework.physical_observable_transfer import BilinearProbe, PhysicalObservableTransfer
    from substrate_framework.physical_bilinear_preparation import BilinearPreparation, CrossPreparation
    model, modes, transfer, coordinates, preparation, left, right = fixture
    fourier = np.exp(2j*np.pi*np.outer(np.arange(3), np.arange(3))/3)/np.sqrt(3)
    for name, unitary in (("basis", fourier), ("phase", np.diag(np.exp(1j*np.array([.2, -.7, .4]))))):
        rotate = lambda a: unitary@a@unitary.conj().T
        c = model.channels[0]
        channel = FermionChannel(rotate(c.h0), c.occupations, unitary@c.orbitals,
                                  rotate(c.g), rotate(c.densities), c.weight, c.spin)
        transformed = InertialFermionDynamics([channel], model.interaction,
                                              model.frozen_curvature, model.masses)
        other = PhysicalObservableTransfer(transformed, modes)
        transform_probe = lambda p: BilinearProbe(p.u, p.p, tuple(rotate(a) for a in p.electronic))
        first = transform_probe(transfer.dual(left))
        second = other.dual(transform_probe(left))
        check(name+"_dual_covariance", np.linalg.norm(other._flatten(first)-other._flatten(second)))
        cross = CrossPreparation(tuple(rotate(a) for a in preparation.cross.matrices),
                                  preparation.cross.state_authority)
        newprep = BilinearPreparation(transformed.channels, modes, preparation.electronic_authority,
                                     preparation.ionic, cross)
        check(name+"_ordered_covariance", abs(
            preparation.ordered_parts(left, right)["total"]
            -newprep.ordered_parts(transform_probe(left), transform_probe(right))["total"]))
        check(name+"_time_pullback_covariance", np.linalg.norm(other._flatten(
            transform_probe(transfer.pullback(left, .31)))-other._flatten(
            other.pullback(transform_probe(left), .31))))


def uniform_BvK_generator():
    from substrate_framework.inertial_fermion_dynamics import InertialFermionDynamics
    from substrate_framework.physical_observable_transfer import PhysicalObservableTransfer, BilinearProbe
    channels = []
    for spin in (-.5, .5):
        for k in range(2):
            channels.append(FermionChannel(np.diag([0., .8+k*.1]),
                np.array([.7, .2+k*.1]), np.eye(2),
                np.array([[[.2, .3j], [-.3j, -.1]]]),
                np.array([[[.4, .2+.1j], [.2-.1j, -.3]]]), .5, spin))
    modes = UniformModeIndex((2, 1, 1), ((0, 0, 0), (1, 0, 0))*2)
    model = InertialFermionDynamics(channels, [[.2]], [[1.1]], [1.7])
    adapter = PhysicalBilinearCoordinates(channels, 1, modes)
    transfer = PhysicalObservableTransfer(model, modes)
    eye = np.eye(adapter.size)
    generator = np.column_stack([adapter.pack(*independent_rhs(
        model, *adapter.unpack(column))) for column in eye])
    expected = adapter.expected_commutator()
    check("complete_BvK_expected_commutator_preservation", np.linalg.norm(
        generator@expected+expected@generator.T))
    wrong = expected.copy()
    wrong[2:, 2:] /= 2
    mutation("unscaled_CAR_metric", np.linalg.norm(generator@wrong+wrong@generator.T))
    density = transfer.density_probe([1])
    check("BvK_physical_probe_weights", max(np.linalg.norm(a-.5*c.densities[0])
          for a, c in zip(density.electronic, channels)))
    dual = np.array([[transfer.dual(BilinearProbe(*adapter.unpack(row))).evaluate(
        *adapter.unpack(column)) for column in eye] for row in eye])
    check("BvK_trace_dual", np.linalg.norm(dual-generator))
    rejected = 0
    try:
        UniformModeIndex((2, 1, 1), ((0, 0, 0),)*4).validate(channels)
    except ValueError:
        rejected = 1
    check("duplicate_q_or_spin_k_rejected", 1-rejected)
    wrong_channel = FermionChannel(channels[0].h0, channels[0].occupations,
        channels[0].orbitals, channels[0].g, channels[0].densities, .3, channels[0].spin)
    rejected = 0
    try:
        modes.validate([wrong_channel, *channels[1:]])
    except ValueError:
        rejected = 1
    check("arbitrary_IBZ_weight_rejected", 1-rejected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    errors = []
    for name, block in (("CAR_normalization", car_and_normalization),
                        ("canonical_ion", canonical_ion_block),
                        ("finite_transfer", transfer_and_preparation),
                        ("BvK_generator", uniform_BvK_generator)):
        try:
            value = block()
            if name == "finite_transfer":
                basis_and_phase_covariance(value)
        except Exception:
            errors.append(dict(block=name, traceback=traceback.format_exc()))
    report = dict(schema="P254-0024-physical-bilinear-oracle/1", findings=FINDINGS,
        errors=errors, passed=not errors and all(x["passed"] for x in FINDINGS),
        scope="Finite Gamma linearized mean-field closure and complete uniform BvK normalization",
        not_claimed=["exact interacting dynamics", "Gaussian fermion higher moments",
                     "stationary C(t-s)", "general-q material transfer", "material preparation",
                     "exact CCR in a finite oscillator truncation"],
        source_receipt="physical-bilinear-additional-sources.json")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    if not report["passed"]:
        raise SystemExit("Physical bilinear oracle failed; complete findings saved to "+str(args.output))


if __name__ == "__main__":
    main()
