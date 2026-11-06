#!/usr/bin/env python3
"""New finite Gamma dyad oracle; execute only under Root's registration.

Expected results use literal Fourier sums and independently placed analytic
projectors in the COMPLETE small original PW basis. SymPy differentiates the
placement exponential before any production action is compared. This is a
synthetic reciprocal-table algebra fixture, not Ni/H material calibration or a
replay of the accepted radial/Native oracles. No material files are decoded.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback

import numpy as np
import sympy as sp

from substrate_framework.qe_material_dyads import GammaDyadCoupling
from substrate_framework.qe_material_grid import PhysicalRhoFFT
from substrate_framework.qe_material_hamiltonian import PWFactors, SpeciesProjectors
from substrate_framework.qe_material_ionic_action import GammaIonicAction
from substrate_framework.qe_nc_projectors import centered_projectors, place_projectors
from substrate_framework.qe_nc_radial import DQ, NCRadialTables
from substrate_framework.qe_nc_upf import NCUPF, expand_channels


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RTOL, ATOL = 1e-10, 1e-12


def binding(path):
    path = Path(path).resolve()
    data = path.read_bytes()
    return dict(path=str(path), bytes=str(len(data)), sha256=hashlib.sha256(data).hexdigest())


def encoded(value):
    value = np.asarray(value)
    return dict(real=value.real.tolist(), imag=value.imag.tolist(), shape=list(value.shape))


class Report:
    def __init__(self):
        self.checks = []
        self.failures = []

    def compare(self, name, actual, expected, *, mutation=False):
        actual, expected = np.asarray(actual), np.asarray(expected)
        same_shape = actual.shape == expected.shape
        finite = bool(np.isfinite(actual).all() and np.isfinite(expected).all())
        close = same_shape and finite and bool(np.allclose(actual, expected, rtol=RTOL, atol=ATOL))
        passed = (same_shape and finite and not close) if mutation else close
        error = float(np.max(np.abs(actual - expected), initial=0)) if same_shape and finite else None
        self.checks.append(dict(name=name, kind="mutation_exposed" if mutation else "positive_comparison",
                                passed=passed, max_abs_error=error, rtol=RTOL, atol=ATOL,
                                actual=encoded(actual), expected=encoded(expected)))
        if not passed:
            self.failures.append(name)
            raise AssertionError(name + ": failed prerequisite; dependent comparisons stopped")
        return passed

    def require(self, name, condition):
        passed = bool(condition)
        self.checks.append(dict(name=name, kind="assertion", passed=passed))
        if not passed:
            self.failures.append(name)
            raise AssertionError(name + ": failed prerequisite; dependent comparisons stopped")

    def rejects(self, name, action):
        try:
            action()
        except ValueError:
            self.require(name, True)
        else:
            self.require(name, False)


def analytic_radial(symbol, q, volume):
    """Independent explicit polynomials, reproduced by cubic table interpolation."""
    scale = 4*np.pi / np.sqrt(volume)
    if symbol == "Ni":
        return scale*np.column_stack((.31 + .07*q*q, -.19 + .04*q*q, q*(.13 + .02*q*q)))
    return scale*(.23 + .03*q*q)[:, None]


def analytic_local(symbol, q, volume):
    a, b, z, g0 = (.17, .021, .7, -.11) if symbol == "Ni" else (-.13, .017, .4, .09)
    out = (4*np.pi/volume)*(a + b*q*q)
    nonzero = q != 0
    out[nonzero] -= 8*np.pi*z/volume*np.exp(-q[nonzero]**2/4)/q[nonzero]**2
    out[~nonzero] = 4*np.pi*g0/volume
    return out


def reciprocal_tables(symbol, volume, qmax):
    """Real existing NCRadialTables with deliberately prescribed finite coefficients.

    The NCUPF dataclass here describes synthetic s/s/p or s channel metadata;
    no parsed UPF identity or quadrature accuracy is asserted by this fixture.
    """
    ell = np.array([0, 0, 1] if symbol == "Ni" else [0], dtype=int)
    dij = np.array([[.83, -.21, 0], [-.21, -.37, 0], [0, 0, .49]]) if symbol == "Ni" else np.array([[-.61]])
    channel_radial, channel_l, channel_lm, expanded_d = expand_channels(ell, dij)
    z = .7 if symbol == "Ni" else .4
    r = np.array([0., .5, 1.])
    metadata = f"0026 synthetic reciprocal-table fixture {symbol}".encode()
    upf = NCUPF(path=f"synthetic:{symbol}", original_bytes=metadata,
                sha256=hashlib.sha256(metadata).hexdigest(), element=symbol, z_valence=z,
                r=r, rab=np.full(3, .5), local_ry=np.zeros(3),
                beta_r=np.zeros((3, len(ell))), angular_l=ell,
                cutoff_indices=np.full(len(ell), 3), raw_dij=dij, local_mesh=3, beta_mesh=3,
                channel_radial=channel_radial, channel_l=channel_l, channel_lm=channel_lm,
                expanded_d=expanded_d)
    q = np.arange(int(qmax/DQ) + 5)*DQ
    a, b, g0 = (.17, .021, -.11) if symbol == "Ni" else (-.13, .017, .09)
    return NCRadialTables(upf, volume, qmax, analytic_radial(symbol, q, volume),
                         4*np.pi/volume*(a + b*q*q), 4*np.pi*g0/volume)


def literal_centered(symbol, momenta, volume):
    """Explicit s/s/p channels, no production harmonics or interpolation."""
    p = np.asarray(momenta)
    q = np.sqrt(np.sum(p*p, axis=1))
    if symbol == "H":
        return (4*np.pi/np.sqrt(volume))*(.23 + .03*q*q)[:, None]/np.sqrt(4*np.pi)
    s1 = (4*np.pi/np.sqrt(volume))*(.31 + .07*q*q)/np.sqrt(4*np.pi)
    s2 = (4*np.pi/np.sqrt(volume))*(-.19 + .04*q*q)/np.sqrt(4*np.pi)
    # q from the p radial factor cancels the unit-vector denominator, also at G=0.
    pref = -1j*(4*np.pi/np.sqrt(volume))*np.sqrt(3/(4*np.pi))*(.13 + .02*q*q)
    return np.column_stack((s1, s2, pref*p[:, 2], -pref*p[:, 0], -pref*p[:, 1]))


def fixture():
    shape = (7, 6, 5)
    miller = np.array([[0, 1, -1, 0, 0, 1, -1, 1],
                       [0, 0, 1, -1, 0, 1, -1, -1],
                       [0, 1, 0, 1, -1, -1, 1, 1]], dtype=int)
    # All literal PW differences, in a nontrivial order, are the local support.
    differences = {tuple(a-b) for a in miller.T for b in miller.T}
    local_miller = np.array(sorted(differences, reverse=True), dtype=int).T
    reciprocal = np.array([[.83, .11, -.07], [0, .97, .13], [0, 0, 1.21]])
    volume = (2*np.pi)**3 / np.linalg.det(reciprocal)
    positions = np.array([[.27, 1.07, -.43], [-.31, .53, .91], [.62, -.77, .19]])
    species = ("Ni", "H", "Ni")  # Species groups are NOT atom-major contiguous.
    p, g = (reciprocal @ miller).T, (reciprocal @ local_miller).T
    qmax = max(np.linalg.norm(p, axis=1).max(), np.linalg.norm(g, axis=1).max()) + .1
    tables = {name: reciprocal_tables(name, volume, qmax) for name in ("Ni", "H")}
    fft = PhysicalRhoFFT.create(shape, miller, reciprocal, volume)
    ionic = GammaIonicAction(shape, miller, local_miller, reciprocal, positions, species, tables)
    projectors = []
    for name in ("Ni", "H"):
        atoms = tuple(i+1 for i, symbol in enumerate(species) if symbol == name)
        centered = centered_projectors(tables[name], p)
        beta = np.asfortranarray(np.column_stack([
            place_projectors(centered, p, positions[:, atom-1]) for atom in atoms]))
        projectors.append(SpeciesProjectors(name, atoms, beta, tables[name].upf.expanded_d))
    factors = PWFactors(fft, np.sum(p*p, axis=1), tuple(projectors), np.zeros(3), 1,
                        ({"fixture": "explicit synthetic original PW coefficients; not material"},))
    coupling = GammaDyadCoupling(factors, ionic, block_size=2)
    return coupling, local_miller, p, g, positions, species


def independent_operators(report, coupling, local_miller, p, g, positions, species):
    """Differentiate the literal placed Hamiltonian, not the action implementation."""
    shape, miller = coupling.factors.fft.physical, coupling.factors.fft.miller
    coordinates = np.indices(shape).reshape(3, -1, order="F") / np.array(shape)[:, None]
    waves = np.exp(2j*np.pi*(coordinates.T @ miller))
    local_waves = np.exp(2j*np.pi*(coordinates.T @ local_miller))
    t = sp.symbols("t0:3", real=True)
    k = sp.symbols("k0:3", real=True)
    phase = sp.exp(-sp.I*sum(k[a]*t[a] for a in range(3)))
    differentiated = [sp.diff(phase, t[a]) for a in range(3)]
    dphase = [sp.lambdify((*k, *t), expression, "numpy") for expression in differentiated]
    report.require("symbolically_differentiated_placement_phase",
                   all(sp.simplify(differentiated[a] + sp.I*k[a]*phase) == 0 for a in range(3)))
    npw, nr, volume = coupling.npw, coupling.nr, coupling.cell_volume_bohr3
    local, bra, ket, whole = [], [], [], []
    for atom, name in enumerate(species):
        tau = positions[:, atom]
        values = analytic_local(name, np.linalg.norm(g, axis=1), volume)
        b0 = literal_centered(name, p, volume)
        # D is independently declared, including a non-diagonal same-l block.
        if name == "Ni":
            d = np.zeros((5, 5))
            d[:2, :2] = [[.83, -.21], [-.21, -.37]]
            d[2:, 2:] = np.eye(3)*.49
        else:
            d = np.array([[-.61]])
        placed = b0*np.exp(-1j*(p @ tau))[:, None]
        original = np.empty((npw, npw), dtype=complex)
        unplaced = np.empty_like(original)
        for i in range(npw):
            for j in range(npw):
                original[i, j] = sum(placed[i, a]*d[a, b]*placed[j, b].conjugate()
                                     for a in range(len(d)) for b in range(len(d)))
                unplaced[i, j] = sum(b0[i, a]*d[a, b]*b0[j, b].conjugate()
                                     for a in range(len(d)) for b in range(len(d)))
        for axis in range(3):
            local_coefficients = values*dphase[axis](*g.T, *tau)
            field = local_waves @ local_coefficients
            lm = waves.conj().T @ (field[:, None]*waves)/nr
            nm = np.empty_like(original)
            for i in range(npw):
                for j in range(npw):
                    delta = p[i]-p[j]
                    nm[i, j] = unplaced[i, j]*dphase[axis](*delta, *tau)
            bm = 1j*original*p[None, :, axis]
            km = -1j*p[:, axis, None]*original
            report.compare(f"primitive_atom{atom}_axis{axis}_both_differentiated_sides", bm + km, nm)
            local.append(lm)
            bra.append(bm)
            ket.append(km)
            whole.append(lm + nm)
    return waves, np.array(whole), np.array(local), np.array(bra), np.array(ket)


def expected_moments(waves, operators, matrix, volume):
    # Original PW-square X exists only inside this deliberately tiny oracle.
    density = np.array([sum(waves[r, i]*matrix[i, j]*waves[r, j].conjugate()
                            for i in range(matrix.shape[0]) for j in range(matrix.shape[1]))
                        for r in range(waves.shape[0])])/volume
    gradient = np.array([sum(operator[i, j]*matrix[j, i]
                             for i in range(matrix.shape[0]) for j in range(matrix.shape[1]))
                         for operator in operators])
    return density, gradient


def exercise(report):
    coupling, local_miller, p, g, positions, species = fixture()
    waves, operators, local, bra, ket = independent_operators(
        report, coupling, local_miller, p, g, positions, species)
    rng = np.random.default_rng(2250026)
    def complex_array(shape):
        return (rng.normal(size=shape) + 1j*rng.normal(size=shape))*.3
    left, right = complex_array((coupling.npw, 3)), complex_array((coupling.npw, 5))
    left[:, 1] += .7*left[:, 0]  # Intentionally nonorthogonal raw frames.
    right[:, 2] -= .4j*right[:, 0]
    coefficient = complex_array((3, 5))
    volume, measure = coupling.cell_volume_bohr3, coupling.grid_measure_bohr3
    report.require("nonorthogonal_raw_frame", abs(np.vdot(left[:, 0], left[:, 1])) > .01)
    report.require("full_original_basis_exceeds_factor_rank", coupling.npw > left.shape[1])
    report.require("all_atoms_and_coordinates", coupling.ndof == 9 and len(operators) == 9)
    report.compare("physical_volume_measure", measure*coupling.nr, volume)
    primitive_matrix = left @ coefficient @ right.conj().T
    primitive_density, _ = expected_moments(waves, operators, primitive_matrix, volume)
    local_expected = expected_moments(waves, local, primitive_matrix, volume)[1]
    nonlocal_expected = expected_moments(waves, bra + ket, primitive_matrix, volume)[1]
    report.compare("primitive_single_transform_local_gradient",
                   coupling._gradient.local(primitive_density), local_expected)
    report.compare("primitive_shared_overlap_nonlocal_gradient",
                   coupling._gradient.nonlocal_block(left @ coefficient, right), nonlocal_expected)
    # Primitive actions: >16 input columns must be streamed to the existing API.
    vectors = complex_array((coupling.npw, 19))
    field, direction = complex_array((coupling.nr,)), complex_array((coupling.ndof,))
    vmat = waves.conj().T @ (field[:, None]*waves)/coupling.nr
    gmat = sum(direction[a]*operators[a] for a in range(coupling.ndof))
    vaction = coupling.potential_action(field, vectors)
    gaction = coupling.ionic_action(direction, vectors)
    report.compare("literal_potential_action_all_rows_19_columns", vaction, vmat @ vectors)
    report.compare("differentiated_ionic_action_all_rows_19_columns", gaction, gmat @ vectors)
    report.compare("potential_complex_linearity", coupling.potential_action(1j*field, vectors), 1j*vaction)
    report.compare("ionic_complex_direction_linearity", coupling.ionic_action(1j*direction, vectors), 1j*gaction)
    # Per-coordinate comparisons prevent cancellation behind one sampled direction.
    for coordinate in range(coupling.ndof):
        e = np.eye(coupling.ndof)[coordinate]
        report.compare(f"literal_coordinate_{coordinate}_action",
                       coupling.ionic_action(e, vectors[:, :3]), operators[coordinate] @ vectors[:, :3])

    cases = [("rectangular_3_by_5", left, right, coefficient),
             ("rectangular_5_by_3", right, left, coefficient.T),
             ("implicit_identity", left, right[:, :3], None),
             ("complex_nonhermitian_square", left, right[:, :3], complex_array((3, 3))),
             ("nineteen_columns_streamed", complex_array((coupling.npw, 19)),
              complex_array((coupling.npw, 19)), None)]
    moment_cases = []
    for name, l, r, c in cases:
        matrix = l @ r.conj().T if c is None else l @ c @ r.conj().T
        expected_density, expected_gradient = expected_moments(waves, operators, matrix, volume)
        density, gradient = coupling.moments(l, r, c)
        report.compare(name + "_literal_density", density, expected_density)
        report.compare(name + "_all_coordinate_gradient", gradient, expected_gradient)
        report.compare(name + "_physical_integrated_trace", measure*np.sum(density), np.trace(matrix))
        moment_cases.append((name, matrix, density, gradient))
    for name, matrix, density, gradient in moment_cases:
        # A full arbitrary complex field and direction test bilinear transpose.
        field = complex_array((coupling.nr,))
        direction = complex_array((coupling.ndof,))
        vmat = waves.conj().T @ (field[:, None]*waves)/coupling.nr
        gmat = sum(direction[a]*operators[a] for a in range(coupling.ndof))
        report.compare(name + "_bilinear_density_potential_transpose",
                       measure*np.dot(field, density), np.trace(vmat @ matrix))
        report.compare(name + "_bilinear_ionic_transpose", np.dot(direction, gradient),
                       np.trace(gmat @ matrix))
        # The genuine Hilbert adjoint uses X†; it is a distinct pairing.
        report.compare(name + "_density_conjugate_adjoint",
                       measure*np.vdot(density, field), np.trace(matrix.conj().T @ vmat))
        report.compare(name + "_ionic_conjugate_adjoint", np.vdot(gradient, direction),
                       np.trace(matrix.conj().T @ gmat))
        report.compare(name + "_wrong_density_adjoint_for_transpose", measure*np.vdot(density, field),
                       np.trace(vmat @ matrix), mutation=True)
        report.compare(name + "_wrong_ionic_adjoint_for_transpose", np.vdot(gradient, direction),
                       np.trace(gmat @ matrix), mutation=True)

    # A saved-span source has a nonzero action in the explicitly empty complement.
    retained = np.zeros((coupling.npw, 3), dtype=complex)
    retained[:3] = [[1, .2j, -.1], [0, .8, .3], [.1j, 0, .6]]
    outside = gmat @ retained
    report.require("nonzero_empty_complement_transition", np.linalg.norm(outside[3:]) > 1e-5)
    report.compare("full_complement_action", coupling.ionic_action(direction, retained), outside)
    projected = outside.copy()
    projected[3:] = 0
    report.compare("mutation_saved_span_output_projection", projected, outside, mutation=True)
    complement_left = outside.copy()
    complement_left[:3] = 0
    complement_x = complement_left @ retained.conj().T
    expected_density, expected_gradient = expected_moments(waves, operators, complement_x, volume)
    density, gradient = coupling.moments(complement_left, retained)
    report.compare("complement_cross_dyad_density", density, expected_density)
    report.compare("complement_cross_dyad_gradient", gradient, expected_gradient)
    report.compare("mutation_complement_cross_dyad_dropped", np.zeros_like(gradient), expected_gradient, mutation=True)

    matrix = left @ coefficient @ right.conj().T
    density, gradient = coupling.moments(left, right, coefficient)
    expected_density, expected_gradient = expected_moments(waves, operators, matrix, volume)
    report.require("every_coordinate_has_complex_signal", np.all(np.abs(expected_gradient.imag) > 1e-8))
    only_bra = expected_moments(waves, local + bra, matrix, volume)[1]
    only_ket = expected_moments(waves, local + ket, matrix, volume)[1]
    report.compare("mutation_missing_projector_ket_side", only_bra, expected_gradient, mutation=True)
    report.compare("mutation_missing_projector_bra_side", only_ket, expected_gradient, mutation=True)
    report.compare("mutation_native_real_part", gradient.real, expected_gradient, mutation=True)
    report.compare("mutation_native_factor_two", 2*gradient, expected_gradient, mutation=True)
    report.compare("mutation_force_minus_sign", -gradient, expected_gradient, mutation=True)
    report.compare("mutation_density_missing_volume", volume*density, expected_density, mutation=True)
    report.compare("mutation_density_unitary_fft_normalization", density/coupling.nr, expected_density, mutation=True)
    conjugated = coupling.moments(left, right, coefficient.conj())
    report.compare("mutation_coefficient_conjugation", conjugated[1], expected_gradient, mutation=True)
    report.compare("mutation_conjugated_direction", coupling.ionic_action(direction.conj(), vectors),
                   gmat @ vectors, mutation=True)
    dropped = gradient.copy()
    dropped[2::3] = 0
    report.compare("mutation_omitted_z_coordinates", dropped, expected_gradient, mutation=True)
    # Explicit dyad adjoint identity; not a silent Hermitization of X.
    adj_density, adj_gradient = coupling.moments(right, left, coefficient.conj().T)
    report.compare("density_explicit_adjoint", adj_density, density.conj())
    report.compare("gradient_explicit_adjoint", adj_gradient, gradient.conj())
    report.compare("moments_complex_linearity", coupling.moments(left, right, 1j*coefficient)[1], 1j*gradient)
    single_block = GammaDyadCoupling(coupling.factors, coupling.ionic, block_size=1)
    report.compare("block_partition_density", single_block.moments(left, right, coefficient)[0], density)
    report.compare("block_partition_gradient", single_block.moments(left, right, coefficient)[1], gradient)
    report.rejects("reject_rectangular_implicit_identity", lambda: coupling.moments(left, right))
    report.rejects("reject_bad_coefficient_shape", lambda: coupling.moments(left, right, coefficient.T))
    report.rejects("reject_truncated_pw_rows", lambda: coupling.moments(left[:-1], right, coefficient))
    report.rejects("reject_non_gamma_factors", lambda: GammaDyadCoupling(
        replace(coupling.factors, k_bohr_inv=np.array([.1, 0, 0])), coupling.ionic))
    report.rejects("reject_wrong_volume", lambda: GammaDyadCoupling(
        replace(coupling.factors, fft=replace(coupling.factors.fft, omega=volume*2)), coupling.ionic))
    displaced = GammaIonicAction(coupling.ionic.shape, coupling.factors.fft.miller, local_miller,
                                coupling.factors.fft.reciprocal_bohr_inv,
                                positions + np.array([[.07], [0], [0]]), species, coupling.ionic.tables)
    report.rejects("reject_mismatched_geometry", lambda: GammaDyadCoupling(coupling.factors, displaced))
    empty = np.empty((coupling.npw, 0), dtype=complex)
    report.compare("empty_dyad_density", coupling.moments(empty, empty)[0], np.zeros(coupling.nr))
    report.compare("empty_dyad_gradient", coupling.moments(empty, empty)[1], np.zeros(coupling.ndof))
    report.require("empty_potential_action", coupling.potential_action(field, empty).shape == empty.shape)
    report.require("empty_ionic_action", coupling.ionic_action(direction, empty).shape == empty.shape)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", required=True, type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    report = Report()
    exception = None
    try:
        exercise(report)
    except Exception:
        exception = traceback.format_exc()
        report.failures.append("uncaught_oracle_exception")
    paths = {Path(__file__).resolve(), HERE / "implementation-registration.json",
             HERE / "reference-preparation-grounding.json", HERE / "action-dynamics-grounding.json",
             HERE / "additional-source-access-001.json"}
    for name, module in tuple(sys.modules.items()):
        if name.startswith("substrate_framework.") and getattr(module, "__file__", None):
            path = Path(module.__file__).resolve()
            if path.is_file():
                paths.add(path)
    result = dict(schema="P254-0026-material-dyad-oracle-v1", passed=not report.failures,
                  checks=report.checks, failures=report.failures, exception=exception,
                  elapsed_seconds=time.monotonic()-started,
                  source_bindings=[binding(path) for path in sorted(paths)],
                  runtime=dict(python=sys.version, numpy=np.__version__, sympy=sp.__version__),
                  scope="New finite original-basis algebra only; ordinary NC/no NLCC/full-complex Gamma. No material state, Native rerun or propagation claim.",
                  independence="Literal non-FFT Fourier sums; analytic reciprocal coefficients and explicit placed-projector sums; SymPy differentiation of the original placement exponential; all coordinates before composed trace pairings.")
    args.result.parent.mkdir(parents=True, exist_ok=True)
    with args.result.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(result=str(args.result.resolve()), passed=result["passed"],
                          checks=len(report.checks), failures=report.failures)))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
