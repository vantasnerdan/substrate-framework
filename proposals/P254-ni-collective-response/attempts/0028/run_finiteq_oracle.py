#!/usr/bin/env python3
"""Independent small Fourier/projector oracle; never reads material WFCs.

Run from repository root with PYTHONPATH=src. Reuses A27's finite XML/WFC
fixture, not its operator results. Dense PW matrices are ORACLE ONLY.
All reported numbers are generated at execution, not stored expectations.
"""
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace as NS
import json
import sys

import numpy as np
from scipy.linalg import norm
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'0027'))
from bloch_oracle_fixture import fixture
from substrate_framework.qe_bloch_grid import build_bloch_grid, UNITARY_POLICY
from substrate_framework.qe_bloch_edge import BlochEdge
from substrate_framework.qe_bloch_finiteq import BlochFiniteQBridge
from substrate_framework.qe_bloch_finiteq_ionic import BlochFiniteQIonic, finiteq_radial_tables
from substrate_framework.qe_nc_radial import build_radial_tables


class TwoProjectorNC:
    """Defined analytic two-s-channel model, including off-diagonal species D."""
    def __init__(self, omega):
        self.omega_bohr3, self.qmax_bohr_inv = omega, 20.
        self.upf = NS(nh=2, angular_l=np.array([0, 0]),
                      channel_radial=np.array([0, 1]), channel_lm=np.array([0, 0]),
                      channel_l=np.array([0, 0]),
                      expanded_d=np.array([[1.3, .41], [.41, -.7]]))

    def beta(self, radius):
        r2 = np.asarray(radius)**2
        return np.column_stack((np.exp(-.2*r2), .7*np.exp(-.37*r2)))

    def vloc(self, radius):
        return -.7*np.exp(-.17*np.asarray(radius)**2)


def literal_beta(p, tau):
    """Direct analytic model, no production projector/radial/FFT routine."""
    columns = np.empty((len(p), 2), complex)
    for row, momentum in enumerate(p):
        radius2 = sum(float(v)**2 for v in momentum)
        for channel, (scale, decay) in enumerate(((1., .2), (.7, .37))):
            columns[row, channel] = (scale*np.exp(-decay*radius2)/np.sqrt(4*np.pi)
                                     * np.exp(-1j*np.dot(momentum, tau)))
    return columns


def literal_local(ms, mt, fold, density_rows, values):
    coefficients = {tuple(g): v for g, v in zip(density_rows.T, values)}
    out = np.empty((mt.shape[1], ms.shape[1]), complex)
    for t, gt in enumerate(mt.T):
        for s, gs in enumerate(ms.T):
            out[t, s] = coefficients.get(tuple(gt-gs+fold), 0.)
    return out


def literal_ionic(inputs, edge, direction, *, epsilon=None, lost=None,
                  wrong_k=False, diagonal_d=False):
    """Differentiate explicit displaced phases, or return displaced bare vertex.

    epsilon=None is the literal analytic derivative. Non-None gives the
    analytic continuation of the undifferentiated finite-q bare operator;
    this permits a real epsilon central difference for COMPLEX directions.
    """
    ref = inputs.reference
    ms, mt = edge.source_basis.miller, edge.target_basis.miller
    ps = (ref.reciprocal_bohr_inv @ ms+edge.source_basis.k_cart_bohr_inv[:, None]).T
    kt = edge.source_basis.k_cart_bohr_inv if wrong_k else edge.target_basis.k_cart_bohr_inv
    pt = (ref.reciprocal_bohr_inv @ mt+kt[:, None]).T
    gq = (ref.reciprocal_bohr_inv @ inputs.density.miller+edge.q_cart_bohr_inv[:, None]).T
    local = np.zeros(len(gq), complex)
    nonlocal_part = np.zeros((len(pt), len(ps)), complex)
    d = np.array([[1.3, .41], [.41, -.7]])
    if diagonal_d:
        d = np.diag(np.diag(d))
    for atom in range(ref.nat):
        tau, u = ref.positions_bohr[:, atom], direction[atom]
        for gi, g in enumerate(gq):
            base = -.7*np.exp(-.17*np.dot(g, g))*np.exp(-1j*np.dot(g, tau))
            local[gi] += (base*(-1j*np.dot(g, u)) if epsilon is None else
                          base*np.exp(-1j*epsilon*np.dot(g, u)))
        bs, bt = literal_beta(ps, tau), literal_beta(pt, tau)
        for t, ptarget in enumerate(pt):
            for s, psource in enumerate(ps):
                base = sum(bt[t, i]*d[i, j]*bs[s, j].conjugate()
                           for i in range(2) for j in range(2))
                if epsilon is None:
                    factor = 0j
                    if lost != 'source':
                        factor += 1j*np.dot(psource, u)
                    if lost != 'target':
                        factor -= 1j*np.dot(ptarget, u)
                    nonlocal_part[t, s] += base*factor
                else:
                    nonlocal_part[t, s] += base*np.exp(-1j*epsilon*np.dot(ptarget-psource, u))
    return literal_local(ms, mt, edge.reciprocal_fold, inputs.density.miller, local)+nonlocal_part


def relative(actual, expected):
    return float(norm(np.asarray(actual)-np.asarray(expected))/max(1., norm(expected)))


def close(label, actual, expected, records, tol=3e-11):
    error = relative(actual, expected)
    records[label] = error
    if error > tol:
        raise AssertionError(f'{label}: relative error {error}')


def rejected(label, mutant, expected, records):
    difference = relative(mutant, expected)
    records[label] = difference
    if difference < 1e-6:
        raise AssertionError(f'mutation was not exposed: {label}')


def symbolic_oracle():
    e, ps, pt, tau = sp.symbols('epsilon ps pt tau', real=True)
    u, b, d, c = sp.symbols('u b d c', complex=True)
    literal = b*d*c*sp.exp(-sp.I*(pt-ps)*(tau+e*u))
    both = sp.I*ps*u*literal-sp.I*pt*u*literal
    assert sp.simplify(sp.diff(literal, e)-both) == 0
    g, q, v = sp.symbols('G q v', real=True)
    local = v*sp.exp(-sp.I*(g+q)*(tau+e*u))
    assert sp.simplify(sp.diff(local, e)+sp.I*(g+q)*u*local) == 0
    assert sp.simplify(sp.diff(literal, e)-sp.I*ps*u*literal) != 0
    assert sp.simplify(sp.diff(literal, e)+sp.I*pt*u*literal) != 0
    return 'exact local and both-projector displacement identities; missing-side mutants nonzero'


def case(inputs, bindings, spin, source_index, q_index, records):
    ref = inputs.reference
    grid = build_bloch_grid(ref, bindings, symmetry_policy=UNITARY_POLICY)
    target_index = tuple((np.array(source_index)+q_index) % 4)
    source = grid.load_frame(next(s for s in grid.slots if s.ispin == spin and s.mesh_index == source_index))
    target = grid.load_frame(next(s for s in grid.slots if s.ispin == spin and s.mesh_index == target_index))
    edge = BlochEdge.create(ref, source.basis, target.basis, np.array(q_index)/4)
    bridge = BlochFiniteQBridge(inputs, edge, block_size=2)
    ionic = BlochFiniteQIonic(inputs, edge, block_size=2)
    ms, mt = source.basis.miller, target.basis.miller
    rng = np.random.default_rng(2818+spin)
    random = lambda shape: rng.normal(size=shape)+1j*rng.normal(size=shape)
    right, left, coefficient = random((len(ms.T), 3)), random((len(mt.T), 2)), random((2, 3))
    direction = random((ref.nat, 3))
    weight = 1/64
    # Literal physical coordinates in the production F-order convention.
    x = np.indices(ref.fft_grid).reshape(3, -1, order='F').T/np.array(ref.fft_grid)
    ug = np.exp(2j*np.pi*(x @ inputs.density.miller))
    us, ut = np.exp(2j*np.pi*(x @ ms)), np.exp(2j*np.pi*(x @ mt))
    values = random((inputs.density.miller.shape[1],))
    field = ug @ values
    a = literal_local(ms, mt, edge.reciprocal_fold, inputs.density.miller, values)
    close('number_forward', bridge.local_apply(field, right), a @ right, records)
    close('number_adjoint', bridge.local_adjoint_apply(field, left), a.conj().T @ left, records)
    reverse = BlochFiniteQBridge(inputs, edge.reverse(ref), block_size=1)
    close('reverse_number_adjoint', reverse.local_apply(field.conj(), left), a.conj().T @ left, records)
    close('reverse_fold', reverse.edge.reciprocal_fold, -edge.reciprocal_fold, records)
    dyad = left @ coefficient @ right.conj().T
    literal_density = np.zeros(len(x), complex)
    for t in range(len(mt.T)):
        for s in range(len(ms.T)):
            literal_density += dyad[t, s]*np.exp(2j*np.pi*(x @ (mt[:, t]-ms[:, s]+edge.reciprocal_fold)))
    omega = abs(np.linalg.det(ref.cell_bohr))
    literal_density *= weight/omega
    density = bridge.density_from_dyads(left, right, coefficient, weight=weight)
    close('literal_density', density[:, spin-1], literal_density, records)
    close('empty_other_spin', density[:, 2-spin], np.zeros(len(x)), records)
    close('physical_density_pair', bridge.density_pair(density, field), weight*np.vdot(dyad, a), records)
    ncoeff = ug.conj().T @ literal_density/len(x)
    gq = (ref.reciprocal_bohr_inv @ inputs.density.miller+edge.q_cart_bohr_inv[:, None]).T
    literal_gradient = np.array([ug @ (1j*gq[:, axis]*ncoeff) for axis in range(3)])
    close('density_qplusG_gradient', bridge.density_gradient(density)[:, :, spin-1], literal_gradient, records)
    prepared = ionic.prepare(direction)
    exact = literal_ionic(inputs, edge, direction)
    close('ionic_forward', prepared.apply(right), exact @ right, records)
    close('ionic_adjoint', prepared.adjoint_apply(left), exact.conj().T @ left, records)
    moment_density, gradient = ionic.moments(left, right, coefficient, weight=weight)
    close('moment_density', moment_density, density, records)
    literal_moments = []
    for atom in range(ref.nat):
        for axis in range(3):
            unit = np.zeros((ref.nat, 3)); unit[atom, axis] = 1.
            literal_moments.append(weight*np.trace(literal_ionic(inputs, edge, unit).conj().T @ dyad))
    close('all_positive_gradients', gradient, literal_moments, records)
    close('complex_gradient_pair', np.vdot(direction.ravel(), gradient),
          weight*np.trace(exact.conj().T @ dyad), records)
    # Opposite factor-width branch and identity coefficient are both exercised.
    l2, r2, c2 = random((len(mt.T), 4)), random((len(ms.T), 2)), random((4, 2))
    n2, g2 = ionic.moments(l2, r2, c2, weight=weight)
    close('smaller_right_density', n2, bridge.density_from_dyads(l2, r2, c2, weight=weight), records)
    close('smaller_right_gradient', np.vdot(direction.ravel(), g2),
          weight*np.trace(exact.conj().T @ (l2 @ c2 @ r2.conj().T)), records)
    ni, gi = ionic.moments(left, right[:, :2], weight=weight)
    close('identity_coefficient_density', ni,
          bridge.density_from_dyads(left, right[:, :2], np.eye(2), weight=weight), records)
    close('identity_coefficient_gradient', np.vdot(direction.ravel(), gi),
          weight*np.trace(exact.conj().T @ (left @ right[:, :2].conj().T)), records)
    refinements = []
    for epsilon in (.02, .01, .005):
        finite_difference = (literal_ionic(inputs, edge, direction, epsilon=epsilon)
                             - literal_ionic(inputs, edge, direction, epsilon=-epsilon))/(2*epsilon)
        refinements.append(relative(finite_difference, exact))
    if not all(a > 3*b for a, b in zip(refinements, refinements[1:])):
        raise AssertionError(f'central-difference refinement failed: {refinements}')
    records['central_difference_errors'] = refinements
    for side in ('source', 'target'):
        rejected('lost_'+side+'_derivative', literal_ionic(inputs, edge, direction, lost=side), exact, records)
    rejected('wrong_endpoint_k', literal_ionic(inputs, edge, direction, wrong_k=True), exact, records)
    rejected('diagonalized_species_D', literal_ionic(inputs, edge, direction, diagonal_d=True), exact, records)
    if np.any(edge.reciprocal_fold):
        rejected('wrong_fold_sign', literal_local(ms, mt, -edge.reciprocal_fold,
                                                 inputs.density.miller, values), a, records)
        rejected('density_wrong_fold_sign', literal_density*np.exp(-4j*np.pi*(x @ edge.reciprocal_fold)),
                 literal_density, records)
    rejected('lost_physical_volume', density[:, spin-1]*omega, literal_density, records)
    rejected('orbital_k_density_gradient', literal_gradient+1j*source.basis.k_cart_bohr_inv[:, None]
             * literal_density[None, :], literal_gradient, records)
    records['dimensions'] = [edge.source_dimension, edge.target_dimension]
    records['fold'] = edge.reciprocal_fold.tolist()
    records['source_symmetry_fold'] = source.basis.slot.reciprocal_shift.tolist()
    records['target_symmetry_fold'] = target.basis.slot.reciprocal_shift.tolist()
    return edge, exact


def table_range_oracle(inputs, edge):
    # Numerical radial extension, not a different pseudopotential. The kernel
    # itself is already accepted; this tests only the NEW range/identity route.
    r = np.linspace(.01, 6., 121)
    upf = NS(r=r, rab=np.full(len(r), r[1]-r[0]), beta_mesh=len(r), local_mesh=len(r),
             beta_r=(r*np.exp(-r*r))[:, None], angular_l=np.array([0]), nbeta=1,
             local_ry=-2/np.sqrt(1+r*r), z_valence=1.)
    original = build_radial_tables(upf, abs(np.linalg.det(inputs.reference.cell_bohr)), .1)
    narrowed = NS(reference=inputs.reference, density=inputs.density, tables={'X': original})
    extended = finiteq_radial_tables(narrowed, edge)['X']
    required = np.max(np.linalg.norm(inputs.reference.reciprocal_bohr_inv @ inputs.density.miller
                                     + edge.q_cart_bohr_inv[:, None], axis=0))
    assert extended.upf is original.upf and extended.omega_bohr3 == original.omega_bohr3
    assert extended.qmax_bohr_inv >= required and extended is not original
    assert np.isfinite(extended.vloc(np.array([required]))).all()
    np.testing.assert_array_equal(extended.vloc(np.array([.05])), original.vloc(np.array([.05])))
    return dict(required_qmax=float(required), extended_qmax=extended.qmax_bohr_inv,
                same_upf_identity=True, prior_range_unchanged=True)


def main():
    report = {'scope': 'synthetic independent small-Fourier/projector oracle, NOT material evaluation',
              'symbolic': symbolic_oracle(), 'cases': {}}
    with TemporaryDirectory(prefix='p254-finiteq-') as directory:
        inputs, bindings, _, _ = fixture(directory)
        omega = abs(np.linalg.det(inputs.reference.cell_bohr))
        inputs.tables['X'] = TwoProjectorNC(omega)
        outputs = []
        for ngrid in (9, 13):
            adjusted = NS(reference=replace(inputs.reference, fft_grid=(ngrid,)*3,
                                            fft_smooth=(ngrid,)*3, fft_box=(ngrid,)*3),
                          density=inputs.density, tables=inputs.tables, sources=inputs.sources)
            for spin, source_index, q_index in ((1, (0, 0, 1), (0, 0, 1)),
                                                (2, (1, 0, 0), (1, 1, 0))):
                record = {}
                edge, exact = case(adjusted, bindings, spin, source_index, q_index, record)
                report['cases'][f'grid{ngrid}-spin{spin}'] = record
                if spin == 1:
                    assert edge.source_dimension != edge.target_dimension
                    np.testing.assert_array_equal(edge.reciprocal_fold, [0, 0, -1])
                    outputs.append(exact)
            if ngrid == 9:
                report['radial_range'] = table_range_oracle(adjusted, edge)
                grid = build_bloch_grid(adjusted.reference, bindings, symmetry_policy=UNITARY_POLICY)
                folded = next(s for s in grid.slots if s.ispin == 1 and np.any(s.reciprocal_shift))
                record = {}
                folded_edge, _ = case(adjusted, bindings, 1, folded.mesh_index, (1, 0, 1), record)
                assert np.any(folded_edge.source_basis.slot.reciprocal_shift)
                report['cases']['distinct-source-star-fold'] = record
                for name, corrupted in (
                        ('wrong_edge_fold', replace(folded_edge, reciprocal_fold=folded_edge.reciprocal_fold+1)),
                        ('source_star_fold_confusion', replace(
                            folded_edge, reciprocal_fold=folded_edge.reciprocal_fold
                            + folded_edge.source_basis.slot.reciprocal_shift))):
                    try:
                        BlochFiniteQBridge(adjusted, corrupted)
                    except ValueError:
                        report[name] = 'rejected'
                    else:
                        raise AssertionError(name+' was accepted')
        close('fft_refinement_same_operator', outputs[1], outputs[0], report)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
