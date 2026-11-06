#!/usr/bin/env python3
"""Independent small Fourier/spherical operator oracle; never a material run.

Independent angular quadrature differentiates exp(-i r p.n), and literal row
pairs assemble rectangular local/projector operators. Production derivative
algorithms are not used to generate the dense reference. Root executes later.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace as NS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0027'), str(HERE)]

import numpy as np
import sympy as sp
from bloch_oracle_fixture import fixture as geometry_fixture
from run_field_derivative_oracle import fixture as radial_fixture, literal_integral
from run_finiteq_oracle import literal_local
from run_physical_spatial_material import (binding, implementation_bindings,
    authenticate_registration, compare, relative, enforce_runtime)
from substrate_framework.qe_bloch_grid import build_bloch_grid, UNITARY_POLICY
from substrate_framework.qe_bloch_edge import BlochEdge
from substrate_framework.qe_material_hamiltonian import MaterialHamiltonian
from substrate_framework.qe_material_grid import PhysicalRhoFFT
from substrate_framework.material_mixed_field_ionic import UPFPreparedVertices
from substrate_framework.material_joint_field_scalar import JointFieldCoefficients
from substrate_framework.material_physical_operators import (
    build_physical_factors, physical_hamiltonian, physical_coupling, physical_finiteq,
    PhysicalFieldVertices, PhysicalFrozenAction, spatial_coefficient_handoff)


def expose(records, name, mutant, expected, floor=1e-9):
    value = relative(mutant, expected)
    records[name] = value
    if not np.isfinite(value) or value <= floor:
        raise AssertionError('mutation not exposed: '+name+' '+str(value))


def symbolic():
    t, ps, pt, tau, u = sp.symbols('t ps pt tau u', real=True)
    b, c, d = sp.symbols('b c d', commutative=False)
    phase = sp.exp(-sp.I*(pt-ps)*(tau+t*u))
    expected = -sp.I*(pt-ps)*u*b*d*c*phase
    assert sp.simplify(sp.diff(b*d*c*phase, t)-expected) == 0
    assert sp.simplify(sp.diff(b*d*c*phase, t)+sp.I*pt*u*b*d*c*phase) != 0
    h, dv, g, dg, x, a, ell, C, dc, contact, mixed = sp.symbols('h dv g dg x a ell C dc contact mixed')
    old = h+dv+x*(g+dg)+(C+dc)*x*x/2+a*ell+contact*a*a/2+x*mixed*a
    new = (h+dv)+x*(g+dg)+(C+dc)*x*x/2+a*ell+contact*a*a/2+x*mixed*a
    assert sp.expand(old-new) == 0
    assert sp.diff(new, a) == ell+contact*a+x*mixed
    assert sp.diff(new, x) == g+dg+(C+dc)*x+a*mixed
    assert sp.expand(old-(new+dv+x*dg)) != 0
    return 'exact rectangular phase and once-only joint scalar derivative identities'


def model(directory):
    inputs, bindings, _, _ = geometry_fixture(directory)
    # Explicit synthetic geometry, not a symmetry claim or modification of Ni/H.
    ref = replace(inputs.reference, positions_bohr=np.array([[.23, 1.17], [-.31, .83], [.19, -1.09]]))
    table = radial_fixture()
    table.omega_bohr3 = abs(float(np.linalg.det(ref.cell_bohr)))
    table.qmax_bohr_inv = 20.
    table.vloc = lambda r: -.7*np.exp(-.17*np.asarray(r)**2)
    def forbidden_beta(radius):
        raise AssertionError('physical builder allocated the legacy interpolated beta bank')
    table.beta = forbidden_beta
    nh = len(table.upf.channel_radial)
    table.upf.nh = nh
    rng = np.random.default_rng(2542841)
    raw = rng.normal(size=(nh, nh))+1j*rng.normal(size=(nh, nh))
    table.upf.expanded_d = (raw+raw.conj().T)/7+np.diag(np.linspace(.2, 1.1, nh))
    inputs = NS(reference=ref, tables={'X': table}, density=inputs.density, sources=inputs.sources)
    grid = build_bloch_grid(ref, bindings, symmetry_policy=UNITARY_POLICY)
    fft = PhysicalRhoFFT.create(ref.fft_smooth, inputs.density.miller, ref.reciprocal_bohr_inv, table.omega_bohr3)
    x = np.indices(ref.fft_smooth).reshape(3, -1, order='F').T/np.asarray(ref.fft_smooth)
    potential = .13+.17*np.cos(2*np.pi*x[:, 0])-.09*np.sin(2*np.pi*x[:, 1])
    background = NS(inputs=inputs, fft=fft, vrs=np.column_stack((potential, -.3*potential)))
    return inputs, grid, background


def momenta(basis):
    return (basis.reciprocal_bohr_inv @ basis.miller+basis.k_cart_bohr_inv[:, None]).T


def integral(table, p, direction=(0., 0., 0.)):
    return literal_integral(table, p, np.asarray(direction), np.asarray(direction))


def dense_ionic(inputs, edge, u, centers, *, step=None, omit=None, phase_sign=-1):
    """Literal row-pair Fourier derivative or displaced operator continuation."""
    ref, table = inputs.reference, inputs.tables['X']
    ps, pt = momenta(edge.source_basis), momenta(edge.target_basis)
    gq = (ref.reciprocal_bohr_inv @ inputs.density.miller+edge.q_cart_bohr_inv[:, None]).T
    source, target = centers
    local = np.zeros(len(gq), complex)
    nonlocal_part = np.zeros((len(pt), len(ps)), complex)
    for atom in range(ref.nat):
        tau, ua = ref.positions_bohr[:, atom], u[atom]
        for gi, g in enumerate(gq):
            value = -.7*np.exp(-.17*np.dot(g, g))*np.exp(phase_sign*1j*np.dot(g, tau))
            local[gi] += value*(-1j*np.dot(g, ua) if step is None else np.exp(-1j*step*np.dot(g, ua)))
        # Dense ONLY here: no production projector/ionic action used as oracle.
        for ti, ptarget in enumerate(pt):
            for si, psource in enumerate(ps):
                value = sum(target[ti, i]*table.upf.expanded_d[i, j]*source[si, j].conjugate()
                            for i in range(table.upf.nh) for j in range(table.upf.nh))
                value *= np.exp(phase_sign*1j*np.dot(ptarget-psource, tau))
                if step is None:
                    factor = (0 if omit == 'source' else 1j*np.dot(psource, ua))
                    factor -= (0 if omit == 'target' else 1j*np.dot(ptarget, ua))
                else:
                    factor = np.exp(-1j*step*np.dot(ptarget-psource, ua))
                nonlocal_part[ti, si] += value*factor
    return literal_local(edge.source_basis.miller, edge.target_basis.miller,
        edge.reciprocal_fold, inputs.density.miller, local)+nonlocal_part


def diagonal_case(inputs, grid, background, mesh, report, mutants):
    frame = grid.load_frame(next(s for s in grid.slots if s.ispin == 1 and s.mesh_index == mesh))
    basis, table = frame.basis, inputs.tables['X']
    h = physical_hamiltonian(background, basis, 1, block_size=2)
    factors = h.factors
    assert factors.sources[:-1] == inputs.sources
    assert factors.sources[-1]['identity'] == basis.identity
    np.testing.assert_array_equal(factors.fft.miller, basis.miller)
    np.testing.assert_array_equal(factors.k_bohr_inv, basis.k_cart_bohr_inv)
    p = momenta(basis)
    e = np.array([.23, -.19, .17])
    center, first, _, _ = integral(table, p, e)
    x = np.indices(inputs.reference.fft_smooth).reshape(3, -1, order='F').T/np.array(inputs.reference.fft_smooth)
    plane = np.exp(2j*np.pi*x @ basis.miller)
    dense = np.diag(np.sum(p*p, axis=1)).astype(complex)
    dense += plane.conj().T @ (background.vrs[:, 0, None]*plane)/len(x)
    linear = np.diag(2*p @ e).astype(complex)
    one_leg = np.zeros_like(dense)
    for slot, atom in enumerate(factors.species[0].atom_indices):
        phase = np.exp(-1j*p @ inputs.reference.positions_bohr[:, atom-1])[:, None]
        beta, ba = center*phase, first*phase
        sl = slice(slot*table.upf.nh, (slot+1)*table.upf.nh)
        compare(report, 'literal_UPF_beta_atom'+str(atom), factors.species[0].beta[:, sl], beta)
        dense += beta @ table.upf.expanded_d @ beta.conj().T
        left = ba @ table.upf.expanded_d @ beta.conj().T
        right = beta @ table.upf.expanded_d @ ba.conj().T
        linear += left+right
        one_leg += right
    rng = np.random.default_rng(2542828+sum(mesh))
    psi = rng.normal(size=(len(p), 4))+1j*rng.normal(size=(len(p), 4))
    d = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    compare(report, 'literal_dense_H0', h.apply(psi), dense @ psi)
    vertices = PhysicalFieldVertices(h, projector_block_size=3)
    action = vertices.prepare(e)
    compare(report, 'literal_angular_field_derivative', action.apply(psi), linear @ psi)
    compare(report, 'field_adjoint', action.adjoint_apply(psi), linear.conj().T @ psi)
    values = vertices.moments(psi, d)
    expected_field = np.trace(d @ (psi.conj().T @ linear @ psi))
    compare(report, 'complex_original_D_field_trace', e @ values['field_gradient'], expected_field)
    compare(report, 'current_minus_one_half', e @ values['electronic_current_hartree'], -expected_field/2)
    edge = BlochEdge.create(inputs.reference, basis, basis, (0., 0., 0.))
    coupling = physical_coupling(inputs, basis, factors, block_size=2)
    u = np.array([[.19+.07j, -.13, .11j], [-.23j, .17+.03j, .29]])
    exact = dense_ionic(inputs, edge, u, (center, center))
    compare(report, 'Gamma_Fourier_derivative', coupling.ionic_action(u.ravel(), psi), exact @ psi)
    _, gradient = coupling.moments(psi, psi, d)
    compare(report, 'Gamma_complex_D_gradient', u.ravel() @ gradient, np.trace(d @ (psi.conj().T @ exact @ psi)))
    # Independently defined legacy bank; comparison is not used as dense oracle.
    legacy_species = tuple(replace(s, beta=np.asfortranarray(.93*s.beta)) for s in factors.species)
    legacy_factors = replace(factors, species=legacy_species)
    legacy = MaterialHamiltonian(background, legacy_factors, 1, 2)
    correction = UPFPreparedVertices(legacy).evaluate(psi, d, u, e)
    legacy_coupling = physical_coupling(inputs, basis, legacy_factors, block_size=2)
    frozen = PhysicalFrozenAction(h, coupling, direction=u.ravel(), shift=e,
                                  fft_block_size=2, projector_block_size=3)
    expected = (dense+exact+linear) @ psi
    compare(report, 'fused_physical_spatial_field', frozen.apply(psi), expected)
    compare(report, 'fused_adjoint', frozen.adjoint_apply(psi), (dense+exact+linear).conj().T @ psi)
    compare(report, 'A25_plus_UPFPreparedVertices_identity', frozen.apply(psi),
            legacy.apply(psi)+legacy_coupling.ionic_action(u.ravel(), psi)+correction['action'])
    zero = vertices.evaluate(psi, d, u, np.zeros(3))
    compare(report, 'zero_field_no_correction', zero['action'], np.zeros_like(psi))
    compare(report, 'zero_field_no_DeltaV0', zero['delta_energy'], 0.)
    compare(report, 'zero_field_no_Delta_g', zero['delta_gradient'], np.zeros(3*inputs.reference.nat))
    # Field action never acquires state-dependent contact or mixed operators.
    scaled = vertices.evaluate(psi, 1.7j*d, 2*u, e)
    compare(report, 'field_action_independent_of_D_u', scaled['action'], linear @ psi)
    compare(report, 'field_reduction_complex_linearity', scaled['field_gradient'], 1.7j*values['field_gradient'])
    small = vertices.prepare(1e-12*e)
    compare(report, 'no_small_shift_cutoff', small.apply(psi)/1e-12, linear @ psi)
    expose(mutants, 'omit_right_field_leg', (linear-one_leg) @ psi, linear @ psi)
    expose(mutants, 'reverse_current_sign', expected_field/2, -expected_field/2)
    expose(mutants, 'discard_complex_D', np.trace(d.real @ (psi.conj().T @ linear @ psi)), expected_field)
    expose(mutants, 'discard_offdiagonal_D', np.trace(np.diag(np.diag(d)) @ (psi.conj().T @ linear @ psi)), expected_field)
    expose(mutants, 'double_zero_field_correction', expected+h.apply(psi)-legacy.apply(psi), expected)
    expose(mutants, 'reverse_positive_gradient', -u.ravel() @ gradient, np.trace(d @ (psi.conj().T @ exact @ psi)))
    expose(mutants, 'wrong_atom_phase', dense_ionic(inputs, edge, u, (center, center), phase_sign=1), exact)
    # Independent finite field difference from angular integrals, fixed atom phase.
    def shifted_dense(t):
        b0 = integral(table, p+t*e)[0]
        result = np.diag(np.sum((p+t*e)**2, axis=1)).astype(complex)
        for tau in inputs.reference.positions_bohr.T:
            b = b0*np.exp(-1j*p @ tau)[:, None]
            result += b @ table.upf.expanded_d @ b.conj().T
        return result
    errors = []
    for step in (.002, .001, .0005):
        finite = (shifted_dense(step)-shifted_dense(-step))/(2*step)
        errors.append(relative(finite, linear))
    if errors[-1] > 2e-7 or not errors[-1] < errors[0]/8:
        raise AssertionError('independent field finite difference refinement '+repr(errors))
    report['field_finite_difference_errors'] = errors
    return frame, h, psi, d, center


def rectangular_case(inputs, source, target, report, mutants):
    sf, sh, right, d, scenter = source
    tf, th, left, _, tcenter = target
    edge = BlochEdge.create(inputs.reference, sf.basis, tf.basis, (0., 0., .25))
    if edge.source_dimension == edge.target_dimension or not np.any(edge.reciprocal_fold):
        raise AssertionError('oracle must retain unequal rows and nontrivial fold')
    owner = physical_finiteq(inputs, edge, sh.factors, th.factors, block_size=2)
    u = np.array([[.11+.07j, -.13j, .17], [-.19, .23+.03j, -.29j]])
    exact = dense_ionic(inputs, edge, u, (scenter, tcenter))
    prepared = owner.prepare(u)
    compare(report, 'rectangular_forward', prepared.apply(right), exact @ right)
    compare(report, 'rectangular_adjoint', prepared.adjoint_apply(left), exact.conj().T @ left)
    coefficient = d[:3, :]
    left = left[:, :3]
    density, gradient = owner.moments(left, right, coefficient, weight=1/64)
    dyad = left @ coefficient @ right.conj().T
    compare(report, 'rectangular_complex_D_pair', np.vdot(u.ravel(), gradient), np.trace(exact.conj().T @ dyad)/64)
    actual_grad = []
    for atom in range(inputs.reference.nat):
        for axis in range(3):
            unit = np.zeros((inputs.reference.nat, 3)); unit[atom, axis] = 1
            g = dense_ionic(inputs, edge, unit, (scenter, tcenter))
            actual_grad.append(np.trace(g.conj().T @ dyad)/64)
    compare(report, 'rectangular_all_moments', gradient, actual_grad)
    x = np.indices(inputs.reference.fft_smooth).reshape(3, -1, order='F').T/np.array(inputs.reference.fft_smooth)
    direct_density = np.zeros(len(x), complex)
    for ti, mt in enumerate(tf.basis.miller.T):
        for si, ms in enumerate(sf.basis.miller.T):
            direct_density += dyad[ti, si]*np.exp(2j*np.pi*x @ (mt-ms+edge.reciprocal_fold))
    direct_density /= 64*sh.factors.fft.omega
    compare(report, 'rectangular_literal_Fourier_density', density[:, 0], direct_density)
    expose(mutants, 'wrong_rectangular_fold', density[:, 0]*np.exp(-4j*np.pi*x @ edge.reciprocal_fold), direct_density)
    for side in ('source', 'target'):
        expose(mutants, 'omit_'+side+'_ionic_leg', dense_ionic(inputs, edge, u, (scenter, tcenter), omit=side), exact)
    expose(mutants, 'rectangular_wrong_atom_phase', dense_ionic(inputs, edge, u, (scenter, tcenter), phase_sign=1), exact)
    errors = []
    for step in (.002, .001, .0005):
        finite = (dense_ionic(inputs, edge, u, (scenter, tcenter), step=step)
                  -dense_ionic(inputs, edge, u, (scenter, tcenter), step=-step))/(2*step)
        errors.append(relative(finite, exact))
    if errors[-1] > 2e-7 or not errors[-1] < errors[0]/8:
        raise AssertionError('independent displacement refinement '+repr(errors))
    report.update(finite_difference_errors=errors, dimensions=[edge.source_dimension, edge.target_dimension],
                  q=edge.q_fractional.tolist(), fold=edge.reciprocal_fold.tolist())


def scalar_case(report, mutants):
    rng = np.random.default_rng(2542870)
    c = rng.normal(size=(6, 6)); c = c+c.T
    dc = rng.normal(size=(6, 6)); dc = (dc+dc.T)/13
    contact = np.array([[.7, .13, -.09], [.13, .9, .17], [-.09, .17, 1.1]])
    mixed = rng.normal(size=(6, 3))/7
    coefficients = JointFieldCoefficients(dc, contact, mixed, 'synthetic independent scalar oracle')
    q = np.stack((c, 1.2*c, .7*c))
    corrected, residual = spatial_coefficient_handoff(q, coefficients)
    u, a = rng.normal(size=6), np.array([.19, -.17, .23])
    for k in range(3):
        expected = u @ q[k] @ u/2+u @ dc @ u/2+a @ contact @ a/2+u @ mixed @ a
        actual = u @ corrected[k] @ u/2+a @ residual.contact @ a/2+u @ residual.mixed @ a
        compare(report, 'once_only_Cq_'+str(k), actual, expected)
        expose(mutants, 'double_curvature_q'+str(k), actual+u @ dc @ u/2, expected)
    compare(report, 'residual_DeltaC_zero', residual.delta_curvature, np.zeros_like(dc))
    expose(mutants, 'double_contact', actual+a @ contact @ a/2, expected)
    expose(mutants, 'double_mixed', actual+u @ mixed @ a, expected)
    expose(mutants, 'omit_contact', actual-a @ contact @ a/2, expected)
    expose(mutants, 'omit_mixed', actual-u @ mixed @ a, expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    enforce_runtime()
    registration, _ = authenticate_registration()
    report = dict(schema='P254-physical-spatial-oracle-v1', registration=registration,
        implementation_sources=implementation_bindings(), symbolic=symbolic(),
        scope='independent small synthetic construction only; no material or trajectory evidence',
        diagonal={}, rectangular={}, scalar={}, mutations={})
    with TemporaryDirectory(prefix='p254-physical-spatial-oracle-') as directory:
        inputs, grid, background = model(directory)
        endpoints = []
        for mesh in ((0, 0, 0), (0, 0, 1), (0, 0, 2)):
            key = '-'.join(map(str, mesh))
            result, mutants = {}, {}
            item = diagonal_case(inputs, grid, background, mesh, result, mutants)
            report['diagonal'][key] = result
            report['mutations'][key] = mutants
            if mesh != (0, 0, 0):
                endpoints.append(item)
        report['mutations']['rectangular'] = {}
        rectangular_case(inputs, *endpoints, report['rectangular'], report['mutations']['rectangular'])
    report['mutations']['scalar'] = {}
    scalar_case(report['scalar'], report['mutations']['scalar'])
    report['status'] = 'passed'
    text = json.dumps(report, indent=2, allow_nan=False)+'\n'
    if len(text.encode()) > 4*1024**2:
        raise ValueError('oracle report byte cap exceeded')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        stream.write(text)
    print(json.dumps(dict(status='passed', output=binding(args.output))))


if __name__ == '__main__':
    main()
