#!/usr/bin/env python3
"""Independent literal radial/angular and dense trace oracle; synthetic only.

No material, DFT, original filling or Gram computations. Execute explicitly.
The literal integrals differentiate exp(-i r p.n), not production jets.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace as NS
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0025'), str(HERE.parent/'0024')]

import numpy as np
from run_field_derivative_oracle import fixture, literal_integral
from run_full_field_moments_material import implementation_bindings
from material_run_contract import binding
from substrate_framework.qe_material_hamiltonian import MaterialHamiltonian
from substrate_framework.material_physical_field_moments import (
    full_column_moments, transport_moments, aggregate_slots, require_scalar_channels, MOMENT_KEYS)

AXES = np.eye(3)
PAIRS = tuple((i, j) for i in range(3) for j in range(i, 3))


def model(p, positions, tables, active_d):
    beta = literal_integral(tables, p, AXES[0], AXES[1])[0]
    # Smooth rotation-scalar stand-in for a DIFFERENT zero-field radial family.
    legacy = beta*(1+.07*np.sum(p*p, axis=1))[:, None]
    bank = np.concatenate([legacy*np.exp(-1j*p@tau)[:, None] for tau in positions], axis=1)
    fft = NS(reciprocal_bohr_inv=np.eye(3), miller=p.T, physical=(2, 2, 2))
    reference = NS(positions_bohr=np.asarray(positions).T)
    inputs = NS(reference=reference, tables={'fixture': tables}, sources=())
    species = NS(name='fixture', atom_indices=tuple(range(1, len(positions)+1)),
                 beta=bank, active_d=active_d)
    factors = NS(fft=fft, k_bohr_inv=np.zeros(3), species=(species,), sources=({},))
    return MaterialHamiltonian(NS(inputs=inputs, fft=fft), factors, 1), beta, legacy


def dense_reference(p, positions, tables, active_d, psi, density):
    def mean(a):
        return np.trace(density @ (psi.conj().T @ a @ psi))
    gradient = np.array([mean(np.diag(2*p[:, i])) for i in range(3)])
    contact = 2*mean(np.eye(len(p)))*np.eye(3, dtype=complex)
    kinetic_g, kinetic_c = gradient.copy(), contact.copy()
    delta_e = 0j
    delta_g = np.zeros((len(positions), 3), complex)
    delta_c = np.zeros((len(positions), 3, 3), complex)
    mixed_cross, phase_bad_gradient = np.zeros((3, 3), complex), np.zeros(3, complex)
    base = literal_integral(tables, p, AXES[0], AXES[1])[0]
    legacy = base*(1+.07*np.sum(p*p, axis=1))[:, None]
    for atom, tau in enumerate(positions):
        phase = np.exp(-1j*p@tau)[:, None]
        b, old = base*phase, legacy*phase
        delta_e += mean(b@active_d@b.conj().T-old@active_d@old.conj().T)
        for i, j in PAIRS:
            values = literal_integral(tables, p, AXES[i], AXES[j])
            b, bi, bj, bij = [v*phase for v in values]
            l = bi@active_d@b.conj().T+b@active_d@bi.conj().T
            cross = bi@active_d@bj.conj().T+bj@active_d@bi.conj().T
            c = bij@active_d@b.conj().T+cross+b@active_d@bij.conj().T
            if i == j:
                gradient[i] += mean(l)
                wrong_b, wrong_bi = values[0]*phase.conj(), values[1]*phase.conj()
                phase_bad_gradient[i] += mean(wrong_bi@active_d@wrong_b.conj().T
                                              +wrong_b@active_d@wrong_bi.conj().T-l)
            contact[i, j] += mean(c)
            mixed_cross[i, j] += mean(cross)
            contact[j, i] = contact[i, j]
            mixed_cross[j, i] = mixed_cross[i, j]
            for sign, v in ((1, b), (-1, old)):
                vi, vj = -1j*p[:, i, None]*v, -1j*p[:, j, None]*v
                vij = -p[:, i, None]*p[:, j, None]*v
                if i == j:
                    delta_g[atom, i] += sign*mean(vi@active_d@v.conj().T+v@active_d@vi.conj().T)
                dc = mean(vij@active_d@v.conj().T+vi@active_d@vj.conj().T
                          +vj@active_d@vi.conj().T+v@active_d@vij.conj().T)
                delta_c[atom, i, j] += sign*dc
                delta_c[atom, j, i] = delta_c[atom, i, j]
    result = dict(number=mean(np.eye(len(p))), gradient=gradient, current_hartree=-gradient/2,
                  contact=contact, kinetic_gradient=kinetic_g, kinetic_contact=kinetic_c,
                  delta_nonlocal_energy_ry=delta_e, delta_nonlocal_positive_gradient_ry_bohr=delta_g,
                  delta_nonlocal_force_ry_bohr=-delta_g, delta_nonlocal_onsite_ry_bohr2=delta_c)
    return result, mixed_cross, phase_bad_gradient


def run():
    tables = fixture()
    nh = len(tables.upf.channel_l)
    active_d = np.diag(.5+.2*tables.upf.channel_l)
    # Mix two distinct p radial channels with identical coefficients for all m.
    for a, b in ((1, 9), (2, 10), (3, 11)):
        active_d[a, b] = active_d[b, a] = .17
    p = np.array([[0., 0., 0.], [.2, -.3, .4], [1.1, .7, -.6], [-.31, 0., .5], [.4, .9, .2]])
    positions = np.array([[.21, -.18, .33], [-.43, .51, .12]])
    rng = np.random.default_rng(22528004)
    psi = rng.normal(size=(len(p), 3))+1j*rng.normal(size=(len(p), 3))
    # Deliberately independent non-Hermitian D: imaginary trace order cannot hide.
    density = rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3))
    h, beta, legacy = model(p, positions, tables, active_d)
    actual = full_column_moments(h, psi, density, row_block=2)
    expected, cross, wrong_phase = dense_reference(p, positions, tables, active_d, psi, density)
    errors = {k: float(np.max(np.abs(actual[k]-expected[k]))) for k in MOMENT_KEYS}
    assert max(errors.values()) < 2e-10, errors
    assert actual['work']['local_fft_count'] == actual['work']['pw_square_matrices'] == 0
    assert actual['work']['wide_projector_columns'] == 29*nh*len(positions)
    # Raw-column rephasing must transform D independently in opposite legs.
    u = np.diag(np.exp(1j*np.array([.3, -.7, 1.1])))
    changed = full_column_moments(h, psi@u, u.conj().T@density@u, row_block=3)
    assert max(np.max(np.abs(changed[k]-actual[k])) for k in MOMENT_KEYS) < 2e-10
    wrong_conjugation = full_column_moments(h, psi, density.T, row_block=2)
    wrong_rephase = full_column_moments(h, psi@u, density, row_block=2)
    # Independent dense differentiation of the literal field scalar and ionic delta.
    def scalar(a, displacement):
        b = literal_integral(tables, p+a, AXES[0], AXES[1])[0]
        operator = np.diag(np.sum((p+a)**2, axis=1)).astype(complex)
        correction = np.zeros_like(operator)
        for atom, tau in enumerate(positions+displacement):
            phase = np.exp(-1j*p@tau)[:, None]
            v, v0, old = b*phase, beta*phase, legacy*phase
            operator += v@active_d@v.conj().T
            correction += v0@active_d@v0.conj().T-old@active_d@old.conj().T
        return (np.trace(density@psi.conj().T@operator@psi),
                np.trace(density@psi.conj().T@correction@psi))
    refinements = []
    for step in (.04, .02, .01):
        zero = np.zeros_like(positions)
        e, f = AXES[0], AXES[1]
        g = (scalar(step*e, zero)[0]-scalar(-step*e, zero)[0])/(2*step)
        c = sum(si*sj*scalar(step*(si*e+sj*f), zero)[0]
                for si in (-1, 1) for sj in (-1, 1))/(4*step**2)
        du = np.zeros_like(positions); du[0, 0] = step
        dg = (scalar(np.zeros(3), du)[1]-scalar(np.zeros(3), -du)[1])/(2*step)
        dc = (scalar(np.zeros(3), du)[1]+scalar(np.zeros(3), -du)[1]
              -2*scalar(np.zeros(3), zero)[1])/step**2
        refinements.append([step, abs(g-actual['gradient'][0]), abs(c-actual['contact'][0, 1]),
                            abs(dg-actual['delta_nonlocal_positive_gradient_ry_bohr'][0, 0]),
                            abs(dc-actual['delta_nonlocal_onsite_ry_bohr2'][0, 0, 0])])
    assert np.all(np.asarray(refinements)[-1, 1:] < np.asarray(refinements)[0, 1:]/10), refinements
    # Rotation, atom permutation and translation: literal radial operators at
    # transformed geometry, not the production transport formula as reference.
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., -1.]])
    permutation = np.array([1, 0])
    translation = np.array([.13, -.21, .17])
    pt = p@rotation.T
    taut = positions[permutation]@rotation.T+translation
    ht, _, _ = model(pt, taut, tables, active_d)
    target = full_column_moments(ht, np.exp(-1j*pt@translation)[:, None]*psi,
                                density, row_block=2)
    transported = transport_moments(actual, rotation, permutation)
    covariance_errors = {k: float(np.max(np.abs(target[k]-transported[k]))) for k in MOMENT_KEYS}
    assert max(covariance_errors.values()) < 2e-10, covariance_errors
    bad_d = active_d.copy(); bad_d[1, 1] += .2
    try:
        require_scalar_channels(tables.upf, bad_d)
    except ValueError:
        pass
    else:
        raise AssertionError('anisotropic channel mutation was admitted')
    # Full weight algebra through the actual aggregator with two independent
    # spin sources, each carrying two different unitary slots. No spin flip.
    slots = []
    for spin in (1, 2):
        for index, (r, weight) in enumerate(((np.eye(3, dtype=int), .25), (rotation.astype(int), .75))):
            op = NS(effective_time_reversal=False, rotation=r, fractional_translation=np.zeros(3),
                    equivalent_atoms=np.array([1, 2]))
            slots.append(NS(ispin=spin, source_ik=1, mesh_index=(index, 0, 0), weight=weight,
                            operation=op, atom_image_offsets=np.zeros((3, 2))))
    reference = NS(fft_smooth=(2, 2, 2), nat=2, atom_names=('X', 'X'), cell_bohr=4*np.eye(3),
                   positions_bohr=positions.T,
                   select=lambda ik, spin: (NS(weight=.3 if spin == 1 else .7), None, None))
    def require_slot(slot):
        assert any(slot is item for item in slots)
    grid = NS(slots=slots, reference=reference, require_slot=require_slot)
    second = {k: 1.7*np.asarray(actual[k]) for k in MOMENT_KEYS}
    second['normalization'] = actual['normalization']
    sources = {(1, 1): actual, (2, 1): second}
    summed = aggregate_slots(grid, sources, geometry_policy='declared-geometric')
    # Independently evaluated target, without the atom permutation used above.
    hrot, _, _ = model(pt, positions@rotation.T, tables, active_d)
    rotated = full_column_moments(hrot, psi, density, row_block=2)
    manual = {k: 2.7*(.25*np.asarray(actual[k])+.75*np.asarray(rotated[k])) for k in MOMENT_KEYS}
    weight_errors = {k: float(np.max(np.abs(summed['moments'][k]-manual[k]))) for k in MOMENT_KEYS}
    assert max(weight_errors.values()) < 3e-10, weight_errors
    for k in MOMENT_KEYS:
        assert np.max(np.abs(summed['original_source_weighted_diagnostic'][k]-1.49*np.asarray(actual[k]))) < 3e-10
    try:
        aggregate_slots(grid, {(1, 1): actual}, geometry_policy='declared-geometric')
    except ValueError:
        pass
    else:
        raise AssertionError('incomplete source coverage admitted')
    mutations = dict(
        omit_kinetic_gradient=float(np.linalg.norm(actual['kinetic_gradient'])),
        omit_kinetic_contact=float(np.linalg.norm(actual['kinetic_contact'])),
        omit_nonlocal_gradient=float(np.linalg.norm(actual['gradient']-actual['kinetic_gradient'])),
        omit_mixed_contact_legs=float(abs(cross[0, 1])),
        reverse_current_sign=float(np.linalg.norm(actual['gradient'])),
        conjugate_trace_order=float(np.linalg.norm(wrong_conjugation['contact']-actual['contact'])),
        omit_D_raw_column_rephase=float(np.linalg.norm(wrong_rephase['contact']-actual['contact'])),
        reverse_atom_phase=float(np.linalg.norm(wrong_phase)),
        source_weight_applied_twice=float(np.linalg.norm(manual['contact']
            -(.3+.7*1.7)*(.25*actual['contact']+.75*rotated['contact']))),
        replace_full_slot_weight_with_original_source=float(np.linalg.norm(manual['gradient']
            -summed['original_source_weighted_diagnostic']['gradient'])),
        discard_imaginary_trace=float(np.linalg.norm(actual['contact'].imag)),
        equate_A25_to_physical_zero_field=float(abs(actual['delta_nonlocal_energy_ry'])))
    assert min(mutations.values()) > 1e-7, mutations
    return dict(status='passed', literal_dense_errors=errors,
                dense_finite_difference_refinement=refinements, covariance_errors=covariance_errors,
                weight_errors=weight_errors, mutation_exposures=mutations,
                anisotropic_channels_rejected=True, incomplete_sources_rejected=True,
                scope='synthetic independent literal radial/dense oracle; no actual material or field state')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run()
    result.update(registration=binding(HERE/'full-field-moments-registration.json'),
                  implementation_sources=implementation_bindings())
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=result['status'], output=binding(args.output))))


if __name__ == '__main__':
    main()
