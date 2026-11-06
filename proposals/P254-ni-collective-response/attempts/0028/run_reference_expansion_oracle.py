#!/usr/bin/env python3
"""Independent literal full-PW, scalar-differentiation and mutation oracles.

No material execution. Root runs this once source settles. All dense matrices
below are deliberately tiny independent oracle objects, never production APIs.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace as NS
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3]/'src'))


def main():
    import numpy as np
    import hashlib
    def binding(path):
        path = Path(path).resolve()
        data = path.read_bytes()
        return dict(path=str(path), bytes=str(len(data)), sha256=hashlib.sha256(data).hexdigest())
    implementation_paths = [Path(__file__), HERE/'reference-expansion-registration.json',
        HERE.parents[3]/'src/substrate_framework/material_reference_transport.py',
        HERE.parents[3]/'src/substrate_framework/material_frozen_state_curvature.py']
    implementation_sources = [binding(path) for path in implementation_paths]
    import sympy as sp
    from substrate_framework.material_reference_state import SourceState
    from substrate_framework.material_reference_transport import (
        slot_transport, raw_state_transport, SourceObservables, aggregate_complete_slots)
    from substrate_framework.material_frozen_state_curvature import (
        frozen_onsite_curvature, local_density_curvature, conditional_minus_mv,
        common_reference_correction)
    from substrate_framework.qe_material_grid import PhysicalRhoFFT
    from substrate_framework.qe_bloch_grid import BlochGrid
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(22528)
    results, mutations = {}, {}

    def close(name, actual, expected, tol=2e-10):
        error = float(np.max(np.abs(np.asarray(actual)-expected)))
        results[name] = error
        if error > tol:
            raise AssertionError((name, error, tol))

    def killed(name, wrong, expected, minimum=1e-7):
        gap = float(np.max(np.abs(np.asarray(wrong)-expected)))
        mutations[name] = gap
        if gap < minimum:
            raise AssertionError(('mutation survived', name, gap))

    def rejected(name, action):
        try:
            action()
        except ValueError:
            mutations[name] = 'rejected'
        else:
            raise AssertionError(('mutation accepted', name))

    # Independent differentiated scalar proves the complete four-term kernel.
    u, v, pa, pb, qa, qb = sp.symbols('u v pa pb qa qb', real=True)
    kernel = sp.exp(-sp.I*((pa-qa)*u+(pb-qb)*v))
    second = sp.diff(kernel, u, v).subs({u: 0, v: 0})
    assert sp.expand(second - (-pa*pb-qa*qb+pa*qb+pb*qa)) == 0
    results['sympy_four_nonlocal_terms'] = True

    shape, omega = (4, 4, 4), 7.3
    m = np.array(list(__import__('itertools').product((-1, 0, 1), repeat=3))).T
    reciprocal = 2*np.pi*np.eye(3)
    fft = PhysicalRhoFFT.create(shape, m, reciprocal, omega)
    coords = np.indices(shape).reshape(3, -1, order='F')/4
    fourier = np.exp(2j*np.pi*(coords.T @ m))  # literal transform, no FFT API
    psi = rng.normal(size=(len(m.T), 4)) + 1j*rng.normal(size=(len(m.T), 4))
    a = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
    coefficient = a @ a.conj().T/100
    gram = psi.conj().T @ psi
    state = SourceState(gram, np.arange(4)/5, coefficient, {}, 'source')
    rho = psi @ coefficient @ psi.conj().T
    density = np.einsum('ri,ij,rj->r', fourier, rho, fourier.conj())/omega
    # This rotation is not self-inverse: transposes and atom inverse maps fail.
    s = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    f = np.array([.25, 0., 0.])
    points = [np.array([.125, .375, .25])]
    for _ in range(3):
        points.append((s.T @ points[-1]-f) % 1)
    positions = np.array(points).T
    atoms = np.array([1, 2, 3, 0])
    op = NS(rotation=s, fractional_translation=f, equivalent_atoms=atoms+1,
            effective_time_reversal=False, index=2)
    identity = NS(rotation=np.eye(3, dtype=int), fractional_translation=np.zeros(3),
                  equivalent_atoms=np.arange(4)+1, effective_time_reversal=False, index=1)
    ref = NS(fft_smooth=shape, cell_bohr=np.eye(3), nat=4, positions_bohr=positions,
             atom_names=('X',)*4, nbnd_by_spin=(4, 4))
    slots = []
    for spin in (1, 2):
        for index, operation in enumerate((identity, op)):
            image = np.rint(operation.rotation.T @ positions
                            - operation.fractional_translation[:, None]
                            - positions[:, operation.equivalent_atoms-1]).astype(int)
            slots.append(NS(ispin=spin, mesh_index=(index, 0, 0), source_ik=1,
                operation=operation, reciprocal_shift=np.zeros(3, int),
                atom_image_offsets=image, source_binding={'sha256': 'source'},
                weight=.5, mesh_slot=index))
    grid = BlochGrid(ref, tuple(slots), 'oracle')
    transform = slot_transport(grid, slots[1])
    ks = np.array([.125, -.25, .0])
    shift = np.array([1, 0, 0])
    target_psi = psi*np.exp(-2j*np.pi*((m.T+ks) @ f))[:, None]
    target_m = s @ m-shift[:, None]
    target_fourier = np.exp(2j*np.pi*(coords.T @ target_m))
    target_rho = target_psi @ coefficient @ target_psi.conj().T
    target_density = np.einsum('ri,ij,rj->r', target_fourier, target_rho, target_fourier.conj())/omega
    close('literal_phase_fft_density', transform.density(density), target_density)
    close('unitary_gram', target_psi.conj().T @ target_psi, gram)
    close('unitary_D_unchanged', raw_state_transport(state).density_coefficients, coefficient)
    killed('wrong_translation_sign', density[np.ravel_multi_index(tuple(((s.T @ (coords*4).astype(int)
          + (f*4).astype(int)[:, None]) % 4)), shape, order='F')], target_density)
    killed('wrong_fft_order', target_density.reshape(shape, order='F').reshape(-1, order='C'), target_density)
    conjugated = raw_state_transport(state, antiunitary=True)
    close('antiunitary_D_and_Gram', psi.conj() @ conjugated.density_coefficients @ psi.T, rho.conj())
    close('antiunitary_gram', conjugated.gram, gram.conj())
    killed('unconjugated_antiunitary_Gram', gram, conjugated.gram)
    killed('unconjugated_antiunitary_D', psi.conj() @ coefficient @ psi.T, rho.conj())
    badop = NS(**vars(op)); badop.effective_time_reversal = True
    badslot = NS(**vars(slots[1])); badslot.operation = badop
    rejected('physical_antiunitary_slot', lambda: slot_transport(BlochGrid(ref, (badslot,), 'oracle'), badslot))

    # A local scalar energy and its independent gradient at each actual atom.
    g = reciprocal @ m
    radial = np.exp(-np.sum(g*g, axis=0)/90)/omega
    def local_gradient(n):
        result = []
        for position in positions.T:
            phase = radial*np.exp(-1j*(g.T @ position))
            fields = fourier @ (-1j*g*phase).T
            result.append(omega/fft.nnr*np.einsum('r,ri->i', n, fields))
        return np.array(result).reshape(-1)
    gradient = local_gradient(density)
    direct_gradient = local_gradient(target_density)
    close('atom_vector_gradient', transform.gradient(gradient), direct_gradient)
    killed('missing_atom_permutation', (gradient.reshape(-1, 3) @ s.T).reshape(-1), direct_gradient)
    killed('inverse_vector_rotation', (gradient.reshape(-1, 3)[atoms] @ s).reshape(-1), direct_gradient)
    ledger = dict(spin=np.array([q.ispin for q in slots]), source_ik=np.ones(4, int),
        mesh_index=np.array([q.mesh_index for q in slots]), source_sha256=np.array(['source']*4),
        measures=np.full(4, .5), populations=np.tile(state.populations, (4, 1)))
    observables = [SourceObservables(slots[i], state, density, gradient) for i in (0, 2)]
    aggregate = aggregate_complete_slots(grid, ledger, observables, gradient_policy='declared-geometric')
    close('full_spin_mesh_density', aggregate['density_bohr3'], np.tile(((density+target_density)/2)[:, None], (1, 2)))
    close('full_spin_mesh_gradient', aggregate['positive_gradient_ry_bohr'], gradient+direct_gradient)
    killed('missing_mesh_weight', 2*aggregate['positive_gradient_ry_bohr'], gradient+direct_gradient)
    rejected('incomplete_source_spin_coverage', lambda: aggregate_complete_slots(grid, ledger, observables[:1], gradient_policy='declared-geometric'))
    rejected('duplicate_source_slots', lambda: aggregate_complete_slots(grid, ledger, observables*2, gradient_policy='declared-geometric'))
    rejected('unproved_default_covariance', lambda: aggregate_complete_slots(grid, ledger, observables))
    ref.select = lambda source_ik, spin: (NS(weight=1.), None, None)
    conditional_blocks = rng.normal(size=(4, 3, 3))
    original_mv_blocks = rng.normal(size=(4, 3, 3))
    conditional_target = np.array([s @ conditional_blocks[b] @ s.T for b in atoms])
    curvature_observables = [SourceObservables(slots[i], state, density, gradient,
        conditional_blocks, original_mv_blocks) for i in (0, 2)]
    curvature_aggregate = aggregate_complete_slots(grid, ledger, curvature_observables,
                                                   gradient_policy='declared-geometric')
    exact_delta = conditional_blocks + conditional_target - 2*original_mv_blocks
    close('original_XML_MV_weight_not_star_transport',
          curvature_aggregate['nonlocal_correction_ry_bohr2'], exact_delta)
    killed('silently_star_transporting_original_MV',
           conditional_blocks + conditional_target - original_mv_blocks
           - np.array([s @ original_mv_blocks[b] @ s.T for b in atoms]), exact_delta)
    direct_aggregate = aggregate_complete_slots(grid, ledger, curvature_observables,
        direct_observables=lambda slot, state: dict(positive_gradient_ry_bohr=direct_gradient,
                                                    conditional_nonlocal_ry_bohr2=conditional_target))
    close('direct_slot_fallback_curvature', direct_aggregate['nonlocal_correction_ry_bohr2'], exact_delta)

    # Curvature tiny factors: two atoms, genuinely complex arbitrary dyad,
    # finite k, and literal independent full PW matrices and scalar differences.
    nat = 2
    p = (g+np.array([.23, -.17, .11])[:, None]).T
    beta = (rng.normal(size=(len(m.T), 4))+1j*rng.normal(size=(len(m.T), 4)))/9
    active_d = np.array([[.7, .2], [.2, -.4]])
    mv = np.array([1.02, .55, -.03, 0.])
    ionic = NS(p=p, g=g.T, vloc={}, basis=NS(slot=NS(source_binding={'sha256': 'source'}, raw_occupations=mv)),
        inputs=NS(reference=NS(nat=nat, atom_names=('X',)*nat, positions_bohr=positions[:, :nat]),
                  density=NS(miller=m), tables={'X': NS(vloc=lambda mod: np.exp(-mod**2/90)/omega)}))
    factors = NS(fft=fft, species=(NS(beta=beta, active_d=active_d, atom_indices=(1, 2)),))
    coupling = NS(ionic=ionic, factors=factors, npw=psi.shape[0], nr=fft.nnr, block_size=3)
    # Non-Hermitian coefficients catch illicit conjugation and real-part repair.
    complex_d = coefficient + .03j*np.eye(4)
    dense_rho = psi @ complex_d @ psi.conj().T
    result = frozen_onsite_curvature(psi, complex_d, coupling)
    expected_local = np.empty((nat, 3, 3), complex)
    expected_nl = np.empty_like(expected_local)
    expected_mv, expected_conditional = np.empty_like(expected_nl), np.empty_like(expected_nl)
    mv_rho = (psi*mv) @ psi.conj().T
    removed_second, removed_cross = np.empty_like(expected_local), np.empty_like(expected_local)
    def local_matrix(atom, displacement):
        phase = radial*np.exp(-1j*(g.T @ (positions[:, atom]+displacement)))
        field = fourier @ phase
        return fourier.conj().T @ (field[:, None]*fourier)/fft.nnr
    def nl_matrix(atom, displacement):
        b = np.exp(-1j*(p @ displacement))[:, None]*beta[:, 2*atom:2*atom+2]
        return b @ active_d @ b.conj().T
    for atom in range(nat):
        b = beta[:, 2*atom:2*atom+2]
        for i in range(3):
            for j in range(3):
                b_i, b_j = -1j*p[:, i, None]*b, -1j*p[:, j, None]*b
                b_ij = -p[:, i, None]*p[:, j, None]*b
                legs = (b_ij @ active_d @ b.conj().T, b @ active_d @ b_ij.conj().T,
                        b_i @ active_d @ b_j.conj().T, b_j @ active_d @ b_i.conj().T)
                expected_mv[atom, i, j] = np.trace(mv_rho @ sum(legs))
                expected_conditional[atom, i, j] = np.trace(rho @ sum(legs))
                expected_nl[atom, i, j] = np.trace(dense_rho @ sum(legs))
                removed_second[atom, i, j] = np.trace(dense_rho @ (sum(legs)-legs[1]))
                removed_cross[atom, i, j] = np.trace(dense_rho @ (legs[0]+legs[1]))
                phase = -g[i]*g[j]*radial*np.exp(-1j*(g.T @ positions[:, atom]))
                h = fourier.conj().T @ ((fourier @ phase)[:, None]*fourier)/fft.nnr
                expected_local[atom, i, j] = np.trace(dense_rho @ h)
    close('literal_full_PW_four_nonlocal_legs', result.nonlocal_term, expected_nl)
    close('literal_full_PW_local', result.local, expected_local)
    killed('missing_second_nonlocal_leg', removed_second, expected_nl)
    killed('missing_nonlocal_cross_terms', removed_cross, expected_nl)
    killed('local_second_derivative_sign', -result.local, expected_local)
    killed('local_volume_measure', result.local/omega, expected_local)
    killed('real_part_projection', result.total.real, result.total)
    # Mixed centered scalar derivative is independent of the analytic kernels.
    direction_a, direction_b = np.array([.4, -.2, .3]), np.array([-.1, .35, .2])
    step = 2e-5
    for atom in range(nat):
        def energy(displacement):
            return np.trace(dense_rho @ (local_matrix(atom, displacement)+nl_matrix(atom, displacement)))
        finite = sum(sa*sb*energy(step*(sa*direction_a+sb*direction_b))
                     for sa in (-1, 1) for sb in (-1, 1))/(4*step*step)
        close('differentiated_scalar_atom_'+str(atom), direction_a @ result.total[atom] @ direction_b, finite, 2e-4)
    physical = conditional_minus_mv(psi, state, coupling, mv, conditional_density=density)
    literal_delta = frozen_onsite_curvature(psi, coefficient-np.diag(mv), coupling)
    close('conditional_minus_numerical_MV', physical['orbital_total_correction_ry_bohr2'], literal_delta.total)
    close('independent_original_MV_nonlocal', physical['original_mv_nonlocal_ry_bohr2'], expected_mv)
    close('independent_conditional_nonlocal', physical['conditional_nonlocal_ry_bohr2'], expected_conditional)
    killed('MV_as_conditional_preparation', np.zeros_like(literal_delta.total), literal_delta.total)
    clipped_delta = frozen_onsite_curvature(psi, coefficient-np.diag(np.clip(mv, 0, 1)), coupling)
    killed('clipping_raw_MV_numerical_coefficients', clipped_delta.total, literal_delta.total)
    # The actual saved local reference is deliberately NOT the orbital-MV rho.
    saved = np.column_stack((physical['mv_density_bohr3']+.01*np.cos(2*np.pi*coords[0]), np.zeros(fft.nnr)))
    conditional = np.column_stack((density, np.zeros(fft.nnr)))
    common = common_reference_correction(conditional, saved, physical['nonlocal_correction_ry_bohr2'], coupling)
    close('actual_saved_density_local_conversion', common.local,
          local_density_curvature(density-saved.sum(axis=1), coupling))
    killed('substituting_orbital_MV_for_saved_rho', literal_delta.total, common.total)
    if [binding(path) for path in implementation_paths] != implementation_sources:
        raise ValueError('oracle implementation changed during execution')
    record = dict(status='passed', oracle_residuals=results, killed_mutations=mutations,
                  implementation_sources=implementation_sources,
                  material_execution=False, scope='small independent algebra only; no full material claim')
    with args.output.open('x') as stream:
        json.dump(record, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    print(json.dumps(record, sort_keys=True))


if __name__ == '__main__':
    main()
