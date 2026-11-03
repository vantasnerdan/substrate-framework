"""Independent finite local/NC action oracle; no material files or Native calls.

Parent once-only command, after the live producer wave has settled::

    python proposals/P254-ni-collective-response/attempts/0029/run_reduced_local_projection_oracle.py

The reference builds explicit plane-wave evaluation matrices and full finite
multiplication/action/adjoint matrices. It does NOT use an FFT, production
feature/profile geometry, production local fields, or the pair-contraction
formula to generate the reference. Both even-grid Nyquist aliases and odd-grid
nonalias support are covered, with noncubic F-order grids, rectangular endpoint
PW bases, complex nonorthogonal active columns, arbitrary atom positions,
complex radial factors and NC projectors/D, both explicit spins, reciprocal
folds and independently closed opposite physical supports. This is a finite
mathematical qualification, NOT material accuracy or performance evidence.
"""
from pathlib import Path
from types import SimpleNamespace
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[3]/'src'), str(HERE)]

import numpy as np
import sympy as sp

from reduced_local_projection import (project_edge, project_local_pairs,
                                      prepare_local_projection)


def relative(left, right):
    return float(np.linalg.norm(left-right)/max(1., np.linalg.norm(left), np.linalg.norm(right)))


def symbolic_alias_pairing():
    """Small exact DFT/alias and physical density normalization identities."""
    values = sp.symbols('n0:4', complex=True)
    ma, mb = sp.symbols('ma mb', complex=True)
    transform = lambda g: sum(values[j]*sp.I**(-j*g) for j in range(4))/4
    assert sp.simplify(transform(2)-transform(-2)) == 0
    assert sp.simplify((ma+mb)*transform(2)
                       -ma*transform(2)-mb*transform(-2)) == 0
    conjugate_transform = sum(sp.conjugate(values[j])*sp.I**(2*j) for j in range(4))/4
    assert sp.simplify(sp.conjugate(transform(2))-conjugate_transform) == 0
    omega = sp.symbols('Omega', positive=True)
    t, s, p, phase_g = sp.symbols('t s p phase_g', complex=True)
    density = t*sp.conjugate(s)*sp.conjugate(p)/omega
    feature_integrand = omega*sp.conjugate(phase_g/sp.sqrt(omega))*density
    action_integrand = sp.conjugate(t)*p*(phase_g/sp.sqrt(omega))*s
    assert sp.simplify(sp.conjugate(feature_integrand)-action_integrand) == 0


def evaluation(points, millers):
    """Explicit e^(+2 pi i G.x), entirely independent of any FFT mapper."""
    return np.exp(2j*np.pi*(points @ np.asarray(millers).T))


def closed_support(original, q, qm):
    shift = np.rint(q+qm).astype(int)
    items = [tuple(row) for row in original]
    for row in original:
        opposite = tuple(-row-shift)
        if opposite not in items:
            items.append(opposite)
    return np.asarray(items, dtype=int)


def finite_case(rng, shape, q, qm, ks, spin, alias):
    z = lambda shape: rng.normal(size=shape)+1j*rng.normal(size=shape)
    shape = tuple(shape)
    nnr = int(np.prod(shape))
    # The list is explicitly generated in F order, independently of reshape.
    integers = np.asarray([(x, y, z_) for z_ in range(shape[2])
                           for y in range(shape[1]) for x in range(shape[0])])
    points = integers/np.asarray(shape)
    all_pw = integers-np.asarray(shape)//2
    source_g = all_pw[::2]
    target_g = all_pw[1::2]
    nsource, ntarget, m, nat, nh = len(source_g), len(target_g), 3, 2, 2
    source, target = z((nsource, m))/4, z((ntarget, m))/4
    omega = 7.3
    fold = np.floor(ks+q+.5).astype(int)
    kt = ks+q-fold
    phase = np.exp(-2j*np.pi*(points @ fold))
    if alias:
        original = np.asarray([(0, 0, 0), (-2, 0, 0), (2, 0, 0), (6, 0, 0),
                               (1, -1, 0), (0, 0, 1), (0, 0, -1), (0, 2, 1)])
    else:
        original = np.asarray([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
    support_g = closed_support(original, q, qm)
    opposite_g = closed_support(original, qm, q)
    # This tests physical q/-q closure without inferring a projected opposite
    # edge or invoking a source-star symmetry in the implementation.
    for row in support_g:
        assert any(np.allclose(row+q, -other-qm, rtol=0, atol=1e-14)
                   for other in opposite_g)
    reciprocal = np.asarray([[1.2, .1, -.2], [.0, .9, .15], [.05, -.1, 1.1]])
    gq = (support_g+q) @ reciprocal.T
    positions = np.asarray([[.13, -.24], [.47, .08], [-.32, .61]])
    radial = z((len(support_g),))*.2
    slots = lambda rows: np.ravel_multi_index(tuple((rows % shape).T), shape, order='F')
    mapping = lambda rows: SimpleNamespace(physical=shape, nnr=nnr, omega=omega, nl=slots(rows))
    support = mapping(support_g)
    support.wavevectors = gq.T
    source_map, target_map = mapping(source_g), mapping(target_g)
    bridge = SimpleNamespace(source_fft=source_map, target_fft=target_map,
                             density_fft=support, multiplication_phase=phase)
    d = z((nh, nh))*.2  # Deliberately no Hermitian assumption in the algebra.
    bs, bt = z((nsource, nat*nh))*.1, z((ntarget, nat*nh))*.1
    ss = SimpleNamespace(name='A', atom_indices=(1, 2), active_d=d, beta=bs)
    ts = SimpleNamespace(name='A', atom_indices=(1, 2), active_d=d, beta=bt)
    ps, pt = (source_g+ks) @ reciprocal.T, (target_g+kt) @ reciprocal.T
    ionic = SimpleNamespace(bridge=bridge, edge=SimpleNamespace(ispin=spin),
        ndof=3*nat, nat=nat, gq=gq, vloc={'A': radial}, species=((ss, ts),),
        source_factors=SimpleNamespace(species=(ss,)),
        target_factors=SimpleNamespace(species=(ts,)), ps=ps, pt=pt,
        inputs=SimpleNamespace(reference=SimpleNamespace(
            atom_names=('A', 'A'), positions_bohr=positions)))
    base_features = np.asarray([(0, 0, 0), (1, 0, 0), (-1, 0, 0),
                                (0, 0, 1), (0, 0, -1)])
    features = closed_support(base_features, q, qm)
    es, et, eg = evaluation(points, source_g), evaluation(points, target_g), evaluation(points, support_g)
    sw, tw = es @ source, et @ target
    density_grid = (tw[:, :, None]*sw[:, None, :].conj()*phase.conj()[:, None, None]/omega)
    profiles = evaluation(points, features)/np.sqrt(omega)
    density = np.zeros((2*len(features), m, m), complex)
    density[(spin-1)*len(features):spin*len(features)] = np.einsum(
        'rg,rab->gab', profiles.conj(), density_grid)*omega/nnr
    potential = np.zeros((m, m, 2*len(features)), complex)
    for j in range(len(features)):
        full_action = et.conj().T @ ((phase*profiles[:, j])[:, None]*es)/nnr
        potential[:, :, (spin-1)*len(features)+j] = target.conj().T @ full_action @ source
    local_forward = np.empty((3*nat, m, m), complex)
    local_reverse = np.empty_like(local_forward)
    full_forward = np.empty_like(local_forward)
    full_reverse = np.empty_like(local_forward)
    fields = np.empty((nnr, 3*nat), complex)
    coefficients = np.empty((len(support_g), 3*nat), complex)
    for atom in range(nat):
        for axis in range(3):
            c = 3*atom+axis
            # Physical derivative of the scalar ionic potential, generated
            # directly rather than from M or a projected gradient formula.
            coefficients[:, c] = (-1j*gq[:, axis]*radial*np.exp(-1j*(gq @ positions[:, atom])))
            field = eg @ coefficients[:, c]
            fields[:, c] = field
            local_action = et.conj().T @ ((phase*field)[:, None]*es)/nnr
            local_adjoint = es.conj().T @ ((phase*field).conj()[:, None]*et)/nnr
            local_forward[c] = target.conj().T @ local_action @ source
            local_reverse[c] = (source.conj().T @ local_adjoint @ target).T
            bsa = bs[:, atom*nh:(atom+1)*nh]
            bta = bt[:, atom*nh:(atom+1)*nh]
            nc = 1j*(bta @ d @ (bsa.conj().T*ps[:, axis])
                      -(pt[:, axis, None]*bta) @ d @ bsa.conj().T)
            # Separately assembled adjoint action; no projected matrix repair.
            nc_adj = -1j*((ps[:, axis, None]*bsa) @ d.conj().T @ bta.conj().T
                           -bsa @ d.conj().T @ (bta.conj().T*pt[:, axis]))
            full_forward[c] = target.conj().T @ (local_action+nc) @ source
            full_reverse[c] = (source.conj().T @ (local_adjoint+nc_adj) @ target).T
    expected_local = dict(ionic=local_forward, ionic_reverse=local_reverse,
                          density=density, potential=potential)
    expected_full = dict(ionic=full_forward, ionic_reverse=full_reverse,
                         density=density, potential=potential)
    context = SimpleNamespace(points=points, es=es, et=et, phase=phase, profiles=profiles,
        source=source, target=target, density_grid=density_grid, coefficients=coefficients,
        eg=eg, fields=fields, nnr=nnr, omega=omega, slots=support.nl,
        features=features, support_g=support_g, spin=spin, shape=shape)
    return ionic, source, target, features, expected_local, expected_full, context


def mutations(context, expected):
    """Plausible wrong physics/conventions must be visibly rejected."""
    c = context
    project_fields = lambda fields: np.asarray([
        c.target.conj().T @ (c.et.conj().T @ ((c.phase*fields[:, j])[:, None]*c.es)/c.nnr)
        @ c.source for j in range(fields.shape[1])])
    # Last-write-wins loses distinct physical Nyquist/extended-support terms.
    last = {}
    for i, slot in enumerate(c.slots):
        last[int(slot)] = i
    keep = np.asarray(list(last.values()))
    wrong_alias = project_fields(c.eg[:, keep] @ c.coefficients[keep])
    retained = np.asarray([tuple(row) in set(map(tuple, c.features)) for row in c.support_g])
    wrong_cutoff = project_fields(c.eg[:, retained] @ c.coefficients[retained])
    wrong_fold_density = np.zeros_like(expected['density'])
    active = slice((c.spin-1)*len(c.features), c.spin*len(c.features))
    # Replacing conj(p) by p in the directed density is a reciprocal-fold bug.
    wrong_fold_density[active] = np.einsum('rg,rab->gab', c.profiles.conj(),
        c.density_grid*c.phase[:, None, None]**2)*c.omega/c.nnr
    wrong_order_points = np.indices(c.shape).reshape(3, -1).T/np.asarray(c.shape)
    wrong_order_fields = evaluation(wrong_order_points, c.support_g) @ c.coefficients
    return dict(alias_overwrite=relative(wrong_alias, expected['ionic']),
        ionic_feature_cutoff=relative(wrong_cutoff, expected['ionic']),
        fold_sign=relative(wrong_fold_density, expected['density']),
        f_order=relative(project_fields(wrong_order_fields), expected['ionic']),
        density_missing_volume=relative(c.omega*expected['density'], expected['density']),
        extra_spin_factor=relative(2*expected['density'], expected['density']),
        potential_missing_conjugation=relative(expected['density'].transpose(1, 2, 0), expected['potential']),
        reverse_replaced_by_forward=relative(expected['ionic'], expected['ionic_reverse']))


def oracle():
    symbolic_alias_pairing()
    rng = np.random.default_rng(2903051)
    cases = [
        ((4, 3, 2), np.asarray([.25, -.25, .25]), np.asarray([-.25, .25, -.25]),
         np.asarray([.4, -.4, .4]), 1, True),
        ((4, 3, 2), np.asarray([-.5, 0., -.5]), np.asarray([-.5, 0., -.5]),
         np.asarray([-.3, .1, -.3]), 2, True),
        ((3, 5, 3), np.zeros(3), np.zeros(3), np.asarray([.1, -.2, .3]), 2, False),
    ]
    worst = {key: 0. for key in ('ionic', 'ionic_reverse', 'density', 'potential')}
    worst_local = dict(worst)
    mutation_report = []
    aliases = []
    for args in cases:
        ionic, source, target, gs, expected_local, expected_full, context = finite_case(rng, *args)
        prepared = prepare_local_projection(ionic)
        ionic.local_gradient = prepared.original_gradient
        aliases.append(prepared.alias_count)
        for batch in (1, 4, 8):  # includes cross-band boundaries and ragged tails
            local, _ = project_local_pairs(ionic, source, target, gs, prepared=prepared, batch_size=batch)
            actual, _ = project_edge(ionic, source, target, gs, prepared=prepared, batch_size=batch)
            for key in worst:
                worst_local[key] = max(worst_local[key], relative(local[key], expected_local[key]))
                worst[key] = max(worst[key], relative(actual[key], expected_full[key]))
        # Cached overlaps must give the same FULL physical matrix, not merely
        # agree with another use of the production projector-overlap helper.
        for name, factors, columns, momenta in (
                ('source_overlaps', ionic.source_factors, source, ionic.ps),
                ('target_overlaps', ionic.target_factors, target, ionic.pt)):
            species = factors.species[0]
            base = np.stack([species.beta[:, j*2:(j+1)*2].conj().T @ columns for j in range(2)])
            derivative = np.stack([np.stack([
                species.beta[:, j*2:(j+1)*2].conj().T @ (momenta[:, axis, None]*columns)
                for j in range(2)]) for axis in range(3)])
            setattr(ionic, name, {'A': (base, derivative)})
        cached, _ = project_edge(ionic, source, target, gs, prepared=prepared, batch_size=4)
        for key in worst:
            worst[key] = max(worst[key], relative(cached[key], expected_full[key]))
        if args[-1]:
            mutation_report.append(mutations(context, expected_local))
    assert max(worst.values()) < 3e-12, worst
    assert max(worst_local.values()) < 3e-12, worst_local
    # Each alias/fold case must independently expose every stated mutation.
    for report in mutation_report:
        assert min(report.values()) > 1e-5, report
    return dict(full_action_relative_defects=worst, local_action_relative_defects=worst_local,
                additive_alias_counts=aliases, mutation_relative_defects=mutation_report,
                physical_material_execution=False, performance_measured=False)


if __name__ == '__main__':
    print(json.dumps(oracle(), indent=2, sort_keys=True))
