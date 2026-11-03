"""Independent finite-basis projection/coordinate oracle; not material execution."""
from pathlib import Path
from types import SimpleNamespace
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[3]/'src'), str(HERE)]
import numpy as np
import scipy.linalg as la
import sympy as sp
from reduced_active_basis import prepare_active_basis
from reduced_projection_algebra import q_feature_geometry, profile, projector_overlaps
from reduced_local_projection import project_edge
from reduced_prepared_projection import PreparedIonicEdge, overlap_arrays, restore_overlaps
from substrate_framework.qe_material_grid import PhysicalRhoFFT
from substrate_framework.qe_bloch_finiteq import BlochFiniteQBridge
from substrate_framework.qe_bloch_finiteq_ionic import BlochFiniteQIonic


def relative(a, b):
    return float(la.norm(a-b)/max(1., la.norm(a), la.norm(b)))


def oracle():
    rng = np.random.default_rng(290031)
    z = lambda shape: rng.normal(size=shape)+1j*rng.normal(size=shape)
    # A fully occupied window is exactly Pauli-inert even with nonzero vertices.
    f0, f1, coupling = sp.symbols('f0 f1 coupling', real=True)
    rho_symbolic = sp.diag(f0, f1)
    vertex_symbolic = sp.Matrix([[0, coupling], [coupling, 0]])
    commutator = vertex_symbolic*rho_symbolic-rho_symbolic*vertex_symbolic
    assert sp.simplify(commutator[0, 1]-coupling*(f1-f0)) == 0
    assert commutator.subs({f0:1, f1:1}) == sp.zeros(2)
    assert commutator.subs({f0:1, f1:0, coupling:1}) != sp.zeros(2)
    # Unequal saved gaps make the old nearest12 recipe select only filled rows.
    contrast_labels = np.r_[np.linspace(-.12, -.01, 16), np.linspace(1.5, 3., 8)]
    contrast_populations = np.r_[np.ones(16), np.zeros(8)]
    contrast_basis = prepare_active_basis(np.eye(24, dtype=complex), np.diag(contrast_populations).astype(complex),
        contrast_labels, (0., 0.), 12, populations=contrast_populations)
    contrast8 = prepare_active_basis(np.eye(24, dtype=complex), np.diag(contrast_populations).astype(complex),
        contrast_labels, (0., 0.), 8, populations=contrast_populations)
    assert np.array_equal(contrast8.selected, contrast_basis.selected[:8])
    assert set(contrast_populations[contrast_basis.selected]) == {0., 1.}
    nearest = np.argsort(np.abs(contrast_labels))[:12]
    old_rho = np.diag(contrast_populations[nearest])
    new_rho = contrast_basis.rho0
    probe = np.ones((12, 12), complex)-np.eye(12)
    old_drive = float(la.norm(probe @ old_rho-old_rho @ probe))
    contrast_drive = float(la.norm(probe @ new_rho-new_rho @ probe))
    assert old_drive == 0 and contrast_drive > 1
    # Independent two-term NC derivative algebra before numerical composition.
    bt, bs, d, pt, ps = sp.symbols('bt bs d pt ps', commutative=False)
    assert sp.expand(sp.I*(bt*d*bs*ps-pt*bt*d*bs)) == sp.I*bt*d*bs*ps-sp.I*pt*bt*d*bs
    n, m = 16, 12
    unitary = la.qr(z((23, n)), mode='economic')[0]
    a = z((n, n)); s = np.eye(n)+a.conj().T @ a/40
    root = la.sqrtm(s).astype(complex)
    inverse = la.inv(root)
    psi = unitary @ root
    populations = np.linspace(.95, .05, n)
    d_raw = inverse @ np.diag(populations) @ inverse
    labels = np.linspace(-1., 1., n)
    basis = prepare_active_basis(psi, d_raw, labels, (-.05, .05), m, populations=populations)
    selected = basis.selected
    expected_rho = unitary[:, selected].conj().T @ (psi @ d_raw @ psi.conj().T) @ unitary[:, selected]
    basis_defect = relative(basis.rho0, expected_rho)
    principal_mutation = relative(d_raw[np.ix_(selected, selected)], expected_rho)
    assert basis_defect < 1e-11 and principal_mutation > 1e-3
    nested = prepare_active_basis(psi, d_raw, labels, (-.05, .05), 8)
    assert np.array_equal(nested.selected, selected[:8])
    closure = 0.
    for iq in range(64):
        index = np.asarray(np.unravel_index(iq, (4, 4, 4)))
        j = int(np.ravel_multi_index(tuple(-index % 4), (4, 4, 4)))
        q = index/4.; q -= np.floor(q+.5)
        qm = np.asarray(np.unravel_index(j, (4, 4, 4)))/4.; qm -= np.floor(qm+.5)
        g, p = q_feature_geometry(iq, 27); gm, pm = q_feature_geometry(j, 27)
        closure = max(closure, float(la.norm(gm[p[:len(g)]]+qm+g+q)))
        assert np.array_equal(pm[p], np.arange(len(p)))
        small, _ = q_feature_geometry(iq, 7)
        assert set(map(tuple, small)) <= set(map(tuple, g))
    assert closure == 0
    # A complete finite PW/grid example exercises the REAL existing physical
    # local/NC action implementation as witness against the batched contraction.
    shape = (5, 5, 5); volume = 7.3
    miller = np.indices(shape).reshape(3, -1)-2
    fft = PhysicalRhoFFT.create(shape, miller, np.eye(3), volume)
    edge = SimpleNamespace(source_dimension=125, target_dimension=125, ispin=1)
    bridge = BlochFiniteQBridge.__new__(BlochFiniteQBridge)
    bridge.edge = edge; bridge.block_size = 4; bridge.ispin = 1
    bridge.source_fft = fft; bridge.target_fft = fft
    bridge.density_fft = fft.at_q(np.array([.25, -.25, 0.]))
    bridge.grid_measure_bohr3 = volume/125
    grid = np.indices(shape)
    bridge.multiplication_phase = np.exp(-2j*np.pi*(grid[0]/5-grid[1]/5)).reshape(-1, order='F')
    nat, width, nh = 2, 3, 2
    source = la.qr(z((125, width)), mode='economic')[0]
    target = la.qr(z((125, width)), mode='economic')[0]
    dij = z((nh, nh)); dij += dij.conj().T
    source_species = SimpleNamespace(name='A', atom_indices=(1, 2), active_d=dij, beta=z((125, nat*nh))/10)
    target_species = SimpleNamespace(name='A', atom_indices=(1, 2), active_d=dij, beta=z((125, nat*nh))/10)
    ionic = BlochFiniteQIonic.__new__(BlochFiniteQIonic)
    ionic.bridge = bridge; ionic.edge = edge; ionic.nat = nat; ionic.ndof = 3*nat
    ionic.source_factors = SimpleNamespace(species=(source_species,))
    ionic.target_factors = SimpleNamespace(species=(target_species,))
    ionic.species = ((source_species, target_species),)
    ionic.ps = miller.T+.1; ionic.pt = miller.T-.2
    ionic.gq = bridge.density_fft.wavevectors.T
    ionic.vloc = {'A': np.exp(-np.linalg.norm(ionic.gq, axis=1)**2/4)}
    ionic.inputs = SimpleNamespace(reference=SimpleNamespace(atom_names=('A', 'A'), positions_bohr=rng.normal(size=(3, nat))))
    gs, _ = q_feature_geometry(0, 7)
    projected, timings = project_edge(ionic, source, target, gs)
    forward_defect = reverse_defect = 0.
    for coordinate in range(3*nat):
        u = np.zeros((nat, 3)); u.reshape(-1)[coordinate] = 1
        action = ionic.prepare(u)
        forward = target.conj().T @ action.apply(source)
        reverse = (source.conj().T @ action.adjoint_apply(target)).T
        forward_defect = max(forward_defect, relative(forward, projected['ionic'][coordinate]))
        reverse_defect = max(reverse_defect, relative(reverse, projected['ionic_reverse'][coordinate]))
    potential_defect = 0.
    for j, g in enumerate(gs):
        v = profile(shape, g, volume)
        direct = target.conj().T @ bridge.local_apply(v, source)
        potential_defect = max(potential_defect, relative(direct, projected['potential'][:, :, j]))
    assert max(forward_defect, reverse_defect, potential_defect) < 1e-11
    source_overlaps = projector_overlaps(ionic.source_factors, source, ionic.ps)
    target_overlaps = projector_overlaps(ionic.target_factors, target, ionic.pt)
    cached_source = restore_overlaps(overlap_arrays(source_overlaps), sorted(source_overlaps))
    cached_target = restore_overlaps(overlap_arrays(target_overlaps), sorted(target_overlaps))
    prepared_context = PreparedIonicEdge(ionic.inputs, edge, bridge, nat, 3*nat, ionic.gq,
        ionic.vloc, ionic.species, cached_source, cached_target, {})
    cached_projected, _ = project_edge(prepared_context, source, target, gs)
    cached_overlap_defect = max(relative(projected[name], cached_projected[name]) for name in projected)
    assert cached_overlap_defect < 1e-13
    x = z((width, width))
    density = bridge.density_from_dyads(target, source, x, weight=1.)
    coeff = np.einsum('rab,ab->r', projected['density'], x)
    density_defect = max(abs(coeff[j]-bridge.grid_measure_bohr3*np.vdot(profile(shape, g, volume), density[:, 0])) for j, g in enumerate(gs))
    assert density_defect < 1e-11
    wrong_transpose = relative(projected['potential'], projected['density'].transpose(1, 2, 0))
    assert wrong_transpose > 1e-3
    return dict(schema='p254-reduced-projection-finite-oracle-v1',
        basis_coordinate_defect=basis_defect, principal_raw_d_mutation=principal_mutation,
        occupation_contrast_drive_norm=contrast_drive, nearest_all_filled_mutation_drive_norm=old_drive,
        cached_overlap_projection_defect=cached_overlap_defect,
        feature_physical_closure_defect=closure, ionic_forward_defect=forward_defect,
        ionic_independent_reverse_defect=reverse_defect, potential_density_pairing_defect=potential_defect,
        density_linear_coefficient_defect=float(density_defect),
        missing_conjugation_mutation=wrong_transpose, physical_material_execution=False)


if __name__ == '__main__':
    print(json.dumps(oracle(), indent=2, allow_nan=False))
