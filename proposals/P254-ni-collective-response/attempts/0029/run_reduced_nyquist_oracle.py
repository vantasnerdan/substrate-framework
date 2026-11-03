"""Independent full-support finite Nyquist action/adjoint coverage; no material run."""
from pathlib import Path
from types import SimpleNamespace
import json
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[3]/'src'), str(HERE)]
import numpy as np
import scipy.linalg as la
from substrate_framework.qe_material_grid import PhysicalRhoFFT
from substrate_framework.qe_bloch_finiteq import BlochFiniteQBridge
from substrate_framework.qe_bloch_finiteq_ionic import BlochFiniteQIonic
from reduced_paired_local_support import create_support, paired_millers
from reduced_projection_algebra import q_feature_geometry
from reduced_local_projection import project_edge


def fixture(q, qm, rng, original, *, paired):
    shape, omega = (7, 7, 7), 11.7
    orbital = np.indices(shape).reshape(3, -1)-3
    fft = PhysicalRhoFFT.create(shape, orbital, np.eye(3), omega)
    density_fft = PhysicalRhoFFT.create(shape, original, np.eye(3), omega, q)
    support = create_support(density_fft, q, qm) if paired else density_fft
    ks = np.array([.25, 0., 0.])
    kt = ks+q; kt -= np.floor(kt+.5)
    fold = np.rint(kt-ks-q).astype(int)
    edge = SimpleNamespace(source_dimension=343, target_dimension=343, ispin=1,
                           q_fractional=q, reciprocal_fold=fold)
    bridge = BlochFiniteQBridge.__new__(BlochFiniteQBridge)
    bridge.edge = edge; bridge.block_size = 4; bridge.ispin = 1
    bridge.source_fft = fft; bridge.target_fft = fft; bridge.density_fft = support
    bridge.grid_measure_bohr3 = omega/343
    grid = np.indices(shape)
    bridge.multiplication_phase = np.exp(-2j*np.pi*sum(fold[a]*grid[a]/shape[a] for a in range(3))).reshape(-1, order='F')
    ionic = BlochFiniteQIonic.__new__(BlochFiniteQIonic)
    ionic.bridge = bridge; ionic.edge = edge; ionic.nat = 1; ionic.ndof = 3
    ionic.species = (); ionic.source_factors = SimpleNamespace(species=())
    ionic.target_factors = SimpleNamespace(species=())
    ionic.ps = orbital.T+ks; ionic.pt = orbital.T+kt
    ionic.gq = support.wavevectors.T
    ionic.vloc = {'A': np.exp(-np.sum(ionic.gq**2, axis=1)/3)}
    ionic.inputs = SimpleNamespace(reference=SimpleNamespace(atom_names=('A',), positions_bohr=np.array([[.13], [.27], [-.19]])))
    return ionic, ks, kt


def reverse_fixture(forward, qm, q, original, *, paired):
    reverse, _, _ = fixture(qm, q, None, original, paired=paired)
    ks, kt = forward.ps[0]-forward.bridge.source_fft.miller[:, 0], forward.pt[0]-forward.bridge.target_fft.miller[:, 0]
    fold = np.rint(ks-kt-qm).astype(int)
    reverse.edge.reciprocal_fold = fold
    grid = np.indices(reverse.bridge.source_fft.physical)
    reverse.bridge.multiplication_phase = np.exp(-2j*np.pi*sum(fold[a]*grid[a]/7 for a in range(3))).reshape(-1, order='F')
    reverse.ps, reverse.pt = forward.pt, forward.ps
    return reverse


def oracle():
    rng = np.random.default_rng(291101)
    original = np.indices((3, 3, 3)).reshape(3, -1)-1
    source = la.qr(rng.normal(size=(343, 3))+1j*rng.normal(size=(343, 3)), mode='economic')[0]
    target = la.qr(rng.normal(size=(343, 3))+1j*rng.normal(size=(343, 3)), mode='economic')[0]
    cases = [(np.array([-.5, 0., 0.]), np.array([-.5, 0., 0.])),
             (np.array([-.5, .25, 0.]), np.array([-.5, -.25, 0.])),
             (np.array([-.5, -.5, -.5]), np.array([-.5, -.5, -.5]))]
    reports = []
    for q, qm in cases:
        projected = []
        for paired in (False, True):
            forward, _, _ = fixture(q, qm, rng, original, paired=paired)
            reverse = reverse_fixture(forward, qm, q, original, paired=paired)
            gs, _ = q_feature_geometry(0, 7)
            pf, _ = project_edge(forward, source, target, gs)
            pr, _ = project_edge(reverse, target, source, gs)
            defect = la.norm(pr['ionic']-pf['ionic'].conj().transpose(0, 2, 1))
            action_defect = 0.
            for axis in range(3):
                u = np.zeros((1, 3)); u[0, axis] = 1
                af = forward.prepare(u); ar = reverse.prepare(u)
                direct = target.conj().T @ af.apply(source)
                opposite = source.conj().T @ ar.apply(target)
                action_defect = max(action_defect, la.norm(opposite-direct.conj().T),
                    la.norm(pf['ionic'][axis]-direct),
                    la.norm(pr['ionic'][axis]-opposite))
            projected.append((float(defect), float(action_defect)))
        assert projected[1][0] < 1e-11 and projected[1][1] < 1e-11
        assert projected[0][0] > 1e-4
        closure = paired_millers(original, q, qm)
        opposite = {tuple(g) for g in paired_millers(original, qm, q).T}
        shift = np.rint(q+qm).astype(int)
        assert all(tuple(-g-shift) in opposite for g in closure.T)
        reports.append(dict(q=q.tolist(), qminus=qm.tolist(), original_unpaired_mutation=projected[0][0],
                            paired_projection_defect=projected[1][0], independent_action_defect=projected[1][1],
                            original_support=original.shape[1], paired_support=closure.shape[1]))
    # Distinct physical derivative/radial modes may share a discrete FFT bin:
    # the physical extension must sum them, never overwrite one mode.
    fft = PhysicalRhoFFT.create((4, 4, 4), np.array([[-2, -1, 0, 1], [0]*4, [0]*4]), np.eye(3), 3.)
    support = create_support(fft, [-.5, 0, 0], [-.5, 0, 0])
    coefficients = np.arange(1, len(support.nl)+1).astype(complex)
    field = support.from_g(coefficients)
    grid = np.indices((4, 4, 4))
    direct = sum(c*np.exp(2j*np.pi*sum(g[a]*grid[a]/4 for a in range(3))) for c, g in zip(coefficients, support.miller.T))
    alias_defect = float(la.norm(field-direct.reshape(-1, order='F')))
    assert alias_defect < 1e-11
    return dict(schema='p254-physical-nyquist-support-oracle-v1', cases=reports,
                additive_alias_defect=alias_defect, physical_material_execution=False)


if __name__ == '__main__':
    print(json.dumps(oracle(), indent=2, allow_nan=False))
