#!/usr/bin/env python3
"""Sensitive finite dense/Fourier oracle. No material calculation is performed."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace as NS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'src'))
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'

import numpy as np
from threadpoolctl import threadpool_limits, threadpool_info
from substrate_framework.qe_material_grid import PhysicalRhoFFT
from substrate_framework.qe_material_hamiltonian import MaterialHamiltonian, PWFactors, SpeciesProjectors
from substrate_framework.qe_bloch_ionic_action import BlochIonicAction
from substrate_framework.qe_bloch_qzero import BlochQZeroBridge
from substrate_framework.material_frozen_action import FrozenMaterialAction
from substrate_framework.qe_material_ionic_action import GammaIonicAction


def binding(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            digest.update(block)
    return dict(path=str(path), sha256=digest.hexdigest(), bytes=path.stat().st_size)


def source_bindings():
    registration = json.loads((HERE/'fused-action-registration.json').read_text())
    paths = [ROOT/path for path in registration['source_authority']]
    for path in paths:
        if binding(path)['sha256'] != registration['source_authority'][str(path.relative_to(ROOT))]:
            raise ValueError('preregistered legacy source changed: '+str(path))
    paths += [ROOT/'src/substrate_framework'/name for name in
              ('material_frozen_action.py', 'qe_bloch_dyads.py')]
    paths += [HERE/name for name in ('fused-action-registration.json',
              'run_fused_action_oracle.py', 'run_fused_action_material.py',
              'run_streamed_material_evolution.py')]
    return [binding(path) for path in paths]


def publish(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def relative(actual, expected):
    return float(np.linalg.norm(actual-expected)/max(1., np.linalg.norm(expected)))


class FiniteLocalTable:
    """Specified analytic finite-model local Fourier coefficients, not UPF data."""
    def __init__(self, amplitude):
        self.amplitude = amplitude

    def vloc(self, q):
        return self.amplitude*np.exp(-.21*q*q)


def fixture():
    rng = np.random.default_rng(22528061)
    shape = (5, 3, 3)
    coordinates = np.indices(shape).reshape(3, -1, order='F').T
    density_miller = coordinates.T.copy()
    for axis, size in enumerate(shape):
        density_miller[axis, density_miller[axis] > size//2] -= size
    rows = np.array([0, 1, 4, 5, 8, 12, 17, 23, 29, 34, 41])
    miller = density_miller[:, rows]
    reciprocal = np.array([[.7, .1, 0.], [0., .8, -.1], [.1, 0., .9]])
    k = np.array([.13, -.08, .21])
    fft = PhysicalRhoFFT.create(shape, miller, reciprocal, 19.)
    density_fft = PhysicalRhoFFT.create(shape, density_miller, reciprocal, 19.)
    p = (reciprocal@ miller+k[:, None]).T
    positions = np.array([[.2, -.3, .1], [.6, .15, -.2], [-.4, .25, .5]])
    reference = NS(nat=3, atom_names=('A', 'B', 'A'), positions_bohr=positions.T)
    inputs = NS(reference=reference, tables={'A': FiniteLocalTable(.6), 'B': FiniteLocalTable(-.4)}, sources=())
    species = []
    for name, atoms, nh in (('A', (1, 3), 2), ('B', (2,), 3)):
        beta = .2*(rng.normal(size=(len(rows), len(atoms)*nh))+1j*rng.normal(size=(len(rows), len(atoms)*nh)))
        active_d = rng.normal(size=(nh, nh))+1j*rng.normal(size=(nh, nh))
        species.append(SpeciesProjectors(name, atoms, np.asfortranarray(beta), active_d))
    factors = PWFactors(fft, np.sum(p*p, axis=1), tuple(species), k, None, ({},))
    saved = rng.normal(size=fft.nnr)+.1j*rng.normal(size=fft.nnr)
    h0 = MaterialHamiltonian(NS(fft=fft, inputs=inputs, vrs=np.column_stack((saved, -saved))), factors, 1)
    # Finite algebra fixture directly supplies the complete checked-factor data
    # boundary. It does not claim a saved-material/BlochBasis admission. All
    # actual local/ionic preparation and application methods remain production.
    bridge = object.__new__(BlochQZeroBridge)
    bridge.basis = NS(miller=miller)
    bridge.orbital_fft, bridge.density_fft = fft, density_fft
    ionic = object.__new__(BlochIonicAction)
    ionic.bridge, ionic.inputs, ionic.factors = bridge, inputs, factors
    ionic.p, ionic.g = p, density_fft.wavevectors.T
    ionic.atoms, ionic.vloc = {}, {}
    coupling = NS(factors=factors, ionic=ionic, ndof=9)
    frame = rng.normal(size=(len(rows), 67))+1j*rng.normal(size=(len(rows), 67))
    raw_d = rng.normal(size=(67, 67))+1j*rng.normal(size=(67, 67))
    stage = .3*(rng.normal(size=fft.nnr)+1j*rng.normal(size=fft.nnr))
    direction = rng.normal(size=(3, 3))+.4j*rng.normal(size=(3, 3))
    # Independent literal Fourier matrices: no production transform is used.
    fractions = coordinates/np.asarray(shape)
    e = np.exp(2j*np.pi*fractions@ miller)
    ed = np.exp(2j*np.pi*fractions@ density_miller)
    return h0, coupling, frame, raw_d, stage, direction, e, ed


def dense_parts(h0, coupling, stage, direction, e, ed):
    f, ionic = h0.factors, coupling.ionic
    parts = dict(kinetic=np.diag(f.kinetic_ry).astype(complex),
        saved_local=e.conj().T@(h0.background.vrs[:, 0, None]*e)/len(e),
        stage_local=e.conj().T@(stage[:, None]*e)/len(e))
    for key in ('species_nonlocal', 'ionic_bra', 'ionic_ket', 'ionic_local'):
        parts[key] = np.zeros((len(ionic.p), len(ionic.p)), complex)
    local_coefficients = np.zeros(len(ionic.g), complex)
    for s in f.species:
        nh = len(s.active_d)
        for slot, atom_id in enumerate(s.atom_indices):
            atom = atom_id-1
            b = s.beta[:, slot*nh:(slot+1)*nh]
            n = b@s.active_d@b.conj().T
            parts['species_nonlocal'] += n
            # Differentiate both exponential phase factors independently.
            parts['ionic_bra'] += n*(1j*(ionic.p@direction[atom]))[None, :]
            parts['ionic_ket'] += (-1j*(ionic.p@direction[atom]))[:, None]*n
            table = ionic.inputs.tables[s.name]
            local_coefficients += (-1j*(ionic.g@direction[atom])
                *table.vloc(np.linalg.norm(ionic.g, axis=1))
                *np.exp(-1j*ionic.g@ionic.inputs.reference.positions_bohr[:, atom]))
    local = ed@local_coefficients
    parts['ionic_local'] = e.conj().T@(local[:, None]*e)/len(e)
    return parts


def run(output):
    sources = source_bindings()
    h0, coupling, frame, raw_d, field, u, e, ed = fixture()
    errors, sensitivity, retained = {}, {}, dict(frame=frame, raw_D=raw_d, stage_field=field,
                                                direction=u, fourier=e, density_fourier=ed)
    def compare(name, actual, expected):
        errors[name] = relative(actual, expected)
        retained[name+'_actual'] = np.asarray(actual)
        retained[name+'_expected'] = np.asarray(expected)
        retained[name+'_difference'] = np.asarray(actual-expected)
    baseline = dense_parts(h0, coupling, np.zeros_like(field), np.zeros_like(u), e, ed)
    dense_h0 = sum(baseline.values())
    compare('h0', FrozenMaterialAction(h0, fft_block_size=3, projector_block_size=7).apply(frame), dense_h0@frame)
    compare('h0_adjoint', FrozenMaterialAction(h0).adjoint_apply(frame), dense_h0.conj().T@frame)
    for label, direction in (('real', u.real), ('complex', u)):
        parts = dense_parts(h0, coupling, field, direction, e, ed)
        dense = sum(parts.values())
        frozen = FrozenMaterialAction(h0, coupling, potential=field, direction=direction.ravel(),
                                      fft_block_size=3, projector_block_size=7)
        actual = frozen.apply(frame)
        compare(label, actual, dense@frame)
        compare(label+'_adjoint', frozen.adjoint_apply(frame), dense.conj().T@frame)
        compare(label+'_single', frozen.apply(frame[:, 0]), dense@frame[:, 0])
        compare(label+'_panel64', FrozenMaterialAction(h0, coupling, potential=field,
                 direction=direction.ravel()).apply(frame), dense@frame)
        expected_trace = np.trace((frame@raw_d@frame.conj().T)@dense)
        compare(label+'_raw_D_trace', np.trace(raw_d@(frame.conj().T@actual)), expected_trace)
        if label == 'complex':
            for name, piece in parts.items():
                sensitivity['omit_'+name] = relative((dense-piece)@frame, dense@frame)
                retained['dense_'+name] = piece
            for name in ('ionic_bra', 'ionic_ket'):
                sensitivity['reverse_'+name] = relative((dense-2*parts[name])@frame, dense@frame)
            conjugated = sum(dense_parts(h0, coupling, field, direction.conj(), e, ed).values())
            sensitivity['conjugate_direction'] = relative(conjugated@frame, dense@frame)
            sensitivity['hermitize_operator'] = relative((dense+dense.conj().T)@frame/2, dense@frame)
            hermitian_n = (parts['species_nonlocal']+parts['species_nonlocal'].conj().T)/2
            sensitivity['hermitize_active_D'] = relative((dense-parts['species_nonlocal']+hermitian_n)@frame, dense@frame)
            sensitivity['hermitize_raw_D'] = relative(np.trace((raw_d+raw_d.conj().T)/2@(frame.conj().T@actual)), expected_trace)
            # Complex linearity is in the ionic direction, not its conjugate.
            a, b = .3+.7j, -.2+.4j
            def delta(v):
                return FrozenMaterialAction(h0, coupling, direction=v.ravel()).apply(frame)-dense_h0@frame
            compare('direction_linearity', delta(a*u.real+b*u.imag), a*delta(u.real)+b*delta(u.imag))
    zero = FrozenMaterialAction(h0, coupling, direction=np.zeros(9))
    if zero.terms or zero.direction is None:
        raise AssertionError('exact zero direction must retain identity but skip ionic terms')
    compare('zero_direction', zero.apply(frame), dense_h0@frame)
    tiny = FrozenMaterialAction(h0, coupling, direction=(1e-100*u).ravel())
    if len(tiny.terms) != 3 or not any(np.any(term[0] != 0) for term in tiny.terms):
        raise AssertionError('nonzero tiny direction was pruned')
    # The cancellation case also exercises the exactly-zero combined-field path.
    no_local = FrozenMaterialAction(h0, potential=-h0.background.vrs[:, 0])
    if no_local.has_local:
        raise AssertionError('exactly zero local field not recognized')
    compare('zero_combined_local', no_local.apply(frame), (dense_h0-baseline['saved_local'])@frame)
    # Exercise the existing Gamma API's 3D field as well as Bloch's F-order
    # vector, with a genuinely zero-k finite operator and both derivative legs.
    f, inputs = h0.factors, coupling.ionic.inputs
    gamma = GammaIonicAction(f.fft.physical, f.fft.miller,
        coupling.ionic.bridge.density_fft.miller, f.fft.reciprocal_bohr_inv,
        inputs.reference.positions_bohr, inputs.reference.atom_names, inputs.tables)
    for species in f.species:
        nh = len(species.active_d)
        for index, atom in enumerate(species.atom_indices):
            beta = species.beta[:, index*nh:(index+1)*nh]
            gamma.atoms[atom-1] = (beta@species.active_d, beta.conj().T)
    gf = PWFactors(f.fft, np.sum(gamma.p*gamma.p, axis=1), f.species,
                   np.zeros(3), None, f.sources)
    gh = MaterialHamiltonian(h0.background, gf, 1)
    gc = NS(factors=gf, ionic=gamma, ndof=9)
    # Only the independent dense reference needs this common input view.
    gamma.inputs = inputs
    gamma_dense = sum(dense_parts(gh, gc, field, u, e, ed).values())
    gamma_action = FrozenMaterialAction(gh, gc, potential=field, direction=u.ravel(),
                                        fft_block_size=3, projector_block_size=7)
    compare('gamma_3d_field', gamma_action.apply(frame), gamma_dense@frame)
    compare('gamma_adjoint', gamma_action.adjoint_apply(frame), gamma_dense.conj().T@frame)
    passed = max(errors.values()) <= 2e-11 and min(sensitivity.values()) > 1e-8
    np.savez(output/'finite-arrays.npz', **retained)
    if source_bindings() != sources:
        raise ValueError('source changed during oracle')
    result = dict(schema='P254-fused-action-oracle-v1', status='passed' if passed else 'failed',
        scope='finite independent dense/Fourier algebra only; no saved material or trajectory',
        errors=errors, mutation_sensitivity=sensitivity, tolerance=2e-11,
        mutation_minimum=1e-8, source_bindings=sources, arrays=binding(output/'finite-arrays.npz'),
        native_threadpools=threadpool_info(), cpu_affinity=sorted(os.sched_getaffinity(0)),
        limitations=['Finite fixture directly supplies operator factors, not a saved-material admission',
                    'No runtime or full-material claim; material consumer must pass separately'])
    publish(output/'result.json', result)
    if not passed:
        raise AssertionError('fused action oracle failed; complete receipts retained')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    with threadpool_limits(limits=1):
        run(args.output)


if __name__ == '__main__':
    main()
