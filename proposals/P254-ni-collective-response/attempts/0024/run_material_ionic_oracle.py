"""Independent fixed-PW displaced NC potential/action checks, no Native/SCF.

QE 7.5 compute_dvloc and dvqpsi_us_only specify the analytic derivative;
this comparator instead places the FULL potential at displaced real tau,
then differences its actions at epsilon=2e-4 and1e-4 bohr. Species radial
and centered-projector factors are the independently Native-compared NC
sources; no analytic derivative helper is used by DisplacedPotential.
Finite differences are fixed basis/no relaxation, NOT a physical ensemble.
Usage: python run_material_ionic_oracle.py --output JSON
"""
import time
BOOT = time.monotonic()
import argparse
import json
from pathlib import Path
import sys
import traceback
import numpy as np
from material_pilot_inputs import HERE, bound_json, load_material, DIRECTIONS
from material_pilot_io import binding, implementation, save_json
from substrate_framework.qe_nc_projectors import centered_projectors
from substrate_framework.qe_material_ionic_action import saved_band_overlap

REGISTRATION = {'path': str(HERE/'implementation-registration.json'), 'bytes': '10830',
                'sha256': 'c316f6eb60d5d5eb676edd4fdfecc9f30cc983910e37ac1123a8fba15046e5cc'}


class DisplacedPotential:
    """Independent finite potential. No GammaIonicAction methods/fields used.

    Fourier coefficients V_G multiply exp(-i G.tau). Projector potential
    B(tau) D B(tau)^H is placed afresh at each displaced real tau. All37
    background atoms are included; unchanged contributions cached by spin
    probe to avoid repeated work. No derivative expression enters this class.
    """
    def __init__(self, ref, wave, density, tables):
        self.shape = ref.fft_grid
        self.n = int(np.prod(self.shape))
        self.p = (wave.reciprocal_bohr_inv @ wave.miller).T
        self.g = (density.reciprocal_bohr_inv @ density.miller).T
        self.pidx = tuple(wave.miller[i] % self.shape[i] for i in range(3))
        self.gidx = tuple(density.miller[i] % self.shape[i] for i in range(3))
        self.positions = ref.positions_bohr
        self.species = ref.atom_names
        self.v = {s: t.vloc(np.linalg.norm(self.g, axis=1)) for s, t in tables.items()}
        self.b0 = {s: centered_projectors(t, self.p) for s, t in tables.items()}
        self.d = {s: t.upf.expanded_d for s, t in tables.items()}
        self.miller_lookup = {tuple(int(x) for x in density.miller[:, i]): i
                              for i in range(density.ngm)}

    def atom(self, atom, position):
        symbol = self.species[atom]
        coefficients = self.v[symbol] * np.exp(-1j*(self.g @ position))
        b = self.b0[symbol] * np.exp(-1j*(self.p @ position))[:, None]
        return coefficients, b

    def apply(self, coefficients, b, d, psi):
        # Deliberately uses unnormalized forward/inverse pair, rather than
        # production's ortho transforms; physical N factors are explicit.
        potential_grid = np.zeros(self.shape, dtype=np.complex128)
        potential_grid[self.gidx] = coefficients
        potential_grid = np.fft.ifftn(potential_grid)*self.n
        output = np.empty_like(psi)
        for j in range(psi.shape[1]):
            grid = np.zeros(self.shape, dtype=np.complex128)
            grid[self.pidx] = psi[:, j]
            physical_wave = np.fft.ifftn(grid)*self.n
            output[:, j] = (np.fft.fftn(potential_grid*physical_wave)/self.n)[self.pidx]
        if b is not None:
            output += b @ (d @ (b.conj().T @ psi))
        return output

    def background(self, psi):
        coefficients = np.zeros(len(self.g), dtype=np.complex128)
        nonlocal_action = np.zeros_like(psi)
        for atom in range(len(self.species)):
            c, b = self.atom(atom, self.positions[:, atom])
            coefficients += c
            nonlocal_action += b @ (self.d[self.species[atom]] @ (b.conj().T @ psi))
        return self.apply(coefficients, None, None, psi)+nonlocal_action

    def displaced(self, atom, displacement, psi, background):
        c0, b0 = self.atom(atom, self.positions[:, atom])
        c, b = self.atom(atom, self.positions[:, atom]+displacement)
        d = self.d[self.species[atom]]
        # Difference of complete finite potentials, not analytic differentiation.
        return background + self.apply(c-c0, b, d, psi) - b0 @ (d @ (b0.conj().T @ psi))

    def finite_submatrix(self, atom, displacement, miller, rows):
        c, b = self.atom(atom, self.positions[:, atom]+displacement)
        result = b[rows] @ self.d[self.species[atom]] @ b[rows].conj().T
        for i, row in enumerate(rows):
            for j, col in enumerate(rows):
                index = self.miller_lookup.get(tuple(int(x) for x in miller[:, row]-miller[:, col]))
                if index is not None:
                    result[i, j] += c[index]
        return result


def check(receipt, name, actual, expected, rtol=1e-6, atol=1e-8, mutation=False):
    actual, expected = np.asarray(actual), np.asarray(expected)
    error = float(np.max(np.abs(actual-expected), initial=0.0))
    scale = float(np.max(np.abs(expected), initial=0.0))
    allowance = atol+rtol*scale
    finite = bool(np.isfinite(actual).all() and np.isfinite(expected).all())
    passed = finite and (error > allowance if mutation else error <= allowance)
    receipt['findings'].append({'name': name, 'max_absolute_error': error if np.isfinite(error) else None,
                               'reference_scale': scale if np.isfinite(scale) else None,
                               'allowance': allowance if np.isfinite(allowance) else None,
                               'rtol': rtol, 'atol_ry_per_bohr': atol,
                               'mutation_exposed': passed if mutation else None, 'pass': passed})


def run(receipt):
    registration, reg = bound_json(REGISTRATION)
    inventory, access = bound_json(registration['access'])
    receipt['registration'] = reg
    receipt['source_access'] = access
    receipt['implementation_sources'] = implementation()
    receipt['actual_inputs'] = inventory['material_inputs']
    ref, waves, density, action, directions = load_material(inventory)
    independent = DisplacedPotential(ref, waves[0], density, action.tables)
    receipt['selection'] = {'spins': [1, 2], 'saved_bands_one_based': [1, 346],
                            'full_pw_rows': waves[0].igwx, 'background_atoms': ref.nat,
                            'directions': DIRECTIONS, 'epsilon_bohr': [2e-4, 1e-4]}
    for spin, wave in enumerate(waves, 1):
        psi = wave.coefficients[:, [0, 345]]
        background = independent.background(psi)
        adjoint = psi.conj().T.copy()
        for spec, u in zip(DIRECTIONS, directions):
            atom = spec['atom_one_based']-1
            label = 'spin'+str(spin)+'/'+spec['symbol']+str(atom+1)
            operation = action.prepare(u)
            local, bra, ket = operation.components(psi)
            analytic = local+bra+ket
            for epsilon in (2e-4, 1e-4):
                plus = independent.displaced(atom, epsilon*u[atom].real, psi, background)
                minus = independent.displaced(atom, -epsilon*u[atom].real, psi, background)
                derivative = (plus-minus)/(2*epsilon)
                # Both errors retained; finer step is the registered acceptance.
                check(receipt, label+'/displaced/'+str(epsilon), analytic, derivative,
                      rtol=1e-6 if epsilon == 1e-4 else 4e-6,
                      atol=1e-8 if epsilon == 1e-4 else 4e-8)
            for name, mutated in (('drop-bra', local+ket), ('wrong-bra-sign', local-bra+ket),
                                  ('drop-ket', local+bra), ('wrong-ket-sign', local+bra-ket),
                                  ('FFT-normalization', local/action.nnr+bra+ket)):
                check(receipt, label+'/'+name, mutated, derivative, mutation=True)
            wrong_phase_actions = []
            for sign in (1, -1):
                shifted = ref.positions_bohr[:, atom] + np.array([0.137, -0.083, 0.051])
                coefficients, projectors = independent.atom(atom, shifted+sign*1e-4*u[atom].real)
                wrong_phase_actions.append(independent.apply(coefficients, projectors,
                    independent.d[spec['symbol']], psi))
            check(receipt, label+'/erroneous-atom-phase',
                  (wrong_phase_actions[0]-wrong_phase_actions[1])/2e-4, derivative, mutation=True)
            complex_factor = 1+0.37j
            complex_action = action.prepare(complex_factor*u).apply(psi)
            check(receipt, label+'/complex-linear', complex_action, complex_factor*analytic,
                  rtol=2e-11, atol=2e-12)
            check(receipt, label+'/conjugated-direction', complex_factor.conjugate()*analytic,
                  complex_action, mutation=True)
            overlap = saved_band_overlap(adjoint, analytic)
            explicit = np.asarray([[np.vdot(psi[:, i], analytic[:, j]) for j in range(2)] for i in range(2)])
            check(receipt, label+'/explicit-projection', overlap, explicit, rtol=2e-11, atol=2e-12)
            check(receipt, label+'/wrong-projection', psi.T @ analytic, explicit, mutation=True)
            check(receipt, label+'/Gamma-adjoint', overlap, overlap.conj().T, rtol=2e-10, atol=2e-11)
            phases = np.exp(1j*np.array([0.37, -0.83]))
            check(receipt, label+'/phase-covariance', operation.apply(psi*phases),
                  analytic*phases, rtol=2e-11, atol=2e-12)
            unitary = np.array([[1, 1j], [1j, 1]], dtype=np.complex128)/np.sqrt(2)
            rotated = operation.apply(psi @ unitary)
            check(receipt, label+'/basis-covariance', rotated, analytic @ unitary,
                  rtol=2e-11, atol=2e-12)
            check(receipt, label+'/bilinear-basis-covariance',
                  (psi @ unitary).conj().T @ rotated, unitary.conj().T @ explicit @ unitary,
                  rtol=2e-11, atol=2e-12)
            # Direct Fourier convolution on independent PW basis probes checks
            # FFT phase, row indexing and normalization without any FFT reference.
            rows = np.asarray([0, len(action.p)//7, len(action.p)//3, len(action.p)-1])
            basis = np.zeros((len(action.p), len(rows)), dtype=np.complex128)
            basis[rows, np.arange(len(rows))] = 1
            epsilon = 1e-4
            direct = (independent.finite_submatrix(atom, epsilon*u[atom].real, wave.miller, rows)-
                      independent.finite_submatrix(atom, -epsilon*u[atom].real, wave.miller, rows))/(2*epsilon)
            check(receipt, label+'/direct-PW-potential-basis', operation.apply(basis)[rows], direct)
    receipt['status'] = 'pass' if all(x['pass'] for x in receipt['findings']) else 'failed'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipt = {'status': 'failed', 'findings': [],
               'scope': 'new fixed-PW displaced bare ionic derivative; no Native/SCF/H0/ensemble',
               'budget_authority': 'remaining original reader+NC cumulative allowance; Root owns launch'}
    try:
        run(receipt)
    except BaseException as error:
        receipt['failure'] = {'type': type(error).__name__, 'message': str(error),
                              'traceback': traceback.format_exc()}
    finally:
        receipt['elapsed_seconds'] = time.monotonic()-BOOT
        args.output.parent.mkdir(parents=True, exist_ok=True)
        save_json(args.output, receipt, rotate=True)
        print(json.dumps(receipt, allow_nan=False), flush=True)
    return 0 if receipt['status'] == 'pass' else 1


if __name__ == '__main__':
    sys.exit(main())
