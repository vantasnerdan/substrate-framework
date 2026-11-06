"""Independent differentiated-lattice and scalar-energy oracle; synthetic only.

Run explicitly by integration owner. This script has not been executed by its
construction worker. SymPy differentiates scalar pair kernels; the direct sums
use uncentered lattice translations, not production rgen or radial derivatives.
"""
from itertools import product
import json

import numpy as np
import sympy as sp
from scipy.special import erfc

from substrate_framework.qe_ionic_curvature import QEIonicCurvature
from substrate_framework.material_frozen_lattice import FrozenLattice
from substrate_framework.material_force_lattice import q_grid_to_force_lattice, force_lattice_to_q_grid


def main():
    x = sp.symbols('x y z', real=True)
    k = sp.symbols('kx ky kz', real=True)
    alpha_symbol = sp.symbols('alpha', positive=True)
    radius = sp.sqrt(sum(v*v for v in x))
    short_kernel = sp.erfc(sp.sqrt(alpha_symbol)*radius)/radius
    short_hessian = sp.lambdify((*x, alpha_symbol), sp.hessian(short_kernel, x), modules='scipy')
    phase_kernel = sp.exp(sp.I * sum(a*b for a, b in zip(k, x)))
    phase_hessian = sp.lambdify((*x, *k), sp.hessian(phase_kernel, x), modules='numpy')
    cell = np.array([[4.7, .2, .1], [0., 5.2, .3], [0., 0., 4.9]])
    positions = np.array([[.2, 4.1], [.3, 2.6], [.4, 3.8]])
    charges = np.array([1.3, .8])
    alpha, alat = .4, 4.7
    reciprocal = 2*np.pi*np.linalg.inv(cell).T
    volume = abs(np.linalg.det(cell))
    cutoff = 5/np.sqrt(alpha)

    def vectors(n):
        return reciprocal @ np.array(list(product(range(-n, n+1), repeat=3))).T

    def build(n=3, radius=cutoff):
        return QEIonicCurvature(cell_bohr=cell, positions_bohr=positions,
            valence_charges=charges, g_cart_bohr_inv=vectors(n), alat_bohr=alat,
            alpha_bohr_inv2=alpha, real_cutoff_bohr=radius, max_real_vectors=1000)

    # Literal uncentered image positions; active image list fixed for derivatives.
    bonds = []
    for a in range(2):
        for b in range(2):
            for integer in product(range(-4, 5), repeat=3):
                translation = cell @ integer
                r = translation - positions[:, a] + positions[:, b]
                r2 = r@r
                if r2 <= cutoff**2 and r2/alat**2 > 1e-10:
                    bonds.append((a, b, translation, r))

    def oracle(q):
        rp = np.zeros((2, 3, 2, 3), complex)
        rs, gp, gs = (np.zeros_like(rp) for _ in range(3))
        for a, b, translation, r in bonds:
            h = 2*charges[a]*charges[b]*np.asarray(short_hessian(*r, alpha), float)
            rp[a, :, b, :] -= np.exp(1j*q@translation)*h
            rs[a, :, a, :] += h
        for a in range(2):
            for b in range(2):
                delta = positions[:, a]-positions[:, b]
                for g in vectors(3).T:
                    for shift, dest, diagonal in ((q, gp, False), (np.zeros(3), gs, True)):
                        p = g+shift
                        p2 = p@p
                        if p2 <= 1e-8:
                            continue
                        mixed = -np.asarray(phase_hessian(*delta, *p), complex)
                        block = 8*np.pi/volume*charges[a]*charges[b]*np.exp(-p2/(4*alpha))/p2*mixed
                        if diagonal:
                            dest[a, :, a, :] -= block.real
                        else:
                            dest[a, :, b, :] += block
        return gp, gs, rp, rs

    ionic = build()
    q = reciprocal @ np.array([.23, -.17, .11])
    native = ionic.terms(q)
    expected = oracle(q)
    term_names = ('reciprocal_pair', 'reciprocal_self', 'real_pair', 'real_self')
    errors = {name:float(np.max(np.abs(getattr(native, name)-value))) for name, value in zip(term_names, expected)}
    assert max(errors.values()) < 2e-11, errors
    full = native.total
    expected_full = sum(expected)
    mutations = {f'missing_{name}':float(np.linalg.norm(full-getattr(native,name)-expected_full)) for name in term_names}
    mutations['wrong_final_sign'] = float(np.linalg.norm(-full-expected_full))
    mutations['missing_q_phases'] = float(np.linalg.norm(ionic.terms(np.zeros(3)).total-expected_full))
    mutations['conjugated_q_phases'] = float(np.linalg.norm(full.conj()-expected_full))
    assert min(mutations.values()) > 1e-8, mutations

    # Independent finite scalar energy, with fixed finite support, at Gamma.
    def energy(displacement):
        pos = positions + displacement.reshape(2,3).T
        real = 0.
        for a,b,translation,_ in bonds:
            r = np.linalg.norm(translation-pos[:,a]+pos[:,b])
            real += charges[a]*charges[b]*erfc(np.sqrt(alpha)*r)/r
        rec = 0.
        for g in vectors(3).T:
            g2=g@g
            if g2 > 1e-8:
                rho = np.sum(charges*np.exp(1j*g@pos))
                rec += 4*np.pi/volume*np.exp(-g2/(4*alpha))/g2*abs(rho)**2
        # Ewald self/background are position-independent at fixed cell.
        return real+rec

    direction = np.array([.3,-.4,.2,-.1,.5,.7])
    gamma = ionic.matrix(np.zeros(3))
    exact = float(np.real(direction @ gamma @ direction))
    e0 = energy(np.zeros(6))
    finite_difference = []
    for step in (.003, .0015, .00075):
        estimate = (energy(step*direction)-2*e0+energy(-step*direction))/step**2
        finite_difference.append(dict(step_bohr=step, value=float(estimate), error=abs(estimate-exact)))
    assert finite_difference[-1]['error'] < 2e-6, finite_difference
    assert finite_difference[-1]['error'] < finite_difference[0]['error'], finite_difference
    hermitian = np.linalg.norm(ionic.matrix(q)-ionic.matrix(q).conj().T)
    opposite = np.linalg.norm(ionic.matrix(-q)-ionic.matrix(q).conj())
    translation = np.linalg.norm(gamma.reshape(2,3,2,3).sum(axis=2))
    assert max(hermitian, opposite, translation) < 2e-11

    axis = np.array([1.,2.,-1.]); axis /= np.linalg.norm(axis)
    threshold_continuity = []
    analytic_continuity = []
    raw_head_norms = []
    for eps in (3e-3,1.5e-3,7.5e-4):
        smallq=eps*axis
        phase=charges*np.exp(1j*positions.T@smallq)
        head = (8*np.pi/volume*np.exp(-eps**2/(4*alpha))
                *np.einsum('a,b,i,j->aibj',phase,phase.conj(),axis,axis)).reshape(6,6)
        raw=ionic.matrix(smallq)-gamma
        analytic_continuity.append(float(np.linalg.norm(raw-head)))
        raw_head_norms.append(float(np.linalg.norm(raw)))
    for eps in (1e-6,5e-7,2.5e-7):
        threshold_continuity.append(float(np.linalg.norm(ionic.matrix(eps*axis)-gamma)))
    assert analytic_continuity[-1] < .4*analytic_continuity[0]
    assert threshold_continuity[-1] < .4*threshold_continuity[0]
    assert raw_head_norms[-1] > 1e-3  # no false physical Gamma continuity

    refinement = [build(n=n, radius=(5+n-3)/np.sqrt(alpha)).matrix(q) for n in (3,4,5)]
    refinement_deltas=[float(np.linalg.norm(b-a)) for a,b in zip(refinement,refinement[1:])]
    assert refinement_deltas[1] < max(refinement_deltas[0],1e-12)
    # Offsite imported residual deliberately nonzero; no onsite repair allowed.
    saved = gamma.copy(); saved[0,4] += .0123j; saved[4,0] -= .0123j
    frozen = FrozenLattice(saved, ionic)
    assert np.array_equal(frozen.matrix(np.zeros(3)),saved)
    assert abs(frozen.electronic_remainder[0,4]-.0123j) < 1e-15
    grid=frozen.complete_q_grid((3,1,1))
    lattice=q_grid_to_force_lattice(grid)
    roundtrip=float(np.max(np.abs(force_lattice_to_q_grid(lattice)-grid)))
    assert roundtrip < 2e-12
    print(json.dumps(dict(scope='independent synthetic finite-sum oracle, not material evidence',
        term_errors=errors, mutation_distances=mutations, scalar_energy_derivatives=finite_difference,
        hermiticity=float(hermitian), opposite_q_conjugacy=float(opposite), gamma_ionic_translation=float(translation),
        threshold_continuity=threshold_continuity, head_subtracted_continuity=analytic_continuity,
        unremoved_head_norms=raw_head_norms, finite_sum_refinement_deltas=refinement_deltas,
        fourier_roundtrip=roundtrip),indent=2,allow_nan=False))


if __name__ == '__main__':
    main()
