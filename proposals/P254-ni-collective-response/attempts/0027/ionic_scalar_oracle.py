"""Independent small Ewald comparisons: literal pair energy, not production sums.

Blocks: finite energy/constants -> independent energy differentiation ->
QE-force-style comparator -> symmetry/units -> truncation refinement/mutations.
Comparators do not import production kernels, intermediates or force formulas.
Only check_scalar imports the public evaluator to obtain the candidate result.
"""
from itertools import product
from math import pi, sqrt

import numpy as np
from scipy.special import erfc


# Deliberately generous fixed cube, unlike production inverse-axis enumeration.
INTEGERS = np.array(list(product(range(-8, 9), repeat=3)), dtype=float)


def direct_energy(cell, positions, charges, alpha, rc, gc, *, part='all'):
    """Ordered pairs and cos(G.(Ri-Rj)); no structure factor or force input."""
    volume = abs(np.linalg.det(cell))
    lattice = INTEGERS @ cell.T
    reciprocal = 2*pi*np.linalg.solve(cell.T, INTEGERS.T).T
    g2 = np.sum(reciprocal**2, axis=1)
    gv = reciprocal[(g2 > 0) & (g2 <= gc*gc)]
    q2 = np.sum(gv**2, axis=1)
    # Explicit finite coverage assertions belong to the comparator, not a fit.
    max_delta = max(np.linalg.norm(a-b) for a in positions for b in positions)
    if 9*np.linalg.svd(cell, compute_uv=False).min()-max_delta <= rc:
        raise ValueError('independent direct real cube does not cover the requested sphere')
    if 9*np.linalg.svd(2*pi*np.linalg.inv(cell).T, compute_uv=False).min() <= gc:
        raise ValueError('independent direct reciprocal cube is insufficient')
    real, reciprocal_energy = 0., 0.
    for i, ri in enumerate(positions):
        for j, rj in enumerate(positions):
            distances = np.linalg.norm(lattice + rj-ri, axis=1)
            if i != j and np.any(distances == 0):
                raise ValueError('singular distinct ions')
            lengths = distances[(distances > 0) & (distances <= rc)]
            real += charges[i]*charges[j]*np.sum(erfc(sqrt(alpha)*lengths)/lengths)
            reciprocal_energy += (4*pi/volume*charges[i]*charges[j]
                                  * np.sum(np.exp(-q2/(4*alpha))/q2*np.cos(gv@(ri-rj))))
    self_energy = -2*sqrt(alpha/pi)*sum(z*z for z in charges)
    background = -pi*sum(charges)**2/(alpha*volume)
    return {'real': real, 'reciprocal': reciprocal_energy,
            'self': self_energy, 'background': background,
            'all': real+reciprocal_energy+self_energy+background}[part]


def energy_force(cell, positions, charges, alpha, rc, gc, h, *, part='all'):
    """Five-point Cartesian energy derivative; production forces are unseen."""
    force = np.zeros_like(positions)
    for atom in range(len(charges)):
        for axis in range(3):
            energies = []
            for step in (-2, -1, 1, 2):
                displaced = positions.copy()
                displaced[atom, axis] += step*h
                energies.append(direct_energy(cell, displaced, charges, alpha, rc, gc, part=part))
            force[atom, axis] = -(energies[0]-8*energies[1]+8*energies[2]-energies[3])/(12*h)
    return force


def force_ew_style(cell, positions, charges, alpha, rc, gc):
    """Independent all-ordered-atom force route with QE r=R-dtau sign.

    Trigonometric sine/cosine charge arrays replace production complex phases;
    no pair antisymmetrization or final total-force subtraction is applied.
    Numerical alpha/cutoffs are matched explicitly, not Native auto-selected.
    """
    result = np.zeros_like(positions)
    vectors = INTEGERS @ cell.T
    reciprocal = 2*pi*np.linalg.solve(cell.T, INTEGERS.T).T
    lengths2 = np.sum(reciprocal**2, axis=1)
    keep = (lengths2 > 0) & (lengths2 <= gc*gc)
    gs, gs2 = reciprocal[keep], lengths2[keep]
    sine = np.sin(gs@positions.T) @ charges
    cosine = np.cos(gs@positions.T) @ charges
    weight = np.exp(-gs2/(4*alpha))/gs2
    for i, position in enumerate(positions):
        phase = gs@position
        sumnb = np.cos(phase)*sine - np.sin(phase)*cosine
        result[i] = -8*pi/abs(np.linalg.det(cell))*charges[i]*((weight*sumnb)@gs)
        for j in range(len(charges)):
            if i == j:
                continue
            r = vectors - (position-positions[j])
            length = np.linalg.norm(r, axis=1)
            keep = (length > 0) & (length <= rc)
            r, length = r[keep], length[keep]
            factor = (2*charges[i]*charges[j]/length**2
                      * (erfc(sqrt(alpha)*length)/length
                         + sqrt(4*alpha/pi)*np.exp(-alpha*length**2)))
            result[i] -= factor@r
    return result


def check_scalar(check, mutation):
    from substrate_framework.qe_ionic_ewald import EwaldSumContract, periodic_ionic_ewald
    cell = np.array([[5.1,.7,-.2],[.1,4.8,.5],[.3,-.4,5.7]])
    positions = np.array([[.37,.61,.83],[2.14,1.28,3.07],[3.21,3.63,1.52]])
    charges = np.array([1.,3.,2.])  # explicit synthetic valence, not Ni/H
    alpha, rc, gc = .45, 7.13, 4.31

    def evaluate(c=cell, r=positions, z=charges, a=alpha, radius=rc, reciprocal=gc):
        return periodic_ionic_ewald(c, r, z, contract=EwaldSumContract(a, radius, reciprocal),
                                   coordinate_identity='synthetic asymmetric ions A,B,C Cartesian xyz',
                                   source_authority='literal synthetic valence/geometry, no material output', sources=())

    result = evaluate()
    check('finite_energy', result.energy_ry, direct_energy(cell,positions,charges,alpha,rc,gc))
    for part in ('real','reciprocal','self','background'):
        expected = direct_energy(cell,positions,charges,alpha,rc,gc,part=part)
        check('energy_'+part, result.energy_terms_ry[part], expected)
        mutation('omit_energy_'+part, result.energy_ry-result.energy_terms_ry[part], result.energy_ry)
    derivatives = []
    for h in (2e-4, 1e-4):
        # A finite difference must stay on the same local finite-shell branch.
        if result.sum_metadata['real_cutoff_gap_bohr'] <= 4*h:
            raise AssertionError('fixture lies too near a real cutoff crossing')
        derivative = energy_force(cell,positions,charges,alpha,rc,gc,h)
        check('direct_energy_derivative_'+str(h),result.force_ry_bohr,derivative,atol=3e-8)
        derivatives.append(derivative)
    check('derivative_refinement', derivatives[0], derivatives[1], atol=3e-8)
    check('qe_force_style',result.force_ry_bohr,force_ew_style(cell,positions,charges,alpha,rc,gc))
    for part, actual in (('real',result.real_force_ry_bohr),('reciprocal',result.reciprocal_force_ry_bohr)):
        check('derivative_'+part,actual,energy_force(cell,positions,charges,alpha,rc,gc,1e-4,part=part),atol=3e-8)
        mutation('omit_force_'+part,result.force_ry_bohr-actual,derivatives[-1],atol=3e-8)
    mutation('wrong_sign',-result.force_ry_bohr,derivatives[-1])
    mutation('hartree_not_ry',result.force_ry_bohr/2,derivatives[-1])
    mutation('wrong_reciprocal_phase',result.real_force_ry_bohr-result.reciprocal_force_ry_bohr,derivatives[-1])
    mutation('nuclear_not_valence',evaluate(z=np.array([2.,7.,5.])).force_ry_bohr,derivatives[-1])
    bad = result.force_ry_bohr.copy()
    bad[0,0] += .17
    mutation('centering_cannot_repair_derivative',bad-bad.mean(axis=0),derivatives[-1])
    translated = evaluate(r=positions+np.array([.71,-.39,1.12]))
    check('common_translation_energy',translated.energy_ry,result.energy_ry)
    check('common_translation_force',translated.force_ry_bohr,result.force_ry_bohr)
    check('unprocessed_net_force',result.force_ry_bohr.sum(axis=0),np.zeros(3))
    permutation = np.array([2,0,1])
    check('atom_order',evaluate(r=positions[permutation],z=charges[permutation]).force_ry_bohr,result.force_ry_bohr[permutation])
    mutation('atom_order_mismatch',evaluate(r=positions[permutation]).force_ry_bohr,result.force_ry_bohr[permutation])
    shifted = positions.copy()
    shifted[1] += cell[:,1]
    check('single_atom_lattice_periodicity',evaluate(r=shifted).force_ry_bohr,result.force_ry_bohr)
    scale = 1.7
    scaled = evaluate(c=cell*scale,r=positions*scale,a=alpha/scale**2,radius=rc*scale,reciprocal=gc/scale)
    check('length_energy_units',scaled.energy_ry,result.energy_ry/scale)
    check('length_force_units',scaled.force_ry_bohr,result.force_ry_bohr/scale**2)
    mutation('cartesian_cell_transpose',evaluate(c=cell.T).force_ry_bohr,result.force_ry_bohr)
    mutation('wrong_length_force_conversion',scaled.force_ry_bohr*scale,result.force_ry_bohr)
    check('gradient_atom_major_sign',result.gradient_ry_bohr,-derivatives[-1].ravel(),atol=3e-8)
    refined = evaluate(radius=9.17, reciprocal=6.13)
    finer = evaluate(radius=11.19, reciprocal=7.31)
    check('cutoff_refinement_energy',refined.energy_ry,result.energy_ry,atol=2e-4)
    check('cutoff_refinement_force',refined.force_ry_bohr,result.force_ry_bohr,atol=2e-4)
    check('fine_cutoff_energy',finer.energy_ry,refined.energy_ry,atol=2e-8)
    check('fine_cutoff_force',finer.force_ry_bohr,refined.force_ry_bohr,atol=2e-8)
    mutation('grossly_insufficient_reciprocal_cutoff',evaluate(reciprocal=1.1).force_ry_bohr,finer.force_ry_bohr)
    mutation('grossly_insufficient_real_cutoff',evaluate(radius=1.1).force_ry_bohr,finer.force_ry_bohr)
