"""Finite SymPy/numerical oracle; no saved material or electronic work.

Run only by the parent after slices settle. Symbolic coefficients come from
an independently differentiated minimally coupled scalar. Numerical regression
also calls the actual existing mixed_components kinetic component. This checks
algebra, normalization and ordering, not material accuracy or new trajectories.
"""
import json
from hashlib import sha256
from types import SimpleNamespace

import numpy as np
import sympy as sp

import saved_canonical_response as response
from saved_canonical_response_readout import geometry
from mixed_excitation_energy import mixed_components
from substrate_framework.material_modes import AMU_RY


def close(actual, expected):
    np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)


def different(actual, expected):
    if np.allclose(actual, expected, rtol=1e-9, atol=1e-10):
        raise AssertionError('selected wrong convention/input escaped the exposing fixture')


def algebra_block():
    lam, eta = sp.symbols('lambda eta')
    # General complex symbolic directions, not a real/conjugate special case.
    ps = sp.symbols('p0 pplus pminus px')
    aa = sp.symbols('a0 aplus aminus ax')
    z, m, cells = sp.symbols('Z M N', nonzero=True)
    generic_k = [x-z*y for x, y in zip(ps, aa)]
    generic_scalar = cells*(generic_k[0]+lam*generic_k[1]+eta*generic_k[2]
                            +lam*eta*generic_k[3])**2/(2*m)
    generic_derivative = sp.diff(generic_scalar, lam, eta).subs({lam: 0, eta: 0})
    assert sp.expand(generic_derivative-cells*(generic_k[0]*generic_k[3]
                                              +generic_k[1]*generic_k[2])/m) == 0
    n = 7
    mass = (2, 2, 2, 11, 11, 11)
    charges = (3, 1)
    # Rectangular, non-axis-aligned E exposes mode/xyz transposition. Directions
    # have independent imaginary parts; none can be reconstructed by conjugacy.
    e = sp.Matrix([[1, 0, 0], [0, sp.Rational(3, 5), sp.Rational(4, 5)]])
    p = [[sp.Rational(2+3*s+j, 5)+sp.I*sp.Rational(1+s-2*j, 7)
          for j in range(6)] for s in range(4)]
    a = [[sp.Rational(1+2*s+j, 11)+sp.I*sp.Rational(2-s+j, 13)
          for j in range(2)] for s in range(4)]
    k = [[p[s][j]-charges[j//3]*sum(a[s][mode]*e[mode, j % 3] for mode in range(2))
          for j in range(6)] for s in range(4)]
    coefficients = (1, lam, eta, lam*eta)
    scalar_atoms, derivatives = [], []
    for atom in range(2):
        scalar = sp.Rational(n, 2)*sum(
            sum(coefficients[s]*k[s][j] for s in range(4))**2/mass[j]
            for j in range(3*atom, 3*atom+3))
        derivative = sp.diff(scalar, lam, eta).subs({lam: 0, eta: 0})
        product_rule = n*sum((k[0][j]*k[3][j]+k[1][j]*k[2][j])/mass[j]
                            for j in range(3*atom, 3*atom+3))
        assert sp.expand(derivative-product_rule) == 0
        scalar_atoms.append(scalar)
        derivatives.append(complex(derivative.evalf()))
    assert sp.expand(sp.diff(sum(scalar_atoms), lam, eta).subs({lam: 0, eta: 0})
                     -sum(sp.diff(x, lam, eta).subs({lam: 0, eta: 0}) for x in scalar_atoms)) == 0
    number = lambda value: np.array(value, dtype=np.complex128)
    stages = {name: dict(u=number([j+2*s+sp.I*(s-j) for j in range(6)]),
                         p=number(p[s]), coordinate=number(a[s]),
                         momentum=number([s+1j, 2*s-3j]))
              for s, name in enumerate(response.STAGES)}
    kwargs = dict(masses=np.array(mass, float), charges=np.array(charges, float),
                  polarizations=np.array(e, float))
    mean = dict(u_bohr=stages['mean']['u'], p_canonical=stages['mean']['p'],
                field_coordinate=stages['mean']['coordinate'], field_momentum=stages['mean']['momentum'])
    small = {c: SimpleNamespace(**stages[c]) for c in response.DIRECTIONS}
    parsed = response.canonical_stages(mean, small, natoms=2, nfield=2)
    try:
        response.canonical_stages(mean, {'plus': small['plus'], 'mixed': small['mixed']}, natoms=2, nfield=2)
    except ValueError:
        pass
    else:
        raise AssertionError('missing independent minus was accepted')
    mechanical = response.mechanical_stages(parsed, **kwargs)
    observed = response.kinetic_atomwise(mechanical, masses=kwargs['masses'], ncells=n)
    close(observed, derivatives)
    close(mechanical['mixed']['shift'], number(a[3]) @ np.array(e, float))
    close(mechanical['minus']['velocity'].reshape(-1), number(k[2])/mass)
    different(mechanical['minus']['u'], mechanical['plus']['u'].conj())
    # Exercise the real existing11-component implementation, not a copied total.
    physical = {name: SimpleNamespace(u=s['u'].reshape(-1), p=s['p'].reshape(-1),
        shift=s['shift'], momentum=s['momentum'] @ kwargs['polarizations'],
        density=np.zeros(2), potential=np.zeros(2), gradient=np.zeros(6),
        field_gradient_electronic=np.zeros(3)) for name, s in mechanical.items()}
    reference = SimpleNamespace(ncells=n, masses=kwargs['masses'], nbar=np.zeros(2),
        curvature=np.zeros((6, 6)), ionic_force=np.zeros(6), support_force=np.zeros(6))
    components = mixed_components(reference=reference,
        coefficients=SimpleNamespace(delta_curvature=np.zeros((6, 6)),
            contact=np.zeros((3, 3)), mixed=np.zeros((6, 3))),
        field=SimpleNamespace(mass=5.), charges=kwargs['charges'], stages=physical,
        one_body={'mixed': 0j}, grid_measure=1., curvature_plus=np.zeros((6, 6)),
        curvature_minus=np.zeros((6, 6)), gamma=True)
    close(np.sum(observed), components['ionic_kinetic_minimal'])
    xyz = {name: {key: value.copy() for key, value in stage.items()} for name, stage in stages.items()}
    mutate = lambda changed: response.kinetic_atomwise(
        response.mechanical_stages(changed, **kwargs), masses=kwargs['masses'], ncells=n)
    xyz['mixed']['p'][:] = 0
    xyz['mixed']['coordinate'][:] = 0
    different(mutate(xyz), observed)  # omitted mean-times-mixed
    for scaling in ('field', 'matter'):
        altered = {name: {key: value.copy() for key, value in stage.items()} for name, stage in stages.items()}
        for name in response.DIRECTIONS:
            if scaling == 'field':
                altered[name]['coordinate'] /= np.sqrt(n)
            else:
                altered[name]['p'] /= np.sqrt(n)
        different(mutate(altered), observed)
    altered = {name: {key: value.copy() for key, value in stage.items()} for name, stage in stages.items()}
    altered['minus'] = {key: value.conj() for key, value in altered['plus'].items()}
    different(mutate(altered), observed)
    different(observed/2, observed)
    kk = [mechanical[name]['kinetic_momentum'] for name in response.STAGES]
    wrong_dagger = n*np.sum((kk[0].conj()*kk[3]+kk[1].conj()*kk[2])/np.array(mass).reshape(2, 3), axis=1)
    different(wrong_dagger, observed)
    wrong_mass = response.kinetic_atomwise(mechanical, masses=np.full(6, mass[0]), ncells=n)
    different(wrong_mass, observed)
    # Source positions are3xnatoms, canonical vectors are natoms x3.
    positions = np.array([[2., 17.], [3., 19.], [5., 23.]])
    table = response.atom_table(mechanical, observed, atom_ids=[7, 41], species=['Ni', 'H'], positions=positions)
    close(table['atoms'][1]['position_bohr_xyz'], [17., 19., 23.])
    close(table['atoms'][1]['response']['plus']['p_bohr_inverse_xyz'], number(p[1][3:6]))
    close(table['totals']['Ni']['kinetic_mixed_global_ry']+table['totals']['H']['kinetic_mixed_global_ry'], np.sum(observed))
    different(positions.reshape(2, 3), positions.T)
    close(table['totals']['all']['response']['minus']['displacement_weight_bohr2'], np.sum(np.abs(stages['minus']['u'])**2))
    zero = response.displacement_weights(np.zeros((2, 3), complex))
    assert zero['atom_fraction'] is None and zero['total_bohr2'] == 0
    return dict(symbolic='exact differentiated product rule and additive grouping',
                numerical='independent complex coefficients, actual kinetic scalar, unequal masses, rectangular E, original ordering',
                mutations_exposed=['mean-times-mixed omission', 'global field sqrtN division', 'matter rescaling',
                    'conjugate-minus substitution', 'extra half', 'dagger', 'equal masses', 'position reshape'])


def identity_block():
    pack = lambda array: dict(real=np.asarray(array).tolist(), imag=np.zeros_like(array).tolist())
    amu = np.array([58.]*31+[1.]*6)
    specification = {'finite_oracle_only': True}
    model = sha256(json.dumps(specification, sort_keys=True).encode()).hexdigest()
    identity = dict(phase_space_order='u(original atoms,xyz),p(original atoms,xyz)',
        field_order='a_xyz,Pi_xyz', control='periodic',
        matched_joint_specification=specification, matched_joint_identity=model,
        atom_indices=pack(np.arange(1, 38)),
        atom_names=dict(real=['Ni']*31+['H']*6, imag=['']*37),
        positions_bohr=pack(np.arange(111).reshape(3, 37)), valence_charges=pack(np.array([18.]*31+[1.]*6)),
        mass=dict(variant='original', original_masses_amu=pack(amu), effective_masses_amu=pack(amu)),
        preparation=dict(mass_ry=pack(np.repeat(amu*AMU_RY, 3))),
        cell_indices=pack(np.indices((4, 4, 4)).reshape(3, -1).T))
    original = dict(schema='P254-matched-clamped-inputs-v1', identity=identity)
    geometry(original)
    from copy import deepcopy
    for mutation in ('model', 'mass', 'positions', 'IDs', 'basis'):
        wrong = deepcopy(original)
        i = wrong['identity']
        if mutation == 'model':
            i['matched_joint_identity'] = '0'*64
        elif mutation == 'mass':
            i['mass']['variant'] = 'H2'
        elif mutation == 'positions':
            i['positions_bohr'] = pack(np.arange(111).reshape(37, 3))
        elif mutation == 'IDs':
            i['atom_indices']['real'][1] = 1
        else:
            i['field_order'] = 'collective/sqrtN'
        try:
            geometry(wrong)
        except ValueError:
            pass
        else:
            raise AssertionError('wrong '+mutation+' identity was accepted')
    return dict(wrong_inputs_refused=['model', 'mass variant', 'Cartesian position shape', 'original atom IDs', 'field basis'])


def main():
    result = dict(schema='P254-saved-canonical-response-finite-oracle-v1',
                  algebra=algebra_block(), identity=identity_block(),
                  scope='Finite encoded algebra/normalization/ordering only; no material accuracy or saved-data claim')
    result['status'] = 'passed'
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
