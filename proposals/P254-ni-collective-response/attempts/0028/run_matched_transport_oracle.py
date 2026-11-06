#!/usr/bin/env python3
"""Independent SymPy chain rule, Fourier sums and SciPy dense-operator oracles.

Toy arithmetic tolerances are NOT material coefficient-error licenses. Every
check has an asymmetric mutation witness, including separate Gram/D, Bloch
phase, improper polar parity, mixed-index order and original atom pullback.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace as NS

from matched_reference_io import runtime, save_json
from matched_transport_io import transport_registration


def symbolic_tensors():
    import numpy as np
    import sympy as sp
    from substrate_framework.material_reference_transport import SlotTransport
    from substrate_framework.material_matched_transport import tensor_transport
    r = np.array([[0, 1, 0], [0, 0, -1], [1, 0, 0]])  # nonsymmetric, improper
    atoms = np.array([2, 0, 1])
    u = sp.symbols('u:9')
    a = sp.symbols('a:3')
    zero = dict.fromkeys((*u, *a), 0)
    us, av = sp.Matrix(u), sp.Matrix(a)
    gradient = sp.Matrix([sp.Rational(i+1, 7)+sp.I*sp.Rational(i*i+1, 31) for i in range(9)])
    field = sp.Matrix([sp.Rational(2*i+1, 5)+sp.I*sp.Rational(i+1, 13) for i in range(3)])
    local, nonlocal_, mixed = [], [], []
    for atom in range(3):
        base = sp.Matrix(3, 3, lambda i, j: sp.Rational((atom+1)*(2*i+3*j+1), 19)
                         +sp.I*sp.Rational(1+i*i+j+atom, 23))
        local.append(base+base.T)
        nonlocal_.append((atom+2)*(base+base.T)+sp.eye(3))
        mixed.append(base)
    contact = sp.Matrix(3, 3, lambda i, j: sp.Rational(i+j+1, 11)+sp.I*sp.Rational(i*j+1, 17))
    def scalar(blocks, with_linear=False):
        value = sum((us[3*b:3*b+3, 0].T*blocks[b]*us[3*b:3*b+3, 0])[0]/2 for b in range(3))
        if with_linear:
            value += (gradient.T*us)[0]+(field.T*av)[0]+(av.T*contact*av)[0]/2
            value += sum((us[3*b:3*b+3, 0].T*mixed[b]*av)[0] for b in range(3))
        return sp.expand(value)
    fl, fn = scalar(local, True), scalar(nonlocal_)
    def numeric(value):
        return np.array(value, dtype=complex)
    def derivatives(l, n):
        return dict(positive_gradient_ry_bohr=numeric([sp.diff(l+n, x).subs(zero) for x in u]),
            local_onsite_ry_bohr2=numeric([[[sp.diff(l, u[3*b+i], u[3*b+j]).subs(zero)
                                           for j in range(3)] for i in range(3)] for b in range(3)]),
            nonlocal_onsite_ry_bohr2=numeric([[[sp.diff(n, u[3*b+i], u[3*b+j]).subs(zero)
                                              for j in range(3)] for i in range(3)] for b in range(3)]),
            field_gradient=numeric([sp.diff(l, x).subs(zero) for x in a]),
            current_hartree=-numeric([sp.diff(l, x).subs(zero) for x in a])/2,
            contact=numeric([[sp.diff(l, x, y).subs(zero) for y in a] for x in a]),
            mixed_ionic_field=numeric([[[sp.diff(l, u[3*b+i], a[j]).subs(zero)
                                        for j in range(3)] for i in range(3)] for b in range(3)]))
    original = derivatives(fl, fn)
    original.update(density_bohr3=np.arange(27)+1j*np.arange(27)**2/17,
                    number=np.asarray(2+.3j), number_pair=np.asarray([2+.3j, 2+.3j]))
    transform = SlotTransport(None, np.arange(27)[::-1], atoms, r, np.zeros((3, 3)))
    errors, mutations = {}, {}
    for anti in (False, True):
        epsilon = -1 if anti else 1
        substitutions = {}
        for target, source in enumerate(atoms):
            pulled = sp.Matrix(r).T*us[3*target:3*target+3, 0]
            substitutions.update({u[3*source+i]: pulled[i] for i in range(3)})
        pulled_a = epsilon*sp.Matrix(r).T*av
        substitutions.update({a[i]: pulled_a[i] for i in range(3)})
        # Conjugate coefficients only, not real displacement/field variables.
        left = fl.xreplace({sp.I: -sp.I}) if anti else fl
        right = fn.xreplace({sp.I: -sp.I}) if anti else fn
        expected = derivatives(left.subs(substitutions, simultaneous=True),
                               right.subs(substitutions, simultaneous=True))
        got = tensor_transport(original, transform, orbital_conjugation=anti)
        for name, value in expected.items():
            errors[f'{anti}:{name}'] = float(np.max(abs(got[name]-value)))
        expected_density = original['density_bohr3'][::-1]
        if anti:
            expected_density = expected_density.conj()
        errors[f'{anti}:density'] = float(np.max(abs(got['density_bohr3']-expected_density)))
        errors[f'{anti}:number'] = float(abs(got['number']-(original['number'].conj() if anti else original['number'])))
        wrong = {
            'rotation_transpose': SlotTransport(None, transform.fft_pullback, atoms, r.T, transform.geometry_residual_bohr),
            'inverse_atom_map': SlotTransport(None, transform.fft_pullback, np.argsort(atoms), r, transform.geometry_residual_bohr),
            'improper_axial_parity': SlotTransport(None, transform.fft_pullback, atoms, -r, transform.geometry_residual_bohr),
        }
        for name, mapping in wrong.items():
            candidate = tensor_transport(original, mapping, orbital_conjugation=anti)
            mutations[f'{anti}:{name}'] = max(float(np.max(abs(candidate[k]-expected[k]))) for k in expected)
        mutations[f'{anti}:mixed_transpose'] = float(np.max(abs(got['mixed_ionic_field'].swapaxes(1, 2)-expected['mixed_ionic_field'])))
        if anti:
            mutations['missing_field_parity'] = float(np.max(abs(-got['field_gradient']-expected['field_gradient'])))
            mutations['missing_mixed_parity'] = float(np.max(abs(-got['mixed_ionic_field']-expected['mixed_ionic_field'])))
            mutations['missing_complex_conjugation'] = float(np.max(abs(got['contact'].conj()-expected['contact'])))
    mutations['density_pullback_omitted'] = float(np.max(abs(original['density_bohr3']-original['density_bohr3'][::-1])))
    return errors, mutations


def fourier_orbital():
    import numpy as np
    from scipy.linalg import fractional_matrix_power
    from substrate_framework.material_reference_state import SourceState
    from substrate_framework.material_reference_transport import slot_transport, raw_state_transport
    from substrate_framework.qe_real_bloch_grid import transport_coefficients
    rng = np.random.default_rng(28193)
    c = rng.normal(size=(7, 3))+1j*rng.normal(size=(7, 3))
    gram = c.conj().T @ c
    inverse = fractional_matrix_power(gram, -.5)
    d = inverse @ np.diag([.17, .43, .81]) @ inverse
    state = SourceState(gram, np.array([.17, .43, .81]), d, {}, 'independent-toy')
    miller = np.array([[0, 1, 0, 0, -1, 1, 0], [0, 0, 1, 0, 0, 1, -1], [0, 0, 0, 1, 1, 0, -1]])
    r = np.array([[0, 1, 0], [0, 0, -1], [1, 0, 0]])
    f, ks, shift = np.array([1/3, 0, 2/3]), np.array([.25, .125, 0]), np.array([1, 0, -1])
    op = NS(rotation=r, fractional_translation=f, equivalent_atoms=np.array([1]), effective_time_reversal=False)
    reference = NS(fft_smooth=(3, 3, 3), nat=1, atom_names=('X',), cell_bohr=np.eye(3), positions_bohr=np.zeros((3, 1)))
    xyz = np.array(list(np.ndindex((3, 3, 3))))[:, ::-1].T/3
    source_wave = np.exp(2j*np.pi*xyz.T @ miller) @ c
    source_density = np.einsum('pi,pi->p', source_wave.conj(), source_wave @ d)
    errors, mutations = {}, {}
    for anti in (False, True):
        epsilon = -1 if anti else 1
        slot = NS(operation=op, orbital_conjugation=anti, epsilon=epsilon,
                  source_k_fractional=ks, atom_image_offsets=np.zeros((3, 1)))
        grid = NS(reference=reference, require_slot=lambda value: None)
        mapping = slot_transport(grid, slot)
        transported = raw_state_transport(state, antiunitary=anti)
        ct = transport_coefficients(slot, miller, c.copy())
        mt = epsilon*r @ miller-shift[:, None]
        actual_wave = np.exp(2j*np.pi*xyz.T @ mt) @ ct
        pulled_wave = np.exp(2j*np.pi*(r.T @ xyz-f[:, None]).T @ miller) @ c
        if anti:
            pulled_wave = pulled_wave.conj()
        expected_wave = np.exp(-2j*np.pi*(shift @ xyz+epsilon*np.dot(ks, f)))[:, None]*pulled_wave
        errors[f'{anti}:bloch_phase'] = float(np.max(abs(actual_wave-expected_wave)))
        target_density = np.einsum('pi,pi->p', actual_wave.conj(), actual_wave @ transported.density_coefficients)
        expected_density = mapping.density(source_density.conj() if anti else source_density)
        errors[f'{anti}:fourier_density'] = float(np.max(abs(target_density-expected_density)))
        errors[f'{anti}:gram'] = float(np.max(abs(ct.conj().T @ ct-transported.gram)))
        errors[f'{anti}:number'] = float(abs(np.trace(transported.gram @ transported.density_coefficients)-sum(state.populations)))
        bad_wave = np.exp(2j*np.pi*xyz.T @ mt) @ (c.conj() if anti else c)
        mutations[f'{anti}:phase_omitted'] = float(np.max(abs(bad_wave-expected_wave)))
        if anti:
            bad_density = np.einsum('pi,pi->p', actual_wave.conj(), actual_wave @ d)
            mutations['D_not_conjugated'] = float(np.max(abs(bad_density-expected_density)))
            mutations['Gram_not_conjugated'] = float(np.max(abs(gram-transported.gram)))
            mutations['spin_or_population_refill'] = float(abs(np.trace(transported.gram @ np.diag(state.populations))-sum(state.populations)))
    return errors, mutations


def independent_projectors():
    import numpy as np
    import sympy as sp
    from scipy.linalg import norm
    from substrate_framework.material_angular_intertwiners import original_channel_intertwiner
    from substrate_framework.material_matched_transport_certificate import nonlocal_factor_bound
    x, y, z = sp.symbols('x y z', real=True)
    polynomials = [1, z, -x, -y, (2*z*z-x*x-y*y)/2, -sp.sqrt(3)*x*z,
                   -sp.sqrt(3)*y*z, sp.sqrt(3)*(x*x-y*y)/2, sp.sqrt(3)*x*y]
    ell = np.array([0, 1, 1, 1, 2, 2, 2, 2, 2])
    upf = NS(angular_l=np.array([0, 1, 2]), channel_radial=np.repeat([0, 1, 2], [1, 3, 5]),
             channel_l=ell, channel_lm=np.arange(9), nh=9)
    evaluate = [sp.lambdify((x, y, z), expr, 'numpy') for expr in polynomials]
    def beta(p, tau):
        harmonics = np.column_stack([np.broadcast_to(fn(*p.T), len(p)) for fn in evaluate])
        return harmonics*((-1j)**ell)[None, :]*np.exp(-.01*np.sum(p*p, axis=1)-1j*(p @ tau))[:, None]
    r = np.array([[0, 1, 0], [0, 0, -1], [1, 0, 0]])
    m = np.array([[0, 0, 0], [1, 0, 0], [0, -1, 0], [0, 0, 1], [1, 1, -1]])
    ks, f, image = np.array([.25, .125, .375]), np.array([.25, .5, 0]), np.array([1, -1, 2])
    ps = 2*np.pi*(m+ks)
    tau = np.array([.1, .2, -.3])
    target_tau = r @ (tau+f+image)
    bs = beta(ps, tau)
    q = original_channel_intertwiner(upf, r)
    rng = np.random.default_rng(291)
    active = rng.normal(size=(9, 9))+1j*rng.normal(size=(9, 9))
    errors, mutations = {}, {}
    for anti in (False, True):
        epsilon = -1 if anti else 1
        phase = np.exp(-epsilon*2j*np.pi*((m+ks) @ f+ks @ image))
        c = phase[:, None]*(bs.conj() if anti else bs)
        bt = beta(epsilon*ps @ r.T, target_tau)
        ds = active.conj() if anti else active
        result = nonlocal_factor_bound(c, bt, ds, active, q)
        errors[f'{anti}:projector_phase_and_Q'] = float(norm(bt-c @ q))
        actual = norm(bt @ active @ bt.conj().T-c @ ds @ c.conj().T, 2)
        errors[f'{anti}:bound_violation'] = float(max(0., actual-result['nonlocal_operator_bound_ry']))
        mutations[f'{anti}:Q_identity'] = float(norm(bt-c))
        mutations[f'{anti}:Q_transpose'] = float(norm(bt-c @ q.T))
        mutations[f'{anti}:atomic_image_phase_omitted'] = float(norm(bt-(c/np.exp(-epsilon*2j*np.pi*(ks @ image))) @ q))
        mutations[f'{anti}:D_conjugacy_ignored'] = float(actual)
    return errors, mutations


def reciprocal_fold_identity():
    import numpy as np
    from substrate_framework.material_matched_transport import identity_representative
    identity = NS(rotation=np.eye(3, dtype=int), fractional_translation=np.zeros(3))
    rotation = NS(rotation=np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]),
                  fractional_translation=np.zeros(3))
    shift = np.array([0, 0, -1])
    representative = NS(ispin=1, source_ik=7, orbital_conjugation=False,
                        reciprocal_shift=shift, operation=identity)
    alternatives = [
        NS(ispin=1, source_ik=7, orbital_conjugation=True,
           reciprocal_shift=np.zeros(3), operation=identity),
        NS(ispin=1, source_ik=7, orbital_conjugation=False,
           reciprocal_shift=np.zeros(3), operation=rotation),
        NS(ispin=2, source_ik=7, orbital_conjugation=False,
           reciprocal_shift=np.zeros(3), operation=identity)]
    selected = identity_representative(NS(slots=alternatives+[representative]), 1, 7)
    m = np.array([[0, 1, -1], [1, 0, 2], [0, -1, 1]])
    k = np.array([.25, 0., .5])
    xyz = np.array([[.1, .23, .41], [.33, .19, .76]])
    c = np.array([[1+.3j, .2], [-.7j, .9+.2j], [.31, -.2+.4j]])
    d = np.array([[.7, .1+.2j], [.1-.2j, .4]])
    source = np.exp(2j*np.pi*(xyz @ m)) @ c
    folded = np.exp(2j*np.pi*(xyz @ (m-shift[:, None]))) @ c
    expected = np.exp(-2j*np.pi*(xyz @ shift))[:, None]*source
    density = lambda psi: np.einsum('pi,pi->p', psi.conj(), psi @ d)
    errors = dict(selector=float(selected is not representative),
        physical_momenta=float(np.max(abs((m-shift[:, None])+(k+shift)[:, None]-(m+k[:, None])))),
        periodic_gauge=float(np.max(abs(folded-expected))),
        coefficient_density=float(np.max(abs(density(folded)-density(source)))))
    mutations = dict(missing_miller_relabel=float(np.max(abs(m+(k+shift)[:, None]-(m+k[:, None])))),
        missing_k_fold=float(np.max(abs((m-shift[:, None])+k[:, None]-(m+k[:, None])))),
        incorrect_frame_phase=float(np.max(abs(folded-source))))
    for name, slots in (('no_original_identity', alternatives),
                        ('duplicate_identity', [representative, representative])):
        try:
            identity_representative(NS(slots=slots), 1, 7)
        except ValueError:
            errors[name] = 0.
        else:
            errors[name] = 1.
    return errors, mutations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    runtime()
    registration, _ = transport_registration()
    groups = {}
    for name, run in (('symbolic_scalar_chain_rule', symbolic_tensors),
                      ('independent_fourier_orbital', fourier_orbital),
                      ('independent_polynomial_projectors', independent_projectors),
                      ('reciprocal_fold_identity', reciprocal_fold_identity)):
        residuals, mutations = run()
        if not all(value < 2e-10 for value in residuals.values()):
            raise AssertionError((name, 'oracle residual', residuals))
        if not all(value > 1e-6 for value in mutations.values()):
            raise AssertionError((name, 'insensitive mutation', mutations))
        groups[name] = dict(residuals=residuals, mutation_witnesses=mutations)
    result = dict(status='passed', registration=registration, groups=groups,
        arithmetic_tolerance=2e-10, material_tolerance_license=False,
        scope='independent initial tensor/orbital/phase equations, not a material covariance claim')
    print(json.dumps(dict(result=save_json(args.output, result)), sort_keys=True))


if __name__ == '__main__':
    main()
