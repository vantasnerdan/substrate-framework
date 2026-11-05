"""Expose the real NESSi result against independent full-Hamiltonian SciPy expm.

Reads ORIGINAL input plus three actual executable outputs: nominal, altered
inactive branch storage, and intentionally wrong thermal preparation. Does not
supply a fallback trajectory to the native solver. Full retarded/lesser triangles
and every equal-time density are compared, not merely a process exit code.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm


def binding(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    return dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


def read_original(path):
    with Path(path).open() as stream:
        header = stream.readline().split()
        if len(header) != 6 or header[0] != 'P254_GAUSSIAN_IVP_V1':
            raise ValueError('wrong input schema')
        ns, nb, nt, order = map(int, header[1:5])
        dt = float(header[5])
        values = np.fromstring(stream.read(), sep=' ')
    if not np.isfinite(values).all() or len(values)%2:
        raise ValueError('invalid complex input')
    values = values[::2]+1j*values[1::2]
    position = 0
    def matrix(rows, cols):
        nonlocal position
        width = rows*cols
        a = values[position:position+width].reshape(rows, cols).copy()
        position += width
        return a
    rs, rb = matrix(ns, ns), matrix(nb, nb)
    hs, hb, v = [], [], []
    for _ in range(nt+1):
        hs.append(matrix(ns, ns)); hb.append(matrix(nb, nb)); v.append(matrix(ns, nb))
    if position != len(values):
        raise ValueError('trailing input numbers')
    return dict(ns=ns, nb=nb, nt=nt, dt=dt, order=order, rho_s=rs, rho_b=rb,
                hs=hs, hb=hb, v=v)


def read_native(path):
    with Path(path).open() as stream:
        header = stream.readline().split()
        if len(header) != 4 or header[0] != 'P254_NESSI_REALTIME_V1':
            raise ValueError('wrong native output schema')
        ns, nt, dt = int(header[1]), int(header[2]), float(header[3])
        result = dict(ns=ns, nt=nt, dt=dt)
        for key in ('retarded_n_m', 'lesser_m_n', 'sigma_retarded_n_m', 'sigma_lesser_m_n', 'density_n'):
            if stream.readline().strip() != key:
                raise ValueError('missing native component '+key)
            rows = []
            for n in range(nt+1):
                for m in ([None] if key == 'density_n' else range(n+1)):
                    values = stream.readline().split()
                    width = 1 if m is None else 2
                    if len(values) != width+2*ns*ns or int(values[0]) != n or (m is not None and int(values[1]) != m):
                        raise ValueError('wrong native component indices/width')
                    numbers = np.array(values[width:], dtype=float)
                    rows.append((numbers[::2]+1j*numbers[1::2]).reshape(ns, ns))
            result[key] = np.asarray(rows)
        mixed = stream.readline().split()
        if len(mixed) != 2 or mixed[0] != 'inactive_mixed_max':
            raise ValueError('missing inactive-branch measure')
        result['inactive_mixed_max'] = float(mixed[1])
        probe = stream.readline().split()
        if len(probe) != 3 or probe[0] != 'inert_branch_probe' or stream.read().strip():
            raise ValueError('missing actual branch probe values or trailing native output')
        result['inert_branch_probe'] = np.array(probe[1:], dtype=float)
    for key, value in result.items():
        if key not in ('ns', 'nt') and not np.isfinite(value).all():
            raise ValueError('nonfinite native quantity: '+key)
    return result


def full_hamiltonian_reference(x):
    """Full coupled system+bath propagation, unlike the NESSi embedding route."""
    ns, nb, nt, dt = x['ns'], x['nb'], x['nt'], x['dt']
    full_rho = np.zeros((ns+nb, ns+nb), dtype=complex)
    full_rho[:ns,:ns], full_rho[ns:,ns:] = x['rho_s'], x['rho_b']
    evolution = [np.eye(ns+nb, dtype=complex)]
    bath = [np.eye(nb, dtype=complex)]
    for n in range(nt):
        full_h = np.block([[x['hs'][n], x['v'][n]], [x['v'][n].conj().T, x['hb'][n]]])
        evolution.append(expm(-1j*dt*full_h)@evolution[-1])
        bath.append(expm(-1j*dt*x['hb'][n])@bath[-1])
    result = dict(retarded_n_m=[], lesser_m_n=[], sigma_retarded_n_m=[], sigma_lesser_m_n=[], density_n=[])
    for n in range(nt+1):
        un = evolution[n]
        result['density_n'].append((un@full_rho@un.conj().T)[:ns,:ns])
        for m in range(n+1):
            um = evolution[m]
            result['retarded_n_m'].append((-1j*un@um.conj().T)[:ns,:ns])
            result['lesser_m_n'].append((1j*um@full_rho@un.conj().T)[:ns,:ns])
            result['sigma_retarded_n_m'].append(-1j*x['v'][n]@bath[n]@bath[m].conj().T@x['v'][m].conj().T)
            result['sigma_lesser_m_n'].append(1j*x['v'][m]@bath[m]@x['rho_b']@bath[n].conj().T@x['v'][n].conj().T)
    return {key: np.asarray(value) for key, value in result.items()}


def compare(original, nominal, inert, wrong, wrong_input):
    for result in (nominal, inert, wrong):
        if (result['ns'],result['nt'],result['dt']) != (original['ns'],original['nt'],original['dt']):
            raise ValueError('native output grid differs from ORIGINAL input')
    reference = full_hamiltonian_reference(original)
    errors = {key: float(np.max(np.abs(nominal[key]-value))) for key, value in reference.items()}
    branch = {key: float(np.max(np.abs(inert[key]-nominal[key]))) for key in reference}
    original_rho = original['rho_s']
    initial_error = float(np.max(np.abs(nominal['density_n'][0]-original_rho)))
    wrong_error = float(np.max(np.abs(wrong['density_n']-reference['density_n'])))
    wrong_initial = float(np.max(np.abs(wrong['density_n'][0]-original_rho)))
    wrong_ret = float(np.max(np.abs(wrong['retarded_n_m']-nominal['retarded_n_m'])))
    sigma_size = float(np.max(np.abs(nominal['sigma_retarded_n_m'])))
    coherence = float(abs(original_rho[0,1]))
    commutator = float(np.linalg.norm(original['hs'][0]@original_rho-original_rho@original['hs'][0]))
    occupations = np.linalg.eigvalsh(original_rho)
    unchanged_mutation_operators = (
        (wrong_input['ns'],wrong_input['nb'],wrong_input['nt'],wrong_input['dt']) ==
        (original['ns'],original['nb'],original['nt'],original['dt']) and
        np.array_equal(wrong_input['rho_b'],original['rho_b']) and
        all(np.array_equal(wrong_input[key],original[key]) for key in ('hs','hb','v')))
    thermal = np.linalg.solve(np.eye(original['ns'])+expm(2*original['hs'][0]),
                              np.eye(original['ns']))
    thermal_input_error = float(np.max(np.abs(wrong_input['rho_s']-thermal)))
    wrong_state_error = float(np.max(np.abs(wrong['density_n'][0]-wrong_input['rho_s'])))
    predicates = dict(
        all_native_components_match_full_original_expm=max(errors.values()) < 2e-6,
        exact_original_initial_state=initial_error < 2e-13,
        branch_storage_inactive=max(branch.values()) < 2e-13,
        nominal_mixed_zero=nominal['inactive_mixed_max'] == 0,
        inert_probe_mixed_zero=inert['inactive_mixed_max'] == 0,
        nonzero_embedding=sigma_size > .05,
        complex_initial_coherence=abs(original_rho[0,1].imag) > .1,
        nonstationary_original=commutator > .1,
        exact_zero_one_occupancies=bool(np.array_equal(original_rho@original_rho, original_rho) and
            np.trace(original_rho) == 1 and np.array_equal(original_rho,original_rho.conj().T)),
        wrong_thermal_overwrite_rejected=wrong_error > .1 and wrong_initial > .1,
        retarded_state_independence=wrong_ret < 2e-13,
        wrong_input_is_actual_thermal_overwrite=unchanged_mutation_operators and thermal_input_error < 2e-13,
        wrong_native_preserves_actual_mutated_input=wrong_state_error < 2e-13,
        actual_branch_probe_values=np.array_equal(nominal['inert_branch_probe'],[1.,0.]) and
            np.array_equal(inert['inert_branch_probe'],[3.7,7.25]))
    predicates = {name: bool(value) for name, value in predicates.items()}
    return dict(status='passed' if all(predicates.values()) else 'failed', predicates=predicates,
        native_vs_full_expm_max_abs=errors, inert_branch_changes_max_abs=branch,
        initial_density_max_abs_error=initial_error, wrong_thermal_max_density_error=wrong_error,
        wrong_input_thermal_density_error=thermal_input_error,
        wrong_native_mutated_initial_density_error=wrong_state_error,
        wrong_thermal_initial_density_error=wrong_initial, wrong_thermal_retarded_change=wrong_ret,
        embedding_retarded_max_abs=sigma_size, initial_coherence_abs=coherence,
        initial_H_rho_commutator_frobenius=commutator, initial_system_occupations=occupations.tolist(),
        original_initial_density_real=original_rho.real.tolist(), original_initial_density_imag=original_rho.imag.tolist(),
        observed_final_native_density_real=nominal['density_n'][-1].real.tolist(),
        observed_final_native_density_imag=nominal['density_n'][-1].imag.tolist(),
        observed_final_full_expm_density_real=reference['density_n'][-1].real.tolist(),
        observed_final_full_expm_density_imag=reference['density_n'][-1].imag.tolist(),
        evidence_scope='one three-orbital coherent product-Gaussian embedding IVP; finite-step audited comparison, not rigorous enclosure, material rate or Ni/H collision kernel')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ('input', 'wrong-input', 'native-output', 'inert-output', 'wrong-output', 'executable', 'report'):
        parser.add_argument('--'+argument, required=True)
    args = parser.parse_args()
    paths = dict(original_input=args.input, native=args.native_output, inert_branch=args.inert_output,
                 wrong_preparation=args.wrong_output, wrong_input=args.wrong_input,
                 native_executable=args.executable)
    result = compare(read_original(args.input), read_native(args.native_output),
                     read_native(args.inert_output), read_native(args.wrong_output),
                     read_original(args.wrong_input))
    result['schema'] = 'P254-A40-NESSi-original-Gaussian-expm-consumer-v1'
    result['bindings'] = {name: binding(path) for name,path in paths.items()}
    result['oracle_source'] = binding(__file__)
    with Path(args.report).open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(result, indent=2, allow_nan=False))
    raise SystemExit(0 if result['status']=='passed' else 1)


if __name__ == '__main__':
    main()
