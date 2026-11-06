#!/usr/bin/env python3
"""Authenticate original four-case inputs and report initial local pulse dose.

No load_case, wavefunction, response leg, FFT, trajectory or H0/K action is
called.  Only saved parameter archive members and JSON metadata are read.
Usage: python run_local_pulse_energy.py --registration PATH --output PATH
The complete JSON report is accompanied by a compact stdout summary.
"""
import argparse
import json
import math
import os
from pathlib import Path
import resource
import signal
import sys
from time import monotonic

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]
from local_pulse_energy_sources import restore_registration, centered_preparation, saved_gamma_comparison
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'


def timeout(signum, frame):
    raise TimeoutError('initial local pulse readout exceeded explicit wall-time cap')


def typed_array(value):
    return dict(dtype=value.dtype.str, shape=list(value.shape), values=value.tolist())




def case_ledger(case, aggregate, impulse, mass_unit, energy_function, fourier_function):
    import numpy as np
    from substrate_framework.streamed_material_state import authenticate
    names = ('atomic_masses_amu', 'valence_charges', 'atom_indices', 'atom_names')
    with np.load(authenticate(case['archive']), allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in names}
    masses, charges, indices, elements = (arrays[name] for name in names)
    if (masses.ndim != 1 or not masses.size or any(a.shape != masses.shape for a in arrays.values())
            or not np.isfinite(masses).all() or not np.isfinite(charges).all()
            or np.any(masses <= 0) or np.any(charges <= 0)
            or not np.issubdtype(indices.dtype, np.integer)
            or len(np.unique(indices)) != len(indices)):
        raise ValueError('complete finite positive original atom parameters required')
    if case['atom'] != impulse['atom'] or case['axis'] != impulse['axis']:
        raise ValueError('per-case canonical pulse target differs')
    positions = np.flatnonzero(indices == case['atom'])
    if len(positions) != 1 or case['axis'] not in ('x', 'y', 'z'):
        raise ValueError('unique original atom index and Cartesian axis required')
    index = int(positions[0])
    if elements[index] != case['element'] or masses[index] != case['original_mass_amu']:
        raise ValueError('original struck atom element/mass differs from registration')
    mass = float(masses[index]*mass_unit['AMU_RY'])
    charge = float(charges[index])
    dose = energy_function(mass, charge, impulse['amplitude'], 0, 0)
    rows = fourier_function(impulse['mesh'], impulse['cell'], impulse['amplitude'])
    shares = [dict(mesh_index=list(row['mesh_index']), q_fractional=list(row['q_fractional']),
        canonical_momentum_coefficient=[row['coefficient'].real, row['coefficient'].imag],
        quadratic_momentum_weight=row['quadratic_weight'],
        fraction_of_local_quadratic_input=row['quadratic_weight']/impulse['amplitude']**2,
        initial_quadratic_input_ry=row['quadratic_weight']/(2*mass)) for row in rows]
    gamma = next(row for row in shares if row['mesh_index'] == [0, 0, 0])
    return dict(case=case['case'], original_registration_entry=case, original_aggregate=aggregate,
        original_parameter_fields={name: typed_array(value) for name, value in arrays.items()},
        struck_atom=dict(original_atom_index=case['atom'], zero_based_array_index=index,
            element=str(elements[index]), axis=case['axis'], mass_amu=float(masses[index]),
            mass_Ry=mass, original_UPF_valence_charge=charge),
        preparation='Declared centered matched preparation; not equilibrium',
        initial_input_components_ry=dose, initial_local_input_ry=dose['total'],
        initial_local_input_joule=dose['total']*mass_unit['RYDBERG_SI'],
        Fourier=dict(convention='Pq=N^(-1/2) sum_R exp(-2pi i q.R) p_R; inverse has positive phase',
            ncells=impulse['ncells'], mesh=impulse['mesh'], all_q_count=len(shares), shares=shares,
            sum_quadratic_input_ry=math.fsum(row['initial_quadratic_input_ry'] for row in shares),
            Parseval_summation_discrepancy_ry=math.fsum(row['initial_quadratic_input_ry'] for row in shares)-dose['quadratic'],
            one_Gamma_share=gamma,
            scope='All64 initial centered quadratic input shares, including q/-q once each. Not propagated component energy; oneGamma is not the whole local pulse.'))






def matched_comparison(ledgers, amplitude, energy_function):
    reference = next(row for row in ledgers if row['case'] == 'loaded-ni31h6')
    reference_energy = reference['initial_local_input_ry']
    reference_mass = reference['struck_atom']['mass_Ry']
    rows = []
    for row in ledgers:
        mass = row['struck_atom']['mass_Ry']
        matched_amplitude = amplitude*math.sqrt(mass/reference_mass)
        matched_energy = energy_function(mass, row['struck_atom']['original_UPF_valence_charge'], matched_amplitude, 0, 0)['total']
        rows.append(dict(case=row['case'], same_impulse_amplitude=amplitude,
            same_impulse_input_ry=row['initial_local_input_ry'],
            ratio_to_loaded_input=row['initial_local_input_ry']/reference_energy,
            difference_from_loaded_input_ry=row['initial_local_input_ry']-reference_energy,
            amplitude_for_loaded_dose=matched_amplitude, input_at_matched_amplitude_ry=matched_energy,
            matched_dose_summation_discrepancy_ry=matched_energy-reference_energy))
    same_mass = all(row['struck_atom']['mass_Ry'] == reference_mass for row in ledgers)
    return dict(reference_case=reference['case'], reference_local_dose_ry=reference_energy,
        actual_struck_masses_equal=same_mass,
        same_impulse_is_matched_dose=all(row['initial_local_input_ry'] == reference_energy for row in ledgers),
        cases=rows, interpretation='Only the initial input is compared. Equal centered Ni1 masses imply equal dose at equal canonical impulse; UPF charge drops out only because the declared mechanical mean is zero. No material-response ranking follows.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpu', type=int, default=10)
    parser.add_argument('--max-seconds', type=float, default=180.)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if (args.cpu not in range(6, 11) or args.cpu not in os.sched_getaffinity(0)
            or not math.isfinite(args.max_seconds) or args.max_seconds <= 0 or args.max_seconds > 180):
        parser.error('available allocated CPU6-10 and positive wall cap <=180 seconds required')
    os.sched_setaffinity(0, {args.cpu})
    memory_cap = 2 << 30
    _, hard = resource.getrlimit(resource.RLIMIT_AS)
    if hard != resource.RLIM_INFINITY:
        memory_cap = min(memory_cap, hard)
    resource.setrlimit(resource.RLIMIT_AS, (memory_cap, hard))
    signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, args.max_seconds)
    started = monotonic()
    result = dict(schema='P254-initial-local-pulse-energy-v1', status='failed',
        units='Finite local input in global Ry and joule; canonical pulse amplitude in bohr^-1 (hbar=1)',
        scope='INITIAL INPUT ONLY: declared centered matched preparation, all4 actual parameter archives, all64 Fourier shares and separately authenticated saved Gamma seed metadata',
        scientific_limits=['No propagated component redistribution, nuclear input, noise, heat or full D/C closure.',
            'No new trajectory, FFT, orbital/contraction read, saved-band substitute or residual allocation.',
            'Original nonstationary sudden-release product is not equilibrium; original finite Taylor/geometry limitations retained.',
            'Fourier input weights are not independently evolved complex-q energies or evidence of positive-time convergence.'],
        execution=dict(cpu=args.cpu, blas_threads=1, max_seconds=args.max_seconds, address_space_cap_bytes=memory_cap))
    try:
        import numpy as np
        import sympy as sp
        from threadpoolctl import threadpool_limits
        from substrate_framework.streamed_material_state import file_binding
        from substrate_framework.material_modes import AMU_RY, AMU_SI, ELECTRONMASS_SI, RYDBERG_SI
        from local_pulse_energy import pulse_components, fourier_pulse
        from local_pulse_energy_oracle import check_pulse
        with threadpool_limits(limits=1):
            reg, aggregates, cases, proof, adoption, binding = restore_registration(args.registration)
            result.update(registration_binding=binding, original_registration=reg,
                implementation_sources=[file_binding(HERE/name) for name in
                    ('local_pulse_energy.py', 'local_pulse_energy_oracle.py', 'local_pulse_energy_sources.py',
                     'run_local_pulse_energy.py', 'joint_mean_energy_sources.py')]
                    +[file_binding(ROOT/'src/substrate_framework'/name) for name in
                      ('streamed_material_state.py', 'partitioned_joint_slots.py')],
                mass_unit=dict(AMU_RY=AMU_RY, AMU_SI=AMU_SI, ELECTRONMASS_SI=ELECTRONMASS_SI,
                    RYDBERG_SI=RYDBERG_SI, convention='M_Ry=m_amu*AMU_SI/ELECTRONMASS_SI/2; authenticated material_modes.py constants'),
                versions=dict(python=sys.version, numpy=np.__version__, sympy=sp.__version__))
            impulse = reg['impulse']
            if (impulse['mesh'] != [4, 4, 4] or impulse['ncells'] != math.prod(impulse['mesh'])
                    or impulse['atom'] != 1 or impulse['axis'] != 'x'
                    or impulse['cell'] != [0, 0, 0] or impulse['amplitude'] != 1.0):
                raise ValueError('registered original amplitude1/cell0/Ni1x/4^3 request required')
            result['centered_preparation'] = centered_preparation(reg, cases, proof)
            result['oracle'] = check_pulse(pulse_components, fourier_pulse,
                mesh=impulse['mesh'], cell=impulse['cell'], mass_unit=result['mass_unit'])
            if result['oracle']['status'] != 'passed':
                raise ValueError('independent differential-operator/Fourier/mutation oracle did not pass')
            result['cases'] = [case_ledger(case, aggregates[case['case']], impulse, result['mass_unit'],
                                          pulse_components, fourier_pulse) for case in reg['cases']]
            loaded = next(row for row in result['cases'] if row['case'] == 'loaded-ni31h6')
            result['saved_Gamma_comparison'] = saved_gamma_comparison(reg, adoption, proof['initialization'], loaded, fourier_pulse)
            result['matched_dose_comparison'] = matched_comparison(result['cases'], impulse['amplitude'], pulse_components)
            result['actual_reads'] = dict(original_parameter_archives=len(result['cases']),
                archive_members_per_case=['atomic_masses_amu', 'valence_charges', 'atom_indices', 'atom_names'],
                saved_seed='adoption002.new_seed JSON metadata and source-bound request; no leg payloads')
            result['status'] = 'complete'
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        result['execution'].update(elapsed_seconds=monotonic()-started,
            process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result, stream, separators=(',', ':'), allow_nan=False)
        stream.write('\n')
    summary = dict(schema=result['schema'], status=result['status'], output=str(args.output.resolve()),
        scope='initial local canonical-pulse input only', oracle=result.get('oracle', {}).get('status'))
    if result['status'] == 'complete':
        summary.update(local_inputs_ry={row['case']: row['initial_local_input_ry'] for row in result['cases']},
            one_Gamma_inputs_ry={row['case']: row['Fourier']['one_Gamma_share']['initial_quadratic_input_ry'] for row in result['cases']},
            matched_dose=result['matched_dose_comparison']['same_impulse_is_matched_dose'],
            actual_saved_Gamma_metadata=result['saved_Gamma_comparison']['status'])
    else:
        summary['error'] = result.get('error')
    print(json.dumps(summary, separators=(',', ':'), allow_nan=False))
    if result['status'] != 'complete':
        raise SystemExit(2)


if __name__ == '__main__':
    main()
