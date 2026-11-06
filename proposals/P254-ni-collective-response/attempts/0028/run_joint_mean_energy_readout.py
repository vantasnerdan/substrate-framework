#!/usr/bin/env python3
"""Read original complete128+128 scalar contributions into a mean-energy ledger.

Usage: python run_joint_mean_energy_readout.py --registration PATH --qualification PATH --output PATH
No trajectory or full-slot contraction is recomputed; admitted parameters and
original cache/stage values are consumed through the existing restoration APIs.
Outputs global Ry components, their changes and unfitted saved-total discrepancies.
Mean Maxwell energy is neither quantum fluctuation energy nor irreversible heat.
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
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'


def timeout(signum, frame):
    raise TimeoutError('saved mean energy readout exceeded explicit wall-time cap')


def evaluate(model, stage, one_body, evidence, component_function):
    import numpy as np
    from substrate_framework.streamed_material_evolution import real_physical
    raw = component_function(reference=model.reference, coefficients=model.coefficients,
        field=model.field, charges=model.charges, stage=stage, one_body=one_body,
        grid_measure=model.measure)
    values = {name: float(real_physical(value, name, model.reality_tolerance)) for name, value in raw.items()}
    total = math.fsum(values.values())
    discrepancy = total-stage.energy
    return dict(**evidence, components_global_ry=values,
        component_imaginary_parts_ry={name: float(np.imag(value)) for name, value in raw.items()},
        total_global_ry=total, saved_original_total_global_ry=stage.energy,
        discrepancy_ry=discrepancy, absolute_discrepancy_ry=abs(discrepancy),
        consistency_scope='Actual summation discrepancy is reported without fitting, redistribution or an invented floating-error bound; not a trajectory-accuracy certificate.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration', type=Path, required=True)
    parser.add_argument('--qualification', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpu', type=int, default=10)
    parser.add_argument('--max-seconds', type=float, default=180.)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if args.cpu not in os.sched_getaffinity(0) or not math.isfinite(args.max_seconds) or args.max_seconds <= 0:
        parser.error('one available allocated CPU and finite positive --max-seconds required')
    os.sched_setaffinity(0, {args.cpu})
    signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, args.max_seconds)
    started = monotonic()
    result = dict(schema='P254-saved-coupled-mean-energy-readout-v1', status='failed',
        units='global Ry; initial/final are original saved means; changes are final minus initial',
        scope='BASELINE homogeneous coupled mean energy, all original128 initial plus128 final slots',
        scientific_limits=[
            'No finite-local-excitation energy, quantum covariance/fluctuation energy, heat or coherence.',
            'Mean Maxwell energy is not quantum fluctuation energy or irreversible heat.',
            'No full E/D/C, spatial/response completion or independent trajectory-accuracy claim.',
            'Original geometry, nonstationary preparation and finite ionic Taylor approximation retained.',
            'Algebraic oracle and saved-total agreement do not certify the original trajectory.'],
        execution=dict(cpu=args.cpu, blas_threads=1, max_seconds=args.max_seconds),
        actual_scope=dict(one_body_scalars_read=0, saved_stages_restored=0,
            consumption='Original scalar contractions and saved stages, not new trajectory/contraction production.'))
    try:
        from threadpoolctl import threadpool_limits
        from substrate_framework.streamed_material_state import file_binding
        from joint_mean_energy_components import components
        from joint_mean_energy_sources import read_stage
        from joint_mean_energy_authority import restore_current_sources
        # Deliberately runtime-only: the independent module must already exist.
        from joint_mean_energy_oracle import check_components
        with threadpool_limits(limits=1):
            reg, model, saved, provenance = restore_current_sources(
                file_binding(args.registration), file_binding(args.qualification))
            result['provenance'] = provenance
            result['implementation_sources'] = [file_binding(HERE/name) for name in
                ('run_joint_mean_energy_readout.py', 'joint_mean_energy_sources.py',
                 'joint_mean_energy_authority.py', 'joint_mean_energy_components.py',
                 'joint_mean_energy_oracle.py')]
            result['oracle'] = check_components(components)
            if result['oracle'].get('status') != 'passed':
                raise ValueError('independent canonical Hamiltonian/mutation oracle did not pass')
            stages = []
            for index, label in enumerate(('initial', 'final')):
                stage, one_body, evidence = read_stage(reg, model, saved, index)
                result['actual_scope']['one_body_scalars_read'] += evidence['scalar_count']
                result['actual_scope']['saved_stages_restored'] += 1
                ledger = evaluate(model, stage, one_body, evidence, components)
                result[label] = ledger
                stages.append(ledger)
            initial, final = stages
            result['changes'] = dict(components_global_ry={name: final['components_global_ry'][name]-value
                for name, value in initial['components_global_ry'].items()},
                total_global_ry=final['total_global_ry']-initial['total_global_ry'],
                saved_original_total_global_ry=final['saved_original_total_global_ry']-initial['saved_original_total_global_ry'],
                discrepancy_ry=final['discrepancy_ry']-initial['discrepancy_ry'])
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
    if result['status'] != 'complete':
        raise SystemExit(2)


if __name__ == '__main__':
    main()
