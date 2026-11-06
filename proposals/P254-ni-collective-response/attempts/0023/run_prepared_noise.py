#!/usr/bin/env python3
"""Actual isolated quasi-free fluctuation input; no Native work or coupled evolution."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import signal
import sys
import time
import traceback

BOOT = time.monotonic()
import numpy as np
from substrate_framework.pw_prepared_noise import prepare_pw_noise
from substrate_framework.qe_pw_response_frames import build_pw_response_models
from substrate_framework.quasifree_spectral_noise import ordered_transition_correlation, paired_correlations
from noise_pilot_io import (HERE, binding, bound_registration, implementation, load_inputs,
                           coordinate_bindings, runtime_context, save_arrays, save_json, verify_buffer)


def now():
    return datetime.now(timezone.utc).isoformat()


def prepare_point(target, point, model, reciprocal_basis, times, budget, event, records):
    branches, metadata = prepare_pw_noise(model, reciprocal_basis)
    save_json(target/(point+'-preparation.json'), metadata)
    total = {name:np.zeros((len(times),8,8),complex) for name in branches}
    population_variance = np.zeros(2)
    spins = {record['identity']:record['spin'] for record in metadata['source_channels']}
    for branch, blocks in branches.items():
        for index, block in enumerate(blocks):
            budget()
            arguments = {key:value for key,value in block.items() if key != 'identity'}
            correlation = ordered_transition_correlation(**arguments, times=times)
            key = f'{point}-{branch}-{index:02d}'
            arrays = {**arguments,'weight':np.asarray(arguments['weight']),
                      'times_tau_ry':times,'ordered_correlation':correlation}
            receipt = save_arrays(target,key,arrays,{'identity':block['identity'],'point':point,
                'branch':branch,'observable_order':metadata['observables'],
                'preparation':metadata['preparation'],'normalization':metadata['normalization']})
            records[key] = receipt
            total[branch] += correlation
            if branch == 'plus':
                f = block['source_occupations']
                population_variance[spins[block['identity']]] += block['weight']*np.sum(f*(1-f))
            event('preparation_block_committed',point=point,branch=branch,identity=block['identity'],archive=receipt['archive'])
            budget()
    reverse_indices = [int(np.flatnonzero(times == -tau).item()) for tau in times]
    paired = paired_correlations(total['plus'],total['minus'][reverse_indices])
    diagnostics = {}
    if point == 'gamma':
        expected = np.broadcast_to(np.diag(population_variance), (len(times),2,2))
        error = float(np.max(np.abs(total['plus'][:,6:,6:] - expected)))
        number_commutator = float(max(np.max(np.abs(paired['commutator'][:,6:,:])),
                                     np.max(np.abs(paired['commutator'][:,:,6:]))))
        diagnostics = {'Gamma_spin_number_variance':population_variance.tolist(),
                       'number_covariance_maxabs_error':error,
                       'number_commutator_maxabs':number_commutator,
                       'allowance':1e-10,'scope':'Auxiliary finite grandcanonical preparation, not material accuracy.'}
        save_json(target/(point+'-diagnostics.json'),diagnostics)
        if error > 1e-10 or number_commutator > 1e-10:
            raise ValueError('actual Gamma preparation violates independent spin-number fluctuation/conservation identity')
    output = save_arrays(target,point+'-correlations',{'times_tau_ry':times,
        'ordered_plus':total['plus'],'ordered_minus':total['minus'],**paired},
        {'preparation':metadata,'diagnostics':diagnostics})
    event('q_preparation_completed',point=point,archive=output['archive'])
    return {'preparation':metadata,'correlations':output,'diagnostics':diagnostics}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--admission',type=Path,required=True)
    args = parser.parse_args()
    target = args.output_dir.resolve()
    target.mkdir(parents=True,exist_ok=False)
    records, points = {}, {}
    report = {'started_at':now(),'status':'running','preregistration_commit':'3b93dfcc',
              'command':[sys.executable,*sys.argv],
              'scope':'Prepared independent-H0 quantum bilinear two-point input at engineering scope. Not interacting electron-ion equilibrium, physical temperature, Gaussian higher moments, loss, full spatial/domain covariance or material transfer.',
              'recovery':'Each committed per-channel spectral recipe and correlation is standalone and read back exactly. No automatic CLI restart, old amplitude migration or unfinished dense eigensolve recovery is claimed; any restart needs a new explicit Root allocation.'}
    settings = None
    def event(kind,**values):
        record = {'at':now(),'elapsed_seconds':time.monotonic()-BOOT,'event':kind,**values}
        with (target/'telemetry.jsonl').open('a') as stream:
            stream.write(json.dumps(record,allow_nan=False)+'\n')
            stream.flush()
        print(json.dumps(record,allow_nan=False),flush=True)
    def budget():
        elapsed = time.monotonic()-BOOT
        if elapsed >= settings['actual_preparation_cap_seconds']:
            raise TimeoutError('registered per-generation preparation allocation exhausted')
        if not records and elapsed >= settings['actual_first_preparation_due_seconds']:
            raise TimeoutError('no committed useful preparation within60s')
    def interrupted(signum,_):
        raise InterruptedError(f'signal{signum}; committed recipes retained')
    for signum in (signal.SIGTERM,signal.SIGINT,signal.SIGALRM):
        signal.signal(signum,interrupted)
    try:
        registration, report['registration'] = bound_registration()
        settings = dict(registration['allocation'])
        launch_raw = args.admission.read_bytes()
        launch = json.loads(launch_raw)
        settings['actual_preparation_cap_seconds'] = launch['scientific_cap_seconds']
        signal.setitimer(signal.ITIMER_REAL,max(0.001,settings['actual_preparation_cap_seconds']-(time.monotonic()-BOOT)))
        for record in launch['source_and_oracle_bindings']:
            verify_buffer(Path(record['path']).read_bytes(),record)
        report['launch_admission'] = binding(args.admission,launch_raw)
        report['implementation'] = implementation()
        report['runtime_context'] = runtime_context()
        context = report['runtime_context']
        if (context['cpu_affinity'] != ['6'] or any(value != '1' for value in context['thread_variables'].values())
                or context['resources']['memory.max'].get('value') != launch['memory_max_bytes']
                or context['resources']['memory.swap.max'].get('value') != '0' or not context['invocation']):
            raise ValueError('actual user-manager CPU/thread/memory/invocation contract differs')
        save_json(target/'runtime-context.json',context)
        event('actual_context_admitted',context=context)
        save_json(target/'result.json',report)
        bundles, _, inputs = load_inputs()
        report['actual_inputs'] = inputs
        budget()
        frames_started = time.monotonic()
        models = build_pw_response_models(bundles['gamma'],bundles['half'],bundles['minus'])
        report['frame_build_seconds'] = time.monotonic()-frames_started
        report['reference'] = models['gamma'].metadata
        event('frames_constructed',seconds=report['frame_build_seconds'])
        times = np.asarray(settings['times_tau_ry'],float)
        cache = {}
        for point, model in models.items():
            budget()
            save_json(target/(point+'-coordinates.json'),coordinate_bindings(model,cache))
            points[point] = prepare_point(target,point,model,bundles['gamma'].background.array('bg'),
                                         times,budget,event,records)
            save_json(target/'result.json',{**report,'points':points,'committed_recipes':records})
        budget()
        report['status'] = 'complete'
    except Exception as error:
        report.update(status='stopped_requires_root_decision',failure={'type':type(error).__name__,
                      'message':str(error),'traceback':traceback.format_exc()})
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        report.update(points=points,committed_recipes=records,finished_at=now(),
                      elapsed_seconds=time.monotonic()-BOOT)
        save_json(target/'result.json',report)
        if report['status'] == 'complete':
            (target/'.gate_done_prepared_noise').write_text(report['finished_at']+'\n')
        event('pilot_finished',status=report['status'],elapsed_seconds=report['elapsed_seconds'])
    return 0 if report['status']=='complete' else 1


if __name__=='__main__':
    raise SystemExit(main())
