"""One actual Picard update: dressed common kernels and common mean reference."""
from pathlib import Path
import argparse
import json
import os
import resource
from time import perf_counter
import numpy as np
from feedback_inputs import load_feedback, file_binding
from supplied_keldysh import solve_supplied, load_solution
from joint_gamma_adapter import ordered_coordinates
from consume_material_step import consume, json_value


def run_feedback_step(components, mean_result, previous, baseline_catalog, output, *, order=1):
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = perf_counter()
    print(json.dumps(dict(stage='actual-feedback-step-start', pid=os.getpid(),
        affinity=sorted(os.sched_getaffinity(0)),
        cgroup=Path('/proc/self/cgroup').read_text(),
        threads={key: os.environ.get(key) for key in (
            'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS')})), flush=True)
    c, reference, catalog, baseline, old, sigma, pi = load_feedback(
        components, mean_result, previous, baseline_catalog)
    try:
        times = reference.times
        nt = len(times)
        e_ret = np.zeros((nt, nt, c.ns, c.m, c.m), complex)
        e_less = np.zeros_like(e_ret)
        for a in range(nt):
            for b in range(nt):
                greater, lesser = reference.electron(a, b)
                e_less[a, b] = lesser
                if a >= b:
                    e_ret[a, b] = greater-lesser
        if not np.array_equal(e_less[0, 0], old['source_electron_less0']):
            raise ValueError('original prepared electron state changed')
        binary = baseline['solver_binary']['path']
        for slot in range(c.ns):
            path = solve_supplied(times, e_ret[:, :, slot], e_less[:, :, slot],
                sigma['retarded'][:, :, slot], sigma['lesser'][:, :, slot],
                statistics=-1, order=order, binary=binary,
                workdir=out/'fermions'/('slot-%03d' % slot))
            solved = load_solution(path)
            e_ret[:, :, slot], e_less[:, :, slot] = solved['ret'], solved['less']
        d = c.d+c.f
        b_ret = np.zeros((nt, nt, d, d), complex)
        b_less = np.zeros_like(b_ret)
        for a in range(nt):
            for b in range(nt):
                greater, lesser = ordered_coordinates(reference, a, b)
                b_less[a, b] = lesser
                if a >= b:
                    b_ret[a, b] = greater-lesser
        if not np.array_equal(b_less[0, 0], old['source_gamma_less'][0, 0]):
            raise ValueError('original prepared Gamma state changed')
        path = solve_supplied(times, b_ret, b_less,
            pi['retarded'], pi['lesser'], statistics=1, order=order,
            binary=binary, workdir=out/'intrinsic')
        boson = load_solution(path)
        dataset = out/'actual-material-step.npz'
        # Controls and original-bare readout stay literal. The new mean-dependent
        # bare reference is saved separately, not mislabeled as original control.
        np.savez(dataset, times=times, n=c.n,
            electron_ret=e_ret, electron_less=e_less,
            source_electron_less0=old['source_electron_less0'],
            gamma_ret=boson['ret'], gamma_less=boson['less'],
            control_ret=old['control_ret'], control_less=old['control_less'],
            source_gamma_ret=old['source_gamma_ret'],
            source_gamma_less=old['source_gamma_less'],
            current_bare_gamma_ret=b_ret, current_bare_gamma_less=b_less)
        readout = consume(dataset)
        readout_path = out/'readout.json'
        readout_path.write_text(json.dumps(readout, default=json_value, indent=2)+'\n')
        delta = dict(
            electron_ret_maxabs=float(np.max(np.abs(e_ret-old['electron_ret']))),
            electron_less_maxabs=float(np.max(np.abs(e_less-old['electron_less']))),
            gamma_ret_raw_coordinate_maxabs=float(np.max(np.abs(boson['ret']-old['gamma_ret']))),
            gamma_less_raw_coordinate_maxabs=float(np.max(np.abs(boson['less']-old['gamma_less']))),
            interpretation='Raw iterate changes only; mixed Gamma coordinate units are not an accuracy norm or convergence criterion.')
        result = dict(schema='p254-actual-common-feedback-material-step-v1', status='complete',
            input_components=file_binding(components), input_mean=file_binding(mean_result),
            input_previous=file_binding(previous), input_baseline=file_binding(baseline_catalog),
            dataset=file_binding(dataset), readout=file_binding(readout_path),
            source_bindings=[file_binding(Path(__file__)), file_binding(Path(__file__).with_name('feedback_inputs.py'))],
            mean_reference_bindings=reference.bindings,
            times_inverse_ry=times.tolist(), integration_order=order,
            actual_fermion_blocks=c.ns, actual_Gamma_coordinates=d,
            actual_existing_solver_invocations=c.ns+1,
            controls='Original represented-only and original-bare saved arrays retained literally; no control solves.',
            reference_source_kind=catalog['reference_source_kind'],
            pi_convention=catalog['pi_convention'],
            mean_feedback_included=True, kernel_feedback_from_previous_state=True,
            selfconsistent=False, temporal_accuracy_established=False,
            delta_from_previous=delta,
            elapsed_seconds=perf_counter()-started,
            maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            limits=['One Picard update is not a converged fixed point.',
                    'Supplied prefix and prescribed non-Gamma bath only, not full0..30/all-q evolution.',
                    'No physical calibration/XC matching, irreversible loss, heat, nuclear rate or whole-E closure.'])
        result_path = out/'result.json'
        result_path.write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(dict(stage='actual-feedback-step-complete', result=file_binding(result_path),
            readout=readout, delta_from_previous=delta), default=json_value), flush=True)
        return result
    finally:
        reference.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('components', 'mean-result', 'previous', 'baseline-catalog', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--order', type=int, default=1)
    args = parser.parse_args()
    run_feedback_step(args.components, args.mean_result, args.previous,
        args.baseline_catalog, args.output, order=args.order)


if __name__ == '__main__':
    main()
