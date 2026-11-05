"""One complete common-state Picard step from genuine saved means and kernels.

No producer is invoked here. Incomplete ownership or absent saved nodes fails
before a solve. Production allocation remains an explicit caller responsibility.
"""
from pathlib import Path
import argparse
import json
import resource
import sys
from time import perf_counter
import numpy as np
from history_storage import (HERE, AUXILIARY, binding, atomic_json, authenticate,
                             save_array, finite_blocks)
from history_provider import PreviousHistory, SavedHistory, SCHEMA
from saved_mean import SavedMean, REFERENCE, PI_CONVENTION, POSITIONS
from history_readout import consume, json_value
from production_transport import FinalBank, solve_supplied
sys.path.insert(0, str(HERE.parent/'0062'))
from history_banks import OwnedKernels


def solve_electron_slot(mean, kernels, bank, slot, *, binary, order, beta,
                        source_bindings, previous, timeout):
    source_ret, source_less = mean.electron_slot(slot)
    if not np.array_equal(source_less[0, 0], previous.aux['source_electron_less0'][slot]):
        raise ValueError('original electron preparation changed')
    return solve_supplied(mean.times, source_ret, source_less,
        kernels.sigma['retarded'][:, :, slot], kernels.sigma['lesser'][:, :, slot],
        statistics=-1, order=order, beta=beta, binary=binary, bank=bank, slot=slot,
        source_bindings=source_bindings, timeout=timeout)


def solve_gamma_slot(mean, kernels, bank, *, binary, order, beta,
                     source_bindings, previous, timeout):
    source_ret, source_less = mean.gamma_slot()
    if not np.array_equal(source_less[0, 0], previous.aux['source_gamma_less'][0, 0]):
        raise ValueError('original full Gamma preparation changed')
    return solve_supplied(mean.times, source_ret, source_less,
        kernels.pi['retarded'], kernels.pi['lesser'], statistics=1, order=order,
        beta=beta, binary=binary, bank=bank, slot=0,
        source_bindings=source_bindings, timeout=timeout)


def continue_history(kernels_path, mean_path, previous_path, output, *, binary, order,
                     beta=1., timeout=120):
    started = perf_counter()
    out = Path(output).resolve()
    # Enter/validate all inputs first. Partial kernel catalogs are never upgraded,
    # completed by zero, or interpreted as full common kernels.
    with OwnedKernels(kernels_path, require_complete=True) as kernels, \
            SavedMean(mean_path) as mean, PreviousHistory(previous_path) as previous:
        catalog = kernels.catalog
        if (catalog['input_dataset'] != previous.binding or
                mean.spec['input_dataset'] != previous.binding):
            raise ValueError('kernels, common mean and previous dressed iterate differ')
        if (not np.array_equal(mean.times, kernels.times) or
                not np.array_equal(previous.times, kernels.times) or
                len(mean.times) <= order or order not in range(1, 6) or
                not np.isfinite(beta) or beta <= 0):
            raise ValueError('same exact genuine history and admitted integration order required')
        for item in catalog.get('source_bindings', []):
            authenticate(item)
        if not all(finite_blocks(a) for bank in (kernels.sigma, kernels.pi) for a in bank.values()):
            raise ValueError('complete finite common kernels required')
        if not np.array_equal(1j*mean.rho0, previous.aux['source_electron_less0']):
            raise ValueError('saved common mean changed the original electron preparation')
        if not np.array_equal(-1j*mean.w0[np.ix_(POSITIONS, POSITIONS)].T,
                              previous.aux['source_gamma_less'][0, 0]):
            raise ValueError('saved common mean changed the original Gamma preparation')
        # No history arrays copied: the legacy archive is mapped in-place; new
        # controls/preparation are written once, subsequent steps reference them.
        sources = [binding(kernels_path), binding(mean_path), previous.binding,
                   binding(binary), *mean.bindings,
                   *(binding(p) for p in sorted(HERE.glob('*.py'))),
                   binding(HERE.parent/'0062'/'history_banks.py')]
        out.mkdir(parents=True, exist_ok=False)
        atomic_json(out/'invocation.json', dict(schema='p254-history-continuation-invocation-v1',
            source_bindings=sources, times_inverse_ry=mean.times.tolist(), order=order, beta=beta,
            layout='one canonical heterogeneous final bank: electron128x12, separate Gamma1x114',
            sourcekind='single-physical-Hermitian-contour',
            reference_source_kind=REFERENCE, pi_convention=PI_CONVENTION))
        with FinalBank(out/'electron', mean.times, 12, 128,
                statistics=-1, order=order, beta=beta) as bank:
            for slot in range(128):
                solve_electron_slot(mean, kernels, bank, slot, binary=binary, order=order,
                    beta=beta, source_bindings=sources, previous=previous, timeout=timeout)
            bank.finalize()
        with FinalBank(out/'gamma', mean.times, 114, 1,
                statistics=1, order=order, beta=beta) as bank:
            solve_gamma_slot(mean, kernels, bank, binary=binary, order=order, beta=beta,
                source_bindings=sources, previous=previous, timeout=timeout)
            bank.finalize()
        if previous.provider is not None:
            auxiliary = previous.provider.spec['auxiliary']
        else:
            (out/'auxiliary').mkdir()
            auxiliary = {name: save_array(out/'auxiliary'/(name+'.npy'), previous.aux[name])
                         for name in AUXILIARY}
        provider = dict(schema=SCHEMA, status='complete', times_inverse_ry=mean.times.tolist(),
            electron=dict(complete=binding(out/'electron'/'complete.json')),
            gamma=dict(complete=binding(out/'gamma'/'complete.json')), auxiliary=auxiliary,
            sourcekind='single-physical-Hermitian-contour', integration_order=order, beta=beta,
            reference_source_kind=REFERENCE, pi_convention=PI_CONVENTION,
            input_kernels=binding(kernels_path), input_mean=binding(mean_path),
            input_previous=previous.binding, source_bindings=sources,
            prescribed_nonGamma_bath='literal original bare provider; no dressed other-q update',
            interpretation='One complete supplied-history Picard iterate; not a converged fixed point or a temporal/error/conservation certificate.')
        atomic_json(out/'provider.json', provider)
        with SavedHistory(out/'provider.json') as history:
            readout = consume(history, previous)
        (out/'readout.json').write_text(json.dumps(readout, default=json_value, indent=2)+'\n')
        result = dict(schema='p254-common-history-continuation-result-v1', status='complete',
            provider=binding(out/'provider.json'), readout=binding(out/'readout.json'),
            invocation=binding(out/'invocation.json'), elapsed_seconds=perf_counter()-started,
            owner_maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            native_invocations=129, times_inverse_ry=mean.times.tolist(),
            physical_scope='All128 active12-band sectors, full114 Gamma; original prescribed nonGamma bath.',
            delta_from_previous=readout['delta_from_previous'],
            limits=['Raw iterate differences are not truncation error, convergence or conservation closure.',
                    'No projection/refill, inferred independent companions, physical calibration, loss, heat or nuclear rate.',
                    'No full31 claim unless the actual complete saved31-node mean and common kernels were supplied.'])
        atomic_json(out/'result.json', result)
        return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('kernels', 'mean', 'previous', 'output', 'binary'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--order', type=int, required=True)
    p.add_argument('--beta', type=float, default=1.)
    p.add_argument('--timeout', type=float, default=120.)
    a = p.parse_args()
    result = continue_history(a.kernels, a.mean, a.previous, a.output,
        binary=a.binary, order=a.order, beta=a.beta, timeout=a.timeout)
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
