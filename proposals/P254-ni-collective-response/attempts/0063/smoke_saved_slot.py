"""One actual paid step003 slot through the A63 source and A61 native boundary.

This does not promote the old kernel archive to the new full-coverage contract
and does not claim a complete128-slot provider or a full31 physical solve.
"""
from pathlib import Path
import json
import numpy as np
from history_storage import StoredNPZ, authenticate, binding, atomic_json
from saved_mean import SavedMean
from history_provider import PreviousHistory, close_solution
from continue_history import solve_electron_slot
from production_transport import FinalBank, load_solution

PAID = Path('/var/tmp/p254-a58-common-feedback-step-003')
COMPONENTS = Path('/var/tmp/p254-a58-dressed-feedback-002/iteration-catalog.json')
MEAN = Path('/var/tmp/p254-a58-cached-mean-003/reference-result.json')
PREVIOUS = Path('/var/tmp/p254-a58-common-feedback-step-002/actual-material-step.npz')


class PaidKernelBoundary:
    """Read only the real paid kernel bytes, not a synthetic/zero partial bank."""
    def __init__(self, catalog):
        self.spec = json.loads(Path(catalog).read_bytes())
        self.archive = StoredNPZ(authenticate(self.spec['archive']))
        self.sigma = {key: self.archive['sigma_'+key] for key in ('retarded', 'lesser')}

    def close(self):
        self.archive.close()


def exercise(binary, root):
    saved_dir = PAID/'fermions'/'slot-000'
    manifest = json.loads((saved_dir/'manifest.json').read_bytes())
    for name in ('input', 'transport', 'output', 'solution'):
        authenticate(manifest[name])
    source_bindings = [binding(saved_dir/'manifest.json'), binding(PAID/'result.json'),
        binding(COMPONENTS), binding(MEAN), binding(PREVIOUS),
        *(manifest[name] for name in ('input', 'transport', 'output', 'solution'))]
    kernels = PaidKernelBoundary(COMPONENTS)
    try:
        with SavedMean(MEAN) as mean, PreviousHistory(PREVIOUS) as previous, \
                StoredNPZ(manifest['input']['path']) as paid_input:
            if (kernels.spec['input_dataset'] != previous.binding or
                    mean.spec['input_dataset'] != previous.binding or
                    not np.array_equal(mean.times, paid_input['times'])):
                raise AssertionError('paid source/mean/kernel previous-state mismatch')
            source_ret, source_less = mean.electron_slot(0)
            for label, observed, expected in (
                    ('source_ret', source_ret, paid_input['source_ret']),
                    ('source_less', source_less, paid_input['source_less']),
                    ('kernel_ret', kernels.sigma['retarded'][:, :, 0], paid_input['kernel_ret']),
                    ('kernel_less', kernels.sigma['lesser'][:, :, 0], paid_input['kernel_less'])):
                if not np.array_equal(observed, expected):
                    raise AssertionError('saved genuine '+label+' not identical')
            del source_ret, source_less
            order, beta = int(paid_input['order']), float(paid_input['beta'])
            with FinalBank(root/'electron-slot000', mean.times, 12, 1,
                    statistics=-1, order=order, beta=beta) as bank:
                receipt = solve_electron_slot(mean, kernels, bank, 0, binary=binary,
                    order=order, beta=beta, source_bindings=[*source_bindings, *mean.bindings],
                    previous=previous, timeout=180)
                bank.finalize()
        saved = load_solution(root/'electron-slot000')
        try:
            with StoredNPZ(manifest['solution']['path']) as paid_output:
                ret_equal = np.array_equal(saved['ret'][0], paid_output['ret'])
                less_equal = np.array_equal(saved['less'][0], paid_output['less'])
                ret_max = float(np.abs(saved['ret'][0]-paid_output['ret']).max())
                less_max = float(np.abs(saved['less'][0]-paid_output['less']).max())
                if not ret_equal or not less_equal:
                    raise AssertionError('A61 saved paid solution differs')
            output_bindings = {name: binding(root/'electron-slot000'/name)
                               for name in ('bank.json', 'ret.npy', 'less.npy', 'complete.json')}
            findings = dict(schema='p254-a63-actual-saved-slot-findings-v1', status='complete',
                scope='One original step003 physical electron slot000, genuine times0,.5,1/order2. No material reload, fullGamma solve or complete new kernel/provider claim.',
                source_arrays_bitwise_equal_to_paid_input=True,
                kernel_arrays_bitwise_equal_to_paid_input=True,
                retarded_bitwise_equal_to_paid_solution=ret_equal,
                lesser_bitwise_equal_to_paid_solution=less_equal,
                retarded_maxabs_difference=ret_max, lesser_maxabs_difference=less_max,
                shape=list(saved['ret'].shape), times_inverse_ry=saved['times'].tolist(),
                statistics=saved['statistics'], order=saved['order'], beta=saved['beta'],
                original_preparation_unchanged=True,
                source_bindings=source_bindings, mean_bindings=mean.bindings,
                paid_common_kernel_archive=kernels.spec['archive'], output_bindings=output_bindings,
                a61_slot_receipt=binding(root/'electron-slot000'/'receipts'/'slot-000.json'),
                native_receipt=binding(root/'electron-slot000'/'.staging'/'slot-000'/'native.json'),
                native_time=binding(root/'electron-slot000'/'.staging'/'slot-000'/'native.time'),
                released_staging=binding(root/'electron-slot000'/'.staging'/'slot-000'/'released.json'),
                packed_result_bytes=receipt['provenance']['output']['bytes'],
                interpretation='Transport/source identity only. No speedup, convergence, temporal accuracy, conservation, physical ranking or rate claim.')
            atomic_json(root/'actual-findings.json', findings)
            return findings
        finally:
            close_solution(saved)
    finally:
        kernels.close()
