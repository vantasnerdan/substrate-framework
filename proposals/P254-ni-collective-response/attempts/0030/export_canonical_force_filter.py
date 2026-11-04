"""Produce and consume one compact Gamma Ni1/x canonical-force filter.

One NEW necessary approved-model baseline reconstruction, not recovery of the
original dense/adaptive history. No finite-dose or full-PW trajectory is run.
"""
import argparse
import inspect
import json
from pathlib import Path
import resource
from time import perf_counter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original-result', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    from export_full_noise_probe import runtime
    runtime()
    directory = Path(args.output).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    from saved_active_baseline_inputs import prepare_saved_active_baseline
    from canonical_force_filter import build_canonical_force_filter, CanonicalForceFilter
    from clamped_scientific_readout import encode
    from substrate_framework.streamed_material_state import file_binding
    import numpy as np
    import scipy
    from scipy.integrate import DOP853

    started = perf_counter()
    inputs = prepare_saved_active_baseline(args.original_result, work_directory=directory/'work')
    filter_started = perf_counter()
    result = build_canonical_force_filter(inputs.correlations, 0, 30., (inputs.coordinate,))
    filter_seconds = perf_counter()-filter_started
    archive_path = directory/'canonical-force-filter.npz'
    with archive_path.open('xb') as stream:
        np.savez(stream, **result.arrays)
    archive_binding = file_binding(archive_path)
    # Exercise the saved numeric consumer, not only the in-memory producer.
    with np.load(archive_path, allow_pickle=False) as archive:
        restored = CanonicalForceFilter(result.final_time,
            **{name: archive[name] for name in result.arrays}, metadata=result.metadata)
    samples = []
    for time in (0., 1., 10., 30.):
        ionic, field = restored.evaluate(time)
        samples.append(dict(force_insertion_time_inverse_ry=time,
                            ionic_force=ionic, field_momentum_force=field))
    work = dict(inputs.work, filter_seconds=filter_seconds,
                total_elapsed_seconds=perf_counter()-started,
                peak_process_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    packet = dict(schema='P254-canonical-force-filter-material-v1', status='complete',
        source=inputs.source, archive=archive_binding, filter=result.metadata,
        sampled_saved_filter=samples, work=work)
    packet['implementation'] = [file_binding(__file__), file_binding(
        Path(__file__).with_name('canonical_force_filter.py')), file_binding(
        Path(__file__).with_name('saved_active_baseline_inputs.py')),
        file_binding(inspect.getsourcefile(DOP853))]
    packet.update(scipy_version=scipy.__version__,
        observable='Ni1/x orthonormal Gamma Fourier displacement at30 inverseRy',
        forcing='Physical canonical ionic force and field-momentum force, inserted into p-prime/Pi-prime; no additional mass/source/spin/cell factors.',
        history='One NEW necessary unexcited moving-baseline reconstruction from original inputs. Not recovered original adaptive/dense history; no finite-dose replay.',
        scope='Conditional active-model causal forcing filter and full initial dual only. No passive two-time noise/kernel is supplied or assumed independent; no full-PW moving mean, coherence correction, irreversible loss, heat or lifetime.',
        artifact_bytes=archive_binding['bytes'])
    packet_path = directory/'packet.json'
    with packet_path.open('x') as stream:
        json.dump(encode(packet), stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(packet=file_binding(packet_path), archive=archive_binding,
                          filter_seconds=filter_seconds, total_seconds=work['total_elapsed_seconds']),
                     allow_nan=False))


if __name__ == '__main__':
    main()
