#!/usr/bin/env python3
"""Exercise saved boundary admission; never substitute finite evidence for Native.

Creates one genuine finite coupled-system RK4 fixture with actual QWork. Checks
its saved four-stage lineage, same-time evaluation/final distinction, and actual
metadata/byte mutations. Separately admits the ORIGINAL supplied Native t0 seed
and refuses an ORIGINAL supplied Native evaluation node. It does not propagate
Native physics or change producer files. Run only with explicit parent approval.

Usage: python run_mixed_energy_boundary_oracle.py --work-directory NEW_DIRECTORY
  --producer DIRECTORY --producer-owner-sha256 SHA --evaluation SAVED_NODE_JSON
  --producer-source-archive ZIP --producer-source-archive-sha256 SHA
  --cpus 6 7 8 --max-seconds 120

Saves result.json and prints the same report only after all checks complete.
Finite fixtures prove boundary/source-replay machinery, not material energy.
"""
import argparse
from contextlib import ExitStack, contextmanager
import json
import os
from pathlib import Path
import sys
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'

from mixed_energy_saved_boundary import _Authority, _stage, admit_saved_boundary
from mixed_energy_source_archive import bank_producer_sources, read_producer_sources
from mixed_finite_coupled_system import CoupledFixture, FiniteDynamics
from mixed_response_rk4 import advance_step
from mixed_response_storage import COMPONENTS, MappedMixedState
from substrate_framework.partitioned_q_storage import QBudget
from substrate_framework.partitioned_q_work import QWork
from substrate_framework.streamed_material_state import file_binding


def require(condition, name, checks):
    if not condition:
        raise AssertionError(name)
    checks.append(name)


def refuses(action, name, checks):
    try:
        action()
    except (ValueError, FileNotFoundError):
        checks.append(name)
        return
    raise AssertionError(name+' was not refused')


@contextmanager
def changed_json(path, change):
    """Mutation scope is ONLY newly generated finite fixture metadata."""
    original = path.read_bytes()
    record = json.loads(original)
    change(record)
    path.write_text(json.dumps(record, allow_nan=False)+'\n')
    try:
        yield
    finally:
        path.write_bytes(original)


def inspect_finite_step(directory, initial, dt, identity):
    authority = _Authority()
    step_binding = authority.at(directory/'request.json')
    seed = authority.read(initial)
    all_outputs = []
    with ExitStack() as stack:
        base = stack.enter_context(MappedMixedState(initial))
        current, accumulator = base, None
        for stage in range(4):
            outputs = _stage(directory/f'stage-{stage}', step_binding=step_binding,
                stage=stage, base=base, current=current, accumulator=accumulator,
                dt=dt, identity=identity, seed=seed, authority=authority)
            all_outputs.append(outputs)
            if stage < 3:
                current = stack.enter_context(MappedMixedState(outputs['input']))
                accumulator = stack.enter_context(MappedMixedState(outputs['accumulator']))
    return all_outputs


def finite_checks(root, cpus, max_seconds, checks):
    fixture = CoupledFixture()
    budget = QBudget(max_actions=100000, max_seconds=max_seconds, max_bytes=512 << 20,
                     max_memory_bytes=3 << 30, worker_memory_bytes=1 << 30)
    sources = [file_binding(HERE/name) for name in (
        'mixed_finite_coupled_system.py', 'mixed_response_rk4.py',
        'mixed_response_edge_transition.py', 'mixed_response_storage.py',
        'mixed_energy_saved_boundary.py')]
    initial = fixture.initial(root/'initial', budget, sources)
    seed = json.loads(Path(initial['path']).read_bytes())
    dynamics = FiniteDynamics(fixture.system, root/'contexts', budget, sources)
    work = QWork(root/'qwork', budget, cpus, 'saved-boundary-genuine-finite-oracle')
    step = root/'step'
    dt = fixture.stop-fixture.start
    final = advance_step(initial, dt, dynamics=dynamics, work=work,
                         directory=step, identity=seed['identity'], budget=budget)
    outputs = inspect_finite_step(step, initial, dt, seed['identity'])
    require(outputs[-1]['final'] == final, 'genuine_four_stage_source_complete_endpoint', checks)
    evaluation = json.loads(Path(outputs[2]['input']['path']).read_bytes())
    physical = json.loads(Path(final['path']).read_bytes())
    require(evaluation['time'] == physical['time']
            and evaluation['kind'] == 'mixed-rk4-evaluation-node'
            and evaluation['provenance']['physical_endpoint'] is False
            and physical['kind'] == 'mixed-rk4-physical-endpoint'
            and physical['provenance']['physical_endpoint'] is True,
            'same_time_does_not_make_evaluation_node_physical', checks)
    require(not (root/'step-receipts').exists(), 'no_post_energy_driver_receipt_needed', checks)
    with MappedMixedState(final) as state:
        for index in range(len(state._shapes['mixed'])):
            for component in COMPONENTS:
                state.leg(component, index)
            state.release_edge(index)
    checks.append('genuine_endpoint_full_payloads_readable')
    inspect = lambda: inspect_finite_step(step, initial, dt, seed['identity'])
    with changed_json(step/'stage-3'/'complete.json', lambda r: r['transitions'].pop()):
        refuses(inspect, 'reject_partial_source_prefix', checks)
    with changed_json(step/'stage-3'/'complete.json', lambda r: r.update(complete=False)):
        refuses(inspect, 'reject_unfinished_stage3', checks)
    with changed_json(step/'stage-3'/'plan.json', lambda r: r.update(stage=2)):
        refuses(inspect, 'reject_changed_stage_plan_binding', checks)
    with changed_json(Path(final['path']), lambda r: r['small']['mixed']['p'][0].__setitem__(0, 999.)):
        refuses(inspect, 'reject_changed_endpoint_bytes', checks)
    payload = next(entry['path'] for entries in physical['legs'].values()
                   for entry in entries if 'path' in entry)
    payload_path = Path(payload)
    with payload_path.open('r+b') as stream:
        stream.seek(-1, os.SEEK_END)
        original = stream.read(1)
        stream.seek(-1, os.SEEK_END)
        stream.write(bytes([original[0] ^ 1]))
    try:
        def read_final():
            with MappedMixedState(final) as state:
                for index in range(len(state._shapes['mixed'])):
                    for component in COMPONENTS:
                        state.leg(component, index)
                    state.release_edge(index)
        refuses(read_final, 'reject_changed_original_payload_bytes', checks)
    finally:
        with payload_path.open('r+b') as stream:
            stream.seek(-1, os.SEEK_END)
            stream.write(original)
    return dict(initial=initial, final=final, observed_budget=budget.snapshot())


def archive_checks(root, native_admission, initial, archive, checks):
    restored = admit_saved_boundary(native_admission['producer_directory'],
        native_admission['owner_binding']['sha256'], initial, source_archive=archive)
    require(restored['owner'] == native_admission['owner']
            and restored['state'] == native_admission['state']
            and restored['source_archive'] == archive,
            'actual_native_archive_keeps_original_typed_owner_state', checks)
    fixture = root/'source-replay'
    fixture.mkdir()
    original = dict(opaque_integer=2**63+37, value=1)
    source = fixture/'parameters.json'
    source.write_text(json.dumps(original)+'\n')
    source_binding = file_binding(source)
    owner = dict(request=dict(schema='P254-native-mixed-q-material-v1-owner',
                              source_bindings=[source_binding]), started_unix=1)
    owner_path = fixture/'owner.json'
    owner_path.write_text(json.dumps(owner)+'\n')
    owner_binding = file_binding(owner_path)
    retained = bank_producer_sources(owner_binding, fixture/'retained.zip')
    source.write_text(json.dumps(dict(original, value=2))+'\n')
    refuses(lambda: _Authority().read(source_binding),
            'changed_live_source_without_explicit_archive_refused', checks)
    authority = _Authority()
    authority.recover_sources(retained, owner_binding, [source_binding])
    recovered = authority.read(source_binding)
    require(recovered == original and type(recovered['opaque_integer']) is int
            and json.loads(source.read_bytes()) != recovered,
            'historical_source_replay_reads_original_bytes_and_integer_not_changed_live_file', checks)
    wrong_owner = dict(owner_binding, sha256='0'*64)
    refuses(lambda: read_producer_sources(retained, wrong_owner, [source_binding]),
            'archive_for_another_owner_refused', checks)
    tampered = fixture/'tampered.zip'
    with ZipFile(retained['path']) as saved, ZipFile(tampered, 'x') as changed:
        for name in saved.namelist():
            changed.writestr(name, source.read_bytes() if name == 'sources/0000' else saved.read(name))
    refuses(lambda: read_producer_sources(file_binding(tampered), owner_binding, [source_binding]),
            'valid_outer_archive_digest_does_not_hide_changed_original_member', checks)
    return dict(actual_native_archive=archive, fixture_owner=owner_binding,
                retained_fixture=retained, changed_live_source=file_binding(source))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--work-directory', type=Path, required=True)
    parser.add_argument('--producer', type=Path, required=True)
    parser.add_argument('--producer-owner-sha256', required=True)
    parser.add_argument('--producer-source-archive', type=Path, required=True)
    parser.add_argument('--producer-source-archive-sha256', required=True)
    parser.add_argument('--evaluation', type=Path, required=True)
    parser.add_argument('--cpus', type=int, nargs='+', choices=range(6, 11), required=True)
    parser.add_argument('--max-seconds', type=float, required=True)
    args = parser.parse_args()
    if len(set(args.cpus)) != len(args.cpus) or not set(args.cpus) <= os.sched_getaffinity(0):
        raise ValueError('unique allocated CPUs6-10 required')
    root = args.work_directory.resolve()
    producer = args.producer.resolve()
    if root == producer or root.is_relative_to(producer):
        raise ValueError('oracle work must not be inside the original producer')
    root.mkdir(parents=True, exist_ok=False)
    checks = []
    initial = file_binding(producer/'initial'/'checkpoint.json')
    admission = admit_saved_boundary(producer, args.producer_owner_sha256, initial)
    require(admission['boundary_kind'] == 'initial'
            and admission['owner'] == json.loads((producer/'owner.json').read_bytes())
            and admission['state'] == json.loads(Path(initial['path']).read_bytes()),
            'actual_native_t0_original_typed_owner_state', checks)
    evaluation_binding = file_binding(args.evaluation)
    evaluation = json.loads(args.evaluation.read_bytes())
    require(args.evaluation.resolve().is_relative_to(producer)
            and evaluation['kind'] == 'mixed-rk4-evaluation-node',
            'actual_original_native_evaluation_metadata', checks)
    refuses(lambda: admit_saved_boundary(producer, args.producer_owner_sha256, evaluation_binding),
            'reject_actual_native_evaluation_node', checks)
    refuses(lambda: admit_saved_boundary(producer, '0'*64, initial), 'reject_changed_owner_digest', checks)
    refuses(lambda: admit_saved_boundary(producer, args.producer_owner_sha256,
            dict(initial, sha256='0'*64)), 'reject_changed_state_digest', checks)
    finite = finite_checks(root, tuple(args.cpus), args.max_seconds, checks)
    archive = file_binding(args.producer_source_archive.resolve())
    require(archive['sha256'] == args.producer_source_archive_sha256,
            'explicit_native_archive_digest_matches_actual_bytes', checks)
    retained_sources = archive_checks(root, admission, initial, archive, checks)
    report = dict(schema='P254-saved-energy-boundary-oracle-v2', status='passed',
        checks=checks, finite=finite, retained_sources=retained_sources,
        native_initial=initial, native_evaluation=evaluation_binding,
        scope='generic saved RK4/source-replay and actual Native metadata; no Native propagation or material energy result')
    data = json.dumps(report, indent=2, allow_nan=False)+'\n'
    (root/'result.json').write_text(data)
    print(data, end='')


if __name__ == '__main__':
    main()
