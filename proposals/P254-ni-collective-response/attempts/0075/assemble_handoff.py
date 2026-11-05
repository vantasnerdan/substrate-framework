"""Package actual global-joint output, consumer and all-exit receipts without replay."""
from pathlib import Path
import hashlib
import json
import os
HERE = Path(__file__).resolve().parent


def binding(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1048576),b''):
            digest.update(chunk)
    return dict(path=str(path),bytes=path.stat().st_size,sha256=digest.hexdigest())


def load(path):
    return json.loads(Path(path).read_bytes())


def main():
    native = load(HERE/'native-findings.json')
    terminal = load(HERE/'receipts-001/exit.json')
    catalog = load(HERE/'component-001/catalog.json')
    consumer = load(HERE/'saved-path-exercise.json')
    complete = (terminal['exitcode']==0 and native['status']=='complete-global-joint-component' and
        consumer['completed_pairs']==961 and len(catalog['coverage']['completed_pairs'])==961)
    record = dict(schema='p254-a75-global-joint-handoff-v1',
        status='complete-global-joint-component' if complete else 'incomplete-global-joint-component',
        scientific_verdict='numeric_evidence' if complete else 'NUMERICALLY_UNRESOLVED',
        component='joint',q_additive=False,full_solve_eligible=False,
        earned_scope=catalog['scope'],
        construction='Exact A62 lines116–125 for every declared ordered pair: common_pair full and ionic; Sigma=fullSigma-ionicSigma; Pi=fullJointPi once.',
        output_binding=binding(HERE/'component-001/catalog.json'),
        output_arrays={kind:catalog[kind] for kind in ('sigma','pi')},
        coverage=dict(actual_complete_pairs=len(catalog['coverage']['completed_pairs']),required_pairs=961,
            times_inverse_ry=catalog['times_inverse_ry'],slots=128,bands=12,spins=2,gamma_coordinates=114,
            persistence=catalog['coverage']['checkpoint']),
        actual_saved_path_consumer=consumer,
        original_input_dataset=catalog['input_dataset'],source_bindings=catalog['source_bindings'],
        coefficient_bindings=catalog['coefficient_bindings'],physical_source_identity=catalog.get('physical_source_identity'),
        canonical_routes=catalog.get('canonical_routes'),source_descriptor=catalog['source_descriptor'],
        carried_owning_implementation=catalog['owning_implementation'],
        reused_implementation_bindings=catalog['component_implementation_bindings'],
        actual_runtime=native.get('runtime'),
        costs=dict(native_elapsed_seconds=native['native_elapsed_seconds'],native_cpu_seconds=native['native_cpu_seconds'],
            peak_rss_bytes=native['peak_rss_bytes'],launch_elapsed_seconds=terminal['elapsed_seconds'],
            production_elapsed_seconds=catalog['elapsed_seconds'],production_cpu_seconds=catalog['native_cpu_seconds'],
            raw_time_output=(HERE/'receipts-001/time.txt').read_text(),
            raw_manager_terminal=load(HERE/'receipts-001/manager-terminal.json')),
        receipts={name:binding(HERE/name) for name in ('preregistration.json','owned-view.json','native-startup.json',
            'native-findings.json','saved-path-exercise.json','receipts-001/launch.json','receipts-001/exit.json',
            'receipts-001/stdout.txt','receipts-001/stderr.txt','receipts-001/time.txt','receipts-001/manager-terminal.json')},
        owned_implementation_bindings=[binding(HERE/name) for name in ('produce_joint.py','component_reader.py','run_native.py','assemble_handoff.py')],
        failure=native.get('failure'),boundary_loss=catalog['boundary_loss'],
        no_replay=dict(old_suites=False,ee=False,raw_ei=False,cubic=False,q_outer_sum=False),
        not_claimed=['complete six-array all64q functional','selfconsistent simulation','trajectory conservation','ranking','rate','heat','merge','rigorous accuracy certification'])
    with (HERE/'handoff.json').open('x') as stream:
        json.dump(record,stream,indent=2,allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    print(json.dumps(dict(status=record['status'],pairs=record['coverage']['actual_complete_pairs'],handoff=binding(HERE/'handoff.json'))))


if __name__=='__main__':
    main()
