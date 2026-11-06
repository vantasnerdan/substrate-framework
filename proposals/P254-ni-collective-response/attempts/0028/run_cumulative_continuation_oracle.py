#!/usr/bin/env python3
"""Exercise cumulative allocation and genuine interrupted finite RK4 recovery.

This is a finite scientific/ownership oracle, never a material or production grant.
The producer uses the original QWork, RK4, storage and finite coupled provider.
Only this oracle's explicit interruption hook kills its own finite subprocess.
"""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from time import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SUPERVISOR = Path('/home/administrator/projects/catalyst-collective/runs/2026-10-19_collective-response/owner-frontiers/p254-original-source-continuation-20261031')
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(SUPERVISOR)]
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
sys.dont_write_bytecode = True

import numpy as np
import sympy as sp
from continuation_allocation import admit_allocation, create_budget
from mixed_finite_coupled_system import CoupledFixture, FiniteDynamics
from mixed_material_driver_support import (
    ContinuationBudget, exact_publish, owned_directory, owned_tree, source_closure)
from mixed_response_rk4 import advance_step
from mixed_response_storage import COMPONENTS, SMALL_NAMES, MappedMixedState
from substrate_framework.partitioned_q_work import QWork
from substrate_framework.streamed_material_state import authenticate, file_binding


def read(binding):
    return json.loads(authenticate(binding).read_bytes())


def state_values(state):
    result = {c+'-small-'+name: np.asarray(getattr(state.small[c], name)).copy()
              for c in COMPONENTS for name in SMALL_NAMES}
    for c in COMPONENTS:
        for index in range(len(state._shapes[c])):
            result[c+'-edge-'+str(index)] = np.asarray(state.leg(c, index)).copy()
            state.release_edge(index)
    return result


class InterruptedFiniteDynamics(FiniteDynamics):
    """Identical physical provider, plus one explicit finite-only kill point."""
    def prepare(self, state, work):
        prepared = super().prepare(state, work)
        if self.interrupt and state.time == self.final_time:
            root = self.directory.parent
            np.savez(root/'internal-evaluation.npz', **state_values(state))
            stages = [file_binding(root/'steps/000000'/f'stage-{i}'/'complete.json')
                      for i in range(3)]
            receipts = [file_binding(p) for p in sorted((root/'scientific-work').glob(
                '*/generation-*/item-*/complete.json'))]
            record = dict(schema='P254-finite-actual-interruption-v1', finite_only=True,
                actual_utc=datetime.now(timezone.utc).isoformat(),
                prepared=prepared.binding, evaluation=state.binding,
                stages=stages, retained_work_receipts=receipts,
                budget=self.budget.snapshot(), ledger=file_binding(self.budget.ledger))
            exact_publish(root/'interrupted.json', record, None)
            os.kill(os.getpid(), signal.SIGKILL)
        return prepared


def original_budget(directory, owner):
    caps = owner['request']['caps']
    return ContinuationBudget(directory, started_at=owner['started_unix'],
        max_actions=caps['max_actions'], max_seconds=caps['max_seconds'],
        max_bytes=caps['max_output_bytes'], max_memory_bytes=caps['max_memory_bytes'],
        worker_memory_bytes=caps['worker_memory_bytes'],
        existing_bytes=sum(owned_tree(directory).values()))


def owner_for(sources, started, seconds):
    return dict(finite_only=True, request=dict(schema='P254-native-mixed-q-material-v1-owner',
        caps=dict(max_actions=4096, max_seconds=seconds, max_output_bytes=512 << 20,
                  max_memory_bytes=3 << 30, worker_memory_bytes=1 << 30),
        cpus=[6], source_bindings=sources), started_unix=started,
        started_at=datetime.fromtimestamp(started, timezone.utc).isoformat())


def allocation_for(directory, owner, additional=180.):
    maximum = owner['request']['caps']['max_seconds']+additional
    return dict(schema='P254-original-source-cumulative-allocation-v1', authorized=True,
        scope='same-original-source-physical-step-continuation',
        original_owner=file_binding(directory/'owner.json'),
        original_started_unix=owner['started_unix'], original_caps=owner['request']['caps'],
        cpus=owner['request']['cpus'], additional_seconds=additional,
        cumulative_max_seconds=maximum,
        cumulative_deadline_unix=owner['started_unix']+maximum,
        issued_at=datetime.now(timezone.utc).isoformat(),
        decision='Root210-finite-oracle-only; NOT Native production authority')


def finite_producer(directory, interrupt, resume):
    fixture = CoupledFixture()
    with owned_directory(directory, resume) as directory:
        if resume:
            owner = json.loads((directory/'owner.json').read_bytes())
            admission = admit_allocation(directory, file_binding(directory/'owner.json')['sha256'],
                                         file_binding(directory/'allocation.json'))
            budget = create_budget(admission)
            initial = read(file_binding(directory/'seed.json'))['initial']
        else:
            sources = source_closure([file_binding(Path(__file__).resolve()),
                                      file_binding(SUPERVISOR/'continuation_allocation.py')])
            owner = owner_for(sources, time(), 180.)
            exact_publish(directory/'owner.json', owner, None)
            budget = original_budget(directory, owner)
            initial = fixture.initial(directory/'initial', budget, sources)
            exact_publish(directory/'seed.json', dict(initial=initial), None)
        sources = owner['request']['source_bindings']
        before = budget.actions
        work = QWork(directory/'scientific-work', budget, (6,), 'finite-cumulative-original-work')
        dynamics = InterruptedFiniteDynamics(fixture.system, directory/'contexts', budget, sources)
        dynamics.interrupt, dynamics.final_time = interrupt, fixture.stop
        final = advance_step(initial, fixture.stop-fixture.start, dynamics=dynamics, work=work,
            directory=directory/'steps/000000', identity='finite-cumulative-original-step', budget=budget)
        record = dict(schema='P254-finite-cumulative-producer-result-v1', finite_only=True,
            final=final, original_owner=file_binding(directory/'owner.json'),
            budget=budget.snapshot(), actions_before=before, actions_after=budget.actions)
        exact_publish(directory/('resumed.json' if resume else 'complete.json'), record, None)
    print(json.dumps(dict(finite_only=True, final=final, actions=budget.actions)))


def require(checks, name, passed, **details):
    checks.append(dict(name=name, passed=bool(passed), **details))
    if not passed:
        raise AssertionError(name)


def invoke(root, label, arguments, expected):
    argv = [sys.executable, str(Path(__file__).resolve()), '--produce-finite', *arguments]
    run = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90)
    record = dict(argv=argv, exit=run.returncode,
                  stdout=run.stdout.decode(), stderr=run.stderr.decode())
    exact_publish(root/(label+'-command.json'), record, None)
    if run.returncode != expected:
        raise AssertionError(label+': '+json.dumps(record))
    return record


def protocol_checks(root, sources, checks):
    directory = root/'expired-original-clock'
    with owned_directory(directory, False) as directory:
        owner = owner_for(sources, time()-2., 1.)
        exact_publish(directory/'owner.json', owner, None)
        # A protocol-only interrupted-admission fixture, not measured H0 work.
        ledger = b''.join(ContinuationBudget.RECORD.pack(i, 3*i) for i in range(1, 4))
        (directory/'action-admissions.bin').write_bytes(ledger)
        (directory/'action-admissions.lock').touch()
        grant = allocation_for(directory, owner, additional=60.)
        exact_publish(directory/'allocation.json', grant, None)
        expired = False
        try:
            original_budget(directory, owner).check()
        except RuntimeError:
            expired = True
        require(checks, 'original_wall_allocation_really_expired', expired)
        admitted = admit_allocation(directory, file_binding(directory/'owner.json')['sha256'],
                                    file_binding(directory/'allocation.json'))
        budget = create_budget(admitted)
        budget.check()
        snapshot = budget.snapshot()
        require(checks, 'explicit_extension_preserves_original_clock_and_admissions',
            budget.started_at == owner['started_unix'] and budget.actions == 3 and budget.columns == 9
            and snapshot['elapsed_seconds'] >= 2.
            and snapshot['original_allocation']['deadline_unix'] == owner['started_unix']+1.
            and snapshot['cumulative_allocation']['deadline_unix'] == owner['started_unix']+61.,
            snapshot=snapshot)
        require(checks, 'cumulative_duration_exact_SymPy_oracle',
            sp.Rational(str(grant['cumulative_max_seconds'])) ==
            sp.Rational(str(owner['request']['caps']['max_seconds']))+
            sp.Rational(str(grant['additional_seconds'])))
        altered = deepcopy(grant)
        altered['original_caps']['max_actions'] += 1
        exact_publish(directory/'drifted-allocation.json', altered, None)
        refused = False
        try:
            admit_allocation(directory, file_binding(directory/'owner.json')['sha256'],
                              file_binding(directory/'drifted-allocation.json'))
        except ValueError:
            refused = True
        require(checks, 'non_time_cap_drift_refused', refused)
        # Bind a well-formed changed owner/grant, but a wrong actual source SHA.
        changed = deepcopy(owner)
        changed['request']['source_bindings'][0]['sha256'] = '0'*64
        (directory/'owner.json').write_text(json.dumps(changed, indent=2)+'\n')
        changed_grant = allocation_for(directory, changed, additional=60.)
        exact_publish(directory/'source-drift-allocation.json', changed_grant, None)
        refused = False
        try:
            admit_allocation(directory, file_binding(directory/'owner.json')['sha256'],
                              file_binding(directory/'source-drift-allocation.json'))
        except ValueError as error:
            refused = 'authenticated authority/source bytes differ' in str(error)
        require(checks, 'independent_numerical_source_digest_refused', refused)


def run(root, checks):
    root.mkdir(parents=True, exist_ok=False)
    resumed, baseline = root/'interrupted', root/'uninterrupted'
    commands = [invoke(root, 'interrupted', [str(resumed), '--interrupt'], -signal.SIGKILL)]
    partial = read(file_binding(resumed/'interrupted.json'))
    owner_binding = file_binding(resumed/'owner.json')
    owner = read(owner_binding)
    require(checks, 'three_real_completed_RK_stages_before_actual_SIGKILL',
        all(read(binding)['complete'] is True for binding in partial['stages'])
        and [read(binding)['stage'] for binding in partial['stages']] == [0, 1, 2]
        and not (resumed/'steps/000000/stage-3/complete.json').exists(),
        original_interruption=file_binding(resumed/'interrupted.json'))
    exact_publish(resumed/'allocation.json', allocation_for(resumed, owner), None)
    commands.append(invoke(root, 'resumed', [str(resumed), '--resume'], 0))
    commands.append(invoke(root, 'uninterrupted', [str(baseline)], 0))
    result, reference = read(file_binding(resumed/'resumed.json')), read(file_binding(baseline/'complete.json'))
    require(checks, 'original_owner_and_completed_stage_receipts_remain_byte_identical',
        file_binding(resumed/'owner.json') == owner_binding
        and all(file_binding(b['path']) == b for b in partial['stages']))
    require(checks, 'all_pre_interruption_QWork_receipts_reused_unchanged',
        all(file_binding(b['path']) == b for b in partial['retained_work_receipts']),
        retained_receipts=partial['retained_work_receipts'])
    require(checks, 'no_H0_replay_and_no_cumulative_action_reset',
        result['actions_before'] == partial['budget']['hamiltonian_actions'] > 0
        and result['actions_after'] == reference['actions_after']
        and result['actions_after'] > result['actions_before'],
        before=result['actions_before'], after=result['actions_after'],
        uninterrupted_actions=reference['actions_after'])
    metadata = read(result['final'])
    require(checks, 'fourth_stage_returns_genuine_physical_endpoint',
        metadata['kind'] == 'mixed-rk4-physical-endpoint'
        and read(file_binding(resumed/'steps/000000/stage-3/complete.json'))['complete'] is True)
    budget = original_budget(baseline, read(file_binding(baseline/'owner.json')))
    with MappedMixedState(result['final'], budget=budget) as state:
        actual = state_values(state)
    with MappedMixedState(reference['final'], budget=budget) as state:
        expected = state_values(state)
    require(checks, 'every_raw_electronic_and_canonical_array_bitwise_matches_uninterrupted',
        actual.keys() == expected.keys() and all(actual[k].dtype == expected[k].dtype
            and actual[k].shape == expected[k].shape
            and actual[k].tobytes() == expected[k].tobytes() for k in actual),
        compared_arrays=sorted(actual), final=result['final'], uninterrupted_final=reference['final'])
    with np.load(resumed/'internal-evaluation.npz', allow_pickle=False) as internal:
        defect = max(float(np.max(np.abs(actual[k]-internal[k]), initial=0)) for k in actual)
    require(checks, 'same_time_internal_evaluation_is_not_a_physical_endpoint_substitute',
        defect > 0 and read(partial['evaluation'])['time'] == metadata['time'],
        max_abs_array_defect=defect, internal=partial['evaluation'])
    protocol_checks(root, owner['request']['source_bindings'], checks)
    return dict(commands=commands, interrupted=partial, resumed=result, uninterrupted=reference)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--produce-finite', type=Path)
    parser.add_argument('--interrupt', action='store_true')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--work-directory', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.produce_finite is not None:
        finite_producer(args.produce_finite.resolve(), args.interrupt, args.resume)
        return 0
    if args.work_directory is None or args.output is None:
        parser.error('--work-directory and --output are required for the oracle')
    checks = []
    result = dict(schema='P254-finite-cumulative-continuation-oracle-v1', status='FAIL',
        finite_only=True, material_result=False, Native_production_allocation_issued=False,
        source_bindings=source_closure([file_binding(Path(__file__).resolve()),
                                       file_binding(SUPERVISOR/'continuation_allocation.py')]))
    try:
        result.update(run(args.work_directory.resolve(), checks))
        result['status'] = 'PASS'
    except BaseException as error:
        result['error'] = type(error).__name__+': '+str(error)
    finally:
        result.update(checks=checks, terminal_at=datetime.now(timezone.utc).isoformat())
        exact_publish(args.output.resolve(), result, None)
        print(json.dumps(dict(status=result['status'], checks=len(checks),
                              output=str(args.output), error=result.get('error'))))
    return 0 if result['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
