"""Collect completed A77 evidence without launching numerical work."""
from pathlib import Path
import json
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'0062'))
from history_banks import binding, save


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=('001', '002'), default='002')
    args = parser.parse_args()
    receipts = HERE / ('receipts-' + args.case)
    exercise = HERE / ('exercise-' + args.case)
    exit_record = json.loads((receipts/'exit.json').read_bytes())
    findings_path = exercise / 'findings.json'
    findings = json.loads(findings_path.read_bytes()) if findings_path.exists() else None
    if exit_record['returncode'] != 0 or findings is None:
        raise ValueError('no successful saved-path exercise; retain all failure evidence')
    evidence = {}
    for name in ('component_accumulator.py','native_exercise.py','launch_native.py','assemble_handoff.py',
                 'preregistration.json','api.json','native-view.json','view-closed.json',
                 'corrected-invocation-decision.json',
                 'receipts-'+args.case+'/launch.json','receipts-'+args.case+'/stdout','receipts-'+args.case+'/stderr',
                 'receipts-'+args.case+'/time.txt','receipts-'+args.case+'/exit.json','receipts-'+args.case+'/manager-terminal.json',
                 'exercise-'+args.case+'/startup.json','exercise-'+args.case+'/findings.json','exercise-'+args.case+'/partial-publication.json',
                 'exercise-'+args.case+'/EI_catalog-admission.json','exercise-'+args.case+'/cubic_catalog-admission.json',
                 'exercise-'+args.case+'/bank/current.json'):
        evidence[name] = binding(HERE/name)
    manager = json.loads((receipts/'manager-terminal.json').read_bytes())
    values = dict(line.split('=',1) for line in manager['stdout'].splitlines() if '=' in line)
    failed = None
    if args.case == '002':
        failure_manager = json.loads((HERE/'receipts-001/manager-terminal.json').read_bytes())
        failure_values = dict(line.split('=',1) for line in failure_manager['stdout'].splitlines() if '=' in line)
        failed = dict(exit=json.loads((HERE/'receipts-001/exit.json').read_bytes()),manager=failure_values,
            reason='transitive archived producer lineage was incorrectly treated as a current live input; failed before any admission')
        for path in (HERE/'receipts-001').iterdir():
            evidence[str(path.relative_to(HERE))] = binding(path)
        evidence['exercise-001/startup.json'] = binding(HERE/'exercise-001/startup.json')
    def measured_cost(directory):
        lines = (directory/'time.txt').read_text().splitlines()
        def value(label):
            return next(line.strip()[len(label):].strip() for line in lines if line.strip().startswith(label))
        elapsed = 0.0
        for part in value('Elapsed (wall clock) time (h:mm:ss or m:ss):').split(':'):
            elapsed = 60*elapsed + float(part)
        return dict(wall_seconds=elapsed, user_seconds=float(value('User time (seconds):')),
                    system_seconds=float(value('System time (seconds):')),
                    max_rss_kib=int(value('Maximum resident set size (kbytes):')),
                    swaps=int(value('Swaps:')), exit=int(value('Exit status:')))
    measured = measured_cost(receipts)
    prior_measured = measured_cost(HERE/'receipts-001') if failed else None
    record = dict(schema='p254-a77-serial-accumulator-terminal-handoff-v1',
        status=findings['status'], module=evidence['component_accumulator.py'], api=evidence['api.json'],
        input_dataset=findings['input_dataset'], original_B_identity=findings['original_B_identity'],
        input_catalogs=findings['input_catalogs'], saved_path_consumer=findings['saved_path_consumer'],
        coverage=findings['coverage'], pointer_reopen=findings['pointer_reopen'],
        identical_duplicate_imports=[dict(component=x['component'],added_cells=x['added_cells'],
            identical_overlap_cells=x['identical_overlap_cells']) for x in findings['duplicates']],
        complete_solve_rejected=findings['complete_solve_rejected'],
        recovery=dict(previous_generation=findings['previous_generation'],
            retained_state_count=findings['retained_state_count'],retained_family_count=findings['retained_family_count'],
            unchanged_EI_family_shared=findings['unchanged_EI_family_shared'],
            scope='atomic pointer reopened after both actual admissions; previous generation authenticated; no injected crash or fabricated complete mesh'),
        cost=dict(elapsed_seconds=findings['elapsed_seconds'],native_cpu_seconds=findings['native_cpu_seconds'],
            peak_rss_bytes=findings['peak_rss_bytes'],manager=values,all_prior_costs_preserved=True,
            automatic_retries=0,native_invocations=1 if failed is None else 2, prior_failed_invocation=failed,
            gnu_time=measured, prior_gnu_time=prior_measured,
            total_gnu_wall_seconds=measured['wall_seconds']+(prior_measured['wall_seconds'] if prior_measured else 0),
            total_gnu_cpu_seconds=measured['user_seconds']+measured['system_seconds']+
                (prior_measured['user_seconds']+prior_measured['system_seconds'] if prior_measured else 0),
            manager_terminal_note='Successful transient unit was already collected when systemctl show ran; preserve raw unset values, not active-limit or memory-peak evidence. Exact unit launch argv, actual startup cgroup/affinity/threads, systemd-run terminal and GNU time provide runtime evidence.'),
        evidence=evidence, exercised_scope='Settled real q1 all961 rawEI Sigma and cubic111 Pi only; no new emitter or pending joint consumed.',
        implemented_not_exercised=['changed overlapping bytes rejection','changed source rejection',
            'non-q1 exact source-time subset admission','next generation copy of an already-existing changed family',
            'interrupted precommit garbage collection','complete all64/all961 final six-array publication'],
        nonclaims=findings['nonclaims'], verification='Only explicitly authorized real saved-path Native exercise; no suite/build/lint/blanket checks.')
    save(HERE/'complete-handoff.json',record)
    print(json.dumps(dict(status=record['status'],handoff=binding(HERE/'complete-handoff.json'),cost=record['cost'])),flush=True)


if __name__ == '__main__':
    main()
