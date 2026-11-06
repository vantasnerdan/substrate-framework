"""Assemble original persisted evidence after Native terminal; no scientific replay."""
import hashlib
import json
from pathlib import Path
import time
HERE=Path(__file__).resolve().parent
OUT=Path('/var/tmp/p254-a45-moving-reference-001')
RECEIPTS=Path('/var/tmp/p254-a45-native-receipts-001')
def binding(path):
    path=Path(path);digest=hashlib.sha256()
    with path.open('rb') as stream:
        while block:=stream.read(1024*1024):digest.update(block)
    return dict(path=str(path),bytes=path.stat().st_size,sha256=digest.hexdigest())
def load(path):return json.loads(Path(path).read_bytes())
def save(path,value):
    data=(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n').encode()
    with Path(path).open('xb') as stream:stream.write(data)
    return binding(path)
wave=load(OUT/'wave-result.json');terminal=load(RECEIPTS/'actual-native-exit.json')
if wave['status']!='complete' or terminal['exitcode']!=0:
    raise RuntimeError('complete Native wave required; retain failure receipts and repair owning cause')
producer=load(OUT/'producer/result.json');kernel=load(OUT/'moving-EI/result.json')
oracle=load(OUT/'oracle/result.json');recovery=load(OUT/'kernel-recovery.json')
contract=load(HERE/'implementation-contract-001.json')
proof_path=HERE.parent/'0041/owning-electron-ion-functional-oracle-001.json'
proof=load(proof_path)
functional=binding(HERE.parent/'0041/electron_ion_functional.py')
if functional['sha256']!='76d38c86664d922d1519ffc2533e9d7caba3698fef76a657a05a47b9b3d226aa':
    raise RuntimeError('carried Fock/kernel source proof premise changed')
source_bindings=[binding(path) for path in sorted(HERE.iterdir()) if path.is_file()]
inventory=[binding(path) for root in (OUT,RECEIPTS) for path in sorted(root.rglob('*')) if path.is_file()]
payload_receipt=save(OUT/'complete-original-evidence-inventory.json',dict(schema='p254-a45-evidence-inventory-v1',
    files=inventory,combined_original_output_bytes=sum(item['bytes'] for item in inventory),
    source_bindings=source_bindings))
findings=[
 dict(id='actual-moving-electronic-transitions',status='complete',
      proposition='U(t,0) from dU/dt=-i h_baseline(t) U for all128slots/12bands at literal0/1/10/30',
      construction='unchanged SavedMovingBaseline and original affine baseline_hamiltonian; ephemeral degree7 h from8 actual function evaluations per original dense segment; no mean replay, density factorization or new preparation',
      evidence=producer['bundle']),
 dict(id='original-bare-Gaussian-bosons',status='complete',
      proposition='all64fullCartesian S_q includingGammafield from original constant bare mechanical generators; literal rho0/W0 and unstable directions unchanged',
      approximation='NOT electron/cubic/collision feedback or moving selfconsistent boson evolution',
      evidence=producer['bundle']),
 dict(id='actual-moving-reference-EI',status='complete',
      proposition='unchanged A41 E-I consumer computes greater/lesser/retarded Sigma external0/64 andPi externalq0/1/2; all64internalq/128spin-k/12bands/111modes/all16pairs',
      actual_edge_pair_evaluations=kernel['actual_edge_pair_evaluations'],
      analytic_edge_pair_evaluations=kernel['analytic_edge_pair_evaluations'],
      evidence=binding(OUT/'moving-EI/result.json')),
 dict(id='exposing-new-timeorder-oracle',status='passed',
      evidence=binding(OUT/'oracle/result.json'),errors=oracle['errors'],
      scope=oracle['scope']),
 dict(id='actual-persisted-consumption-and-recovery',status='passed',
      evidence=binding(OUT/'kernel-recovery.json'),
      consumer=producer['actual_persisted_consumer']['reader'],
      resumed_U_boundaries=[0,1,10,30],
      recovery=producer['recovery'],
      recovered_kernel_new_pairs=recovery['recovered_new_pairs'],
      recovered_kernel_saved_pairs=recovery['recovered_saved_pairs'],
      all_kernel_arrays_byte_equal=recovery['all_arrays_byte_equal']),
 dict(id='raw-diagnostic-scope',status='diagnostics-not-certification',
      evidence=binding(OUT/'producer/result.json'),
      statement='All U/rho/hermiticity/Ssymplectic/CCR deviations are retained without repair; original local IVP tolerances are not rigorous global error bounds; mean interpolation, integration, unstable amplification and representation/model errors remain uncertified')]
report=dict(schema='P254-A45-complete-original-moving-reference-handoff-v1',status='complete',
    generated_unix=time.time(),findings=findings,blockers=[],failed_scientific_attempts=[],
    original_contract=binding(HERE/'implementation-contract-001.json'),source_bindings=source_bindings,
    original_input_bindings=dict(baseline=contract['settled_original_moving_baseline'],
        material_manifest=contract['settled_material_manifest'],settled_sources=contract['settled_sources']),
    outputs=dict(producer_result=binding(OUT/'producer/result.json'),bundle=producer['bundle'],
        raw_transitions=load(OUT/'producer/bundle.json')['archive'],
        moving_kernel_result=binding(OUT/'moving-EI/result.json'),raw_complex_kernels=kernel['output'],
        new_oracle=binding(OUT/'oracle/result.json'),actual_recovery=binding(OUT/'kernel-recovery.json')),
    all_actual_stage_receipts=wave['stages'],actual_native_launcher=terminal,
    actual_native_argv=load(RECEIPTS/'actual-native-argv.json'),
    actual_manager_startup=load(OUT/'manager-startup.json'),
    actual_manager_terminal=load(RECEIPTS/'actual-manager-terminal.json'),
    native_stdout=binding(RECEIPTS/'native-stdout.txt'),native_stderr=binding(RECEIPTS/'native-stderr.txt'),
    actual_kernel_startup=binding(OUT/'actual-startup.json'),
    costs=dict(cumulative_native_wave_seconds=wave['cumulative_native_seconds'],
        whole_limit_seconds=1800,first_actual_producer_limit_seconds=600,
        child_user_seconds=sum(s['child_user_seconds'] for s in wave['stages']),
        child_system_seconds=sum(s['child_system_seconds'] for s in wave['stages']),
        carried_original_A41_proof_native_seconds=19.008,
        carried_original_A42_baseline_native_seconds=64.333,
        previous_A42_density_join_not_spent_by_this_assignment_native_seconds=452.649),
    output_inventory=payload_receipt,
    output_cap_bytes=402653184,
    raw_U_rho_diagnostics=producer['raw_U_rho_diagnostics'],
    raw_all64_S_diagnostics=producer['raw_S_diagnostics'],
    moving_kernel_raw_diagnostics=kernel['diagnostics'],
    actual_moving_kernel_readout_maxabs=recovery['computed_maxabs'],
    settings=producer['settings'],units=kernel['units'],
    allocation=binding(OUT/'first-producer-allocation.json'),
    unaffected_carried_proofs=dict(original_A41_fock_kernel_proof=binding(proof_path),
        original_A41_admission=binding(HERE.parent/'0041/actual-new-functional-oracle-admission-001.json'),
        original_A41_native_source_grant=binding(HERE.parent/'0041/actual-new-functional-native-grant-001.json'),
        unchanged_functional_source=functional,predicates=proof['predicates'],replayed=False),
    actual_task_labeled_herdr_view=load(OUT/'actual-herdr-view-receipt.json'),
    actual_herdr_view_stdout=binding(OUT/'actual-herdr-view-stdout.json'),
    actual_native_manager_footer=(RECEIPTS/'native-stderr.txt').read_text(),
    provenance_caveats=[
        'A41 historical MATCHING missing text and moving_mean_reconstructed=False are unchanged historical CLI metadata; A45 bundle/producer establishes the now-supplied genuine moving U. No A41 source or CLI replacement was made.',
        'Bare original Gaussian S is the original constant quadratic reference, not feedback from moving electrons/cubic/collisions.',
        'Original Native launched before Root delegated new w1:t6C/w1:p6D. Postcompletion finite view displayed original captured scientific stdout and actual terminal receipt; not a live-through-pane launch. No replay; pane retained.',
        'Postterminal systemctl show returned inactive/dead uninstantiated defaults after manager unload, NOT the runtime resource properties. Actual startup shows InvocationID/CPU7/4GiB/swap0/Restartno/30min; original systemd-run stderr owns terminal success and measured service/CPU costs.'],
    precisely_earned_scope=wave['scope'],
    explicitly_not_claimed=['moving selfconsistent KB','passive/complement loops','full spatial closure',
        'cubic/collision feedback in S','lifetime/rate/heat','wholeE completion','rigorous accuracy certification'])
path=OUT/'complete-original-handoff.json'
result=save(path,report)
if result['bytes']>102400:raise RuntimeError('handoff exceeds100KiB')
combined=sum(item['bytes'] for item in inventory)+payload_receipt['bytes']+result['bytes']
if combined>402653184:raise RuntimeError('combined output cap breached')
print(json.dumps(dict(handoff=result,combined_all_original_output_bytes=combined,
    findings=len(findings),blockers=0)),flush=True)
