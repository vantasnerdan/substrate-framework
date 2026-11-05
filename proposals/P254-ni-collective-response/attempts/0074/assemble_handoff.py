"""Collect already-terminal A74 evidence without another physics execution."""
from pathlib import Path
import json
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'0063'))
from history_storage import binding,atomic_json
OUT = Path('/var/tmp/p254-a74-full31-bare-force-001')


def main():
    production = json.loads((HERE/'production-result.json').read_bytes())
    result = json.loads((OUT/'result.json').read_bytes())
    consumer = json.loads((OUT/'saved-consumer.json').read_bytes())
    native = json.loads((HERE/'receipts-001/native-exit.json').read_bytes())
    launch = json.loads((HERE/'receipts-001/launch.json').read_bytes())
    exit_record = json.loads((HERE/'receipts-001/exit.json').read_bytes())
    inventory = dict(receipts=[binding(p) for p in sorted((HERE/'receipts-001').iterdir())],
        numeric_boundaries=[binding(p) for p in sorted((OUT/'boundaries').iterdir())],
        outputs=[binding(p) for p in sorted(OUT.iterdir()) if p.is_file()],
        source=[binding(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.json')
            and p.name not in ('complete-handoff.json','evidence-inventory.json')])
    atomic_json(HERE/'evidence-inventory.json',inventory)
    handoff = dict(schema='p254-a74-complete-fullwindow-bare-force-handoff-v1',
        status='complete',route_verdict='established',evidence_scope='audited numerical prescribed bare covariance-force input; no mean evolution',
        actual_outputs=production,api=binding(HERE/'production-api.json'),
        all_evidence=binding(HERE/'evidence-inventory.json'),
        findings=dict(actual_complete_nodes=list(range(31)),actual_complete_internal_q=list(range(64)),
            original_source=result['original_input_dataset'],
            genuine_gamma_source='Actual original saved228-phase bare source, raw full114 greater covariance; only ionic111 contracted, no field cubic',
            saved_consumer=dict(force_byte_identical=consumer['force_byte_identical'],
                gamma_byte_identical=consumer['gamma_byte_identical'],
                rho_byte_identical=consumer['rho_byte_identical'],
                covariance_byte_identical=consumer['covariance_byte_identical'],
                nonGamma_contractions=consumer['nonGamma_contractions'],
                Gamma_contractions=consumer['Gamma_contractions'],
                original_owning_cubic_scalar_probe_evaluations=len(consumer['cubic_scalar_probe'])),
            force_norm_t0=consumer['force_norm'][0],force_norm_t30=consumer['force_norm'][-1],
            gamma_norm_t0=consumer['gamma_norm'][0],gamma_norm_t30=consumer['gamma_norm'][-1],
            Phi000_norm=consumer['Phi000_norm'],
            shared_old_new_nodes=result['shared_node_reuse'],reuse=result['reuse'],
            actual_transitions=result['transitions']),
        cost=dict(native_exit=native,owner_exit=exit_record,argv=launch['argv'],
            actual_start_utc=launch['actual_start_utc'],
            GNU_time_user_seconds=375.61,GNU_time_system_seconds=81.45,
            GNU_time_wall_seconds=465.40,systemd_stderr_runtime_seconds=465.414,
            systemd_stderr_CPU_seconds=457.080,peak_rss_kib=994472,
            swaps=0,cap_seconds=1200,cpu=7,memory_max_bytes=3221225472,threads=1,
            notes=['Native elapsed excludes interpreter/import shutdown portions present in GNU time/service wall.',
                'Post-exit systemctl show observes an unloaded transient unit (unset/infinite values), not the live limits; runtime.json records actual enforced cgroup memory/swap and affinity. Raw terminal stderr and GNU time retained.',
                'q0 checkpoint seconds=0 is not a measured q0 duration; Gamma contraction/setup cost is included in total actual native/service cost, not independently timed.',
                'A57 original333.84841468697414s and paid-cache costs remain historical separate allocations; no old cost reset or merged ledger.']),
        input_bindings=result['input_bindings'],implementation_bindings=result['implementation_bindings'],
        view_closure=binding(HERE/'terminal-view-closure.json'),
        unsupported_recovery=result['recovery'],
        verification_scope='One actual persisted complete bare bank was consumed through the public reader; all31 rho/W/Gamma/full forces recovered byte-identically and the owning cubic scalar evaluated using storedPhi000. No separate canonical dressed same31 dataset was available or fabricated; that reader path is implemented, not claimed exercised.',
        limits=result['limits'],historical_source_edits=0,mean_trajectories_produced=0)
    atomic_json(HERE/'complete-handoff.json',handoff)
    print(json.dumps(binding(HERE/'complete-handoff.json')),flush=True)


if __name__ == '__main__':
    main()
