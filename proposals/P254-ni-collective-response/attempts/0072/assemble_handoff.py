"""Document actual terminal receipts, saved physics outputs and retained costs."""
from pathlib import Path
from decimal import Decimal
import json
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'0066'))
from column_common import atomic_json, file_binding


def assemble():
    out = Path('/var/tmp/p254-a72-cubic-full-history-q1-001')
    receipt_dir = HERE/'receipts-001'
    exit_record = json.loads((receipt_dir/'exit.json').read_bytes())
    result = json.loads((out/'result.json').read_bytes()) if (out/'result.json').exists() else None
    manager = json.loads((receipt_dir/'manager-terminal.json').read_bytes())
    consumed = json.loads((out/'consumer.json').read_bytes()) if (out/'consumer.json').exists() else None
    cost_exit = json.loads((receipt_dir/'exit.json').read_bytes(),parse_float=Decimal)
    old = Decimal('9332.827233181074257')
    charged = cost_exit['elapsed_seconds']
    ceiling = Decimal('11000')
    inventory = {str(p.relative_to(HERE)):file_binding(p) for p in sorted(receipt_dir.iterdir()) if p.is_file()}
    output_files = {p.name:file_binding(p) for p in sorted(out.iterdir()) if p.is_file()}
    complete = exit_record['returncode'] == 0 and result is not None and result['completed_pairs'] == 961 and consumed is not None and consumed['pairs_consumed'] == 961
    record = dict(schema='p254-a72-full-cubic-q1-terminal-handoff-v1',
        status='complete-owned-component-q1' if complete else 'incomplete-owned-component-q1',
        scoped_result='Actual externalGamma/full111 cubic Pi, internalq1/partnerq3, all genuine31x31 pairs; original A51 ordered_matrix, independent raw ordered boson legs and independently assembled left/right vertices.',
        achieved=dict(completed_pairs=result['completed_pairs'] if result else None,
            literally_reused_pairs=result['reused_pairs'] if result else 33,
            newly_evaluated_pairs=result['new_pairs'] if result else None,
            saved_consumer=consumed,arrays='complex128 NPY greater/lesser/retarded[31,31,111,111]'),
        frozen_preregistration=file_binding(HERE/'preregistration.json'),
        output_files=output_files,receipts=inventory,
        native=dict(returncode=exit_record['returncode'],whole_launcher_wall_seconds=exit_record['elapsed_seconds'],
            actual_end_utc=exit_record['actual_end_utc'],manager_terminal=manager,
            producer_elapsed_seconds=result['elapsed_seconds'] if result else None,
            peak_rss_kib=result['peak_rss_kib'] if result else None,
            manager_receipt_caveat='Terminal systemctl show queried an unloaded transient unit and returned default/unset properties; these are not active execution constraints. Actual CPU9/thread/memory/swap constraints are in actual-runtime.json, exact launch properties in launch.json. systemd stderr reports 421.979s service runtime and 401.227s CPU; its588.0K memory peak is not the Native Python RSS. GNU time independently records1788820KiB Native maximum RSS and zero swaps.',
            bounded_incomplete_boundary_loss=result['bounded_incomplete_boundary_loss'] if result else 'Read persisted coverage.npy; at most one active pair can be uncommitted; no uncommitted coverage is advertised.'),
        cost=dict(old_all_attempts_retained_seconds=str(old),new_conservatively_charged_launcher_wall_seconds=str(charged),
            cumulative_cubic_seconds=str(old+charged),prospectively_amended_ceiling_seconds=str(ceiling),
            remaining_seconds=str(ceiling-old-charged),native_single_launch_cap_seconds=600,automatic_retries=0,
            exact_charge_basis='Decimal parse of actual exit.json elapsed_seconds; includes launcher/manager overhead conservatively, not only numeric contraction.'),
        carried_unchanged=['Original16 simulations','A51 algebra/Schur/Fock/source qualifications'],
        exercised_verification='One actual full961 saved-path consumer inside the same bounded Native; all three saved arrays scanned, numeric full-history/action norms computed, literal33 paid cells and raw causal prescription consumed. No old tests/build/lint/formatter replay.',
        nonclaims=['complete all64q six-array functional','dressed source','stationarity','rank approximation','field-cubic contribution','selfconsistency','physical rates','heat','full simulation','merge'],
        owned_view=file_binding(HERE/'native-view.json'),owned_view_closed=file_binding(HERE/'view-closed.json'))
    atomic_json(HERE/'complete-handoff.json',record)
    print(json.dumps(dict(status=record['status'],handoff=file_binding(HERE/'complete-handoff.json'),cost=record['cost'])),flush=True)

if __name__ == '__main__':
    assemble()
