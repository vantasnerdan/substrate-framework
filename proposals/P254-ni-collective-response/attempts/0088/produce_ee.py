"""EE-only full31 actual-provider producer using unchanged A87 missing-cell EE.

produce_ee(q, output, *, provider, resume=False, pair_batch=2,
           external_batch=16, static_slots=64) returns an A76 catalog path.
No production or numerical validation has been run for this A88 implementation.
"""
from pathlib import Path
import json
import resource
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent / '0087'), str(HERE.parent / '0076')]
import electron_inputs
from electron_inputs import load_inputs
import numpy as np
from history_banks import authenticate, binding
from component_consumer import consume_component
from component_store import ComponentStore, authenticate_tree, durable_save
from dressed_ee import stream_missing_ee
from functional_bank import source_identity


def _implementation_bindings():
    """Bind selected, import-reachable project Python; not Native-library closure."""
    pending = [sys.modules[__name__], electron_inputs]
    visited, paths = set(), {Path(__file__).resolve()}
    project = HERE.parents[3]
    while pending:
        module = pending.pop()
        if module.__name__ in visited:
            continue
        visited.add(module.__name__)
        path = getattr(module, '__file__', None)
        if not path:
            continue
        path = Path(path).resolve()
        if project not in path.parents or path.suffix != '.py':
            continue
        paths.add(path)
        for value in vars(module).values():
            dependency = (value if isinstance(value, type(sys)) else
                          sys.modules.get(getattr(value, '__module__', '')))
            if dependency is not None and dependency.__name__ not in visited:
                pending.append(dependency)
    return [binding(path) for path in sorted(paths)]


def _provider_payloads(provider):
    """Bind the actual A63 mappings and bank receipts, not just provider metadata."""
    spec = json.loads(Path(provider).read_bytes())
    payloads = dict(provider=binding(provider), auxiliary=spec['auxiliary'])
    if spec.get('represented_control') is not None:
        payloads['represented_control'] = spec['represented_control']
    for name in ('electron', 'gamma'):
        complete = spec[name]['complete']
        final = json.loads(authenticate(complete).read_bytes())
        payloads[name] = dict(complete=complete, bank=final['bank'],
                              ret=final['ret'], less=final['less'], slots=final['slots'])
    authenticate_tree(payloads)
    return payloads


def _new_run(resume, completed, pair_batch, external_batch, static_slots,
             began, cpu_began, *, complete_saved_reuse=False):
    # ComponentStore's legacy checkpoint accounting adds this preparation time
    # to time measured from started_monotonic. It includes input authentication,
    # not a claim that static operators were constructed before the EE stream.
    now = time.monotonic()
    return dict(resume=resume, started_monotonic=now,
        started_process_cpu=cpu_began, started_unix_seconds=time.time(),
        reused_committed_slot_cells=completed,
        numerical_output_cells_evaluated=0, committed_slot_cells=0, commits=0,
        commit_flush_seconds=0., pair_batch=pair_batch, external_batch=external_batch,
        static_slots=static_slots, static_rebuild_seconds=now-began,
        preparation_accounting='legacy static_rebuild_seconds includes measured input authentication/loading; EE operator build time is ee_stream.static_seconds',
        complete_saved_reuse=complete_saved_reuse)


def _finish_run(record, run, began, cpu_began):
    run['elapsed_seconds'] = time.monotonic() - began
    run['process_cpu_seconds'] = time.process_time() - cpu_began
    run['peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    if 'ee_stream' in run:
        telemetry = run['ee_stream']
        run['numerical_output_cells_evaluated'] = telemetry['exchange_output_cells']
        run['uncommitted_direct_output_cells'] = telemetry['direct_output_cells'] - run['committed_slot_cells']
        run['uncommitted_exchange_output_cells'] = telemetry['exchange_output_cells'] - run['committed_slot_cells']
        run['ee_static_operator_build_seconds'] = telemetry['static_seconds']
    record['elapsed_seconds'] = sum(
        previous.get('elapsed_seconds', previous.get('checkpoint_elapsed_seconds', 0.))
        for previous in record['runs'])
    record['elapsed_accounting'] = 'measured through final publication start, excluding trailing publication/consumer; hard-killed runs retain measured last-commit lower bounds; earlier run costs retained; pre-store failures have separate sibling receipts'


def _record(q, identity, routes, bindings, implementation):
    return dict(schema='p254-additive-history-component-v1', component='ee',
        status='partial-owned-component-q%d' % q, full_mesh=False, full_solve_eligible=False,
        input_dataset=identity['input_dataset'], numerical_reference=identity['numerical_reference'],
        preparation_origin=identity['preparation_origin'],
        times_inverse_ry=identity['times_inverse_ry'], selected_original_time_indices=list(range(31)),
        owned_q=[q], completed_q=[], axes=['time', 'time', 'slot', 'band', 'band'],
        shape=[31, 31, 128, 12, 12], slots=128, bands=12, spins=2, ionic_directions=111,
        canonical_routes=routes, source_bindings=bindings,
        preregistration=binding(HERE / 'preregistration.json'),
        implementation_bindings=implementation,
        implementation_scope='selected import-reachable project Python only; not complete Native/shared-library closure',
        input_bindings=None,
        coverage=dict(owned_pairs=[[a, b] for a in range(31) for b in range(31)],
            completed_pairs=[], completed_slot_cells=0, external_targets=list(range(128)),
            all111_ionic_directions=True,
            commit='all three numeric tile slices and cell digests flush/fsync; bitmap flush/fsync; atomic checkpoint plus directory fsync'),
        full_time_square=False, raw_nonstationary=True, theta_diagonal=1,
        scope='Actual dressed provider single-q additive EE Sigma only; not joint/cubic/all64q/closure/full solve',
        recovery=dict(authority='atomic checkpoint packed bitmap, not ahead-of-checkpoint bitmap bits',
            committed_science_replayed=False,
            static_rebuild_cost='original Coulomb density decode, A63 provider authentication/maps, full128 operands/direct bubbles and missing-group static operators are paid again on partial resume; static operators are not checkpointed',
            uncommitted_loss_bound='at most pair_batch*min(static_slots,128) direct output cells plus one pair_batch*min(external_batch,static_slots,128) exchange tile, plus in-flight static/input/operand/bubble work',
            unchanged_retry='never automatic'), runs=[])


def produce_ee(q, output, *, provider, resume=False, pair_batch=2,
               external_batch=16, static_slots=64):
    """Produce all961/all128 cells for one q, or explicitly resume missing cells.

Complete resume authenticates committed cells and republishes/consumes their A76
catalog without Coulomb loading, history construction, callbacks or static rebuild.
The atomic checkpoint remains the only ownership authority after interruption.
    """
    if type(q) is not int or not 0 <= q < 64:
        raise ValueError('original integer outer q 0..63 required')
    if type(resume) is not bool:
        raise ValueError('resume must be boolean')
    if any(type(value) is not int or value <= 0
           for value in (pair_batch, external_batch, static_slots)):
        raise ValueError('positive integer tile dimensions required')
    out = Path(output).resolve()
    if not resume and out.exists():
        raise FileExistsError('single-writer fresh output required; use explicit resume')
    began, cpu_began = time.monotonic(), time.process_time()
    factors = ref = store = run = record = None
    try:
        prereg = json.loads((HERE / 'preregistration.json').read_bytes())
        identity = source_identity(provider, prereg['original_preparation']['path'])
        if (identity['preparation_origin'] != prereg['original_preparation'] or
                identity['numerical_reference'] != prereg['original_numerical_reference']):
            raise ValueError('unchanged original preparation/numerical reference required')
        implementation = _implementation_bindings()
        if resume:
            record = json.loads((out / 'checkpoint.json').read_bytes())
            if (record['schema'] != 'p254-additive-history-component-v1'
                    or record['component'] != 'ee' or record['owned_q'] != [q]
                    or any(record[key] != identity[key] for key in
                           ('input_dataset', 'numerical_reference', 'preparation_origin', 'times_inverse_ry'))):
                raise ValueError('resume EE/q/provider/original identity differs')
            authenticate_tree(record['implementation_bindings'])
            if record['implementation_bindings'] != implementation:
                raise ValueError('selected project Python dependency identity differs')
            authenticate(record['preregistration'])
            inputs = json.loads(authenticate(record['input_bindings']).read_bytes())
            authenticate_tree(inputs)
            if record['coverage']['completed_slot_cells'] == 31 * 31 * 128:
                run = _new_run(True, 31*31*128, pair_batch, external_batch,
                               static_slots, began, cpu_began, complete_saved_reuse=True)
                record['runs'].append(run)
                store = ComponentStore(out, record, resume=True)
                run['recovered_ahead_of_checkpoint_slot_cells'] = store.recovered_uncommitted_cells
                record['failure'] = None
                # No early return: finally records the actual authentication/
                # digest/publication cost, retaining every earlier run.
            else:
                factors, ref, inputs = _prepare(provider, identity, q, inputs)
        else:
            factors, ref, inputs = _prepare(provider, identity, q)
            record = _record(q, identity, inputs['canonical_routes'], ref.bindings, implementation)
        if store is None:
            run = _new_run(resume, record['coverage']['completed_slot_cells'],
                           pair_batch, external_batch, static_slots, began, cpu_began)
            record['runs'].append(run)
            record['failure'] = None
            if record['canonical_routes'] != inputs['canonical_routes']:
                raise ValueError('changed resume canonical routes')
            store = ComponentStore(out, record, resume=resume, inputs=inputs if not resume else None)
            run['recovered_ahead_of_checkpoint_slot_cells'] = store.recovered_uncommitted_cells
            store.checkpoint()
            stream_missing_ee(factors, ref.electron, ref.times, q, store,
                pair_batch=pair_batch, external_batch=external_batch, static_slots=static_slots)
            if not store.seen.all():
                raise RuntimeError('owning primitive did not complete all961/all128 cells')
            record['failure'] = None
    except BaseException as exc:
        if store is not None:
            run['failure'] = dict(type=type(exc).__name__, message=str(exc), traceback=traceback.format_exc())
            record['failure'] = run['failure']
        else:
            # Never create/adopt a scientific output or overwrite an existing
            # checkpoint merely to report a failed preparation/lock attempt.
            receipt = out.with_name(out.name + '.preparation-failure-%d.json' % time.time_ns())
            try:
                receipt.parent.mkdir(parents=True, exist_ok=True)
                durable_save(receipt, dict(
                    schema='p254-a88-ee-preparation-failure-v1', q=q,
                    output=str(out), provider=str(provider), resume=resume,
                    failure=dict(type=type(exc).__name__, message=str(exc), traceback=traceback.format_exc()),
                    elapsed_seconds=time.monotonic()-began,
                    process_cpu_seconds=time.process_time()-cpu_began,
                    store_constructor_returned=False,
                    recovery='no output/checkpoint is adopted by this receipt; interrupted constructor files may exist; only a durable checkpoint permits explicit resume'))
            except OSError as receipt_error:
                print('could not persist preparation failure receipt: %s' % receipt_error, file=sys.stderr)
        raise
    finally:
        try:
            if store is not None:
                _finish_run(record, run, began, cpu_began)
                store.publish()
        finally:
            if store is not None:
                store.close()
            if ref is not None:
                ref.close()
    catalog = out / 'catalog.json'
    observation = consume_component(catalog)
    if run['complete_saved_reuse']:
        durable_save(out / 'resume-consumer.json', dict(consumer=observation,
            numerical_reevaluations=0, static_rebuilds=0, complete_saved_reuse=True,
            actual_saved_catalog=binding(catalog)))
    else:
        durable_save(out / 'consumer.json', observation)
    return catalog


def _prepare(provider, identity, q, saved_inputs=None):
    factors, ref = load_inputs(provider)
    try:
        c = factors.c
        if c.bindings['manifest'] != identity['original_B_identity']:
            raise ValueError('original Coulomb coefficient manifest differs from numerical reference')
        payloads = _provider_payloads(provider)
        if (ref.bindings[0] != identity['input_dataset'] or
                payloads['provider'] != identity['input_dataset']):
            raise ValueError('actual mapped provider differs from authenticated ancestry')
        targets = np.asarray(c.target[q], dtype=int).tolist()
        if sorted(targets) != list(range(128)):
            raise ValueError('original canonical routing permutation required')
        routes = dict(source=list(range(128)), target=targets,
            slot_spin=np.asarray(c.slot_spin).tolist(), slot_k=np.asarray(c.slot_k).tolist(),
            opposite_q=int(c.minus[q]), outer_q=q)
        inputs = dict(reference=ref.bindings, provider_payloads=payloads,
            ee_material=c.bindings, coefficient_manifest=c.bindings['manifest'],
            actual_consumed_item_keys=(['slot/%03d' % slot for slot in range(128)] +
                ['edge/%02d/%03d' % (inner_q, slot) for inner_q in range(64) for slot in range(128)]),
            actual_consumed_members=dict(slot=['h0', 'rho0'], edge=['density', 'feature_g', 'feature_partner'],
                background=['reciprocal_bohr_inv', 'cell_bohr']),
            authenticated_archives=c.bindings['authenticated_archives'],
            canonical_routes=routes, numerical_reference=identity['numerical_reference'],
            preparation_origin=identity['preparation_origin'],
            material_inventory_scope='original Coulomb reader consumption only; no mechanical/kernel/ionic member claim; original producer metadata is provenance, not selected runtime implementation')
        authenticate_tree(inputs)
        if saved_inputs is not None and saved_inputs != inputs:
            raise ValueError('actual resume Coulomb/provider input bindings differ')
        return factors, ref, inputs
    except BaseException:
        ref.close()
        raise


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('q', type=int)
    parser.add_argument('output', type=Path)
    parser.add_argument('--provider', type=Path, required=True)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--pair-batch', type=int, default=2)
    parser.add_argument('--external-batch', type=int, default=16)
    parser.add_argument('--static-slots', type=int, default=64)
    args = parser.parse_args()
    catalog = produce_ee(args.q, args.output, provider=args.provider,
        resume=args.resume, pair_batch=args.pair_batch,
        external_batch=args.external_batch, static_slots=args.static_slots)
    print(json.dumps(dict(catalog=str(catalog), component='ee', q=args.q)), flush=True)


if __name__ == '__main__':
    main()
