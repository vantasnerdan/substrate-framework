"""Produce one actual G1 EE or raw mixed-EI q component, never joint/cubic.

API: produce_component(component,q,output,*,provider,resume=False,
                      pair_batch=2,external_batch=2,static_slots=64).
All31 saved nodes/all961 ordered pairs/all128 outputs are the only ownership
scope. Interruption leaves an explicit partial catalog and durable cell resume;
a hard kill needs only checkpoint.json, not a newly published catalog.
"""
from pathlib import Path
import json
import resource
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE)] + [str(HERE.parent / n) for n in (
    '0076', '0063', '0062', '0060', '0056', '0055', '0054', '0053',
    '0052', '0051', '0050', '0048', '0046', '0045', '0044', '0043',
    '0042', '0041', '0030', '0029', '0028', '0027', '0025', '0024', '0064')]
import numpy as np
from history_banks import ORDERS, authenticate, binding
from history_sources import open_reference
import history_provider  # Freeze the dynamically selected provider adapter too.
from material_electron_ion_inputs import load_material, edge_pair
from coulomb_material_inputs import load_coulomb_inputs
from bare_coulomb_factors import BareCoulombFactors
from electron_ion_functional import evaluate_components
from component_consumer import consume_component
from component_emitter import numeric_identity
from component_store import ComponentStore, authenticate_tree, durable_save
from dressed_ee import stream_missing_ee


def _implementation_bindings():
    # Snapshot the owning dependency graph, not unrelated peer modules which an
    # embedding process may also have imported. Module/function references cover
    # the existing loaders' imported project primitives transitively.
    roots = ('produce_component', 'component_store', 'dressed_ee',
             'history_banks', 'history_sources', 'history_provider',
             'material_electron_ion_inputs', 'coulomb_material_inputs',
             'bare_coulomb_factors', 'electron_ion_functional',
             'component_consumer', 'component_emitter')
    pending = [sys.modules[n] for n in roots if n in sys.modules]
    # __main__ is this file when invoked by the CLI.
    paths = {Path(__file__).resolve()}
    visited = set()
    project = HERE.parents[3]
    while pending:
        module = pending.pop()
        name = module.__name__
        if name in visited:
            continue
        visited.add(name)
        path = getattr(module, '__file__', None)
        if not path:
            continue
        path = Path(path).resolve()
        if project not in path.parents or path.suffix != '.py':
            continue
        paths.add(path)
        for value in vars(module).values():
            dep = value if isinstance(value, type(sys)) else sys.modules.get(getattr(value, '__module__', ''))
            if dep is not None and dep.__name__ not in visited:
                pending.append(dep)
    return [binding(p) for p in sorted(paths)]


def _provider_payloads(provider):
    """Bind actual solved arrays, not just the provider's JSON receipt."""
    spec = json.loads(Path(provider).read_bytes())
    payloads = dict(provider=binding(provider), auxiliary=spec['auxiliary'])
    if spec.get('represented_control') is not None:
        payloads['represented_control'] = spec['represented_control']
    for name in ('electron', 'gamma'):
        complete = spec[name]['complete']
        final = json.loads(authenticate(complete).read_bytes())
        payloads[name] = dict(complete=complete, bank=final['bank'],
                              ret=final['ret'], less=final['less'])
    authenticate_tree(payloads)
    return payloads


def _new_run(resume, completed, pair_batch, external_batch, static_slots):
    return dict(resume=bool(resume), started_monotonic=time.monotonic(),
        started_process_cpu=time.process_time(), started_unix_seconds=time.time(),
        reused_committed_slot_cells=completed,
        numerical_output_cells_evaluated=0, committed_slot_cells=0, commits=0,
        commit_flush_seconds=0., pair_batch=pair_batch, external_batch=external_batch,
        static_slots=static_slots, static_rebuild_seconds=0., complete_saved_reuse=False)


def _mixed_ei(c, ref, q, data, store, pair_batch, external_batch):
    """A62 owning EI loop, unchanged full legs, only missing target outputs."""
    run = store.record['runs'][-1]
    run['ei_stream'] = dict(electron_callback_calls=0, boson_callback_calls=0,
                            evaluated_slot_cells=0, contraction_seconds=0.)
    telemetry = run['ei_stream']
    inverse = np.empty(128, dtype=int)
    inverse[np.asarray(c.target[q], dtype=int)] = np.arange(128)
    for a in range(31):
        for b in range(31):
            missing = np.flatnonzero(~store.seen[a, b])
            if not len(missing):
                continue
            # Callbacks retain all128 slots and both orientations, even if only
            # one missing external cell survives. No reverse leg is inferred.
            forward, reverse = ref.electron(a, b), ref.electron(b, a)
            boson = ref.boson(q, a, b)
            telemetry['electron_callback_calls'] += 2
            telemetry['boson_callback_calls'] += 1
            for start in range(0, len(missing), external_batch):
                external = missing[start:start + external_batch]
                tick = time.monotonic()
                values = {order: np.empty((len(external), 1, 12, 12), dtype=np.complex128)
                          for order in ORDERS}
                for j, target in enumerate(external):
                    source = int(inverse[target])
                    pair = edge_pair(data, source, int(target), forward, reverse, boson)
                    result = evaluate_components(pair, n=c.n, t=ref.times[a], s=ref.times[b],
                                                 components=('sigma',))
                    for order in ORDERS:
                        values[order][j, 0] = result['sigma_' + order]
                telemetry['contraction_seconds'] += time.monotonic() - tick
                telemetry['evaluated_slot_cells'] += len(external)
                run['numerical_output_cells_evaluated'] += len(external)
                store.commit(((a, b),), external, values)


def produce_component(component, q, output, *, provider, resume=False,
                      pair_batch=2, external_batch=2, static_slots=64):
    """Return the actual authenticated A76 catalog path on complete production.

No automatic retry. Resume authenticates source/code/committed cells and skips
all saved scientific output cells. Complete resume authenticates the saved
catalog with A76.consume_component before returning, without constructing a
material/history evaluator, rebuilding static data, or evaluating a callback.
    """
    if component not in ('ee', 'mixed_ei'):
        raise ValueError('only ee or mixed_ei components are emitted here')
    if type(q) is not int or not 0 <= q < 64:
        raise ValueError('original integer outer q 0..63 required')
    if type(resume) is not bool:
        raise ValueError('resume must be boolean')
    if any(type(v) is not int or v <= 0 for v in (pair_batch, external_batch, static_slots)):
        raise ValueError('positive integer tile dimensions required')
    out = Path(output).resolve()
    prereg = json.loads((HERE / 'preregistration.json').read_bytes())
    actual = binding(provider)
    if actual != prereg['actual_G1']:
        raise ValueError('actual preregistered G1 provider binding required; no bare/M2 substitution')
    if resume:
        record = json.loads((out / 'checkpoint.json').read_bytes())
        if (record['schema'] != 'p254-additive-history-component-v1'
                or record['component'] != component or record['owned_q'] != [q]
                or record['input_dataset'] != actual
                or record['times_inverse_ry'] != list(np.arange(31, dtype=float))):
            raise ValueError('resume component/q/provider/grid identity differs')
        authenticate_tree(record['implementation_bindings'])
        authenticate(record['preregistration'])
        inputs = json.loads(authenticate(record['input_bindings']).read_bytes())
        authenticate_tree(inputs)
        if record['coverage']['completed_slot_cells'] == 31 * 31 * 128:
            # No callbacks, material decode, static rebuild or native work here.
            catalog = out / 'catalog.json'
            if catalog.exists():
                saved = json.loads(catalog.read_bytes())
                if (record.get('saved_catalog') == binding(catalog)
                        and saved['input_dataset'] == actual and saved['component'] == component
                        and saved['owned_q'] == [q]
                        and saved['status'] == 'complete-owned-component-q%d' % q
                        and saved['checkpoint_bitmap_base64'] == record['checkpoint_bitmap_base64']):
                    observation = consume_component(catalog)
                    durable_save(out / 'resume-consumer.json', dict(
                        consumer=observation, numerical_reevaluations=0, static_rebuilds=0,
                        actual_saved_catalog=binding(catalog), complete_saved_reuse=True))
                    return catalog
            # Final checkpoint can precede catalog publication after a kill.
            # Reopen/authenticate saved cells only, without a numerical producer.
            run = _new_run(True, 31*31*128, pair_batch, external_batch, static_slots)
            run['complete_saved_reuse'] = True
            record['runs'].append(run)
            store = ComponentStore(out, record, resume=True)
            try:
                catalog = store.publish()
            finally:
                store.close()
            durable_save(out / 'resume-consumer.json', dict(consumer=consume_component(catalog),
                numerical_reevaluations=0, static_rebuilds=0, complete_saved_reuse=True))
            return catalog
    else:
        if out.exists():
            raise FileExistsError('single-writer fresh output required; use explicit resume')
        record = None

    began, cpu_began = time.monotonic(), time.process_time()
    ref = store = None
    try:
        c, manifest, coefficient = load_material()
        ref = open_reference(c, coefficient, provider=provider, times=np.arange(31, dtype=float))
        if not np.array_equal(ref.times, np.arange(31, dtype=float)):
            raise ValueError('complete exact31 dressed provider required')
        targets = np.asarray(c.target[q], dtype=int).tolist()
        if sorted(targets) != list(range(128)):
            raise ValueError('original canonical routing permutation required')
        routes = dict(source=list(range(128)), target=targets,
            slot_spin=np.asarray(c.slot_spin).tolist(), slot_k=np.asarray(c.slot_k).tolist(),
            opposite_q=int(c.minus[q]), outer_q=q)
        data = factors = None
        if component == 'ee':
            factors = BareCoulombFactors(load_coulomb_inputs())
            if (not np.array_equal(factors.c.target, c.target)
                    or not np.array_equal(factors.c.slot_spin, c.slot_spin)):
                raise ValueError('EE/material canonical routing disagrees')
        else:
            data = c.q(q)
        if record is None:
            keys = ['mechanical'] + ['slot/%03d' % s for s in range(128)]
            if component == 'mixed_ei':
                keys += ['kernel/%02d' % q] + ['edge/%02d/%03d' % (q, s) for s in range(128)]
            archives = {manifest['items'][key]['archive']['path']: manifest['items'][key]['archive'] for key in keys}
            inputs = dict(reference=ref.bindings, provider_payloads=_provider_payloads(provider),
                coefficient_manifest=coefficient, actual_consumed_item_keys=keys,
                authenticated_archives=list(archives.values()), rho0=numeric_identity(c.rho0),
                numerical_reference=prereg['original_numerical_reference'],
                preparation_origin=prereg['original_preparation'])
            if factors is not None:
                inputs['ee_material'] = factors.c.bindings
            else:
                inputs.update(ionic=numeric_identity(data.ionic),
                              ionic_reverse=numeric_identity(data.ionic_reverse))
            # The original non-Gamma bath is still full111, never a G1->U fit.
            qs = sorted({q, int(c.minus[q])})
            if component == 'mixed_ei':
                inputs['mechanical'] = {str(k): dict(generator=numeric_identity(c.mechanical.generators[k]),
                    W0=numeric_identity(c.mechanical.w0[k])) for k in qs}
            authenticate_tree(inputs)
            record = dict(schema='p254-additive-history-component-v1', component=component,
                status='partial-owned-component-q%d' % q, full_mesh=False, full_solve_eligible=False,
                input_dataset=actual, numerical_reference=prereg['original_numerical_reference'],
                preparation_origin=prereg['original_preparation'],
                times_inverse_ry=ref.times.tolist(), selected_original_time_indices=list(range(31)),
                owned_q=[q], completed_q=[], axes=['time', 'time', 'slot', 'band', 'band'],
                shape=[31, 31, 128, 12, 12], slots=128, bands=12, spins=2, ionic_directions=111,
                canonical_routes=routes, source_bindings=ref.bindings,
                preregistration=binding(HERE / 'preregistration.json'),
                implementation_bindings=_implementation_bindings(),
                coverage=dict(owned_pairs=[[a, b] for a in range(31) for b in range(31)],
                    completed_pairs=[], completed_slot_cells=0, external_targets=list(range(128)),
                    all111_ionic_directions=True,
                    commit='all three numeric tile slices and cell digests flush/fsync; bitmap flush/fsync; atomic checkpoint plus directory fsync'),
                full_time_square=False, raw_nonstationary=True, theta_diagonal=1,
                scope='Actual dressed G1 single-q additive Sigma only; not joint/cubic/all64q/closure/full solve',
                recovery=dict(authority='atomic checkpoint packed bitmap, not ahead-of-checkpoint bitmap bits',
                    committed_science_replayed=False,
                    static_rebuild_cost='material decode, provider authentication, full128 operands/direct bubbles and missing-group static operators are paid again',
                    uncommitted_loss_bound='EE: at most pair_batch*min(static_slots,128) direct output cells plus one pair_batch*min(external_batch,static_slots,128) exchange tile; EI: one min(external_batch,128) output tile',
                    unchanged_retry='never automatic'), runs=[])
            # ComponentStore owns exclusive creation; inputs are published only
            # after that creation and before the initial numerical work.
            record['input_bindings'] = None
        else:
            if record['canonical_routes'] != routes:
                raise ValueError('changed resume canonical routes')
        del manifest
        run = _new_run(resume, record['coverage']['completed_slot_cells'],
                       pair_batch, external_batch, static_slots)
        run['static_rebuild_seconds'] = time.monotonic() - began
        record['runs'].append(run)
        store = ComponentStore(out, record, resume=resume, inputs=inputs if not resume else None)
        run['recovered_ahead_of_checkpoint_slot_cells'] = store.recovered_uncommitted_cells
        store.checkpoint()
        if component == 'ee':
            stream_missing_ee(factors, ref.electron, ref.times, q, store,
                pair_batch=pair_batch, external_batch=external_batch, static_slots=static_slots)
            run['numerical_output_cells_evaluated'] = run['ee_stream']['exchange_output_cells']
        else:
            _mixed_ei(c, ref, q, data, store, pair_batch, external_batch)
        if not store.seen.all():
            raise RuntimeError('owning primitive did not complete all961/all128 cells')
        record['failure'] = None
    except BaseException as exc:
        if store is not None:
            record['failure'] = dict(type=type(exc).__name__, message=str(exc), traceback=traceback.format_exc())
        raise
    finally:
        if store is not None:
            try:
                run['elapsed_seconds'] = time.monotonic() - began
                run['process_cpu_seconds'] = time.process_time() - cpu_began
                run['peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
                if 'ee_stream' in run:
                    run['numerical_output_cells_evaluated'] = run['ee_stream']['exchange_output_cells']
                    run['uncommitted_direct_output_cells'] = run['ee_stream']['direct_output_cells'] - run['committed_slot_cells']
                    run['uncommitted_exchange_output_cells'] = run['ee_stream']['exchange_output_cells'] - run['committed_slot_cells']
                    run['ee_static_operator_build_seconds'] = run['ee_stream']['static_seconds']
                record['elapsed_seconds'] = sum(
                    r.get('elapsed_seconds', r.get('checkpoint_elapsed_seconds', 0.))
                    for r in record['runs'])
                record['elapsed_accounting'] = 'completed/caught runs measured; hard-killed runs include measured last-commit lower bounds only'
                store.publish()
            finally:
                store.close()
                if ref is not None:
                    ref.close()
                    ref = None
        if ref is not None:
            ref.close()
    catalog = out / 'catalog.json'
    durable_save(out / 'consumer.json', consume_component(catalog))
    return catalog


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('component', choices=('ee', 'mixed_ei'))
    parser.add_argument('q', type=int)
    parser.add_argument('output', type=Path)
    parser.add_argument('--provider', type=Path, required=True)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--pair-batch', type=int, default=2)
    parser.add_argument('--external-batch', type=int, default=2)
    parser.add_argument('--static-slots', type=int, default=64)
    args = parser.parse_args()
    catalog = produce_component(args.component, args.q, args.output, provider=args.provider,
        resume=args.resume, pair_batch=args.pair_batch, external_batch=args.external_batch,
        static_slots=args.static_slots)
    print(json.dumps(dict(catalog=str(catalog), component=args.component, q=args.q)), flush=True)


if __name__ == '__main__':
    main()
