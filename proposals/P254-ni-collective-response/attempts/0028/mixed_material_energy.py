"""Authenticated full-source mixed derivative of the eleven-term mean scalar.

Public API: mixed_energy_readout(state_binding, *, dynamics, work, directory,
budget) -> authenticated JSON report binding. Uses NativeMixedDynamics and its
actual prepared n/g/ell/K context, plus the genuine saved mean at the requested
time. Only the mixed bare H0 trace is newly contracted: four borrowed Gamma
dyads per source, original raw D/weights once, arbitrary full-PW column actions
counted against the same QBudget and executed through the real QWork barrier.

Callers must set dynamics.energy_mean_readout to an explicit authenticated
saved-mean ledger with current-source authority and separate historical proof.
Unsupported saved-mean times return status=unavailable;
authentication/resource/worker failures retain a durable failed attempt and
raise. A completed result is replayable without repeating physical actions.
This reports complex d2E/(d lambda d eta), with no half or real projection. It
is NOT finite-dose energy, an energy destination, noise power, heat, a conserved
quantity certificate, or closure of the full spatial/scientific controls.
"""
import json
import os
from pathlib import Path
import sys
from time import monotonic, process_time
from types import SimpleNamespace
from uuid import uuid4

import numpy as np

from substrate_framework.partitioned_q_storage import QBudget, durable_publish, sync_directory
from substrate_framework.streamed_material_state import authenticate, file_binding
import joint_mean_energy_sources
import joint_mean_energy_authority
import mixed_excitation_energy
import mixed_material_energy_sources
from mixed_frame_response import DensityDyad
from substrate_framework.partitioned_q_work import QWork
from mixed_material_dynamics import NativeMixedDynamics
from mixed_material_edge import NativeEdgeOperators
from mixed_material_moments import _one_body
from mixed_material_energy_sources import MeanSourceUnavailable, ledger_binding, saved_mean
from mixed_response_storage import COMPONENTS, MappedMixedState

SCHEMA = 'P254-native-mixed-material-energy-readout-v1'
TERMINAL_BYTES = 16384
TERMS = ('bare_electronic', 'hartree', 'electron_ionic', 'ionic_reference_linear',
         'ionic_C', 'ionic_kinetic_minimal', 'delta_C', 'electron_field',
         'field_contact', 'ion_field_mixed', 'maxwell')


def _json(value):
    return json.loads(json.dumps(value, allow_nan=False))


def _publish_exact(path, value, budget):
    if path.exists():
        binding = file_binding(path)
        if json.loads(authenticate(binding).read_text()) != _json(value):
            raise ValueError('existing mixed energy request has different source/input authority')
        return binding
    return durable_publish(path, value, budget)


def _terminal(path, value):
    """Use the precharged terminal reserve even after a budget/deadline failure."""
    data = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode()
    if len(data) > TERMINAL_BYTES:
        raise RuntimeError('mixed energy terminal receipt reservation exceeded')
    temporary = path.with_name(path.name+'.pending-'+uuid4().hex)
    with temporary.open('xb') as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    sync_directory(path.parent)
    return file_binding(path)


def _sources(dynamics, work, budget):
    modules = (mixed_excitation_energy, mixed_material_energy_sources,
        joint_mean_energy_sources, joint_mean_energy_authority, sys.modules[type(work).__module__],
        sys.modules[type(budget).__module__])
    sources = dict(dynamics.sources, energy_adapter=file_binding(__file__),
        **{m.__name__: file_binding(m.__file__) for m in modules})
    for binding in sources.values():
        authenticate(binding)
    return sources


def _mixed_bare_item(index, state, dynamics, stage):
    """One diagonal endpoint, four existing physical-eta density dyads only."""
    psi = None
    incoming_plus = dynamics.target_minus[index]
    incoming_minus = dynamics.target_plus[index]
    try:
        x = state.leg('mixed', index)
        p = state.leg('plus', incoming_plus)
        m = state.leg('minus', incoming_minus)
        if not any(np.any(a) for a in (x, p, m)):
            return dict(mixed_bare_terms=np.zeros(4, dtype=np.complex128))
        psi = dynamics.baseline.frame(index, state.time)
        # Same physical eta=lambda* bra-direction swap as density_dyads_item.
        # Incoming norms keep their ORIGINAL source D and quadrature weight.
        dyads = (DensityDyad(x, dynamics.D[index], psi, index, dynamics.weights[index]),
            DensityDyad(psi, dynamics.D[index], x, index, dynamics.weights[index]),
            DensityDyad(p, dynamics.D[incoming_plus], p, incoming_plus,
                        dynamics.weights[incoming_plus]),
            DensityDyad(m, dynamics.D[incoming_minus], m, incoming_minus,
                        dynamics.weights[incoming_minus]))
        operators = NativeEdgeOperators(dynamics, index, stage)

        def bare_action(vectors):
            dynamics.budget.action(vectors)
            return operators.bare(index).apply(vectors)

        values = np.asarray([_one_body(d, bare_action, dynamics.block_size)
                             if np.any(d.left) and np.any(d.right) else 0j
                             for d in dyads], dtype=np.complex128)
        return dict(mixed_bare_terms=values)
    finally:
        for source in {index, incoming_plus, incoming_minus}:
            state.release_edge(source)
        if isinstance(psi, np.memmap):
            psi._mmap.close()


def _stages(mean, state, prepared, dynamics):
    polarizations = dynamics.case.field.polarizations
    fields = ('density', 'potential', 'u', 'p', 'gradient', 'field_gradient_electronic')
    stages = dict(mean=SimpleNamespace(**{name: getattr(mean, name) for name in fields},
        shift=mean.shift, momentum=mean.momentum @ polarizations))
    for c in COMPONENTS:
        small, response = state.small[c], dynamics.histories[c].response
        if response.gamma:
            shift = small.coordinate @ polarizations
            momentum = small.momentum @ polarizations
        else:
            if np.any(small.coordinate != 0) or np.any(small.momentum != 0):
                raise ValueError('nonGamma first direction cannot carry global Maxwell modes')
            shift, momentum = np.zeros(3, complex), np.zeros(3, complex)
        stages[c] = SimpleNamespace(**prepared.moments[c], potential=prepared.potential[c],
            u=small.u, p=small.p, shift=shift, momentum=momentum)
    return stages


def _complex(value):
    value = complex(value)
    if not np.isfinite(value):
        raise ValueError('nonfinite actual mixed energy cannot be committed')
    return [value.real, value.imag]


def mixed_energy_readout(state_binding, *, dynamics, work, directory, budget):
    """Consume raw mixed state and genuine read-time mean; return durable binding.

    directory is exclusively owned by this state/source request, not a trajectory
    directory. State payloads are read-only and never retired here. Completed
    per-source scalar QWork receipts survive a failure and are adopted on replay.
    Existing prepared moments/potentials are reused by NativeMixedDynamics.prepare.
    No RHS, unused bare mean/first trace, or mean kernel is recomputed by readout.
    """
    if not isinstance(dynamics, NativeMixedDynamics):
        raise TypeError('genuine NativeMixedDynamics required')
    if not isinstance(work, QWork) or not isinstance(budget, QBudget):
        raise TypeError('genuine QWork and original qualified QBudget required')
    if work.budget is not budget or dynamics.budget is not budget:
        raise ValueError('readout must retain the original cumulative qualified QBudget')
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    budget.reserve(TERMINAL_BYTES)
    terminal = directory/('attempt-'+uuid4().hex+'.json')
    started, cpu_started = monotonic(), process_time()
    before = budget.snapshot()
    request_binding = None
    prepared = None
    try:
        with MappedMixedState(state_binding, budget=budget) as state:
            if state.identity != dynamics.case.identity:
                raise ValueError('raw mixed checkpoint is not bound to the original saved model')
            if state.kind == 'mixed-rk4-quadrature-accumulator':
                raise ValueError('RK4 quadrature accumulator is not a physical read-time state')
            dynamics._validate_state(state)
            baseline_context = dynamics._identity(state)
            sources = _sources(dynamics, work, budget)
            mean_ledger = ledger_binding(dynamics)
            request = dict(schema=SCHEMA, input=state.binding,
                native_context=baseline_context, mean_ledger=mean_ledger,
                sources=sources, work_identity=work.identity,
                scalar='complex mixed d2E/dlambda_deta, physical eta=lambda*, no half')
            request_binding = _publish_exact(directory/'request.json', request, budget)
            result_path = directory/'result.json'
            if result_path.exists():
                binding = file_binding(result_path)
                prior = json.loads(authenticate(binding).read_text())
                if (prior.get('schema') != SCHEMA or prior.get('status') not in ('complete', 'unavailable')
                        or prior.get('request') != request_binding):
                    raise ValueError('existing mixed energy result authority differs')
                _terminal(terminal, dict(schema=SCHEMA, status='reused', result=binding,
                    request=request_binding, budget=budget.snapshot()))
                return binding
            try:
                mean, measure, mean_evidence = saved_mean(state, dynamics, mean_ledger)
            except MeanSourceUnavailable as error:
                report = dict(schema=SCHEMA, status='unavailable', complete=False,
                    request=request_binding, input=state.binding, reason=str(error),
                    requested_time_ry_inverse=error.requested_time,
                    available_mean_sources=error.available, mean_ledger=mean_ledger,
                    physical_actions_performed=False, budget=budget.snapshot(),
                    scope='No fabricated canonical means, endpoint substitution or interpolation; no scalar result claimed.')
                binding = durable_publish(result_path, report, budget)
                _terminal(terminal, dict(schema=SCHEMA, status='unavailable', result=binding,
                    request=request_binding, budget=budget.snapshot()))
                return binding
            budget.check_processes(os.getpid(), [])
            prepared = dynamics.prepare(state, work)
            authenticate(prepared.binding)
            if (prepared.state_binding != state.binding
                    or prepared.record['identity'] != baseline_context
                    or prepared.time != state.time):
                raise ValueError('prepared moments/kernel context is not the actual requested raw state')
            # prepare() has completed the authentic all-source moment/K barrier.
            # This second barrier performs bare H0 ONLY, not density/g/ell/K.
            if state._maps:
                raise ValueError('bare H0 workers must inherit an empty raw leg map cache')
            scalar_before = budget.snapshot()
            identity = dict(request=request_binding, input=state.binding,
                prepared=prepared.binding, mean_source=mean_evidence['stage_record'])
            records = work.execute('native-mixed-energy-bare-h0', identity, dynamics.count,
                lambda i: _mixed_bare_item(i, state, dynamics, prepared.stage))
            key, scalar_directory, count = work.last
            one_body, source_terms = 0j, []
            for index, record in enumerate(records):
                if record['index'] != index or set(record['arrays']) != {'mixed_bare_terms'}:
                    raise ValueError('complete original ordered bare-H0 scalar barrier required')
                source = record['arrays']['mixed_bare_terms']
                value = np.load(authenticate(source), allow_pickle=False)
                if value.shape != (4,) or value.dtype != np.complex128 or not np.isfinite(value).all():
                    raise ValueError('four genuine finite complex128 Gamma dyad traces required')
                subtotal = sum(value, 0j)
                one_body += subtotal
                source_terms.append(dict(endpoint=index,
                    original_sources=[index, index, dynamics.target_minus[index], dynamics.target_plus[index]],
                    weighted_dyad_terms_per_cell_ry=[_complex(z) for z in value],
                    subtotal_per_cell_ry=_complex(subtotal), scalar_array=source,
                    receipt=file_binding(scalar_directory/('item-%04d'%index)/'complete.json')))
            stages = _stages(mean, state, prepared, dynamics)
            raw = mixed_excitation_energy.mixed_components(reference=dynamics.case.reference,
                coefficients=dynamics.case.coefficients, field=dynamics.case.field,
                charges=dynamics.case.charges, stages=stages,
                one_body={'mixed': one_body}, grid_measure=measure,
                curvature_plus=dynamics.common_curvature['plus'],
                curvature_minus=dynamics.common_curvature['minus'],
                gamma=dynamics.histories['plus'].response.gamma)
            if set(raw) != set(TERMS):
                raise ValueError('complete eleven-term actual mixed scalar required')
            components = {name: _complex(raw[name]) for name in TERMS}
            total = sum((complex(raw[name]) for name in TERMS), 0j)
            report = dict(schema=SCHEMA, status='complete', complete=True,
                request=request_binding, input=state.binding, time_ry_inverse=state.time,
                model_identity=dynamics.case.identity,
                source_identity=mean_evidence['provenance']['original_model']['source_identity'],
                components_global_ry=components, total_global_ry=_complex(total),
                mixed_bare_one_body_per_cell_ry=_complex(one_body),
                mean=mean_evidence, prepared=prepared.binding,
                raw_provenance=dict(checkpoint_provenance=state.record['provenance'],
                    native_context=baseline_context, moment_arrays=prepared.record['arrays'],
                    moment_work=prepared.record['moment_work'],
                    original_raw_D=[s.coefficient for s in dynamics.baseline.initial.slots],
                    original_source_weights=list(dynamics.weights),
                    common_curvature='independent actual C(+q) and C(-q), not response C+DeltaC',
                    onsite_delta_C='actual fixed coefficient once; exact zero remains zero',
                    field_conversion='all actual mode coordinates and Pi mapped through original polarizations; nonGamma first shift/Pi explicit Cartesian zero'),
                bare_h0=dict(source_count=count, dyads_per_endpoint=4, source_terms=source_terms,
                    work=dict(key=key, directory=str(scalar_directory), count=count),
                    receipts=records,
                    contraction='bounded arbitrary full-PW columns using mixed_material_moments._one_body; original D/weight once; no PW-square density'),
                measurement=dict(invocation_elapsed_seconds=monotonic()-started,
                    parent_cpu_seconds=process_time()-cpu_started,
                    bare_worker_cpu_seconds=sum(r['cpu_seconds'] for r in records),
                    bare_worker_hamiltonian_actions=sum(r['hamiltonian_actions'] for r in records),
                    receipt_scope='stored actual work measurements include adopted earlier completed items; not all necessarily incurred by this replay',
                    budget_before=before, budget_before_bare=scalar_before,
                    budget_after=budget.snapshot()),
                source_bindings=sources,
                scope=dict(quantity='complex mixed derivative d2E/(dlambda deta), no half or real projection',
                    physical_curve='eta=lambda*: bra derivative direction swaps before dagger',
                    not_claimed=['finite-dose input energy', 'energy destination', 'noise energy',
                        'irreversible heat', 'conservation from a short frozen baseline',
                        'full64q/finite-amplitude/numerical/size/physical-control completion'],
                    approximations='original raw D including singular support, every original PW row/raw column/source and original precision; no stationarity, refilling, residual fit or additional DeltaC'))
            binding = durable_publish(result_path, report, budget)
            _terminal(terminal, dict(schema=SCHEMA, status='complete', result=binding,
                request=request_binding, budget=budget.snapshot()))
            return binding
    except BaseException as error:
        _terminal(terminal, dict(schema=SCHEMA, status='failed', complete=False,
            input=state_binding, request=request_binding,
            error=type(error).__name__+': '+str(error)[:4096],
            elapsed_seconds=monotonic()-started, budget=budget.snapshot(),
            recovery='authenticated completed QWork items retained; no scientific input payload retired'))
        raise
    finally:
        if prepared is not None:
            prepared.close()
