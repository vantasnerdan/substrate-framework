#!/usr/bin/env python3
"""Execute a finite-only self-consistent coupled mixed response proof.

Usage: python run_mixed_coupled_response_oracle.py --output REPORT.json
       [--work-directory NEW_DIRECTORY] [--max-seconds REMAINING_SECONDS]

PASS is written only after genuine QWork/RK4/checkpoint/readout calculations,
independent nonlinear full-density SciPy continuations, parameter/time refinement,
SymPy scalar-derived canonical checks, numeric mutations, and retired-stage replay.
This Gamma finite fixture is NOT material evidence. Previously executed separate
frame/UPF programs, not this comparator, cover nonGamma geometry.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import resource
import sys
from time import perf_counter, process_time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'src'))

import numpy as np
import scipy
import sympy
from mixed_finite_coupled_system import (
    CoupledFixture, FiniteDynamics, ENERGIES, ORDERS, small_vector)
from mixed_response_storage import (
    MappedMixedState, COMPONENTS, SMALL_NAMES, payload_path)
from mixed_response_rk4 import advance_step
from substrate_framework.partitioned_q_storage import QBudget
from substrate_framework.partitioned_q_work import QWork
from substrate_framework.streamed_material_state import authenticate, file_binding


class Evidence:
    def __init__(self):
        self.checks = []

    def require(self, name, condition, **details):
        self.checks.append(dict(name=name, passed=bool(condition), **details))
        if not condition:
            raise AssertionError(name)

    def close(self, name, actual, expected, tolerance):
        actual, expected = np.asarray(actual), np.asarray(expected)
        defect = float(np.max(np.abs(actual-expected), initial=0))
        self.require(name, actual.shape == expected.shape
            and np.isfinite(actual).all() and np.isfinite(expected).all()
            and defect <= tolerance, max_abs_defect=defect, tolerance=tolerance)
        return defect


def source_bindings():
    # Bind actual loaded repository modules, including producer/storage/fixture
    # and their transitive imports; never substitute a predeclared text list.
    paths = {Path(__file__).resolve()}
    for module in tuple(sys.modules.values()):
        value = getattr(module, '__file__', None)
        if value:
            path = Path(value).resolve()
            if path.suffix == '.py' and path.is_relative_to(ROOT):
                paths.add(path)
    return [file_binding(p) for p in sorted(paths)]


def norm(value):
    return float(np.max(np.abs(value), initial=0))


def parameter_derivative(function, h, cache):
    def value(z):
        if z not in cache:
            cache[z] = function(z)
        return cache[z]
    # Wirtinger d_lambda d_eta = (d_xx+d_yy)/4 on eta=lambda*.
    return (value(complex(h))+value(complex(-h))+value(complex(0, h))
            +value(complex(0, -h))-4*value(0j))/(4*h*h)


def independent_continuations(fixture, evidence):
    cache = {}
    hs = (.04, .02, .01)
    densities = [parameter_derivative(fixture.integrate, h, cache) for h in hs]
    energies = []
    for h in hs:
        energies.append(parameter_derivative(
            lambda z: fixture.symbolic.energy(fixture, cache[z]), h, {}))
    rich_density = [(4*densities[i+1]-densities[i])/3 for i in range(2)]
    rich_energy = [(4*energies[i+1]-energies[i])/3 for i in range(2)]
    d0, d1 = norm(densities[1]-densities[0]), norm(densities[2]-densities[1])
    e0, e1 = norm(energies[1]-energies[0]), norm(energies[2]-energies[1])
    evidence.require('physical_parameter_density_refinement', d0 > 0 and d1 < .4*d0,
                     coarse_difference=d0, fine_difference=d1)
    evidence.require('physical_parameter_scalar_refinement', e0 > 0 and e1 < .4*e0,
                     coarse_difference=e0, fine_difference=e1)
    evidence.close('density_Richardson_stability', rich_density[1], rich_density[0], 3e-7)
    evidence.close('eleven_scalar_Richardson_stability', rich_energy[1], rich_energy[0], 3e-7)
    return rich_density[1], rich_energy[1], dict(
        steps=list(hs), full_density_continuations=len(cache),
        density_raw_differences=[d0, d1], scalar_raw_differences=[e0, e1],
        density_Richardson_uncertainty=norm(rich_density[1]-rich_density[0]),
        scalar_Richardson_uncertainty=norm(rich_energy[1]-rich_energy[0]))


def evolve(fixture, initial, count, root, budget, sources, mutation='none'):
    root.mkdir(parents=True, exist_ok=False)
    work = QWork(root/'qwork', budget, (6, 7, 8),
                 'finite-coupled-'+mutation)
    dynamics = FiniteDynamics(fixture.system, root/'contexts', budget, sources, mutation)
    binding = initial
    dt = (fixture.stop-fixture.start)/count
    calls = []
    for step in range(count):
        directory = root/f'step-{step:03d}'
        identity = f'finite-coupled-{mutation}-{count}-{step}'
        before = binding
        binding = advance_step(before, dt, dynamics=dynamics, work=work,
            directory=directory, identity=identity, budget=budget)
        calls.append(dict(initial=before, directory=directory, identity=identity, dt=dt))
    return binding, dynamics, work, calls


def canonical_check(fixture, initial, root, budget, sources, evidence):
    work = QWork(root/'canonical-qwork', budget, (6, 7, 8), 'finite-canonical')
    dynamics = FiniteDynamics(fixture.system, root/'canonical-contexts', budget, sources)
    with MappedMixedState(initial, budget=budget) as state:
        evidence.require('genuine_one_mode_Gamma_checkpoint',
            all(state.small[c].coordinate.shape == (1,)
                and state.small[c].momentum.shape == (1,) for c in COMPONENTS))
        prepared = dynamics.prepare(state, work)
        actual = np.concatenate([prepared.small_rhs['mixed']['u'],
            prepared.small_rhs['mixed']['p'], prepared.small_rhs['mixed']['coordinate'][:1],
            prepared.small_rhs['mixed']['momentum'][:1]])
        def direct_force(z):
            rho, small = fixture.unpack_density(fixture.continuation_initial(z))
            return fixture.symbolic.canonical(small, fixture.symbolic.moments(rho))
        cache = {}
        coarse = parameter_derivative(direct_force, .02, cache)
        fine = parameter_derivative(direct_force, .01, cache)
        rich = (4*fine-coarse)/3
        evidence.close('canonical_mixed_force_same_SymPy_scalar', actual, rich, 2e-8)
        restored = dynamics.restore(prepared.binding, state, work)
        for c in COMPONENTS:
            for name in SMALL_NAMES:
                evidence.close('durable_initial_context_'+c+'_'+name,
                    restored.small_rhs[c][name], prepared.small_rhs[c][name], 0.)
    return fixture.symbolic.force_derivatives


def replay_and_restore(binding, dynamics, work, calls, budget, evidence):
    # The producer's normal successful cutover has already retired old stages.
    # Repeat the exact operation with a brand-new dynamics object. No moments
    # may be recomputed from retired input; complete producer journals decide.
    first = calls[0]
    dynamics = FiniteDynamics(dynamics.s, dynamics.directory, budget,
                             dynamics.sources, dynamics.mutation)
    retired = []
    for path in first['directory'].rglob('*.retired.json'):
        record = json.loads(path.read_text())
        actual_paths = [actual for entry in record.get('payloads', [])
                        if (actual := payload_path(entry)) is not None]
        if actual_paths and any(not Path(actual).exists() for actual in actual_paths):
            retired.append(file_binding(path))
    evidence.require('old_actual_input_payloads_retired_before_replay', bool(retired),
                     retirement_receipts=retired)
    replay = advance_step(first['initial'], first['dt'], dynamics=dynamics,
        work=work, directory=first['directory'], identity=first['identity'], budget=budget)
    evidence.require('actual_advance_step_completed_replay', replay == binding,
                     final_checkpoint=binding)
    restored_count = 0
    for path in dynamics.directory.rglob('context.json'):
        record = json.loads(path.read_text())
        old = record['state']
        checkpoint = json.loads(authenticate(old).read_text())
        retired_legs = [(component, index, entry)
                        for component, entries in checkpoint['legs'].items()
                        for index, entry in enumerate(entries)
                        if (actual := payload_path(entry)) is not None
                        and not Path(actual).exists()]
        if not retired_legs:
            continue
        with MappedMixedState(old, budget=budget) as state:
            context = dynamics.restore(file_binding(path), state, work)
            evidence.require('retired_stage_context_restore_'+str(restored_count),
                all(np.isfinite(a).all() for values in context.small_rhs.values() for a in values.values()),
                context=file_binding(path), retired_checkpoint=old)
            for component, index, entry in retired_legs:
                try:
                    state.leg(component, index)
                except FileNotFoundError:
                    refused = True
                else:
                    refused = False
                evidence.require(
                    'retired_stage_is_not_zero_substituted_'+str(restored_count)
                    +'_'+component+'_'+str(index), refused,
                    component=component, index=index, retired_payload=entry)
            restored_count += 1
    evidence.require('actual_retired_stage_restore_exercised', restored_count > 0,
                     restored_contexts=restored_count)


def receipt_evidence(root):
    result = []
    for path in root.rglob('complete.json'):
        record = json.loads(path.read_text())
        if 'process_peak_rss_bytes' in record:
            result.append(dict(binding=file_binding(path), pid=record['pid'],
                cpu=record['cpu'], affinity=record['affinity'],
                cpu_seconds=record['cpu_seconds'], wall_seconds=record['elapsed_seconds'],
                peak_rss_bytes=record['process_peak_rss_bytes'],
                hamiltonian_actions=record['hamiltonian_actions']))
    return result


def run(args, result, evidence):
    fixture = CoupledFixture()
    budget = QBudget(max_actions=100000, max_seconds=args.max_seconds,
        max_bytes=512 << 20, max_memory_bytes=3 << 30, worker_memory_bytes=1 << 30)
    sources = source_bindings()
    result.update(source_bindings=sources, library_versions=dict(
        numpy=np.__version__, scipy=scipy.__version__, sympy=sympy.__version__))
    root = args.work_directory.resolve()
    root.mkdir(parents=True, exist_ok=False)
    result['work_directory'] = str(root)
    initial = fixture.initial(root/'initial', budget, sources)
    result['initial_checkpoint'] = initial
    result['canonical_SymPy_derivatives'] = canonical_check(
        fixture, initial, root, budget, sources, evidence)
    initial_density, initial_energy = fixture.readout(initial, budget)
    s = fixture.system
    _, u0, p0, a0, pi0, rho0, _ = s.mean_values(s.mean.sol(fixture.stop))
    baseline = fixture.pack_density(rho0, np.concatenate([u0, p0, a0, pi0]))
    evidence.close('independent_full_density_reproduces_exact_moving_mean',
                   fixture.integrate(0j), baseline, 2e-10)
    evidence.require('genuinely_nonstationary_moving_mean',
        norm(baseline-fixture.continuation_initial(0j)) > 1e-4)
    density_cache, energy_cache = {}, {}
    direct_density = lambda z: fixture.continuation_initial(z)
    direct_energy = lambda z: fixture.symbolic.energy(fixture, direct_density(z))
    for name, actual, function, cache in (
            ('initial_full_mixed_density', initial_density, direct_density, density_cache),
            ('initial_eleven_mixed_energy', initial_energy, direct_energy, energy_cache)):
        coarse = parameter_derivative(function, .02, cache)
        fine = parameter_derivative(function, .01, cache)
        evidence.close(name, actual, (4*fine-coarse)/3, 3e-8)
    independent_density, independent_energy, continuation = independent_continuations(fixture, evidence)
    result['independent_parameter_refinement'] = continuation
    result['time_refinement'] = []
    errors = []
    fine_binding = None
    for count in (1, 2, 4):
        binding, dynamics, work, calls = evolve(fixture, initial, count,
            root/f'correct-{count}', budget, sources)
        actual_density, actual_energy = fixture.readout(binding, budget)
        error = norm(actual_density-independent_density)
        errors.append(error)
        result['time_refinement'].append(dict(steps=count,
            dt=(fixture.stop-fixture.start)/count, full_density_and_small_max_error=error,
            energy_components_max_error=norm(actual_energy-independent_energy),
            final_checkpoint=binding))
        if count == 1:
            replay_and_restore(binding, dynamics, work, calls, budget, evidence)
        if count == 4:
            fine_binding, fine_density, fine_energy = binding, actual_density, actual_energy
    evidence.require('actual_coupled_RK4_step_refinement',
        errors[0] > 0 and errors[1] < .6*errors[0] and errors[2] < .45*errors[1],
        max_errors=errors, expected_order='fourth-order global RK4; independent parameter error floor')
    evidence.close('full_self_consistent_mixed_density_and_small', fine_density, independent_density, 6e-7)
    for j, name in enumerate(ENERGIES):
        evidence.close('independent_full_scalar_mixed_'+name,
                       fine_energy[j], independent_energy[j], 6e-7)
    mutant_results = []
    for mutation in ('missing_mixed_source', 'half_mixed_cross_source',
                     'wrong_source_weight', 'wrong_canonical_N', 'wrong_canonical_charge'):
        mutated, _, _, _ = evolve(fixture, initial, 1, root/mutation,
                                 budget, sources, mutation)
        density, energy = fixture.readout(mutated, budget)
        error = norm(density-independent_density)
        evidence.require('numeric_consumer_rejects_'+mutation,
            error > max(2e-6, 5*errors[0]), max_abs_defect=error,
            correct_coarse_defect=errors[0], actual_mutant_checkpoint=mutated)
        mutant_results.append(dict(mutation=mutation,
            full_density_and_small_defect=error, energy_defect=norm(energy-independent_energy),
            checkpoint=mutated))
    _, stationary_energy = fixture.readout(fine_binding, budget, omit_mean_mixed=True)
    stationary_defect = norm(stationary_energy-independent_energy)
    evidence.require('numeric_consumer_rejects_omitted_mean_times_mixed',
        stationary_defect > max(2e-6, 10*norm(fine_energy-independent_energy)),
        max_abs_defect=stationary_defect)
    result['numeric_mutations'] = mutant_results
    receipts = receipt_evidence(root)
    evidence.require('actual_QWork_numeric_workers_exercised',
        bool(receipts) and any(r['hamiltonian_actions'] > 0 for r in receipts)
        and all(r['affinity'] == [r['cpu']] and r['cpu'] in (6, 7, 8) for r in receipts),
        actual_worker_receipts=len(receipts))
    result.update(final_checkpoint=fine_binding,
        mixed_energy_components={name: [float(value.real), float(value.imag)]
                                 for name, value in zip(ENERGIES, fine_energy)},
        worker_measurements=receipts, budget=budget.snapshot(),
        final_source_bindings=source_bindings())
    # Every consumed repository source must still be the actual original bytes.
    for source in sources:
        authenticate(source)
    result['status'] = 'PASS'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--work-directory', type=Path)
    parser.add_argument('--max-seconds', type=float, default=180.)
    args = parser.parse_args()
    if args.work_directory is None:
        args.work_directory = args.output.with_name(args.output.stem+'-work')
    if args.max_seconds <= 0:
        parser.error('--max-seconds must be positive')
    evidence = Evidence()
    started, cpu_started = perf_counter(), process_time()
    children_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = dict(schema='P254-finite-coupled-mixed-response-oracle-v1', status='FAIL',
        finite_only=True, material_result=False, gamma_only=True,
        scope='finite nonstationary exact moving mean; fixed original raw nonorthogonal frames and singular D',
        nonGamma_boundary='not exercised here; separate previously executed frame/UPF programs cover geometry',
        canonical_photon_modes=1,
        full_objective_status='E/D/C open; no material energy construction or ranking inferred',
        started_utc=datetime.now(timezone.utc).isoformat())
    exit_code = 0
    try:
        run(args, result, evidence)
    except BaseException as error:
        result.update(status='FAIL', error=type(error).__name__+': '+str(error),
                      traceback=traceback.format_exc())
        exit_code = 1
    finally:
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        result.update(checks=evidence.checks,
            measurements=dict(wall_seconds=perf_counter()-started,
                parent_cpu_seconds=process_time()-cpu_started,
                children_cpu_seconds=(children.ru_utime+children.ru_stime
                    -children_before.ru_utime-children_before.ru_stime),
                parent_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                children_max_process_peak_rss_bytes=children.ru_maxrss*1024),
            ended_utc=datetime.now(timezone.utc).isoformat())
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('x') as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write('\n')
        print(json.dumps(dict(status=result['status'], output=str(args.output),
                              checks=len(evidence.checks))))
    return exit_code


if __name__ == '__main__':
    raise SystemExit(main())
