"""Two recoverable material preparation boundaries inside the existing runtime."""
from pathlib import Path
import time
import numpy as np

from material_run_contract import bound_json, implementation, verify_consumers
from .material_bloch_output import load_saved_bloch_coverage
from .material_preparation_inputs import (KEYS, validate_grant, load_population_inputs,
                                         load_factor_terminal, population_metadata)
from .material_preparation_gamma import gamma_metadata, instantiate_gamma
from .material_preparation_output import (population_arrays, load_saved_preparation_block,
    validate_population_arrays, validate_selected_arrays)
from .material_preparation_restore import restore_population_specification


def _recover(run, grant, key):
    exists = any((run.target / (key + suffix)).exists() for suffix in ('.json', '.npz'))
    admitted = grant['resume_receipts'].get(key)
    if not exists:
        if admitted is not None:
            raise ValueError('admitted completed preparation boundary is missing')
        return None
    if (not grant['resume'] or admitted is None
            or Path(admitted['path']).resolve() != run.target / (key + '.json')):
        raise ValueError('explicit immutable receipt admission required for existing boundary')
    receipt = bound_json(admitted)
    arrays = load_saved_preparation_block(run.target, key, run.identity, receipt)
    return arrays, receipt


def produce(run, grant, paths):
    validate_grant(grant)
    run.check()
    started = time.monotonic()
    cached = _recover(run, grant, KEYS[0])
    if cached is None:
        coverage, specification = load_population_inputs(grant['inputs'], grant['settings'])
        arrays = population_arrays(specification)
        metadata = population_metadata(grant['inputs'], grant['settings'], coverage, specification)
        validate_population_arrays(arrays, metadata, coverage)
        verify_consumers(grant, implementation(paths))
        receipt = run.commit(KEYS[0], arrays, metadata)
    else:
        arrays, receipt = cached
        coverage = load_saved_bloch_coverage(grant['inputs']['bloch_terminal'])
        validate_population_arrays(arrays, receipt['metadata'], coverage)
        specification = restore_population_specification(coverage, arrays, receipt['metadata'])
        metadata = population_metadata(grant['inputs'], grant['settings'], coverage, specification)
        if receipt['metadata'] != metadata:
            raise ValueError('restored population recipe/policy/source metadata differs')
        run.record(KEYS[0], receipt, 'completed-population-ledger-reuse-no-refill')
    run.summary['timing_seconds']['coverage_intake_population_archive'] = time.monotonic() - started
    run.summary['coverage_identity'] = coverage.coverage_identity
    run.summary['neutral_mean_number'] = specification.filling.weighted_number
    run.summary['neutral_target'] = specification.filling.target_mean_number
    run.check()
    # Exercise the actual saved block reader even after the commit's exact IO readback.
    ledger = load_saved_preparation_block(run.target, KEYS[0], run.identity, receipt)
    validate_population_arrays(ledger, metadata, coverage)
    verify_consumers(grant, implementation(paths))
    started = time.monotonic()
    terminal = load_factor_terminal(grant['inputs'], coverage)
    record, factor_receipt, gamma = gamma_metadata(grant['inputs'], grant['settings'], coverage, terminal)
    run.check()
    cached = _recover(run, grant, KEYS[1])
    if cached is None:
        selected, admissibility = instantiate_gamma(specification, record, factor_receipt, gamma)
        selected['pauli_roundoff_allowance'] = np.array([admissibility['pauli_roundoff_allowance']])
        selected['gram_min_eigenvalue'] = np.array([admissibility['gram_min_eigenvalue']])
        validate_selected_arrays(selected, ledger, coverage)
        verify_consumers(grant, implementation(paths))
        selected_receipt = run.commit(KEYS[1], selected, gamma)
    else:
        selected, selected_receipt = cached
        if selected_receipt['metadata'] != gamma:
            raise ValueError('restored Gamma source/frame/recipe differs')
        validate_selected_arrays(selected, ledger, coverage)
        run.record(KEYS[1], selected_receipt, 'completed-Gamma-D-and-observables-reuse-no-WFC-or-Gram-recompute')
    run.summary['timing_seconds']['Gamma_join_instantiation_archive'] = time.monotonic() - started
    started = time.monotonic()
    saved = load_saved_preparation_block(run.target, KEYS[1], run.identity, selected_receipt)
    validate_selected_arrays(saved, ledger, coverage)
    run.summary['timing_seconds']['selected_saved_reader'] = time.monotonic() - started
    run.summary['selected_number_pair'] = [[float(x.real), float(x.imag)] for x in saved['selected_number_pair']]
    run.summary['scope'] = 'complete neutral population metadata; ONE actual spin1 Gamma raw-D/Pauli/number/probe consumer'
    sources = implementation(paths)
    verify_consumers(grant, sources)
    run.summary['final_implementation_sources'] = sources
    run.check()
