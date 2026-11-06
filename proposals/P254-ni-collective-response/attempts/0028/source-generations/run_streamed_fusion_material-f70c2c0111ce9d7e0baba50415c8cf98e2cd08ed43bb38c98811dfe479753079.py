#!/usr/bin/env python3
"""Actual cutover smoke; carry original scientific proofs, never replay them."""
import argparse
from pathlib import Path
from types import SimpleNamespace

from run_streamed_material_evolution import (
    ROOT, HERE, np, explicit_binding, bound_json, preload_parsers, restore_inputs,
    build_evolution_model, publish, authenticate, DiskMean)
from substrate_framework.streamed_material_state import file_binding, save_array
from threadpoolctl import threadpool_limits


def historical(entry):
    """Preserve an executed epoch; a changed source is not relabeled current."""
    path = Path(entry['path'])
    current = file_binding(path)
    if (current['sha256'], int(current['bytes'])) == (entry['sha256'], int(entry['bytes'])):
        retained = current
    else:
        retained = file_binding(HERE/'source-generations'/
                                (path.stem+'-'+entry['sha256']+path.suffix))
        if (retained['sha256'], int(retained['bytes'])) != (entry['sha256'], int(entry['bytes'])):
            raise ValueError('original executed source bytes unavailable')
    return dict(original=entry, retained=retained)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('fused-report', 'prior-oracle'):
        parser.add_argument('--'+name, type=Path, required=True)
        parser.add_argument('--'+name+'-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    fused_binding = explicit_binding(args.fused_report, args.fused_report_sha256)
    prior_binding = explicit_binding(args.prior_oracle, args.prior_oracle_sha256)
    fused, prior = bound_json(fused_binding), bound_json(prior_binding)
    if (fused['schema'] != 'P254-fused-action-material-v1' or fused['status'] != 'passed'
            or prior['schema'] != 'P254-streamed-evolution-oracle-v1' or prior['status'] != 'passed'):
        raise ValueError('executed original numerical and material authorities required')
    fusion_oracle = bound_json(fused['oracle'])
    if fusion_oracle['status'] != 'passed':
        raise ValueError('original independent fusion oracle failed')
    preserved = [historical(entry) for report in (prior, fused, fusion_oracle)
                 for entry in report['source_bindings']]
    current_paths = [ROOT/'src/substrate_framework'/name for name in (
        'streamed_material_state.py', 'streamed_material_evolution.py',
        'streamed_sector_executor.py', 'material_frozen_action.py')]
    current_paths += [Path(__file__), HERE/'run_streamed_material_evolution.py',
                      HERE/'run_streamed_sector_pilot.py']
    current_sources = [file_binding(path) for path in current_paths]
    inputs = fused['inputs']
    for entry in inputs.values():
        authenticate(entry)
    for entry in (fused['probe_columns'], *fused['probe_outputs'].values()):
        authenticate(entry)
    preload_parsers()
    arrays, receipt, manifest, saved, original, grid, ordered, pairs = restore_inputs(
        inputs['aggregate_receipt'], inputs['manifest'])
    ions = bound_json(inputs['initial_ions'])
    state = DiskMean(inputs['checkpoint']['path'], verify=False)
    selected = [slot for slot in state.slots if slot.identity == fused['sector']['identity']]
    if len(selected) != 1 or selected[0].metadata() != fused['sector']:
        raise ValueError('original selected sector identity differs')
    slot = selected[0]
    settings = SimpleNamespace(block_size=8, reality_tolerance=1e-10, max_h0_actions=2)
    model, background = build_evolution_model(settings, arrays, manifest, saved,
        original, grid, ordered, inputs['aggregate_receipt'], ions, state.identity)
    model.validate(state)
    stage = model.stage(arrays['conditional_density_bohr3'],
                        arrays['positive_gradient_ry_bohr'].reshape(-1), state.u)
    probes = np.load(fused['probe_columns']['path'], mmap_mode='r', allow_pickle=False)
    args.output.mkdir(parents=True, exist_ok=False)
    errors, outputs, refusals = {}, {}, []
    with model.operators(slot) as (coupling, h0):
        if h0.material is None:
            raise AssertionError('actual material factory did not select the fused path')
        full_action = model.frozen_action(coupling, h0, slot, stage)
        for name, action in (('h0', h0), ('actual_initial_stage', full_action)):
            actual = action(probes)
            with np.load(fused['probe_outputs'][name]['path'], allow_pickle=False) as expected:
                error = float(np.linalg.norm(actual-expected['legacy']) /
                              max(1., np.linalg.norm(expected['legacy'])))
            errors[name] = error
            outputs[name] = save_array(args.output/(name+'-actual.npy'), actual)
            if error > 2e-11:
                raise AssertionError('actual cutover differs from retained full-PW reference')
        before = model.counts
        for name, action in (('bare', h0), ('frozen', full_action)):
            try:
                action(probes)
            except RuntimeError as error:
                if 'H0 action cap exhausted' not in str(error):
                    raise
                refusals.append(name)
            else:
                raise AssertionError('shared H0 work cap was bypassed')
        if (before != model.counts or before['h0_calls'] != 2
                or before['h0_vector_columns'] != 2*probes.shape[1]):
            raise AssertionError('work reservation doubled or consumed after cap refusal')
    for entry in current_sources:
        authenticate(entry)
    publish(args.output/'result.json', dict(
        schema='P254-streamed-evolution-oracle-v1', status='passed',
        kind='carried original scientific proofs plus actual material cutover smoke',
        prior_scientific_oracle=prior_binding, fused_material=fused_binding,
        independent_fusion_oracle=fused['oracle'], original_executed_sources=preserved,
        source_bindings=current_sources, inputs=inputs, sector=slot.metadata(),
        retained_reference_probes=fused['probe_outputs'], actual_probe_arrays=outputs,
        relative_errors=errors, work_counts=model.counts, action_cap_refusals=refusals,
        background_receipt=background,
        scope='same frozen full-PW action and shared work accounting; no time trajectory',
        original_proofs_replayed=False))


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
