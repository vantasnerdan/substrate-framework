#!/usr/bin/env python3
"""Bounded full-PW original-source covariance assessment, not trajectory reuse.

Examples are in operator-covariance-registration.json. --source selects original
spin:ik, --target selects spin:i,j,k. --source-orbit selects all targets of those
explicit sources. --all-slots admits the complete original128-slot inventory;
--slot-offset/--slot-count split it into bounded batches without a new model.
Every selected target is compared with its ORIGINAL identity representative.
Physical mode imports only explicitly SHA-bound, settled factory bytes supplied
by Root, after that factory's independent oracle has settled. No fixture model,
projector fitting, occupation refilling, geometry repair or antiunitary shortcut.
"""
import argparse
import importlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
from time import perf_counter
from types import SimpleNamespace as NS

from run_operator_covariance_oracle import (
    HERE, ROOT, np, binding, authenticate, explicit, publish, sources, sha256)
from run_streamed_material_evolution import (
    restore_inputs, receipt_arrays, bound_json, preload_parsers, slot_identity,
    restore_accepted, load_wavefunctions, transport_bloch_frame,
    source_population_indices, restore_source)
from substrate_framework.qe_bloch_frame import admit_source
from substrate_framework.qe_material_hamiltonian import build_pw_factors_for_basis, MaterialHamiltonian
from substrate_framework.material_reference_transport import slot_transport
from substrate_framework.material_operator_covariance import (
    UnitaryPWMap, fourier_witness, covariance_bound, action_witness)
from substrate_framework.streamed_material_evolution import real_physical
from scipy.linalg import cho_factor, cho_solve
from threadpoolctl import threadpool_limits, threadpool_info


def array_binding(value):
    value = np.asarray(value)
    digest = sha256()
    if value.flags.c_contiguous:
        digest.update(memoryview(value).cast('B'))
    else:
        # Hash C-order payload in small row slabs, not a second full raw frame.
        for start in range(0, value.shape[0], 64):
            digest.update(memoryview(np.ascontiguousarray(value[start:start+64])).cast('B'))
    return dict(shape=list(value.shape), dtype=value.dtype.str,
                sha256=digest.hexdigest(),
                bytes=int(value.nbytes), encoding='C-order raw array bytes; shape/dtype separately bound')


def memory():
    result = dict(peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    for line in Path('/proc/self/status').read_text().splitlines():
        name, _, value = line.partition(':')
        if name in ('VmRSS', 'VmHWM', 'VmSize'):
            result[name+'_kib'] = int(value.split()[0])
    return result


def loaded_sources():
    paths = {Path(__file__).resolve(), HERE/'operator-covariance-registration.json'}
    for module in tuple(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path:
            path = Path(path).resolve()
            if path.suffix == '.py' and path.is_relative_to(ROOT):
                paths.add(path)
    return [binding(path) for path in sorted(paths)]


def source_key(text):
    try:
        spin, ik = (int(part) for part in text.split(':'))
    except (ValueError, TypeError):
        raise argparse.ArgumentTypeError('source must be original spin:ik')
    if spin not in (1, 2) or not 1 <= ik <= 10:
        raise argparse.ArgumentTypeError('source must be original spin1/2 and ik1..10')
    return spin, ik


def select_slots(args, ordered):
    by_id = {slot_identity(s): s for s in ordered}
    if args.all_slots:
        if args.source or args.target or args.source_orbit:
            raise ValueError('--all-slots cannot be mixed with source/target selection')
        stop = None if args.slot_count is None else args.slot_offset+args.slot_count
        selected = ordered[args.slot_offset:stop]
    else:
        if not args.source or bool(args.target) == bool(args.source_orbit) or args.slot_offset or args.slot_count is not None:
            raise ValueError('explicit --source and exactly --target or --source-orbit required')
        keys = set(args.source)
        if len(keys) != len(args.source):
            raise ValueError('duplicate source selection')
        if args.source_orbit:
            selected = [s for s in ordered if (s.ispin, s.source_ik) in keys]
        else:
            if len(set(args.target)) != len(args.target) or any(t not in by_id for t in args.target):
                raise ValueError('target identities must be distinct original full-grid slots')
            selected = [by_id[t] for t in args.target]
            if {(s.ispin, s.source_ik) for s in selected} != keys:
                raise ValueError('every target must match exactly the explicit original source selection')
    if not selected or len(selected) > args.max_slots:
        raise ValueError('nonempty selected-slot count must fit explicit --max-slots')
    if any(s.operation.effective_time_reversal is not False for s in selected):
        raise ValueError('no physical antiunitary positive-time reuse')
    return selected


def original_state(record, saved, grid, representative, source, frame):
    population, _ = source_population_indices(grid, saved.population_ledger, source.ispin, source.ik)
    state_binding = record['state']
    if record['reuse_selected_gamma']:
        if state_binding is not None or (source.ispin, source.ik) != (1, 1) or representative.mesh_index != (0, 0, 0):
            raise ValueError('selected Gamma state authority differs')
        arrays = {name: saved.selected_gamma[name] for name in
                  ('gram', 'populations', 'density_coefficients', 'density_eigenvalues')}
    else:
        arrays, receipt = receipt_arrays(state_binding)
        metadata = receipt['metadata']
        if (metadata['source']['sha256'] != source.source['sha256']
                or metadata['original_basis_identity'] != frame.basis.identity
                or metadata['original_frame_identity'] != frame.gauge_identity):
            raise ValueError('source state has different original raw basis/gauge')
    if not np.array_equal(arrays['populations'], population):
        raise ValueError('original all-slot populations differ')
    return restore_source(arrays, source_identity=source.source['sha256'])


def complement_probes(psi, count, seed):
    """Probe-only projection; never alters the original frame, Gram or D."""
    if not 1 <= count <= psi.shape[0]-psi.shape[1]:
        raise ValueError('positive complement probe count must fit ambient complement')
    rng = np.random.default_rng(seed)
    z = rng.normal(size=(psi.shape[0], count))+1j*rng.normal(size=(psi.shape[0], count))
    gram = psi.conj().T@psi
    factor = cho_factor(gram, lower=True)
    for _ in range(2):
        z -= psi@cho_solve(factor, psi.conj().T@z)
    z, triangular = np.linalg.qr(z, mode='reduced')
    if np.any(np.abs(np.diag(triangular)) < 1e-12):
        raise ValueError('rank-deficient complement probes; no silent replacement')
    return z, dict(seed=seed, columns=count, raw_frame_orthogonality_frobenius=float(np.linalg.norm(psi.conj().T@z)),
                   probe_gram_defect_frobenius=float(np.linalg.norm(z.conj().T@z-np.eye(count))),
                   raw_frame_gram=array_binding(gram), qualification='empirical full-PW complement probes only')


def operators(background, factors, ispin, frozen_vrs, block_size):
    h0 = MaterialHamiltonian(background, factors, ispin, block_size=block_size)
    stage_background = NS(fft=background.fft, inputs=background.inputs, vrs=frozen_vrs)
    frozen = MaterialHamiltonian(stage_background, factors, ispin, block_size=block_size)
    return {'h0': h0, 'initial_frozen': frozen}


def validate_factors(factors, basis):
    if (not np.array_equal(factors.fft.miller, basis.miller)
            or not np.array_equal(factors.k_bohr_inv, basis.k_cart_bohr_inv)
            or not np.array_equal(factors.fft.reciprocal_bohr_inv, basis.reciprocal_bohr_inv)
            or factors.kinetic_ry.shape != (basis.miller.shape[1],)):
        raise ValueError('factory changed the checked original full-PW basis')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--oracle', required=True, type=Path)
    parser.add_argument('--oracle-sha256', required=True)
    parser.add_argument('--projector-model', required=True, choices=('legacy-a25', 'physical-upf'))
    parser.add_argument('--physical-module')
    parser.add_argument('--physical-sha256')
    parser.add_argument('--source', action='append', type=source_key)
    parser.add_argument('--target', action='append')
    parser.add_argument('--source-orbit', action='store_true')
    parser.add_argument('--all-slots', action='store_true')
    parser.add_argument('--slot-offset', type=int, default=0)
    parser.add_argument('--slot-count', type=int)
    parser.add_argument('--max-slots', required=True, type=int)
    parser.add_argument('--complement-probes', required=True, type=int)
    parser.add_argument('--block-size', type=int, default=8)
    parser.add_argument('--intertwiners', type=Path, help='optional NPZ with declared angular Q under keys opINDEX_SPECIES; never fitted')
    parser.add_argument('--intertwiners-sha256')
    parser.add_argument('--cpu', required=True, type=int, choices=range(6, 11))
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if (not 1 <= args.max_slots <= 128 or args.complement_probes < 1 or args.block_size < 1
            or args.slot_offset < 0 or args.slot_count is not None and args.slot_count < 1):
        parser.error('positive bounded slot/probe/block allocation required')
    if bool(args.intertwiners) != bool(args.intertwiners_sha256):
        parser.error('declared intertwiners require both path and digest')
    if args.projector_model == 'physical-upf' and (not args.physical_module or not args.physical_sha256):
        parser.error('physical mode requires explicit settled module and digest')
    if args.projector_model == 'legacy-a25' and (args.physical_module or args.physical_sha256):
        parser.error('legacy mode does not consume a physical factory')
    os.sched_setaffinity(0, {args.cpu})
    start = perf_counter()
    registration = json.loads((HERE/'operator-covariance-registration.json').read_text())
    for path, digest in registration['source_authority'].items():
        explicit(ROOT/path, digest)
    oracle_binding = explicit(args.oracle, args.oracle_sha256)
    oracle = bound_json(oracle_binding)
    if (oracle.get('schema') != 'P254-operator-covariance-oracle-v1' or oracle.get('status') != 'passed'
            or not oracle.get('checks') or not all(oracle['checks'].values())
            or oracle.get('source_bindings') != sources()):
        raise ValueError('executed sensitive independent oracle for exact current covariance source bytes required')
    for entry in oracle['source_bindings']:
        authenticate(entry)
    factory, physical_binding = build_pw_factors_for_basis, None
    if args.projector_model == 'physical-upf':
        if args.physical_module != 'substrate_framework.material_physical_operators':
            raise ValueError('only the declared physical-UPF factory interface is admitted')
        specification = importlib.util.find_spec(args.physical_module)
        if specification is None or specification.origin is None:
            raise ValueError('settled physical factory module unavailable')
        physical_binding = explicit(Path(specification.origin), args.physical_sha256)
        factory = importlib.import_module(args.physical_module).build_physical_factors
    intertwiners, intertwiner_binding = {}, None
    if args.intertwiners:
        intertwiner_binding = explicit(args.intertwiners, args.intertwiners_sha256)
        with np.load(args.intertwiners, allow_pickle=False) as archive:
            intertwiners = {name: archive[name] for name in archive.files}
    aggregate_binding = binding(HERE/'full-common-reference-004/full-common-reference.json')
    manifest_binding = binding(HERE/'full-reference-inputs-004.json')
    ions_binding = binding(HERE/'initial-mean-ions-001.json')
    ions = bound_json(ions_binding)
    args.output.mkdir(parents=True, exist_ok=False)
    publish(args.output/'started.json', dict(oracle=oracle_binding, projector_model=args.projector_model,
        physical_factory=physical_binding, aggregate=aggregate_binding, manifest=manifest_binding,
        initial_ions=ions_binding, intertwiners=intertwiner_binding, max_slots=args.max_slots,
        complement_probes=args.complement_probes, cpu_affinity=sorted(os.sched_getaffinity(0)), memory=memory()))
    preload_parsers()
    with threadpool_limits(limits=1):
        arrays, aggregate, manifest, saved, original, grid, ordered, pairs = restore_inputs(aggregate_binding, manifest_binding)
        selected = select_slots(args, ordered)
        background, background_receipt = restore_accepted(original)
        if not np.array_equal(arrays['saved_density_bohr3'], background.rho_updown):
            raise ValueError('actual saved background differs from complete aggregate')
        displacement = np.asarray(ions['u_bohr'])
        if displacement.shape != (3*grid.reference.nat,) or np.any(displacement):
            raise ValueError('this initial frozen witness requires actual zero displacement; no moving-stage term omitted')
        density = real_physical(arrays['conditional_density_bohr3'], 'actual initial density', 1e-10)
        kernel_raw = background.apply_kernel(np.asarray(density-arrays['saved_density_bohr3'], dtype=complex))
        kernel = real_physical(kernel_raw, 'actual initial K(delta n)', 1e-10)
        if kernel.shape != background.vrs.shape:
            raise ValueError('initial frozen kernel has wrong spin/grid shape')
        frozen_vrs = background.vrs+kernel
        provenance = loaded_sources()
        groups = {}
        for slot in selected:
            groups.setdefault((slot.ispin, slot.source_ik), []).append(slot)
        records, source_records = [], []
        record_by_key = {(r['spin'], r['source_ik']): r for r in manifest['sources']}
        used_intertwiners = set()
        for key, targets in groups.items():
            source_started = perf_counter()
            representatives = [s for _, s in pairs[key] if np.array_equal(s.operation.rotation, np.eye(3, dtype=int))
                               and not np.any(s.operation.fractional_translation) and not np.any(s.reciprocal_shift)]
            if len(representatives) != 1:
                raise ValueError('unique unchanged original identity representative required')
            rep = representatives[0]
            source = load_wavefunctions(rep.source_binding)
            admission = admit_source(grid.reference, rep, source)
            source_frame = transport_bloch_frame(grid.reference, rep, source, source_admission=admission)
            state = original_state(record_by_key[key], saved, grid, rep, source, source_frame)
            psi = source.coefficients
            d_binding = array_binding(state.density_coefficients)
            probes, probe_report = complement_probes(psi, args.complement_probes, 2540028+100*key[0]+key[1])
            sf = factory(background.inputs, source_frame.basis)
            validate_factors(sf, source_frame.basis)
            source_ops = operators(background, sf, key[0], frozen_vrs, args.block_size)
            # Cache source actions once per original source, not once per target.
            cached = {name: operator.apply(psi) for name, operator in source_ops.items()}
            cached_probes = {name: operator.apply(probes) for name, operator in source_ops.items()}
            source_records.append(dict(spin=key[0], source_ik=key[1], original_slot=slot_identity(rep),
                original_wavefunction=dict(source.source), original_basis=source_frame.basis.identity,
                original_gauge=source_frame.gauge_identity, raw_frame=array_binding(psi), original_D=d_binding,
                original_state_receipt=record_by_key[key]['state'], selected_gamma_preparation=manifest['preparation_terminal'] if record_by_key[key]['reuse_selected_gamma'] else None,
                admission=dict(admission), complement=probe_report, source_action_cache_seconds=perf_counter()-source_started,
                source_factors_provenance=list(sf.sources), raw_gram_saved_defect_frobenius=float(np.linalg.norm(psi.conj().T@psi-state.gram))))
            for slot in targets:
                slot_started = perf_counter()
                target_frame = transport_bloch_frame(grid.reference, slot, source, source_admission=admission)
                transform = slot_transport(grid, slot)
                phase = np.ones(source.igwx, complex)
                if np.any(slot.operation.fractional_translation):
                    phase = np.exp(-2j*np.pi*((source.miller.T+slot.source_k_fractional)@slot.operation.fractional_translation))
                mapping = UnitaryPWMap(target_frame.basis.row_map, phase, source_frame.basis.identity, target_frame.basis.identity)
                frame_defect = 0.
                for first in range(0, psi.shape[1], args.block_size):
                    last = first+args.block_size
                    frame_defect = max(frame_defect, float(np.max(np.abs(target_frame.coefficients[:, first:last]-mapping.apply(psi[:, first:last])))))
                tf = sf if target_frame.basis.identity == source_frame.basis.identity else factory(background.inputs, target_frame.basis)
                validate_factors(tf, target_frame.basis)
                target_ops = operators(background, tf, key[0], frozen_vrs, args.block_size)
                fw = fourier_witness(sf.fft, tf.fft, mapping, rotation=slot.operation.rotation,
                    translation=slot.operation.fractional_translation, reciprocal_shift=slot.reciprocal_shift,
                    source_k_fractional=slot.source_k_fractional, fft_pullback=transform.fft_pullback)
                declared = {}
                for species in tf.species:
                    name = 'op'+str(slot.operation.index)+'_'+species.name
                    if name in intertwiners:
                        declared[species.name] = intertwiners[name]; used_intertwiners.add(name)
                cases = {}
                for name in ('h0', 'initial_frozen'):
                    potential = background.vrs[:, key[0]-1] if name == 'h0' else frozen_vrs[:, key[0]-1]
                    report = covariance_bound(sf, tf, mapping, potential, potential,
                        fft_pullback=transform.fft_pullback, fourier=fw,
                        atom_pullback=transform.atom_pullback, angular_intertwiners=declared)
                    original_action = action_witness(source_ops[name].apply, target_ops[name].apply, mapping, psi,
                        block_size=args.block_size, cached_source_action=cached[name])
                    complement_action = action_witness(source_ops[name].apply, target_ops[name].apply, mapping, probes,
                        block_size=args.block_size, cached_source_action=cached_probes[name])
                    cases[name] = dict(bound=report, all_original_columns=original_action, complement=complement_action,
                        largest_empirical_witness_minus_analytic_bound_ry=max(original_action['empirical_max_relative_ry'], complement_action['empirical_max_relative_ry'])-report['total_bound_ry'],
                        numerical_bound_zero=report['total_bound_ry'] == 0., automatic_transport_authorized=False)
                indices = np.indices(sf.fft.physical).reshape(3, -1, order='F').T/np.asarray(sf.fft.physical)
                multiplier = np.exp(-2j*np.pi*(indices@slot.reciprocal_shift+slot.source_k_fractional@slot.operation.fractional_translation))
                fft_error = tf.fft.from_g(mapping.apply(probes))-multiplier[:, None]*sf.fft.from_g(probes)[transform.fft_pullback]
                fft_probe = np.linalg.norm(fft_error, axis=0)/np.sqrt(sf.fft.nnr)
                record = dict(spin=slot.ispin, source_ik=slot.source_ik, source_slot=slot_identity(rep),
                    target_slot=slot_identity(slot), mesh_index=list(slot.mesh_index), weight=slot.weight,
                    operation=slot.operation.index, source_basis=source_frame.basis.identity,
                    target_basis=target_frame.basis.identity, target_gauge=target_frame.gauge_identity,
                    unitary=True, orbital_conjugation=False, row_map=array_binding(mapping.row_pullback),
                    phase=array_binding(mapping.phase), fft_pullback=array_binding(transform.fft_pullback),
                    atom_pullback=transform.atom_pullback.tolist(), cartesian_rotation=transform.cartesian_rotation.tolist(),
                    geometry_residual_bohr=transform.geometry_residual_bohr.tolist(),
                    geometry_max_bohr=float(np.max(np.abs(transform.geometry_residual_bohr))),
                    frame_transport_max_residual=frame_defect, fft_probe_residuals=fft_probe.tolist(),
                    background_diagnostics=transform.diagnostics(background),
                    kernel_local_covariance_defect_ry=float(np.max(np.abs(kernel[:, slot.ispin-1]-transform.density(kernel)[:, slot.ispin-1]))),
                    target_factors_provenance=list(tf.sources), cases=cases,
                    seconds=perf_counter()-slot_started, memory=memory())
                publish(args.output/('slot-'+str(len(records))+'.json'), record)
                records.append(record)
                del target_frame, tf, target_ops, fft_error
            if array_binding(state.density_coefficients) != d_binding:
                raise ValueError('original fixed D changed')
            del cached, cached_probes, source_ops, sf, source_frame, source, state, probes, psi
        if set(intertwiners) != used_intertwiners:
            raise ValueError('unused declared intertwiner entries; supply exactly the selected operation/species scope')
        pools = threadpool_info()
    for entry in provenance:
        authenticate(entry)
    result = dict(schema='P254-operator-covariance-material-v1', status='assessed-not-auto-authorized',
        projector_model=args.projector_model, physical_factory=physical_binding, oracle=oracle_binding,
        inputs=dict(aggregate=aggregate_binding, manifest=manifest_binding, initial_ions=ions_binding,
                    original=original, background=background_receipt, intertwiners=intertwiner_binding),
        implementation_sources=provenance, source_records=source_records, slots=records,
        selected_slots=len(records), complete_original_slot_coverage=len(records) == 128,
        actual_initial_kernel=array_binding(kernel), kernel_discarded_imaginary_max=float(np.max(np.abs(np.imag(kernel_raw)))),
        initial_displacement=array_binding(displacement), initial_external_field='zero; no field operator in this H0/current-local-stage model',
        seconds=perf_counter()-start, memory=memory(), threadpools=pools,
        cpu_affinity=sorted(os.sched_getaffinity(0)),
        conditions=dict(exact_mathematical='Exact unitary map, Fourier embedding, zero kinetic/local/projector residuals and D-conjugacy defects imply exact covariance of this operator only.',
            approximate_current_stage='The reported analytic sum bounds the declared current-stage operator in exact arithmetic. Floating evaluation is not interval enclosed. Nonzero bounds remain nonzero; no material tolerance fitted or automatic reuse authorized.',
            frozen_propagation='For Hermitian stages only: norm(exp(-it Ht)-U exp(-it Hs) U*) <= abs(t)*operator_bound; apply your explicitly allocated physical error budget, not probe maxima.',
            future_stages='Not licensed. Reassess actual current K(delta n), displacement/field transform conditions and every additional moving/field term. Initial density symmetry is insufficient.',
            antiunitary='No same-positive-time reuse: antiunitary conjugation flips i; spin/time reversal requires a separate admitted model.',
            empirical='Original346-column and complement actions are empirical witnesses, never upper certificates.',
            geometry='Saved rounded atom positions unchanged; geometry defects retained, not fitted away.'),
        whole_model_reduction_claim=False, trajectory_executed=False)
    publish(args.output/'result.json', result)


if __name__ == '__main__':
    main()
