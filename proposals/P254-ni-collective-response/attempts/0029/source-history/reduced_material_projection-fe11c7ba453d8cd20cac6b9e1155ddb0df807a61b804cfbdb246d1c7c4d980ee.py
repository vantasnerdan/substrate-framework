"""Export current-physical active-band coefficients from authenticated t=0 states.

Run with the repository's existing numerical environment; no Native cache miss
is regenerated. Numerical payloads must remain outside the repository. Completed
items survive interruption; the manifest becomes complete only after all128
sectors, all64 q, kernels and explicit nested passive splits are published.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import resource
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0028')]
import numpy as np
from substrate_framework.streamed_material_state import DiskMean, authenticate, file_binding
from substrate_framework.material_control_reference import bound_arrays, bound_json
from substrate_framework.material_physical_operators import (
    physical_hamiltonian, physical_finiteq, PhysicalFieldVertices, PhysicalFrozenAction)
from substrate_framework.material_joint_response import RealBlochEdge
from substrate_framework.qe_bloch_finiteq import BlochFiniteQBridge
from matched_joint_case import load_case, slot_identity
from reduced_active_basis import ActiveBasis, SELECTION_RECIPE, prepare_active_basis
from reduced_projection_algebra import (q_feature_geometry, project_edge, project_kernel,
    projection_witness, projector_overlaps)
from reduced_projection_cache import ProjectionCache, load_projected_coefficients
from reduced_prepared_projection import (basis_arrays, restore_basis, overlap_arrays,
    restore_overlaps, prepare_ionic_edge)
from reduced_preparation_storage import MappedPreparations
from reduced_paired_local_support import attach_paired_local_support
from reduced_projection_reuse import admit_gamma_pilot, adopt_items


def indexing():
    indices = np.indices((4, 4, 4)).reshape(3, -1).T
    target = np.empty((64, 128), int)
    for iq, q in enumerate(indices):
        for s in range(128):
            k = indices[s % 64]
            target[iq, s] = (s//64)*64+np.ravel_multi_index(tuple((k+q) % 4), (4, 4, 4))
    return dict(k_indices=indices.tolist(), q_indices=indices.tolist(),
                slot_spin=np.repeat([1, 2], 64).tolist(), slot_k=np.tile(np.arange(64), 2).tolist(),
                target=target.tolist(), minus=np.ravel_multi_index(tuple((-indices % 4).T), (4, 4, 4)).tolist())


def selected(text, size):
    if text == 'all':
        return list(range(size))
    out = set()
    for part in text.split(','):
        if '-' in part:
            a, b = map(int, part.split('-'))
            out.update(range(a, b+1))
        else:
            out.add(int(part))
    if not out or min(out) < 0 or max(out) >= size:
        raise ValueError('selection outside explicit indexing')
    return sorted(out)


def mechanical_arrays(case):
    a = case.arrays; d = len(case.reference.masses)
    return dict(force_lattice=a['common_force_lattice_ry_bohr2'], masses=case.reference.masses,
        charges=a['valence_charges'], delta_curvature=np.zeros((d, d)), contact=a['contact'],
        mixed=a['mixed_ionic_field'].reshape(d, 3),
        force=a['ionic_force_ry_bohr'].reshape(d)-a['positive_gradient_ry_bohr'].reshape(d),
        field_gradient=a['field_gradient'], polarizations=case.field.polarizations,
        field_mass=np.asarray(case.field.mass), ionic_ordered=case.ionic_preparation.local.ordered,
        field_ordered=case.field_preparation['ordered_covariance'], trap_omega=case.ionic_preparation.omega_ry)


class MaterialExporter:
    def __init__(self, args):
        self.started = perf_counter(); self.args = args
        directory = Path(args.cache).resolve()
        if directory == ROOT or ROOT in directory.parents:
            raise ValueError('numerical projection cache must remain outside Git')
        self.case = load_case(args.case)  # restore_control_background, never build Native
        restore_seconds = perf_counter()-self.started
        self.boundary_binding = file_binding(args.boundary)
        self.original = DiskMean(authenticate(self.boundary_binding), verify=False)
        if self.original.time != 0 or np.any(self.original.u) or np.any(self.original.p):
            raise ValueError('literal original t=0 zero-coordinate boundary required')
        self.slots = tuple(sorted(self.case.ordered, key=lambda s: (s.ispin, s.mesh_index)))
        self.records = {s.identity: s for s in self.original.slots}
        if len(self.records) != 128 or set(self.records) != {slot_identity(s) for s in self.slots}:
            raise ValueError('original complete128 explicit sectors required')
        if self.original.evidence.get('case') != args.case:
            raise ValueError('boundary belongs to a different material case')
        # Initializer and exact-rebind receipts have different evidence fields.
        # Actual native frame/D/slot values, not one receipt shape, bind basis().
        result_path = HERE.parent/'0028'/('control-'+args.case+'-001')/'result.json'
        result_binding = file_binding(result_path)
        result = bound_json(result_binding)
        self.interval = result['chemical_potential_interval_ry']
        code = [file_binding(HERE/name) for name in ('reduced_active_basis.py',
            'reduced_projection_algebra.py', 'reduced_projection_cache.py', 'reduced_material_projection.py',
            'reduced_paired_local_support.py', 'reduced_prepared_projection.py',
            'reduced_projection_storage.py', 'reduced_preparation_storage.py', 'reduced_projection_reuse.py',
            'reduced_local_gradient.py')]
        # Full material source inventory is authenticated by load_case and bound
        # again by the exact current specification, boundary, and ledger context.
        material_paths = sorted({str(Path(module.__file__).resolve()) for name, module in tuple(sys.modules.items())
            if (name.startswith('substrate_framework.') or name.startswith('matched_'))
            and getattr(module, '__file__', '').endswith('.py')})
        spec = dict(case=args.case, case_identity=self.case.identity, context=self.case.specification,
            original_boundary=self.boundary_binding, neutral_labels=result_binding,
            implementation=code, physical_implementation=[file_binding(p) for p in material_paths],
            bands=args.bands, features=args.features,
            selection_recipe=SELECTION_RECIPE,
            feature_g=[q_feature_geometry(i, args.features)[0].tolist() for i in range(64)],
            feature_partner=[q_feature_geometry(i, args.features)[1].tolist() for i in range(64)],
            mesh=[4, 4, 4], ndof=len(self.case.reference.masses), indexing=indexing(),
            approximation='fixed symmetric prepared saved-band subspace plus explicit frozen passive initial totals',
            positive_time_symmetry=False, original_source_weight_in_rho=False,
            initial_potential='restored physical K(full original density - saved background density)',
            ionic_support='original full density G union exact canonical opposite-edge partners; original UPF radial extension',
            preparation_storage='per-slot original admission once; immutable active-column mmap and cached beta/p-axis overlaps',
            reverse_convention='independent physical positive-gradient adjoint moments indexed X_ab')
        admitted = None
        if args.reuse_gamma_pilot:
            admitted = admit_gamma_pilot(args.reuse_gamma_pilot, spec, HERE/'source-history')
            spec['explicit_gamma_pilot_reuse'] = admitted[0]
        self.cache = ProjectionCache(directory, spec)
        if admitted is not None:
            adopt_items(self.cache, admitted)
        self.cache.manifest['timings']['last_restore'] = dict(seconds=restore_seconds)
        self.cache.publish()
        self.lookup = spec['indexing']
        self.last_index = None; self.last_bundle = None
        self.prepared_maps = MappedPreparations(self.cache)
        self.local_cache = {}
        if not self.cache.has('initial'):
            initial_density = np.asarray(self.case.arrays['density_bohr3'])
            saved = np.asarray(self.case.arrays['saved_density_bohr3'])
            t = perf_counter()
            potential = self.case.background.apply_kernel(initial_density-saved)
            self.cache.put('initial', dict(full_density=initial_density, saved_density=saved,
                initial_potential=potential, full_positive_gradient=self.case.arrays['positive_gradient_ry_bohr'],
                full_field_gradient=self.case.arrays['field_gradient']),
                metadata=dict(aggregate=self.case.aggregate_binding, original_full_initial_totals=True),
                timings=dict(initial_kernel_seconds=perf_counter()-t))
        if not self.cache.has('mechanical'):
            self.cache.put('mechanical', mechanical_arrays(self.case), metadata=dict(
                authority='matched_clamped_case.build_clamped original periodic scalar; DeltaC exactly zero'))
        self.initial_potential = self.cache.get('initial')['initial_potential']
        self.check_budget()

    def check_budget(self):
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
        if perf_counter()-self.started > self.args.seconds or rss > self.args.maxrss_mib:
            raise TimeoutError('projection stage budget reached; completed cache preserved')

    def _admit_basis(self, index):
        slot = self.slots[index]; record = self.records[slot_identity(slot)]
        t = perf_counter()
        authenticate(record.frame); authenticate(record.coefficient)
        psi, raw_d = self.original.arrays(record)
        frame = self.case.grid.load_frame(slot)
        if (frame.basis.identity != record.basis_identity or frame.gauge_identity != record.gauge_identity
                or frame.coefficients.shape != record.shape or slot.source_binding['sha256'] != record.source_identity
                or record.component != slot.ispin-1 or record.weight != slot.weight
                or not np.array_equal(psi, frame.coefficients)):
            raise ValueError('original raw frame values/basis/gauge/source/measure differ')
        # Independent source bytes establish basis geometry; do not retain WFC buffers.
        basis = replace(frame.basis, reciprocal_bohr_inv=frame.basis.reciprocal_bohr_inv.copy(),
                        source_miller=frame.basis.source_miller.copy())
        del frame
        state_binding = self.case.source_states[(slot.ispin, slot.source_ik)]
        state = bound_json(state_binding)
        source_arrays = bound_arrays(state['archive'])
        gram = source_arrays['gram'].conj() if slot.orbital_conjugation else source_arrays['gram']
        expected_d = source_arrays['density_coefficients'].conj() if slot.orbital_conjugation else source_arrays['density_coefficients']
        if not np.array_equal(raw_d, expected_d):
            raise ValueError('actual original D differs from source-state transport')
        active = prepare_active_basis(psi, raw_d, slot.eigenvalues_ry, self.interval, self.args.bands,
            gram=gram, populations=source_arrays['populations'])
        metadata = dict(slot=index, identity=slot_identity(slot), frame=record.frame,
            raw_d=record.coefficient, source=dict(slot.source_binding), source_state=state_binding,
            source_archive=state['archive'], basis_identity=basis.identity, frame_identity=record.gauge_identity,
            orbital_conjugation=slot.orbital_conjugation, operation=slot.operation.index,
            reciprocal_shift=slot.reciprocal_shift.tolist(), selected=active.selected.tolist(),
            source_admission=dict(basis.source_admission),
            diagnostics=active.diagnostics)
        return basis, active, metadata, perf_counter()-t

    def _export_slot(self, index):
        """Admit raw data and build physical factors once, not once per edge."""
        self.check_budget()
        basis, active, metadata, basis_seconds = self._admit_basis(index)
        t = perf_counter()
        h0 = physical_hamiltonian(self.case.background, basis, basis.slot.ispin)
        key = 'slot/%03d' % index
        if not self.cache.has(key):
            vertices = PhysicalFieldVertices(h0)
            frozen = PhysicalFrozenAction(h0, potential=self.initial_potential[:, basis.slot.ispin-1])
            hp = active.columns.conj().T @ frozen.apply(active.columns)
            fields, reverse = [], []
            for axis in range(3):
                action = vertices.prepare(np.eye(3)[axis])
                fields.append(active.columns.conj().T @ action.apply(active.columns))
                reverse.append((active.columns.conj().T @ action.adjoint_apply(active.columns)).T)
            diagnostic = dict(metadata['diagnostics'],
                h0_rho0_commutator=float(np.linalg.norm(hp @ active.rho0-active.rho0 @ hp)),
                h0_hermiticity=float(np.linalg.norm(hp-hp.conj().T)))
            self.cache.put(key, dict(h0=hp, rho0=active.rho0, field=np.asarray(fields),
                field_reverse=np.asarray(reverse), selected_labels=active.labels, selected=active.selected),
                metadata=dict(metadata, diagnostics=diagnostic),
                timings=dict(basis_seconds=basis_seconds, operators_seconds=perf_counter()-t,
                             h0_actions=1, field_forward_actions=3, field_adjoint_actions=3))
        momenta = (basis.reciprocal_bohr_inv @ basis.miller+basis.k_cart_bohr_inv[:, None]).T
        overlap_start = perf_counter()
        overlaps = projector_overlaps(h0.factors, active.columns, momenta)
        prepared = dict(basis_arrays(basis), active_columns=active.columns, rho0=active.rho0,
                        labels=active.labels, selected=active.selected, **overlap_arrays(overlaps))
        self.cache.put('preparation/%03d' % index, prepared,
            metadata=dict(metadata, overlap_names=sorted(overlaps)),
            timings=dict(native_admissions=1, full_factor_builds=1,
                         overlap_seconds=perf_counter()-overlap_start,
                         total_preparation_seconds=basis_seconds+perf_counter()-t))
        self.check_budget()

    def basis(self, index):
        key = 'preparation/%03d' % index
        if not self.cache.has(key):
            self._export_slot(index)
        arrays, metadata = self.prepared_maps(index)
        basis = restore_basis(arrays, metadata, self.slots[index], self.case.grid)
        active = ActiveBasis(arrays['active_columns'], arrays['rho0'], arrays['labels'],
                             arrays['selected'], metadata['diagnostics'])
        return basis, active, metadata, 0.

    def source(self, index):
        if self.last_index != index:
            self.last_bundle = None
            self.last_bundle = self.basis(index)
            self.last_index = index
        return self.last_bundle

    def edge(self, iq, index, witness=False):
        key = 'edge/%02d/%03d' % (iq, index)
        if self.cache.has(key):
            return
        basis, active, source_metadata, _ = self.source(index)
        target_index = self.lookup['target'][iq][index]
        t = perf_counter()
        target_preparation_created = (target_index != index and
            not self.cache.has('preparation/%03d' % target_index))
        if target_index == index:
            target_basis, target_active, target_metadata = basis, active, source_metadata
        else:
            target_basis, target_active, target_metadata, _ = self.basis(target_index)
        source_arrays, _ = self.prepared_maps(index)
        target_arrays, _ = self.prepared_maps(target_index)
        source_overlaps = restore_overlaps(source_arrays, source_metadata['overlap_names'])
        target_overlaps = restore_overlaps(target_arrays, target_metadata['overlap_names'])
        q = np.asarray(self.lookup['q_indices'][iq], float)/4
        q -= np.floor(q+.5)
        edge = RealBlochEdge.create(self.case.control.reference, basis, target_basis, q)
        minus_q = np.asarray(self.lookup['q_indices'][self.lookup['minus'][iq]], float)/4
        minus_q -= np.floor(minus_q+.5)
        coupling = prepare_ionic_edge(self.case.background.inputs, edge, source_overlaps,
                                     target_overlaps, minus_q, local_cache=self.local_cache)
        setup_seconds = perf_counter()-t
        gs, partner = q_feature_geometry(iq, self.args.features)
        arrays, timings = project_edge(coupling, active.columns, target_active.columns, gs)
        arrays.update(feature_g=gs, feature_partner=partner)
        metadata = dict(source=index, target=target_index, q_index=iq, q_fractional=q.tolist(),
                        counterpart_q=self.lookup['minus'][iq], counterpart_source=target_index,
                        reciprocal_fold=edge.reciprocal_fold.tolist(),
                        source_basis=basis.identity, target_basis=target_basis.identity,
                        source_binding=self.records[slot_identity(self.slots[index])].frame,
                        target_binding=self.records[slot_identity(self.slots[target_index])].frame,
                        source_d=self.records[slot_identity(self.slots[index])].coefficient,
                        target_d=self.records[slot_identity(self.slots[target_index])].coefficient,
                        feature_g=gs.tolist(), feature_partner=partner.tolist(),
                        local_support=coupling.support_metadata,
                        target_preparation_created=target_preparation_created,
                        target_full_factor_builds=int(target_preparation_created),
                        target_native_admissions=int(target_preparation_created))
        if witness:
            witness_start = perf_counter()
            h0 = physical_hamiltonian(self.case.background, basis, basis.slot.ispin)
            if target_index == index:
                target_factors = h0.factors
            else:
                from substrate_framework.material_physical_operators import build_physical_factors
                target_factors = build_physical_factors(self.case.background.inputs, target_basis)
            physical = physical_finiteq(self.case.background.inputs, edge, h0.factors, target_factors)
            attach_paired_local_support(physical, minus_q)
            fields = self.cache.get('slot/%03d' % index)['field']
            metadata['physical_witness'] = projection_witness(physical, active.columns, target_active.columns,
                arrays, PhysicalFieldVertices(h0), fields, gs)
            metadata['physical_witness']['full_factor_rebuilds'] = 1+int(target_index != index)
            timings['witness_seconds'] = perf_counter()-witness_start
        self.cache.put(key, arrays, metadata=metadata, timings=dict(timings, setup_seconds=setup_seconds))
        self.check_budget()

    def kernel(self, iq):
        key = 'kernel/%02d' % iq
        if self.cache.has(key):
            return
        q = np.asarray(self.lookup['q_indices'][iq], float)/4
        q -= np.floor(q+.5)
        gs, partner = q_feature_geometry(iq, self.args.features)
        kernel, timings = project_kernel(self.case.background,
            self.case.control.reference.reciprocal_bohr_inv @ q, gs)
        self.cache.put(key, dict(kernel=kernel, feature_g=gs, feature_partner=partner),
            metadata=dict(q_fractional=q.tolist()), timings=timings)
        self.check_budget()

    def finish(self):
        gamma_required = ['slot/%03d' % s for s in range(128)]
        gamma_required += ['edge/00/%03d' % s for s in range(128)]+['kernel/00']
        required = gamma_required+['edge/%02d/%03d' % (q, s) for q in range(64) for s in range(128)]
        required += ['kernel/%02d' % q for q in range(64)]
        if not all(k in self.cache.manifest['items'] for k in gamma_required):
            return
        # Explicit frozen passive fields are derived from the original totals,
        # not silently regenerated/refilled, for both nested electronic models.
        for m in (8, 12):
            if m > self.args.bands or self.cache.has('passive/%02d' % m):
                continue
            active_density = np.zeros_like(self.case.arrays['density_bohr3'], dtype=complex)
            active_gradient = np.zeros(len(self.case.reference.masses), complex)
            active_field = np.zeros(3, complex)
            for chunk in range(16):
                chunk_key = 'passive-chunk/%02d/%02d' % (m, chunk)
                if not self.cache.has(chunk_key):
                    density = np.zeros_like(active_density)
                    gradient = np.zeros_like(active_gradient)
                    field = np.zeros_like(active_field)
                    for s in range(8*chunk, 8*(chunk+1)):
                        self.check_budget()
                        b, a, _, _ = self.basis(s)
                        edge = RealBlochEdge.create(self.case.control.reference, b, b, np.zeros(3))
                        bridge = BlochFiniteQBridge(self.case.background.inputs, edge)
                        fixed = self.cache.get('slot/%03d' % s)
                        rho = fixed['rho0'][:m, :m]
                        density += bridge.density_from_dyads(a.columns[:, :m], a.columns[:, :m], rho, weight=1/64)
                        ionic = self.cache.get('edge/00/%03d' % s)['ionic_reverse'][:, :m, :m]
                        gradient += np.einsum('iab,ab->i', ionic, rho)/64
                        field += np.einsum('iab,ab->i', fixed['field_reverse'][:, :m, :m], rho)/64
                    self.cache.put(chunk_key, dict(density=density, gradient=gradient, field=field),
                                   metadata=dict(bands=m, slots=list(range(8*chunk, 8*(chunk+1)))))
                chunk_data = self.cache.get(chunk_key)
                active_density += chunk_data['density']
                active_gradient += chunk_data['gradient']
                active_field += chunk_data['field']
            initial = self.cache.get('initial')
            self.cache.put('passive/%02d' % m, dict(
                density=initial['full_density']-active_density,
                positive_gradient=initial['full_positive_gradient'].reshape(-1)-active_gradient,
                field_gradient=initial['full_field_gradient']-active_field,
                active_density=active_density, active_positive_gradient=active_gradient,
                active_field_gradient=active_field), metadata=dict(bands=m, frozen=True,
                density_split_defect=float(np.linalg.norm((initial['full_density']-active_density)+active_density-initial['full_density']))))
        self.cache.manifest['status'] = ('complete' if all(k in self.cache.manifest['items'] for k in required)
                                         else 'gamma-complete')
        self.cache.publish()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    export = sub.add_parser('export', help='resume arbitrary explicit slot/q projection')
    export.add_argument('--case', required=True)
    export.add_argument('--boundary', required=True)
    export.add_argument('--cache', required=True)
    export.add_argument('--reuse-gamma-pilot', help='explicit unchanged Gamma reuse into a NEW cache; preserve original producer snapshots')
    export.add_argument('--slots', default='0', help='all or comma/ranges in spin-major k-minor indexing')
    export.add_argument('--q', default='0', help='all or comma/ranges in np.indices((4,4,4)) order')
    export.add_argument('--bands', type=int, choices=(8, 12), default=12)
    export.add_argument('--features', type=int, choices=(7, 27), default=27)
    export.add_argument('--seconds', type=float, default=600.)
    export.add_argument('--maxrss-mib', type=float, default=3072.)
    export.add_argument('--witness', action='store_true')
    export.add_argument('--kernels', action='store_true', help='also export each selected q kernel (54 actions for27G)')
    inspect = sub.add_parser('inspect', help='inspect only manifest metadata and actual completed item coverage')
    inspect.add_argument('manifest')
    args = parser.parse_args(argv)
    if args.command == 'inspect':
        manifest = json.loads(Path(args.manifest).read_text())
        print(json.dumps(dict(status=manifest['status'], specification=manifest['specification'],
            completed_items=len(manifest['items']), timings=manifest['timings'], failures=manifest['failures']), indent=2))
        return 0
    exporter = None
    try:
        exporter = MaterialExporter(args)
        for q in selected(args.q, 64):
            for s in selected(args.slots, 128):
                exporter.edge(q, s, witness=args.witness)
        if args.kernels:
            exporter.last_bundle = None; exporter.last_index = None
            exporter.local_cache.clear()
            for q in selected(args.q, 64):
                exporter.kernel(q)
        exporter.finish()
        exporter.cache.manifest.setdefault('runs', []).append(dict(
            requested_slots=args.slots, requested_q=args.q,
            elapsed_seconds=perf_counter()-exporter.started,
            maxrss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
            completed_items=len(exporter.cache.manifest['items'])))
        exporter.cache.publish()
        print(str(exporter.cache.path))
        return 0
    except (TimeoutError, MemoryError) as exc:
        if exporter is not None:
            exporter.cache.manifest['failures'].append(dict(kind=type(exc).__name__, message=str(exc),
                elapsed_seconds=perf_counter()-exporter.started,
                maxrss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
            exporter.cache.publish()
        print(str(exc), file=sys.stderr)
        return 3


if __name__ == '__main__':
    raise SystemExit(main())
