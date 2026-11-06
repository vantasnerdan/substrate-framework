"""Deterministic real finite-input provenance admission/refusal regressions.

Run with PYTHONPATH=src:proposals/P254-ni-collective-response/attempts/0029:proposals/P254-ni-collective-response/attempts/0028.
All files are temporary; no campaign result, source archive or process is used.
"""
from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from clamped_scientific_readout import encode
from oracle_reduced_electronic_dynamics import fixture
from reduced_comparison_inputs import digest, read_finite
from substrate_framework.streamed_material_state import file_binding


class FiniteSourceRecovery(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.source = self.root/'consumer.py'
        self.old_bytes = b'original_value = 17\n'
        self.source.write_bytes(self.old_bytes)
        self.original = file_binding(self.source)
        self.original['bytes'] = str(self.original['bytes'])
        self.snapshot = self.history(self.original)
        self.snapshot.parent.mkdir()
        self.snapshot.write_bytes(self.old_bytes)

        # Reuse the independent dynamics oracle's deterministic moving baseline
        # and ionic preparation, extending the periodic geometry to the required
        # 64 q / 128 explicit spin slots and three global field coordinates.
        c, _, baseline, _ = fixture()
        mechanical = c.mechanical
        lattice = np.zeros((4, 4, 4, 3, 3))
        lattice[0, 0, 0] = np.diag([.7, .9, 1.2])
        lattice[1, 0, 0] = lattice[3, 0, 0] = -.06*np.eye(3)
        self.masses = mechanical.masses.real
        prepared = dict(force_lattice=lattice, masses=self.masses, charges=np.array([1.3]),
            delta_curvature=np.zeros((3, 3)), contact=np.diag([.08, .07, .06]),
            mixed=np.full((3, 3), .011), force=-mechanical.linear[:3].real/np.sqrt(c.n),
            field_gradient=np.array([.012, -.003, .002]), polarizations=np.eye(3),
            field_mass=np.array(1.7), ionic_ordered=mechanical.w,
            field_ordered=np.eye(6)*.5, trap_omega=mechanical.omega.real)
        mechanical_binding = self.arrays('mechanical.npz', **prepared)
        geometry = dict(atom_indices=np.array([1]), atom_names=np.array(['Ni']),
                        cell_bohr=np.eye(3)*5, valence_charges=np.array([1.3]))
        geometry_binding = self.arrays('geometry.npz', **geometry)
        aggregate_binding = self.json('aggregate.json', dict(archive=geometry_binding))
        boundary = self.json('boundary.py', dict(preparation='original raw electronic boundary'))
        self.context = dict(aggregate=aggregate_binding, original_raw=boundary)
        cells = np.indices((4, 4, 4)).reshape(3, -1).T
        targets = [np.ravel_multi_index(tuple(((cells+q) % 4).T), (4, 4, 4)) for q in cells]
        indexing = dict(q_indices=cells.tolist(), k_indices=cells.tolist(),
            slot_spin=np.repeat([1, 2], 64).tolist(), slot_k=np.tile(np.arange(64), 2).tolist(),
            target=[np.r_[row, row+64].tolist() for row in targets],
            minus=np.ravel_multi_index(tuple((-cells % 4).T), (4, 4, 4)).tolist())
        spec = dict(case='empty-vacancy-ni31', context=self.context, case_identity=digest(self.context),
            original_boundary=boundary, implementation=[deepcopy(self.original)],
            physical_implementation=[], mesh=[4, 4, 4], indexing=indexing)
        # Coefficients are authenticated, not reconstructed by this read-only
        # consumer. Repeated item identities deliberately share immutable bytes.
        item = dict(archive=self.arrays('coefficient.npz', value=np.array([1.])))
        items = {key: deepcopy(item) for key in ('initial', 'passive/08')}
        items['mechanical'] = dict(archive=mechanical_binding)
        for q in range(64):
            items['kernel/%02d' % q] = deepcopy(item)
            for slot in range(128):
                items['edge/%02d/%03d' % (q, slot)] = deepcopy(item)
        for slot in range(128):
            for prefix in ('slot', 'preparation'):
                items['%s/%03d' % (prefix, slot)] = deepcopy(item)
        self.manifest = dict(schema='p254-physical-activeband-cache-v1', status='complete',
                             specification=spec, items=items)
        x = np.r_[baseline.x[:6], np.zeros(6)]
        self.baseline_x = x
        rho = np.broadcast_to(np.diag(np.linspace(.9, .1, 8)), (128, 8, 8)).copy()
        baseline_binding = self.arrays('baseline.npz', rho=rho, x=x)
        self.amplitude = .2
        finite = dict(rho=np.zeros((2, 512, 512), dtype=complex))
        for q in range(64):
            delta = np.zeros(12 if q == 0 else 6, complex)
            delta[3] = self.amplitude/8
            finite['x_%02d' % q] = delta
        finite_binding = self.arrays('finite.npz', **finite)
        config = dict(coefficient_specification=spec, baseline_only=False, bands=8, features=7,
            geometry_archive=geometry_binding, atom=1, axis='x', times=[0.], correlations=False,
            correlation_reference=0., amplitude=self.amplitude, rtol=1e-8, atol=1e-10)
        self.result = dict(schema='P254-activeband-material-run-v1', status='complete',
            configuration=config, consumer_sources=[deepcopy(self.original)],
            completed=[dict(time_inverse_ry=0., baseline_state=baseline_binding,
                finite_state=finite_binding, energy_change_ry={key: 0. for key in
                    ('mechanical', 'band', 'interaction', 'screening', 'total')})])
        self.result_path = self.root/'result.json'
        self.publish()

    def json(self, name, value):
        path = self.root/name
        path.write_text(json.dumps(encode(value), sort_keys=True))
        return file_binding(path)

    def arrays(self, name, **arrays):
        path = self.root/name
        np.savez(path, **arrays)
        return file_binding(path)

    @staticmethod
    def history(binding):
        path = Path(binding['path'])
        return path.parent/'source-history'/(path.stem+'-'+binding['sha256']+path.suffix)

    def publish(self):
        manifest = self.json('manifest.json', self.manifest)
        self.result['original_input_manifest'] = manifest
        self.result['original_coefficient_manifest_snapshot'] = deepcopy(manifest)
        self.result['identity'] = digest(self.result['configuration'])
        self.json('result.json', self.result)

    def changed(self):
        self.source.write_bytes(b'new_code_value = 19\n')

    def refuse(self):
        with self.assertRaises((ValueError, FileNotFoundError)):
            read_finite(self.result_path)

    def test_exact_recovery_preserves_result_and_physical_identity(self):
        before = read_finite(self.result_path)
        immutable = file_binding(self.result_path)
        self.assertEqual(before['source_recoveries'], [])
        self.changed()
        after = read_finite(self.result_path)
        self.assertEqual(file_binding(self.result_path), immutable)
        self.assertEqual(after['identity'], before['identity'])
        self.assertEqual(after['model_source_identity'], before['model_source_identity'])
        receipt, = after['source_recoveries']  # implementation + consumer deduplicate
        self.assertEqual(receipt['original_binding'], self.original)
        self.assertEqual(receipt['archived_snapshot'], dict(path=str(self.snapshot),
            sha256=self.original['sha256'], bytes=self.original['bytes']))
        self.assertEqual(receipt['recovery_timing'],
                         'retrospective source recovery, not original execution binding')
        expected = (self.amplitude*self.baseline_x[3]/8+self.amplitude**2/2)/self.masses[0]
        for run in (before, after):
            kinetic = run['observations'][0]['finite_kinetic_exchange']
            np.testing.assert_allclose(kinetic['ionic_total_ry'], expected, atol=1e-15, rtol=1e-14)
            by_cell = kinetic['ionic_by_cell_atom_ry'][:, 0]
            np.testing.assert_allclose(by_cell[0], expected, atol=1e-15, rtol=1e-14)
            np.testing.assert_allclose(by_cell[1:], 0., atol=1e-15, rtol=0)


    def test_integer_byte_count_identity_survives_recovery(self):
        for source in (self.manifest['specification']['implementation'][0], self.result['consumer_sources'][0]):
            source['bytes'] = int(source['bytes'])
        self.publish()
        self.changed()
        receipt, = read_finite(self.result_path)['source_recoveries']
        self.assertIsInstance(receipt['original_binding']['bytes'], int)
        self.assertEqual(receipt['original_binding']['bytes'], len(self.old_bytes))

    def test_missing_live_source_uses_exact_history(self):
        self.source.unlink()
        self.assertEqual(read_finite(self.result_path)['source_recoveries'][0]['original_binding'], self.original)

    def test_live_exact_source_does_not_consume_unused_bad_snapshot(self):
        self.snapshot.write_bytes(b'wrong')
        self.assertEqual(read_finite(self.result_path)['source_recoveries'], [])

    def test_same_size_changed_snapshot_refuses(self):
        self.changed()
        self.snapshot.write_bytes(b'x'*len(self.old_bytes))
        self.refuse()

    def test_truncated_snapshot_refuses(self):
        self.changed()
        self.snapshot.write_bytes(self.old_bytes[:-1])
        self.refuse()

    def test_missing_snapshot_refuses(self):
        self.changed()
        self.snapshot.unlink()
        self.refuse()

    def test_wrong_bound_size_refuses_even_with_correct_digest(self):
        self.changed()
        for source in (self.manifest['specification']['implementation'][0], self.result['consumer_sources'][0]):
            source['bytes'] = str(int(source['bytes'])+1)
        self.publish()
        self.refuse()

    def test_wrong_bound_digest_refuses_despite_matching_filename(self):
        self.changed()
        for source in (self.manifest['specification']['implementation'][0], self.result['consumer_sources'][0]):
            source['sha256'] = '0'*64
        self.history(self.result['consumer_sources'][0]).write_bytes(self.old_bytes)
        self.publish()
        self.refuse()

    def test_explicit_relocated_snapshot_retains_original_identity(self):
        self.changed()
        relocated = self.root/'explicit-original.py'
        relocated.write_bytes(self.old_bytes)
        self.snapshot.unlink()
        self.result['source_history'] = [dict(original_binding=deepcopy(self.original),
                                             archived_snapshot=file_binding(relocated))]
        self.publish()
        receipt, = read_finite(self.result_path)['source_recoveries']
        self.assertEqual(receipt['original_binding'], self.original)
        self.assertEqual(receipt['archived_snapshot'], file_binding(relocated))
        self.assertEqual(receipt['recovery_method'], 'explicit original_binding/archived_snapshot pair')

    def test_explicit_wrong_size_refuses_without_canonical_fallback(self):
        self.changed()
        saved = file_binding(self.snapshot)
        saved['bytes'] += 1
        self.result['source_history'] = [dict(original_binding=self.original, archived_snapshot=saved)]
        self.publish()
        self.refuse()

    def test_noncode_consumer_has_no_canonical_recovery(self):
        raw = self.root/'consumer.npz'
        raw.write_bytes(self.old_bytes)
        original = file_binding(raw)
        self.result['consumer_sources'] = [original]
        self.history(original).write_bytes(self.old_bytes)
        self.publish()
        raw.unlink()
        self.refuse()

    def test_raw_boundary_never_recovers_even_with_python_suffix(self):
        self.changed()
        binding = self.manifest['specification']['original_boundary']
        raw = Path(binding['path'])
        self.history(binding).write_bytes(raw.read_bytes())
        raw.write_text('changed raw scientific input')
        self.refuse()

    def test_state_bytes_never_recover(self):
        self.changed()
        binding = self.result['completed'][0]['baseline_state']
        raw = Path(binding['path'])
        self.history(binding).write_bytes(raw.read_bytes())
        raw.write_bytes(b'truncated state')
        self.refuse()

    def test_recovery_does_not_relax_q_spin_or_preparation(self):
        self.changed()
        for key in ('edge/27/127', 'kernel/27', 'slot/127', 'preparation/127'):
            with self.subTest(missing=key):
                item = self.manifest['items'].pop(key)
                self.publish()
                self.refuse()
                self.manifest['items'][key] = item
        self.manifest['specification']['indexing']['slot_spin'][64] = 1
        self.publish()
        self.refuse()


if __name__ == '__main__':
    unittest.main()
