"""Exact-node bare bath cache; update only the genuine raw Gamma contribution.

No mean trajectory is produced. ZIP_STORED histories remain mapped, canonical
histories are read sectorwise, and bare equal-time nodes use saved U/S factors.
"""
from pathlib import Path
import json
import sys
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/n) for n in ('0063','0057')]
from history_storage import StoredNPZ, authenticate, binding
from saved_mean import SavedMean, BARE_SCHEMA, POSITIONS
from history_provider import SavedHistory, SCHEMA as PROVIDER_SCHEMA
from mean_feedback_functional import CompactCubic, SOURCE

SCHEMA = 'p254-fullwindow-bare-covariance-force-cache-v1'


def bare_nodes(mean):
    """A63 SavedMean equal-time operations only, without expanding pair histories."""
    nt = len(mean.times)
    rho = np.empty((nt,128,12,12), complex)
    w = np.empty((nt,114,114), complex)
    for a in range(nt):
        rho[a] = -1j*(1j*(mean.u[a]@mean.rho0@mean.u[a].conj().swapaxes(-1,-2)))
        phase = mean.s[a]@mean.w0@mean.s[a].T
        lesser = -1j*phase[np.ix_(POSITIONS,POSITIONS)].T
        greater = -1j*phase[np.ix_(POSITIONS,POSITIONS)]
        retarded = greater-lesser
        w[a] = 1j*(lesser+retarded)
    return rho,w


def gamma_force(source, covariance):
    qs = source.indices[[0,0,0]]
    return np.asarray([-.5*source.contract(qs,covariance=w[:111,:111]) for w in covariance])


def provider_origin(path, expected):
    """Follow authenticated actual previous inputs; no schema/prose identity alias."""
    seen = set()
    while True:
        item = binding(path)
        if item == expected:
            return
        if item['path'] in seen or Path(path).suffix == '.npz':
            raise ValueError('canonical provider does not descend from this bare source')
        seen.add(item['path'])
        spec = json.loads(Path(path).read_bytes())
        if spec.get('schema') != PROVIDER_SCHEMA or spec.get('status') != 'complete':
            raise ValueError('canonical provider source lineage differs')
        path = authenticate(spec['input_previous'])


def actual_nodes(dataset, bare_binding, times):
    """Return original raw equal-time data and bindings, never missing-time fills."""
    path = Path(dataset).resolve()
    with SavedMean(authenticate(bare_binding)) as bare:
        if not bare.is_bare or not np.array_equal(bare.times,times):
            raise ValueError('genuine exact-node original bare source required')
        if path.suffix != '.npz':
            spec = json.loads(path.read_bytes())
            if spec.get('schema') == BARE_SCHEMA:
                if binding(path) != bare_binding:
                    raise ValueError('different prescribed bare source identity')
                rho,w = bare_nodes(bare)
                return rho,w,bare.bindings
            provider_origin(path,bare_binding)
            history = SavedHistory(path)
            times_actual, aux = history.times, history.aux
        else:
            history = StoredNPZ(path)
            times_actual = history['times']
            aux = {k:history[k] for k in ('source_electron_less0','source_gamma_less','source_gamma_ret','n')}
        try:
            nt = len(times)
            if not np.array_equal(times_actual,times):
                raise ValueError('prescribed nonGamma exact nodes changed')
            if (int(aux['n']) != int(bare.archive['n']) or
                    not np.array_equal(aux['source_electron_less0'],1j*bare.rho0) or
                    aux['source_gamma_less'].shape != (nt,nt,114,114) or
                    aux['source_gamma_ret'].shape != (nt,nt,114,114) or
                    not np.array_equal(aux['source_gamma_less'][0,0],
                        -1j*bare.w0[np.ix_(POSITIONS,POSITIONS)].T)):
                raise ValueError('original raw source preparation/normalization differs')
            for a in range(nt):
                phase = bare.s[a]@bare.w0@bare.s[a].T
                expected_less = -1j*phase[np.ix_(POSITIONS,POSITIONS)].T
                expected_ret = -1j*phase[np.ix_(POSITIONS,POSITIONS)]-expected_less
                if (not np.array_equal(aux['source_gamma_less'][a,a],expected_less) or
                        not np.array_equal(aux['source_gamma_ret'][a,a],expected_ret)):
                    raise ValueError('raw original prescribed Gamma source nodes changed')
            rho = np.empty((nt,128,12,12),complex)
            if path.suffix == '.npz':
                el,gr,gl = history['electron_less'],history['gamma_ret'],history['gamma_less']
                if el.shape != (nt,nt,128,12,12):
                    raise ValueError('full original128 electron density required')
                for slot in range(128):
                    for a in range(nt):
                        rho[a,slot] = -1j*el[a,a,slot]
            else:
                for slot in range(128):
                    _,el = history.electron_slot(slot)
                    for a in range(nt):
                        rho[a,slot] = -1j*el[a,a]
                gr,gl = history.gamma_slot()
            if gr.shape != (nt,nt,114,114) or gl.shape != gr.shape:
                raise ValueError('full original114 Gamma covariance required')
            w = np.asarray([1j*(gl[a,a]+gr[a,a]) for a in range(nt)])
            if not np.isfinite(rho).all() or not np.isfinite(w).all():
                raise ValueError('nonfinite raw equal-time input')
            return rho,w,[binding(path),*getattr(history,'bindings',[]),bare_binding]
        finally:
            history.close()


class FullWindowForceCache:
    def __init__(self,result_path,coefficient_binding,*,historical_source_snapshots=()):
        self.path = Path(result_path).resolve()
        self.spec = json.loads(self.path.read_bytes())
        s = self.spec
        if (s['schema'] != SCHEMA or s['status'] != 'complete' or
                s['coefficient_manifest'] != coefficient_binding or
                s['completed_original_internal_q'] != list(range(64))):
            raise ValueError('complete all31/all64 source/material-bound force bank required')
        # Numeric/material/physical inputs always authenticate their actual paths.
        # Only original implementation bytes may be restored from explicit,
        # byte-identical archives; no digest alias or new-source substitution.
        actual = [binding(self.path),s['archive'],s['original_input_dataset'],
            s['coefficient_manifest'],s['compact_source'],s['compact_source_arrays'],
            *s['input_bindings']]
        for item in actual:
            authenticate(item)
        snapshots = {}
        for record in historical_source_snapshots:
            original = record['original_binding']
            archived = record['archived_snapshot']
            if (original not in s['implementation_bindings'] or original['path'] in snapshots or
                    any(original[k] != archived[k] for k in ('bytes','sha256'))):
                raise ValueError('historical snapshot must exactly preserve a cache implementation')
            authenticate(archived)
            snapshots[original['path']] = archived
        self.source_authentication = []
        restored = []
        for original in s['implementation_bindings']:
            archived = snapshots.get(original['path'])
            selected = original if archived is None else archived
            authenticate(selected)
            restored.append(selected)
            self.source_authentication.append(dict(original_binding=original,
                authenticated_binding=selected,restored_from_snapshot=archived is not None))
        self.bindings = [*actual,*restored]
        self.current_reader_binding = binding(Path(__file__))
        self._archive = StoredNPZ(s['archive']['path'])
        try:
            self.times = self._archive['times']
            self.phi000 = self._archive['Phi000']
            self.full_force = self._archive['bare_full_force']
            self.old_gamma_force = self._archive['bare_gamma_force']
            if (not np.array_equal(self.times,np.arange(31,dtype=float)) or
                    self.phi000.shape != (111,111,111) or
                    self.full_force.shape != (31,111) or self.old_gamma_force.shape != (31,111)):
                raise ValueError('complete full-window cache axes differ')
        except BaseException:
            self.close()
            raise

    def forcing_for(self,dataset,source):
        """A63 bare/canonical provider or original stored NPZ; no bath rebuild."""
        if (source.source_binding != self.spec['compact_source'] or
                source.source_record['arrays'] != self.spec['compact_source_arrays']):
            raise ValueError('cubic source changed: prescribed bath is not reusable')
        rho,w,inputs = actual_nodes(dataset,self.spec['original_input_dataset'],self.times)
        new_gamma = gamma_force(source,w)
        full = self.full_force+(new_gamma-self.old_gamma_force)
        return dict(times=self.times,rho=rho,ordered_gamma_covariance=w,
            cubic_covariance_force=full,gamma_force=new_gamma,input_dataset=binding(dataset),
            cache_bindings=self.bindings,input_bindings=inputs,
            nonGamma_contractions=0,Gamma_contractions=len(self.times))

    def close(self):
        self._archive.close()

    def __enter__(self):
        return self

    def __exit__(self,*exc):
        self.close()
