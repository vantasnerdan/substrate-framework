"""Authenticate existing Ni/H inputs and construct actual Gaussian two-time legs.

No mean trajectory/preparation/DFT is regenerated. FrozenReference explicitly
uses fixed ORIGINAL h0 and clamped Cartesian fluctuation generators, not rho(t)
as a surrogate for U(t,s). A saved transition bundle can instead supply moving
U(t,0) and S_q(t,0), with all original axes and byte-bound provenance.
"""
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0029')]
from substrate_framework.streamed_material_state import authenticate, file_binding
from substrate_framework.material_clamped_array import ClampedArray
from reduced_coefficients import ReducedCoefficients, QCoefficients
from reduced_projection_algebra import q_feature_geometry
from reduced_projection_storage import ImmutableItems
from electron_ion_functional import EdgePair

MANIFEST_SHA = 'f851d2df2b2f1136d22b606f059cf6cd9bf6df7c4681340153e12c8569dcd479'
MANIFEST_BYTES = 59561233
DEFAULT_MANIFEST = '/var/tmp/p254-activeband-loaded-fullq-union-001/manifest.json'


class ReadOnlyQBanks:
    """Decode ONE original q bank in RAM, never assemble/copy a disk bank.

    Same original immutable NPZ members and QCoefficients axes as A29's loader.
    Full12/27 only: no projection, subset, weight, conjugation or physical rule
    is introduced. This read-only path avoids MappedQBanks' absent-cache writer.
    """
    def __init__(self, reader, specification):
        self.reader, self.spec = reader, specification
        self.index, self.value = None, None

    def __call__(self, q):
        q = int(q)
        if not 0 <= q < 64:
            raise IndexError(q)
        if self.index == q:
            return self.value
        self.index, self.value = None, None
        gs, partner = q_feature_geometry(q,27)
        names = ('ionic','ionic_reverse','density','potential')
        arrays = None
        for source in range(128):
            key = 'edge/%02d/%03d' % (q,source)
            raw = self.reader.get(key,(*names,'feature_g','feature_partner'))
            metadata = self.reader.items[key]['metadata']
            if (metadata['source'] != source or metadata['q_index'] != q
                    or metadata['target'] != self.spec['indexing']['target'][q][source]
                    or not np.array_equal(raw['feature_g'],gs)
                    or not np.array_equal(raw['feature_partner'],partner)):
                raise ValueError('original directed edge/feature map differs: '+key)
            if arrays is None:
                arrays = {name: np.empty((128,*raw[name].shape),dtype=raw[name].dtype) for name in names}
            for name in names:
                arrays[name][source] = raw[name]
        kernel = self.reader.get('kernel/%02d' % q,('kernel','feature_g','feature_partner'))
        if not np.array_equal(kernel['feature_g'],gs) or not np.array_equal(kernel['feature_partner'],partner):
            raise ValueError('kernel q partner geometry differs')
        arrays.update(kernel=kernel['kernel'],feature_g=gs,feature_partner=partner,
                      feature_spin=np.repeat([1,2],len(gs)))
        self.value = QCoefficients(**arrays)
        self.index = q
        return self.value


def load_material(path=DEFAULT_MANIFEST):
    """Full manifest, every original archive, original Gaussian and all128/64 maps.

    Hash authentication is streaming, no persistent material cache copy. The
    read-only provider owns one original full q decode, discarded before next q.
    Original exact JSON integer metadata stays Python int throughout.
    """
    binding = file_binding(path)
    if binding['sha256'] != MANIFEST_SHA or binding['bytes'] != MANIFEST_BYTES:
        raise ValueError('not the pinned complete full12band27feature material manifest')
    manifest = json.loads(Path(path).read_bytes())
    spec = manifest['specification']
    if (manifest['status'] != 'complete' or manifest['failures'] or len(manifest['items']) != 8548
            or (spec['case'],spec['bands'],spec['features'],spec['ndof'],spec['mesh'])
               != ('loaded-ni31h6',12,27,111,[4,4,4])
            or spec['original_source_weight_in_rho'] is not False):
        raise ValueError('full material/state/normalization contract differs')
    reader = ImmutableItems(manifest['items'])
    for item in manifest['items'].values():
        reader.authenticate(item['archive'])
    authenticate(spec['original_boundary'])
    for source in spec['implementation']+spec['physical_implementation']:
        authenticate(source)
    if 'passive/12' not in manifest['items']:
        raise ValueError('original active/passive split is absent')
    fixed = [reader.get('slot/%03d' % s,('h0','rho0','field','field_reverse')) for s in range(128)]
    mechanical_arrays = reader.get('mechanical')
    mechanical = ClampedArray(**mechanical_arrays)
    indexing = spec['indexing']
    c = ReducedCoefficients(mesh=tuple(spec['mesh']),
        **{name: np.asarray(indexing[name]) for name in
           ('k_indices','q_indices','slot_spin','slot_k','target','minus')},
        **{name: np.asarray([record[name] for record in fixed]) for name in
           ('h0','rho0','field','field_reverse')},
        mechanical=mechanical,q_provider=ReadOnlyQBanks(reader,spec),
        field_polarizations=mechanical_arrays['polarizations'])
    if (c.ns,c.nq,c.d,c.m,c.nspin,c.n) != (128,64,111,12,2,64):
        raise ValueError('original full axis inventory differs')
    return c, manifest, binding


class FrozenReference:
    """Conditional, nonstationary original Gaussian on fixed physical reference.

    Electron mean h0 is frozen; original coherent rho0 is NOT stationary-refilled.
    Cartesian D is the connected clamped ionic/field product preparation evolved
    by the original raw quadratic generator (unstable directions retained).
    This is a reference-kernel evaluation, NOT a moving/selfconsistent solution.
    """
    label = 'conditional-frozen-physical-reference-original-coherent-Gaussian'

    def __init__(self, coefficients, times):
        self.c = coefficients
        self.times = np.asarray(times, float)
        self.u = np.asarray([[np.eye(coefficients.m) if t == 0 else expm(-1j*t*h)
                             for h in coefficients.h0] for t in self.times])
        self.bindings = []
        self._q = None
        self._transitions = None

    def transitions(self, q):
        if self._q != q:
            qm = int(self.c.minus[q])
            def sequence(index):
                a = self.c.mechanical.generators[index]
                return np.asarray([np.eye(len(a)) if t == 0 else expm(t*a) for t in self.times])
            sq = sequence(q)
            sm = sq if q == qm else sequence(qm)
            self._q, self._transitions = q, (sq,sm)
        return self._transitions

    def electron(self, a, b):
        """All explicit spin/k blocks G(t_a,t_b), no sector weights."""
        ua, ub = self.u[a], self.u[b]
        end = ub.conj().swapaxes(-1,-2)
        rho = ua@self.c.rho0@end
        empty = ua@(np.eye(self.c.m)-self.c.rho0)@end
        return -1j*empty, 1j*rho

    def boson(self, q, a, b):
        sq, sm = self.transitions(q)
        qm, d = int(self.c.minus[q]), self.c.d
        # Independent ordered reverse leg; transpose is not conjugation.
        w = sq[a]@self.c.mechanical.w0[q]@sm[b].T
        reverse = sm[b]@self.c.mechanical.w0[qm]@sq[a].T
        return -1j*w[:d,:d], -1j*reverse[:d,:d].T


class SavedTransitions(FrozenReference):
    """Read actual retained transitions, never reconstruct U from density alone.

    JSON schema p254-electron-ion-transitions-v1: coefficient_manifest binding,
    times_inverse_ry, archive binding, physical_history_bindings (nonempty),
    approximation (nonempty). NPZ: U[nt,128,12,12], S_00..S_63[nt,nx_q,nx_q].
    U and S must map literal time0 prepared coordinates, U(0)=S(0)=I. S maps the
    original FULL Cartesian phase space, including Gamma field coordinates.
    Caller owns genuine propagation and its accuracy evidence; this reader does
    not promote a provenance receipt into selfconsistency/accuracy evidence.
    """
    label = 'conditional-saved-moving-transition-reference-original-Gaussian'

    def __init__(self, coefficients, times, bundle, coefficient_binding):
        self.c, self.times = coefficients, np.asarray(times,float)
        record_binding = file_binding(bundle)
        spec = json.loads(Path(bundle).read_bytes())
        if (spec['schema'] != 'p254-electron-ion-transitions-v1'
                or spec['coefficient_manifest'] != coefficient_binding
                or spec['times_inverse_ry'] != self.times.tolist()
                or not spec['physical_history_bindings'] or not spec['approximation']):
            raise ValueError('saved transitions lack exact state/time/material/history binding')
        for binding in spec['physical_history_bindings']:
            authenticate(binding)
        self._archive = np.load(authenticate(spec['archive']),allow_pickle=False)
        self.u = self._archive['U']
        if self.u.shape != (len(times),128,12,12) or not np.isfinite(self.u).all():
            raise ValueError('saved U must retain all128 sectors/all12 bands at every requested time')
        zero = np.flatnonzero(self.times == 0.)
        if len(zero) != 1 or not np.array_equal(self.u[zero[0]],np.broadcast_to(np.eye(12),(128,12,12))):
            raise ValueError('saved transitions must bind literal original time0 identity')
        self.bindings = [record_binding,spec['archive'],*spec['physical_history_bindings']]
        self._q = None
        self._transitions = None

    def transitions(self, q):
        if self._q != q:
            qm = int(self.c.minus[q])
            arrays = []
            for index in (q,qm):
                value = self._archive['S_%02d' % index]
                nx = self.c.mechanical_sizes[index]
                if (value.shape != (len(self.times),nx,nx) or not np.isfinite(value).all()
                        or not np.array_equal(value[np.flatnonzero(self.times == 0.)[0]],np.eye(nx))):
                    raise ValueError('saved full Cartesian transition shape/time0 differs')
                arrays.append(value)
            self._q, self._transitions = q, tuple(arrays)
        return self._transitions

    def close(self):
        self._archive.close()


def edge_pair(data, source, target, forward, backward, boson):
    """Bind BOTH exported physical legs without adjoint inference or rescaling."""
    gg, gl = forward
    bg, bl = backward
    dg, dl = boson
    left = data.ionic[source]
    right = data.ionic_reverse[source].transpose(0,2,1)
    return EdgePair(left,right,left,right,gg[source],gl[source],gg[target],gl[target],
                    bg[source],bl[source],dg,dl)


def reference_diagnostics(reference):
    c = reference.c
    eye = np.eye(c.m)
    return dict(
        rho0_hermiticity_defect=float(np.linalg.norm(c.rho0-c.rho0.conj().swapaxes(-1,-2))),
        rho0_h0_commutator_norm=float(np.linalg.norm(c.h0@c.rho0-c.rho0@c.h0)),
        electronic_transition_unitarity_defect=float(max(np.linalg.norm(u@u.conj().swapaxes(-1,-2)-eye)
            for u in reference.u)),
        mechanical=c.mechanical.coefficient_diagnostics(),
        state_policy='literal projected original coherent Gaussian; no refill/FDT/repair',
        passive_scope='frozen initial passive totals/contact retained; active-passive dynamical loops omitted')
