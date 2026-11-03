"""Numeric active-band contract; physical exporters own all operator conventions.

No equilibrium refill, coefficient repair, implicit spin weight, or q-bank copy is
performed here. A provider may keep a single memory-mapped q record resident.
"""
from dataclasses import dataclass
from typing import Callable, Protocol

import numpy as np


class MechanicalCoefficients(Protocol):
    d: int
    f: int
    n: int
    shape: tuple
    h: object
    j: object
    generators: object
    linear: np.ndarray
    w0: object


@dataclass(frozen=True)
class QCoefficients:
    """One immutable provider-owned q record with a local feature dimension.

    feature_partner[r] indexes the minus-q feature representing -(q+G_r),
    including spin. Screening energy pairs n_minus[feature_partner] with K n_q.
    feature_g contains the spin-independent integer G list; feature_spin labels
    each spin-major feature. Nominal 7/27 families can have larger Nyquist
    closures, so no consumer may infer a q-independent feature count.
    """
    ionic: np.ndarray
    ionic_reverse: np.ndarray
    density: np.ndarray
    potential: np.ndarray
    kernel: np.ndarray
    feature_partner: np.ndarray
    feature_g: np.ndarray | None = None
    feature_spin: np.ndarray | None = None


def _numeric(value, shape, name):
    array = np.asarray(value)
    if array.shape != shape or not np.issubdtype(array.dtype, np.number):
        raise ValueError(f'{name}: expected numeric {shape}, got {array.shape}')
    if not np.all(np.isfinite(array)):
        raise ValueError(f'{name}: nonfinite coefficient')
    return array


@dataclass(kw_only=True)
class ReducedCoefficients:
    mesh: tuple
    k_indices: np.ndarray
    q_indices: np.ndarray
    slot_spin: np.ndarray
    slot_k: np.ndarray
    target: np.ndarray
    minus: np.ndarray
    h0: np.ndarray
    rho0: np.ndarray
    field: np.ndarray
    field_reverse: np.ndarray
    mechanical: MechanicalCoefficients
    q_provider: Callable[[int], QCoefficients]
    field_polarizations: np.ndarray | None = None

    def __post_init__(self):
        self.mesh = tuple(int(v) for v in self.mesh)
        if len(self.mesh) != 3 or min(self.mesh) < 1:
            raise ValueError('three positive periodic mesh dimensions required')
        self.n = self.nk = self.nq = int(np.prod(self.mesh))
        expected = np.indices(self.mesh).reshape(3, -1).T
        self.k_indices = _numeric(self.k_indices, (self.n, 3), 'k_indices')
        self.q_indices = _numeric(self.q_indices, (self.n, 3), 'q_indices')
        if not np.array_equal(self.k_indices, expected) or not np.array_equal(self.q_indices, expected):
            raise ValueError('complete np.indices-order k/q meshes required')
        self.slot_spin = np.asarray(self.slot_spin)
        self.slot_k = np.asarray(self.slot_k)
        self.ns = len(self.slot_spin)
        self.spin_labels = tuple(dict.fromkeys(self.slot_spin.tolist()))
        self.nspin = len(self.spin_labels)
        if not self.nspin or self.ns != self.nspin*self.nk:
            raise ValueError('explicit complete spin/k slots required')
        if self.slot_k.shape != (self.ns,) or not np.issubdtype(self.slot_k.dtype, np.integer):
            raise ValueError('integer slot_k required')
        self.slots_by_spin = np.empty((self.nspin, self.nk), dtype=int)
        for i, spin in enumerate(self.spin_labels):
            slots = np.flatnonzero(self.slot_spin == spin)
            if len(slots) != self.nk or not np.array_equal(self.slot_k[slots], np.arange(self.nk)):
                raise ValueError('slots must be spin-major and k-minor without missing k')
            if not np.array_equal(slots, np.arange(i*self.nk, (i+1)*self.nk)):
                raise ValueError('spin-major slot ordering required')
            self.slots_by_spin[i] = slots
        self.target = _numeric(self.target, (self.nq, self.ns), 'target')
        self.minus = _numeric(self.minus, (self.nq,), 'minus')
        if not np.issubdtype(self.target.dtype, np.integer) or not np.issubdtype(self.minus.dtype, np.integer):
            raise ValueError('integer target/minus maps required')
        minus = np.ravel_multi_index(tuple((-expected % self.mesh).T), self.mesh)
        if not np.array_equal(self.minus, minus):
            raise ValueError('minus must be the exact modulo-mesh inverse')
        for q, vector in enumerate(expected):
            tk = np.ravel_multi_index(tuple(((expected[self.slot_k]+vector) % self.mesh).T), self.mesh)
            wanted = (np.arange(self.ns)//self.nk)*self.nk + tk
            if not np.array_equal(self.target[q], wanted):
                raise ValueError(f'target[{q}] must preserve spin and implement k+q')
        if np.ndim(self.h0) != 3 or np.shape(self.h0)[1] != np.shape(self.h0)[2]:
            raise ValueError('h0 must contain square active-band blocks')
        self.m = np.shape(self.h0)[1]
        self.h0 = _numeric(self.h0, (self.ns,self.m,self.m), 'h0')
        self.rho0 = _numeric(self.rho0, self.h0.shape, 'rho0')
        self.field = _numeric(self.field, (self.ns,3,self.m,self.m), 'field')
        self.field_reverse = _numeric(self.field_reverse, self.field.shape, 'field_reverse')
        self.d, self.f = int(self.mechanical.d), int(self.mechanical.f)
        if int(self.mechanical.n) != self.n or tuple(self.mechanical.shape) != self.mesh:
            raise ValueError('mechanical and electronic lattices must coincide')
        if self.field_polarizations is None:
            if self.f != 3:
                raise ValueError('explicit field_polarizations required unless all Cartesian modes are present')
            self.field_polarizations = np.eye(3)
        self.field_polarizations = _numeric(self.field_polarizations, (self.f,3), 'field_polarizations')
        self.field_modes = np.einsum('fi,siab->sfab', self.field_polarizations, self.field)
        self.field_reverse_modes = np.einsum('fi,siab->sfab', self.field_polarizations, self.field_reverse)
        self.mechanical_sizes = tuple(2*self.d+(2*self.f if q == 0 else 0) for q in range(self.nq))
        for q, size in enumerate(self.mechanical_sizes):
            for name in ('h','j','generators'):
                _numeric(getattr(self.mechanical,name)[q], (size,size), f'mechanical.{name}[{q}]')
        _numeric(self.mechanical.linear, (self.mechanical_sizes[0],), 'mechanical.linear')
        first = self.q_provider(0)
        if not isinstance(first, QCoefficients):
            raise TypeError('q_provider must return QCoefficients')
        self.r0 = np.shape(first.kernel)[0]
        self._validate_q(first)
        self._validated_q = {0}
        partner = np.asarray(first.feature_partner)
        if not np.array_equal(partner[partner],np.arange(len(partner))):
            raise ValueError('Gamma feature_partner must be an involution')
        self._feature_partners = {0:partner}

    def _validate_q(self, value):
        kernel_shape = np.shape(value.kernel)
        if len(kernel_shape) != 2 or kernel_shape[0] != kernel_shape[1] or kernel_shape[0] < 1:
            raise ValueError('nonempty square local-q feature kernel required')
        r = kernel_shape[0]
        for name, shape in (
            ('ionic',(self.ns,self.d,self.m,self.m)),
            ('ionic_reverse',(self.ns,self.d,self.m,self.m)),
            ('density',(self.ns,r,self.m,self.m)),
            ('potential',(self.ns,self.m,self.m,r)),
            ('kernel',(r,r)),
            ('feature_partner',(r,)),
        ):
            _numeric(getattr(value,name), shape, name)
        partner = np.asarray(value.feature_partner)
        if not np.issubdtype(partner.dtype,np.integer) or np.any(partner < 0) or np.any(partner >= r):
            raise ValueError('feature_partner must index the equal-sized minus-q family')
        if len(np.unique(partner)) != r:
            raise ValueError('feature_partner must be a bijection onto the paired feature family')
        if value.feature_g is not None:
            gs = np.asarray(value.feature_g)
            if gs.ndim != 2 or gs.shape[1] != 3 or not np.issubdtype(gs.dtype,np.integer):
                raise ValueError('feature_g must contain integer reciprocal triples')
        if value.feature_spin is not None:
            if np.asarray(value.feature_spin).shape != (r,):
                raise ValueError('feature_spin must label every local-q feature')
        return value

    def q(self, index):
        """Return one provider-owned q record, never copy or retain the q bank."""
        index = int(index)
        if not 0 <= index < self.nq:
            raise IndexError('q index outside complete mesh')
        value = self.q_provider(index)
        if not isinstance(value, QCoefficients):
            raise TypeError('q_provider must return QCoefficients')
        if index not in self._validated_q:
            self._validate_q(value)
            partner = np.asarray(value.feature_partner)
            other = self._feature_partners.get(int(self.minus[index]))
            if other is not None:
                if len(other) != len(partner) or not np.array_equal(other[partner],np.arange(len(partner))):
                    raise ValueError('paired feature_partner maps must be mutual inverses')
            self._feature_partners[index] = partner
            self._validated_q.add(index)
        return value

    def feature_count(self, index):
        return np.shape(self.q(index).kernel)[0]

    @property
    def sqrt_n(self):
        return np.sqrt(self.n)
