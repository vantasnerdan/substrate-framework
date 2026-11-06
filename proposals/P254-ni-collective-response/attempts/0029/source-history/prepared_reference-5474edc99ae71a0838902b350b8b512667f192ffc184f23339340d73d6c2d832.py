"""Declared joint-second-order effective scalar, not new-state full PBE.

E = sum_s w_s Tr(H0_s rho_s) + dv/2 sum(delta_n * K delta_n)
    + u.(G_active + G_held - F_II - Fsupport)
    + u.C.u/2 + sum(p*p/M)/2,
delta_n = n_active + n_held - nbar.  All pairings are bilinear.
The optional additive energy constant is irrelevant and is fixed to zero.
C is the supplied unmassweighted, fixed joint-second-order coefficient, already
including ionic curvature. Only the ionic LINEAR term enters here; adding a
displaced Ewald scalar would double count curvature. W*delta_rho*u*u is omitted
at this Taylor order, not physically zero. No processed force or support is inferred.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Protocol
import numpy as np

from .physical_factor_preparation import FactorPreparation


class DyadCoupling(Protocol):
    npw: int
    nr: int
    ndof: int
    cell_volume_bohr3: float
    grid_measure_bohr3: float

    def moments(self, left, right, coefficient=None): ...
    def potential_action(self, field, vectors): ...
    def ionic_action(self, direction, vectors): ...


def array(value, shape, name, *, real=False):
    raw = np.asarray(value)
    if real and np.iscomplexobj(raw) and np.any(raw.imag != 0):
        raise ValueError(name + ' must be real, not projected to its real part')
    result = np.asarray(raw.real if real else raw, dtype=float if real else complex)
    if result.shape != shape or not np.isfinite(result).all():
        raise ValueError(name + ' must have finite shape ' + str(shape))
    return result


def frozen(value):
    result = np.array(value, copy=True)
    result.flags.writeable = False
    return result


def authority(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('explicit ' + name + ' required')


@dataclass(frozen=True)
class PreparedSector:
    """One explicit independent spin/k mode; no implicit degeneracy factor.

The initial frame is the exact raw frame admitted by FactorPreparation.
Its complement is initially empty, not excluded from subsequent actions.
"""
    preparation: FactorPreparation
    initial_frame: np.ndarray
    coupling: DyadCoupling
    h0: Callable
    density_component: int
    weight: float
    basis_identity: str
    operator_authority: str

    def __post_init__(self):
        if not isinstance(self.preparation, FactorPreparation):
            raise TypeError('existing FactorPreparation required')
        f = self.preparation.functional
        authority(self.operator_authority, 'H0/g/density common-basis authority')
        if self.basis_identity != f.pw_basis_identity:
            raise ValueError('preparation and action basis identities differ')
        psi = array(self.initial_frame, (self.coupling.npw, len(f.gram)), 'initial raw frame')
        if f.full_pw_dimension != self.coupling.npw:
            raise ValueError('full PW dimension differs from preparation')
        if not np.allclose(psi.conj().T @ psi, f.gram, rtol=1e-12, atol=1e-12):
            raise ValueError('initial frame does not have the admitted raw Gram')
        if not callable(self.h0) or not np.isfinite(self.weight) or self.weight <= 0:
            raise ValueError('H0 action and positive explicit mean weight required')
        if not isinstance(self.density_component, int) or self.density_component < 0:
            raise ValueError('nonnegative explicit density component required')
        object.__setattr__(self, 'initial_frame', frozen(psi))

    @property
    def d(self):
        return self.preparation.density_coefficients


@dataclass(frozen=True)
class PreparedReference:
    """Mandatory ionic baseline and identities of ONE conservative effective model.

K must be complex-linear and real/self-transpose in the physical dv pairing;
H0 and the real-field g/V actions must be Hermitian. These are source-bound
operator contracts, not established by sampling or by repairing coefficients.
Fsupport is an explicitly supplied CONSTANT conservative force (-Fsupport.u).
More general time-dependent drives/constraints are not represented here.
ionic_force is the single unprocessed F_II=-dU_II/dR, in Ry/bohr and the
declared coordinate order; ionic_authority binds its scalar/source/sum contract.
Held electronic density/gradient and physical support remain separate inputs.
"""
    nbar: np.ndarray
    ionic_force: np.ndarray
    held_density: np.ndarray
    held_gradient: np.ndarray
    curvature: np.ndarray
    masses: np.ndarray
    support_force: np.ndarray
    kernel: Callable
    ncells: int
    expansion_authority: str
    ionic_authority: str
    coverage_authority: str
    support_authority: str
    coordinate_identity: str
    state_authority: str

    def __post_init__(self):
        for name in ('expansion_authority', 'ionic_authority', 'coverage_authority',
                     'support_authority', 'coordinate_identity', 'state_authority'):
            authority(getattr(self, name), name)
        nbar = np.asarray(self.nbar)
        mass = np.asarray(self.masses)
        if nbar.ndim != 2 or min(nbar.shape) < 1 or mass.ndim != 1 or len(mass) < 1:
            raise ValueError('nonempty grid/component density and full-coordinate masses required')
        for name in ('nbar', 'held_density'):
            object.__setattr__(self, name, frozen(array(getattr(self, name), nbar.shape, name, real=True)))
        for name in ('ionic_force', 'held_gradient', 'masses', 'support_force'):
            object.__setattr__(self, name, frozen(array(getattr(self, name), mass.shape, name, real=True)))
        c = array(self.curvature, (len(mass), len(mass)), 'unmassweighted C', real=True)
        if not np.allclose(c, c.T, atol=1e-12, rtol=1e-12) or np.any(self.masses <= 0):
            raise ValueError('symmetric unmassweighted C and positive Ry masses required; no repair')
        object.__setattr__(self, 'curvature', frozen(c))
        if (not isinstance(self.ncells, int) or isinstance(self.ncells, bool)
                or self.ncells < 1 or not callable(self.kernel)):
            raise ValueError('positive integer uniform-cell N and kernel action required')

    @property
    def ndof(self):
        return len(self.masses)

    def validate_sectors(self, sectors):
        if not sectors:
            raise ValueError('at least one explicit active sector required')
        measure = sectors[0].coupling.grid_measure_bohr3
        volume = sectors[0].coupling.cell_volume_bohr3
        for s in sectors:
            c = s.coupling
            if (c.nr != self.nbar.shape[0] or c.ndof != self.ndof
                    or s.density_component >= self.nbar.shape[1]
                    or c.grid_measure_bohr3 != measure or c.cell_volume_bohr3 != volume):
                raise ValueError('sectors must share full coordinate/grid and physical measure')
            # X=sqrt(N)*delta rho and canonical collective ions use N*w=1.
            # Multiplicities must be separate independent modes, not IBZ weights.
            if s.weight != 1 / self.ncells:
                raise ValueError('explicit uniform independent sectors require weight=1/N')
        if not np.isfinite(measure) or measure <= 0 or not np.isfinite(volume) or volume <= 0:
            raise ValueError('positive physical grid measure and volume required')
        return measure
