"""Finite-rank, full-PW equal-time CAR contraction with explicit state authority.

Raw Psi is never orthonormalized or interpreted as a density. The functional
stores S=Psi^dagger Psi, M_A=Psi^dagger A Psi, and
T_AB=(A^dagger Psi)^dagger(B Psi), retaining the full finite-PW complement.
For a declared gauge-invariant quasifree state rho=Psi D Psi^dagger,
C_AB=Tr(D T_AB)-Tr(D M_A D M_B). D contains density coefficients, not
probabilities in nonorthogonal coordinates. No PW-square object is formed.
No physical ensemble, stationarity, dynamics, spectrum, or ionic state follows.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
import numpy as np


def _matrix(value, name):
    value = np.asarray(value, dtype=np.complex128)
    if value.ndim != 2 or not np.isfinite(value).all():
        raise ValueError(name+' must be a finite matrix')
    return value


def _frozen(value):
    value = np.array(value, dtype=np.complex128, copy=True)
    value.flags.writeable = False
    return value


def _tolerance(value):
    return 128*np.finfo(float).eps*max(1, len(value))*(1+np.max(np.abs(value), initial=0))


def _hermitian(value, name):
    value = _matrix(value, name)
    if value.shape[0] != value.shape[1] or not len(value):
        raise ValueError(name+' must be nonempty square')
    if np.max(np.abs(value-value.conj().T)) > _tolerance(value):
        raise ValueError(name+' is not Hermitian within roundoff allowance')
    return value


def _gram_root(gram):
    gram = _hermitian(gram, 'raw factor Gram')
    values, vectors = np.linalg.eigh(gram)
    if values.min() <= _tolerance(gram):
        raise ValueError('raw factors must be linearly independent at roundoff resolution; no rank truncation')
    return (vectors*np.sqrt(values))@vectors.conj().T, values, vectors


@dataclass(frozen=True)
class ActionFactors:
    """Actual A Psi and A^dagger Psi on identical full finite-PW rows.

The caller binds operator/adjoint and full-row provenance. For a non-Hermitian
operator the adjoint action must be supplied independently, not guessed from
its saved-band overlap. Arrays are consumed without modifying them.
"""
    action: np.ndarray
    adjoint_action: np.ndarray
    operator_identity: str
    adjoint_identity: str
    units: str
    pw_basis_identity: str
    source_authority: str

    def validate(self, shape, basis_identity):
        for name in ('operator_identity', 'adjoint_identity', 'units',
                     'pw_basis_identity', 'source_authority'):
            if not getattr(self, name).strip():
                raise ValueError('explicit '+name+' required')
        if self.pw_basis_identity != basis_identity:
            raise ValueError('operator and factor PW basis identities differ')
        for name in ('action', 'adjoint_action'):
            if _matrix(getattr(self, name), name).shape != shape:
                raise ValueError('action must retain all full-PW rows and source columns')


@dataclass(frozen=True)
class FactorFunctional:
    """State-independent small matrices; no occupations or material numbers.

All saved columns remain present. ordered_products[(A,B)] is directed: it is
not a covariance, is not Hermitized, and does not conjugate B. Serialization
and original-byte receipt checking belong to the producing caller.
"""
    gram: np.ndarray
    overlaps: object
    ordered_products: object
    operator_metadata: object
    factor_identity: str
    pw_basis_identity: str
    full_pw_dimension: int

    @classmethod
    def from_actions(cls, psi, probes, *, factor_identity, pw_basis_identity):
        psi = _matrix(psi, 'raw saved factors')
        if (not psi.shape[1] or psi.shape[1] > psi.shape[0] or not probes
                or not factor_identity.strip() or not pw_basis_identity.strip()):
            raise ValueError('nonempty independent-factor dimensions and explicit identities required')
        for name, probe in probes.items():
            if not isinstance(name, str) or not name.strip():
                raise ValueError('nonempty operator key required')
            probe.validate(psi.shape, pw_basis_identity)
        adjoint = psi.conj().T
        gram = adjoint@psi
        _gram_root(gram)
        overlaps = {name: _frozen(adjoint@probe.action) for name, probe in probes.items()}
        products = {}
        for left, a in probes.items():
            left_adjoint = np.asarray(a.adjoint_action).conj().T
            for right, b in probes.items():
                if (right, left) in products and (
                        a.adjoint_action is a.action and b.adjoint_action is b.action
                        and a.adjoint_identity == a.operator_identity
                        and b.adjoint_identity == b.operator_identity):
                    # For explicitly identical Hermitian action factors, this is
                    # the exact cross-Gram adjoint relation, not Hermitization.
                    products[left, right] = _frozen(products[right, left].conj().T)
                else:
                    products[left, right] = _frozen(left_adjoint@b.action)
        metadata = {name: MappingProxyType({key: getattr(probe, key) for key in
                    ('operator_identity', 'adjoint_identity', 'units', 'pw_basis_identity', 'source_authority')})
                    for name, probe in probes.items()}
        return cls(_frozen(gram), MappingProxyType(overlaps), MappingProxyType(products),
                   MappingProxyType(metadata), factor_identity, pw_basis_identity, psi.shape[0])

    def parameterization(self):
        return {'state': None, 'density_coefficients': None,
                'factor_identity': self.factor_identity, 'pw_basis_identity': self.pw_basis_identity,
                'full_pw_dimension': self.full_pw_dimension, 'factor_count': len(self.gram),
                'rho': 'Psi D Psi^dagger; raw saved Psi unchanged',
                'admissibility': 'D Hermitian; 0 <= S^(1/2) D S^(1/2) <= I; S positive definite',
                'connected_AB': 'Tr(D T_AB)-Tr(D M_A D M_B)',
                'complement': 'full finite-PW orthogonal complement exactly empty in this family',
                'operators': {key: dict(value) for key, value in self.operator_metadata.items()},
                'scope': 'parameterized electronic equal-time functional, not a material preparation'}

    def archive_payload(self):
        """Return numeric arrays and JSON-safe metadata without any preparation.

        Keys use stable insertion-order indices rather than embedding arbitrary
        operator labels. Both directed products are exported, including exact
        cross-Gram adjoint partners. Arrays retain original computed values.
        """
        names = tuple(self.overlaps)
        arrays = {'gram': self.gram}
        layout = {'gram': 'gram', 'overlaps': {}, 'ordered_products': []}
        for i, name in enumerate(names):
            key = 'M_'+str(i).zfill(3)
            arrays[key] = self.overlaps[name]
            layout['overlaps'][name] = key
            for j, right in enumerate(names):
                key = 'T_'+str(i).zfill(3)+'_'+str(j).zfill(3)
                arrays[key] = self.ordered_products[name, right]
                layout['ordered_products'].append({'left': name, 'right': right, 'array': key})
        metadata = self.parameterization()
        metadata['array_layout'] = layout
        return arrays, metadata


def coefficients_from_orthonormal_populations(gram, populations, *, mapping_authority):
    """Explicit symmetric-frame mapping, NOT normalization of the source Psi.

Q=Psi S^-1/2 is orthonormal. A separately authorized diagonal population F
in THIS Q frame gives rho=Q F Q^dagger=Psi D Psi^dagger with
D=S^-1/2 F S^-1/2. Populations are not raw MV values or saved-band labels
unless the caller separately justifies that identification. D is returned;
the input factors are never modified. A different frame specifies a different
state unless F is transported consistently.
"""
    if not mapping_authority.strip():
        raise ValueError('explicit orthonormal-frame population mapping authority required')
    _, values, vectors = _gram_root(gram)
    populations = np.asarray(populations)
    if (np.iscomplexobj(populations) or populations.shape != values.shape
            or not np.isfinite(populations).all() or np.any(populations < 0)
            or np.any(populations > 1)):
        raise ValueError('explicit real probabilities in [0,1] required; no clipping')
    inverse_root = (vectors/np.sqrt(values))@vectors.conj().T
    return (inverse_root*populations)@inverse_root


class FactorPreparation:
    """An explicitly authorized independent quasifree sector, empty complement.

Admissibility uses the nonzero spectrum of rho, that of S^1/2 D S^1/2.
Hermiticity/Pauli allowances are 128 eps rank (1+max abs(matrix)). Raw
coefficients/eigenvalues are retained: no clipping, fitting, or Hermitization.
This family does not encode exact fixed-N states or interacting equilibrium.
"""
    def __init__(self, functional, density_coefficients, state_authority):
        if not state_authority.strip():
            raise ValueError('explicit state authority required')
        d = _hermitian(density_coefficients, 'density coefficients D')
        if d.shape != functional.gram.shape:
            raise ValueError('D and raw Gram dimensions differ')
        root, gram_values, _ = _gram_root(functional.gram)
        represented = root@d@root
        _hermitian(represented, 'represented density')
        spectrum = np.linalg.eigvalsh(represented)
        tolerance = _tolerance(represented)
        if spectrum.min() < -tolerance or spectrum.max() > 1+tolerance:
            raise ValueError('density fails 0 <= rho <= I; coefficients are not probabilities')
        self.functional = functional
        self.density_coefficients = _frozen(d)
        self.state_authority = state_authority
        self.admissibility = {'gram_min_eigenvalue': float(gram_values.min()),
                              'density_eigenvalues': spectrum.tolist(),
                              'pauli_roundoff_allowance': float(tolerance),
                              'coefficient_meaning': 'density coefficients in raw nonorthogonal frame',
                              'complement': 'full finite-PW orthogonal complement exactly empty'}

    def mean(self, name):
        return complex(np.trace(self.density_coefficients@self.functional.overlaps[name]))

    def ordered(self, left, right):
        """Connected <delta A delta B>, no conjugation of B or its coefficients."""
        d, f = self.density_coefficients, self.functional
        return complex(np.trace(d@f.ordered_products[left, right])
                       - np.trace(d@f.overlaps[left]@d@f.overlaps[right]))

    def correlations(self, left, right):
        ab, ba = self.ordered(left, right), self.ordered(right, left)
        return {'mean_A': self.mean(left), 'mean_B': self.mean(right),
                'ordered_AB': ab, 'ordered_BA': ba, 'symmetrized': (ab+ba)/2,
                'commutator': ab-ba, 'state_authority': self.state_authority,
                'admissibility': self.admissibility,
                'units_A': self.functional.operator_metadata[left]['units'],
                'units_B': self.functional.operator_metadata[right]['units']}
