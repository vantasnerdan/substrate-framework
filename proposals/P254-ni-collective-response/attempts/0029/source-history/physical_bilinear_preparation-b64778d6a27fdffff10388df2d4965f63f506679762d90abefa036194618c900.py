"""Selected ordered two-point contractions, never amplitude-noise preparation.

Electrons are independent gauge-invariant quasi-free channels with RAW
probabilities. This does not assert Gaussian higher bilinear moments, a fixed
particle-number ensemble, material temperature, or coupled equilibrium.
Caller-supplied ionic and cross state authority is mandatory for TOTAL output.
Two-point positivity alone is not a representability proof for a joint state.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


def canonical_commutator(nion):
    identity = np.eye(nion)
    zero = np.zeros_like(identity)
    return 1j * np.block([[zero, identity], [-identity, zero]])


@dataclass(frozen=True)
class IonicPreparation:
    """Centered ordered <(u,p)_a (u,p)_b>, hbar=1; explicit state authority."""
    ordered: np.ndarray
    state_authority: str

    def __post_init__(self):
        value = np.asarray(self.ordered, complex)
        if (value.ndim != 2 or value.shape[0] != value.shape[1]
                or value.shape[0] % 2 or not value.shape[0]
                or not np.all(np.isfinite(value)) or not self.state_authority.strip()):
            raise ValueError("finite even square ionic covariance and state authority required")
        if not np.allclose(value, value.conj().T, atol=1e-12, rtol=1e-12):
            raise ValueError("ordered Hermitian-coordinate covariance must be Hermitian")
        if not np.allclose(value-value.T, canonical_commutator(len(value)//2), atol=1e-12):
            raise ValueError("ionic antisymmetric part must be the canonical commutator")
        variances = value.diagonal().real
        if np.any(variances <= 0):
            raise ValueError("canonical ionic coordinates require positive variances")
        # Positive diagonal congruence preserves PSD without mixing u/p units.
        # An absolute eigenvalue check on raw squeezed covariances is unstable.
        scales = np.sqrt(variances)
        normalized = value / scales[:, None] / scales[None, :]
        if np.linalg.eigvalsh(normalized).min() < -1e-12:
            raise ValueError("ionic ordered covariance must be positive semidefinite")
        value = value.copy()
        value.flags.writeable = False
        object.__setattr__(self, "ordered", value)


@dataclass(frozen=True)
class CrossPreparation:
    """One (2*nion,d,d) array per channel: <delta_rho_X * (u,p)_j>.

Each matrix is Hermitian because electronic and ionic initial operators
commute. Cross scaling already refers to X=sqrt(N)*centered bilinears.
The authority must identify an actual joint preparation, not just an arbitrary
PSD matrix. Zero cross is explicit product preparation, never an implicit default.
"""
    matrices: tuple[np.ndarray, ...]
    state_authority: str

    @classmethod
    def product(cls, channels, nion, state_authority):
        return cls(tuple(np.zeros((2*nion, *c.h0.shape), complex) for c in channels),
                   state_authority)


class BilinearPreparation:
    """Contract chosen probes directly without a full electronic covariance."""
    def __init__(self, channels, modes, electronic_authority, ionic=None, cross=None):
        self.channels = tuple(channels)
        modes.validate(self.channels)
        if not electronic_authority.strip():
            raise ValueError("explicit electronic preparation authority required")
        self.modes, self.electronic_authority = modes, electronic_authority
        self.ionic, self.cross = ionic, cross
        if cross is not None:
            if ionic is None or not cross.state_authority.strip():
                raise ValueError("cross preparation requires ionic and joint state authority")
            if len(cross.matrices) != len(self.channels):
                raise ValueError("one cross array per electronic channel required")
            for c, matrix in zip(self.channels, cross.matrices):
                matrix = np.asarray(matrix)
                if (matrix.shape != (len(ionic.ordered), *c.h0.shape)
                        or not np.all(np.isfinite(matrix))
                        or not np.allclose(matrix, matrix.conj().swapaxes(-1, -2), atol=1e-12)):
                    raise ValueError("cross arrays must contain finite Hermitian matrices")

    def ordered_parts(self, left, right):
        """Return <O_left O_right>, with NO conjugation of right coefficients.

For <O_left O_right^dagger>, pass right.adjoint(). Missing pieces remain None;
electronic-only is partial even if both probes are initially electronic, since
pulled probes can have ionic components. TOTAL asserts caller state authority,
not that arbitrary supplied joint moments have been independently certified.
"""
        if len(left.electronic) != len(self.channels) or len(right.electronic) != len(self.channels):
            raise ValueError("one probe matrix per channel required")
        electronic = 0j
        for c, a, b in zip(self.channels, left.electronic, right.electronic):
            if a.shape != c.h0.shape or b.shape != c.h0.shape:
                raise ValueError("probe matrix shape differs from channel")
            rho = c.rho0
            electronic += self.modes.ncells * np.trace(rho @ a @ (np.eye(len(rho))-rho) @ b)
        li, ri = np.concatenate((left.u, left.p)), np.concatenate((right.u, right.p))
        ionic = None if self.ionic is None else li @ self.ionic.ordered @ ri
        cross = None
        if self.cross is not None:
            cross = sum(np.einsum("ij,sji,s->", a, x, ri)
                        + np.einsum("ij,sji,s->", b, x, li)
                        for a, b, x in zip(left.electronic, right.electronic, self.cross.matrices))
        complete = ionic is not None and cross is not None
        return dict(electronic=electronic, ionic=ionic, cross=cross,
                    total=electronic+ionic+cross if complete else None,
                    complete=complete, scope="TOTAL supplied preparation" if complete else "PARTIAL contributions")

    def correlations(self, left, right):
        forward, reverse = self.ordered_parts(left, right), self.ordered_parts(right, left)
        symmetric, commutator = {}, {}
        for key in ("electronic", "ionic", "cross", "total"):
            a, b = forward[key], reverse[key]
            symmetric[key] = None if a is None else (a+b)/2
            commutator[key] = None if a is None else a-b
        return dict(ordered=forward, symmetric=symmetric, commutator=commutator)
