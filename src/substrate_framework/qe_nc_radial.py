"""QE NC radial quadrature and reusable dq=.01 reciprocal-space tables.

Atomic Ry units: q in bohr^-1, omega in bohr^3. Local output is one Fourier
coefficient in Ry, not a radial potential or real-space FFT field. Sources:
QE 7.5 simpsn, sph_bes, beta_mod, vloc_mod; ordinary 3D Coulomb only.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from numpy.typing import NDArray
from scipy.special import erf

from .qe_nc_upf import NCUPF

DQ = 0.01


def qe_simpson_weights(rab: NDArray[np.float64]) -> NDArray[np.float64]:
    """Weights of QE simpson, including its even-mesh endpoint correction."""
    rab = np.asarray(rab, dtype=np.float64)
    if rab.ndim != 1 or rab.size < 3:
        raise ValueError("Simpson needs at least three radial points")
    n = rab.size
    weights = np.ones(n)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    if n % 2 == 0:
        weights[-3] -= 0.25
        weights[-2] += 1.0
        weights[-1] = 1.25
    return weights * rab / 3.0


def qe_spherical_bessel(l: int, x: NDArray[np.float64]) -> NDArray[np.float64]:
    """QE sph_bes s/p/d branches and its small-x polynomial, not sinc beta."""
    x = np.asarray(x, dtype=np.float64)
    if l not in (0, 1, 2) or np.any(x < 0) or not np.isfinite(x).all():
        raise ValueError("finite nonnegative arguments and l=0,1,2 required")
    out = np.empty_like(x)
    small = x <= 0.05
    z = x[small]
    z2 = z*z
    series = 1.0 - z2/(8*(2*l+9))
    series = 1.0 - z2/(6*(2*l+7))*series
    series = 1.0 - z2/(4*(2*l+5))*series
    series = 1.0 - z2/(2*(2*l+3))*series
    out[small] = z**l / math.prod(range(1, 2*l+2, 2)) * series
    z = x[~small]
    if l == 0:
        out[~small] = np.sin(z)/z
    elif l == 1:
        out[~small] = (np.sin(z)/z - np.cos(z))/z
    else:
        out[~small] = ((3/z-z)*np.sin(z) - 3*np.cos(z))/(z*z)
    return out


def _cubic(table: NDArray[np.float64], q: NDArray[np.float64]) -> NDArray[np.float64]:
    q = np.asarray(q, dtype=np.float64)
    if q.ndim != 1 or not np.isfinite(q).all() or np.any(q < 0):
        raise ValueError("q must be a finite nonnegative one-dimensional modulus array")
    scaled = q / DQ
    index = scaled.astype(np.int64)
    if np.any(index + 3 >= table.shape[0]):
        raise ValueError("requested momentum exceeds precomputed radial table")
    p = scaled - index
    if table.ndim == 2:
        p = p[:, None]
    u, v, w = 1-p, 2-p, 3-p
    return (table[index]*u*v*w/6 + table[index+1]*p*v*w/2
            - table[index+2]*p*u*w/2 + table[index+3]*p*u*v/6)


@dataclass(frozen=True)
class NCRadialTables:
    """Build once per species/volume/range; reuse over atoms, spins and k.

    beta_values includes 4*pi/sqrt(omega). local_short_values and local_g0
    include 4*pi/omega. No atom position, angular harmonic, or spin is stored.
    """
    upf: NCUPF
    omega_bohr3: float
    qmax_bohr_inv: float
    beta_values: NDArray[np.float64]  # (nqx, nbeta)
    local_short_values: NDArray[np.float64]  # q=(0,1,...)*DQ
    local_g0: float

    def beta(self, q_bohr_inv: NDArray[np.float64]) -> NDArray[np.float64]:
        """Return (nq,nbeta) radial transforms; angular expansion is separate."""
        return _cubic(self.beta_values, q_bohr_inv)

    def vloc(self, q_bohr_inv: NDArray[np.float64]) -> NDArray[np.float64]:
        """Return centered atomic Fourier coefficients (nq,) in Ry.

        G=0 uses the finite alpha-Z convention, not the q->0 short-range
        table value. Nonzero lattice G uses the analytic erf Coulomb tail.
        Inputs are physical moduli; only actual zero denotes the G=0 slot.
        """
        q = np.asarray(q_bohr_inv, dtype=np.float64)
        values = _cubic(self.local_short_values, q)
        nonzero = q != 0
        q2 = q[nonzero]**2
        values[nonzero] -= (8*np.pi*self.upf.z_valence/self.omega_bohr3
                            * np.exp(-q2/4)/q2)
        values[~nonzero] = self.local_g0
        return values


def _beta_table(upf: NCUPF, q: NDArray[np.float64], omega: float) -> NDArray[np.float64]:
    n = upf.beta_mesh
    r = upf.r[:n]
    weighted = upf.beta_r[:n] * (r*qe_simpson_weights(upf.rab[:n]))[:, None]
    table = np.empty((q.size, upf.nbeta))
    for start in range(0, q.size, 256):
        stop = min(start + 256, q.size)
        qr = q[start:stop, None] * r
        for l in np.unique(upf.angular_l):
            columns = np.flatnonzero(upf.angular_l == l)
            table[start:stop, columns] = qe_spherical_bessel(int(l), qr) @ weighted[:, columns]
    return table * (4*np.pi/np.sqrt(omega))


def _local_table(upf: NCUPF, q: NDArray[np.float64], omega: float
                 ) -> tuple[NDArray[np.float64], float]:
    n = upf.local_mesh
    r, v = upf.r[:n], upf.local_ry[:n]
    weights = qe_simpson_weights(upf.rab[:n])
    short = r * (r*v + 2*upf.z_valence*erf(r)) * weights
    values = np.empty(q.size)
    for start in range(0, q.size, 256):
        stop = min(start + 256, q.size)
        values[start:stop] = np.sinc(q[start:stop, None]*r/np.pi) @ short
    pref = 4*np.pi/omega
    g0 = pref * np.dot(r*(r*v + 2*upf.z_valence), weights)
    return pref*values, float(g0)


def build_radial_tables(upf: NCUPF, omega_bohr3: float,
                        qmax_bohr_inv: float) -> NCRadialTables:
    """Construct both native tables once, with four-point interpolation guard.

    qmax must cover both beta |k+G| and local |G| uses. This function does not
    run implicitly inside any atom projector or coefficient evaluation.
    """
    if not np.isfinite(omega_bohr3) or omega_bohr3 <= 0:
        raise ValueError("positive finite volume in bohr^3 required")
    if not np.isfinite(qmax_bohr_inv) or qmax_bohr_inv < 0:
        raise ValueError("nonnegative finite maximum momentum required")
    q = np.arange(int(qmax_bohr_inv/DQ + 4), dtype=np.float64) * DQ
    beta = _beta_table(upf, q, omega_bohr3)
    local, g0 = _local_table(upf, q, omega_bohr3)
    beta.setflags(write=False)
    local.setflags(write=False)
    return NCRadialTables(upf, float(omega_bohr3), float(qmax_bohr_inv), beta, local, g0)
