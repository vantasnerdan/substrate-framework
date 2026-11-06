"""QE real harmonics, centered NC projectors, and atom placement phases.

Input momenta are Cartesian bohr^-1 and positions Cartesian bohr. The PW
rows remain in caller order. Channels follow init_us_1 radial-then-angular
order. Across atoms QE orders species then atom; no dense PW matrix is made.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .qe_nc_radial import NCRadialTables


def _momenta(value: NDArray[np.float64]) -> NDArray[np.float64]:
    p = np.asarray(value, dtype=np.float64)
    if p.ndim != 2 or p.shape[1] != 3 or not np.isfinite(p).all():
        raise ValueError("finite momenta (npw,3) in bohr^-1 required")
    return p


def qe_real_harmonics(momenta_bohr_inv: NDArray[np.float64], lmax: int
                      ) -> NDArray[np.float64]:
    """Return QE ylmr2 order: m=0, cos(phi), sin(phi), cos(2phi), ... .

    Includes the Condon-Shortley sign (p channels are z,-x,-y), not the
    common x,y,z real-harmonic convention. At zero use QE's theta=pi/2,
    phi=pi/2; all l>0 radial projectors vanish there independently.
    """
    p = _momenta(momenta_bohr_inv)
    if lmax not in (0, 1, 2):
        raise ValueError("selected Ni/H s/p/d harmonics only")
    radius = np.linalg.norm(p, axis=1)
    cost = np.divide(p[:, 2], radius, out=np.zeros_like(radius), where=radius >= 1e-9)
    sent = np.sqrt(np.maximum(0.0, 1-cost*cost))
    phi = np.arctan2(p[:, 1], p[:, 0])
    axis = np.abs(p[:, 0]) <= 1e-9
    phi[axis] = np.where(p[axis, 1] >= 0, np.pi/2, -np.pi/2)
    y = np.empty((len(p), (lmax+1)**2))
    y[:, 0] = 1/np.sqrt(4*np.pi)
    if lmax >= 1:
        c = np.sqrt(3/(4*np.pi))
        y[:, 1] = c*cost
        y[:, 2] = -c*sent*np.cos(phi)
        y[:, 3] = -c*sent*np.sin(phi)
    if lmax >= 2:
        c = np.sqrt(5/(4*np.pi))
        y[:, 4] = c*(3*cost*cost-1)/2
        y[:, 5] = -c*np.sqrt(3)*cost*sent*np.cos(phi)
        y[:, 6] = -c*np.sqrt(3)*cost*sent*np.sin(phi)
        y[:, 7] = c*np.sqrt(3)/2*sent*sent*np.cos(2*phi)
        y[:, 8] = c*np.sqrt(3)/2*sent*sent*np.sin(2*phi)
    return y


def centered_projectors(tables: NCRadialTables,
                         momenta_bohr_inv: NDArray[np.float64]) -> NDArray[np.complex128]:
    """Return B(npw,nh), including (-i)^l, without an atomic phase.

    Reuse this result for every atom of this species with the same PW
    momenta, including both spins when their physical row ordering agrees.
    No radial integration or table construction occurs here.
    """
    p = _momenta(momenta_bohr_inv)
    upf = tables.upf
    radial = tables.beta(np.linalg.norm(p, axis=1))
    harmonics = qe_real_harmonics(p, int(upf.angular_l.max()))
    return (radial[:, upf.channel_radial] * harmonics[:, upf.channel_lm]
            * ((-1j)**upf.channel_l)[None, :])


def place_projectors(centered: NDArray[np.complex128],
                      momenta_bohr_inv: NDArray[np.float64],
                      position_bohr: NDArray[np.float64]) -> NDArray[np.complex128]:
    """Multiply precomputed centered B by exp[-i (k+G).tau]."""
    p = _momenta(momenta_bohr_inv)
    tau = np.asarray(position_bohr, dtype=np.float64)
    centered = np.asarray(centered, dtype=np.complex128)
    if tau.shape != (3,) or not np.isfinite(tau).all():
        raise ValueError("finite Cartesian atom position (3,) in bohr required")
    if centered.ndim != 2 or centered.shape[0] != len(p):
        raise ValueError("centered projectors must have the same physical PW rows")
    return centered * np.exp(-1j*(p @ tau))[:, None]


def atom_projectors(tables: NCRadialTables,
                     momenta_bohr_inv: NDArray[np.float64],
                     position_bohr: NDArray[np.float64], *,
                     centered: NDArray[np.complex128] | None = None
                     ) -> tuple[NDArray[np.complex128], NDArray[np.float64]]:
    """Return atom B(npw,nh), real active D(nh,nh) in Native conventions.

    B D B-dagger is a Ry PW operator. Pass centered= to avoid repeating
    interpolation/harmonics for another atom or spin at identical momenta.
    D is the shared immutable expanded UPF matrix, not a per-atom copy.
    """
    if centered is None:
        centered = centered_projectors(tables, momenta_bohr_inv)
    if centered.shape[1:] != (tables.upf.nh,):
        raise ValueError("centered projector channel count disagrees with species")
    return place_projectors(centered, momenta_bohr_inv, position_bohr), tables.upf.expanded_d
