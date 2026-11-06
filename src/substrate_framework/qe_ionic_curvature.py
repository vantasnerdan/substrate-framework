"""Finite Cartesian QE 7.5 PH d2ionq ion-ion Hessian, in Ry/bohr**2.

This is the ordinary 3D, full signed-G, all-atom branch, not screened DFPT.
Inputs retain the producer's saved G support, alpha and rgen radius. QE's
final dy3=-dy3 is already included in every returned term. No ASR, Hermitian
projection, atomic gauge rotation or mass weighting is applied. The Coulomb
head makes the physical Gamma limit directional; QE excludes |G+q|²<=1e-8.
The threshold is in bohr**-2, unlike rgen's 1e-10 exclusion in alat units.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
from scipy.special import erfc


@dataclass(frozen=True)
class IonicCurvatureTerms:
    reciprocal_pair: np.ndarray
    reciprocal_self: np.ndarray
    real_pair: np.ndarray
    real_self: np.ndarray

    @property
    def total(self):
        return self.reciprocal_pair + self.reciprocal_self + self.real_pair + self.real_self


def qe_alpha(*, total_valence, density_cutoff_ry):
    """Literal d2ionq alpha decrement/energy-tail heuristic, not Hessian error."""
    if not np.isfinite(total_valence) or not np.isfinite(density_cutoff_ry) or density_cutoff_ry <= 0:
        raise ValueError("finite charge and positive density cutoff required")
    alpha = 2.9
    while alpha > 0.1:
        alpha -= 0.1
        upper = 2 * total_valence**2 * np.sqrt(alpha / np.pi) * erfc(np.sqrt(density_cutoff_ry / (4 * alpha)))
        if upper <= 1e-9:
            return alpha
    raise ValueError("QE alpha selection exhausted positive candidates")


class QEIonicCurvature:
    """Prepared exact finite-sum convention; cell/reciprocal vectors are columns.

    positions_bohr: [3,nat]; charges: NC valences, not nuclear charges;
    g_cart_bohr_inv: [3,ngm], all original signed density vectors including zero;
    q passed to terms/matrix is physical Cartesian bohr^-1, never fractional.
    alpha_bohr_inv2 must be chosen from the original PH source. real_cutoff_bohr
    defaults to its 5/sqrt(alpha); explicit overrides are refinement experiments.
    max_real_vectors reproduces rgen's per-pair allocation limit (100 in PH).
    """

    def __init__(self, *, cell_bohr, positions_bohr, valence_charges,
                 g_cart_bohr_inv, alat_bohr, alpha_bohr_inv2,
                 real_cutoff_bohr=None, max_real_vectors=100, g_block=8192):
        self.cell = np.asarray(cell_bohr, dtype=float)
        self.positions = np.asarray(positions_bohr, dtype=float)
        self.charges = np.asarray(valence_charges, dtype=float)
        self.g = np.asarray(g_cart_bohr_inv, dtype=float)
        self.alat = float(alat_bohr)
        self.alpha = float(alpha_bohr_inv2)
        self.real_cutoff = 5 / np.sqrt(self.alpha) if real_cutoff_bohr is None else float(real_cutoff_bohr)
        self.g_block = int(g_block)
        if (self.cell.shape != (3, 3) or self.charges.ndim != 1 or self.charges.size == 0
                or self.positions.shape != (3, self.charges.size)
                or self.g.ndim != 2 or self.g.shape[0] != 3 or self.g.shape[1] == 0
                or not all(np.isfinite(x).all() for x in (self.cell, self.positions, self.charges, self.g))
                or not np.isfinite([self.alat, self.alpha, self.real_cutoff]).all()
                or min(self.alat, self.alpha, self.real_cutoff, self.g_block) <= 0):
            raise ValueError("invalid finite ionic geometry/sum inputs")
        self.volume = abs(float(np.linalg.det(self.cell)))
        if self.volume == 0:
            raise ValueError("nonsingular cell required")
        self.nat = self.charges.size
        self._reciprocal_zero = self._reciprocal_pair(np.zeros(3))
        self._reciprocal_self = np.zeros_like(self._reciprocal_zero)
        # d2ionq fnat uses cos(G.dtau) even if finite support has a residual.
        for a in range(self.nat):
            self._reciprocal_self[a, :, a, :] = -self._reciprocal_zero[a].real.sum(axis=1)
        self._bonds = []
        self._real_self = np.zeros_like(self._reciprocal_zero)
        inverse = np.linalg.inv(self.cell)
        bounds = (np.linalg.norm(inverse, axis=1) * self.real_cutoff).astype(int) + 2
        integers = np.asarray(list(product(*(range(-n, n + 1) for n in bounds))), dtype=int)
        translations = integers @ self.cell.T
        for a in range(self.nat):
            for b in range(self.nat):
                delta = self.positions[:, a] - self.positions[:, b]
                fractional = inverse @ delta
                # Fortran ANINT is ties-away-from-zero, not numpy round-to-even.
                centered = fractional - np.copysign(np.floor(np.abs(fractional) + 0.5), fractional)
                r = translations - self.cell @ centered
                radius2 = np.einsum('ni,ni->n', r, r)
                keep = (radius2 <= self.real_cutoff**2) & (radius2 / self.alat**2 > 1e-10)
                r, radius2 = r[keep], radius2[keep]
                if len(r) > max_real_vectors:
                    raise ValueError("QE rgen per-pair maximum exceeded; declare refinement limit explicitly")
                order = np.argsort(radius2)
                r, radius2 = r[order], radius2[order]
                radius = np.sqrt(radius2)
                ar = np.sqrt(self.alpha) * radius
                gaussian = 2 / np.sqrt(np.pi) * ar * np.exp(-ar**2)
                df = (-erfc(ar) - gaussian) / radius**3
                d2f = (3 * erfc(ar) + gaussian * (3 + 2 * ar**2)) / radius**5
                hessian = 2 * self.charges[a] * self.charges[b] * (
                    d2f[:, None, None] * r[:, :, None] * r[:, None, :]
                    + df[:, None, None] * np.eye(3))
                # Original, not centered delta is required in qrg=q.(r+dtau).
                self._bonds.append((a, b, r + delta, hessian))
                self._real_self[a, :, a, :] += hessian.sum(axis=0)

    def _reciprocal_pair(self, q):
        result = np.zeros((self.nat, 3, self.nat, 3), dtype=complex)
        for start in range(0, self.g.shape[1], self.g_block):
            k = (self.g[:, start:start + self.g_block] + q[:, None]).T
            k2 = np.einsum('ni,ni->n', k, k)
            weight = np.zeros_like(k2)
            keep = k2 > 1e-8
            weight[keep] = 8 * np.pi / self.volume * np.exp(-k2[keep] / (4 * self.alpha)) / k2[keep]
            phase = self.charges[:, None] * np.exp(1j * (self.positions.T @ k.T))
            for i in range(3):
                for j in range(i, 3):
                    block = (phase * (weight * k[:, i] * k[:, j])) @ phase.conj().T
                    result[:, i, :, j] += block
                    if i != j:
                        result[:, j, :, i] += block
        return result

    def terms(self, q_cart_bohr_inv):
        q = np.asarray(q_cart_bohr_inv, dtype=float)
        if q.shape != (3,) or not np.isfinite(q).all():
            raise ValueError("finite Cartesian q in bohr^-1 required")
        reciprocal = self._reciprocal_zero.copy() if not np.any(q) else self._reciprocal_pair(q)
        real = np.zeros_like(reciprocal)
        for a, b, translations, hessian in self._bonds:
            real[a, :, b, :] -= np.einsum('n,nij->ij', np.exp(1j * (translations @ q)), hessian)
        return IonicCurvatureTerms(reciprocal, self._reciprocal_self.copy(), real, self._real_self.copy())

    def matrix(self, q_cart_bohr_inv):
        """Return atom-major [3*nat,3*nat] unmassweighted Cartesian curvature."""
        return self.terms(q_cart_bohr_inv).total.reshape(3 * self.nat, 3 * self.nat)
