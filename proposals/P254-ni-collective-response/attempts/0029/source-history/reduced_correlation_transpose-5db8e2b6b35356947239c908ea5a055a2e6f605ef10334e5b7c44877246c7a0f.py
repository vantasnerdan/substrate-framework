"""Bilinear (never Hermitian) duals for canonical active-density blocks.

Every dual uses sum(z * sigma), not vdot(z, sigma).  The last axis batches
independent probes.  Matrix indices remain (target band, source band).
"""
from __future__ import annotations

import numpy as np


def commutator_tangent(h_target, h_source, rho_target, rho_source, sigma, delta_h):
    """-i(h_t sigma - sigma h_s + delta_h rho_s - rho_t delta_h)."""
    return -1j * (
        np.einsum('sab,sbcp->sacp', h_target, sigma, optimize=True)
        - np.einsum('sabp,sbc->sacp', sigma, h_source, optimize=True)
        + np.einsum('sabp,sbc->sacp', delta_h, rho_source, optimize=True)
        - np.einsum('sab,sbcp->sacp', rho_target, delta_h, optimize=True)
    )


def commutator_transpose(h_target, h_source, rho_target, rho_source, dual):
    """Return duals of sigma and delta_h for the *nonconjugating* pairing.

    For z_ab dot sigma_ab the two matrix duals are
      -i(h_t.T z - z h_s.T), -i(z rho_s.T - rho_t.T z).
    A dagger here would propagate a different, incorrect observable.
    """
    sigma_dual = -1j * (
        np.einsum('sba,sbcp->sacp', h_target, dual, optimize=True)
        - np.einsum('sabp,scb->sacp', dual, h_source, optimize=True)
    )
    potential_dual = -1j * (
        np.einsum('sabp,scb->sacp', dual, rho_source, optimize=True)
        - np.einsum('sba,sbcp->sacp', rho_target, dual, optimize=True)
    )
    return sigma_dual, potential_dual


def wick_edge_pair(left, right_reversed, rho_source, rho_target):
    """Connected ordered quasifree pairing of opposite directed edges.

    left[s,a,b,p] multiplies sigma_(target,a;source,b).  Its paired
    right_reversed[s,c,d,r] is the -q variable whose source is that target.
    Wick gives rho_source[c,b] (I-rho_target)[a,d].  The canonical global
    density coordinates carry neither a spin weight nor a 1/N here.
    """
    hole = np.eye(rho_target.shape[-1])[None, :, :] - rho_target
    return np.einsum('sabp,scdr,scb,sad->pr', left, right_reversed,
                     rho_source, hole, optimize=True)
