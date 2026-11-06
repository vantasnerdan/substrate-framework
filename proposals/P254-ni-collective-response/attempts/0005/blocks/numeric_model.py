"""Evaluate smooth finite curves directly; no native-vertex reconstruction is used."""
import numpy as np
from scipy.linalg import eigh, expm

T0 = np.array([[1, .3+.2j], [.1-.4j, 1.4]], dtype=complex)
M = np.array([[.15+.2j, .3+.1j], [-.2+.25j, -.1-.15j]])
H0 = np.array([[-1.3, .4+.35j], [.4-.35j, 1.7]])
V = np.array([[.6, .2-.4j], [.2+.4j, -.5]])
Q = np.array([[.2, -.3+.1j], [-.3-.1j, -.4]])


def curves(r):
    T = expm(r*M)@T0
    h = H0+np.sin(r)*V+(np.cos(.7*r)-1)*Q
    return T, h, T.conj().T@h@T, T.conj().T@T


def solve(r):
    T, h, H, S = curves(r)
    energies, coefficients = eigh(H, S, driver="gvd")
    return energies, coefficients, T@coefficients


def align(coefficients, reference, S0):
    overlaps = np.diag(reference.conj().T@S0@coefficients)
    assert np.min(np.abs(overlaps)) > .5, "eigenbranch overlap lost"
    phase = overlaps.conj()/np.abs(overlaps)
    return coefficients*phase


def analytic_jets():
    T, h, H, S = curves(0.)
    Tp = M@T
    Hp = Tp.conj().T@h@T+T.conj().T@V@T+T.conj().T@h@Tp
    Sp = Tp.conj().T@T+T.conj().T@Tp
    return Tp, Hp, Sp


def inputs():
    return {k: {"real": a.real.tolist(), "imag": a.imag.tolist()}
            for k, a in {"T0": T0, "M": M, "H0": H0, "V": V, "Q": Q}.items()}
