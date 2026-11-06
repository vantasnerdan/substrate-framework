"""Prepared independent-H0 quasi-free bilinear two-point data, not a bath.

For A_plus[a] = sum_ts g[a,t,s] c_t† c_s and A_minus=A_plus†,
C_ab(tau) = <delta A_minus[a](tau) delta A_plus[b](0)>
is sum_ts conj(g[a,t,s]) g[b,t,s] f_s (1-f_t)
exp(-i (epsilon_t-epsilon_s) tau). Energies and inverse-energy times
must use compatible units (e.g. Ry and t*Ry/hbar). Occupations are raw
Bernoulli probabilities, never k-weighted occupations. Independent sectors
are added by the caller; explicit spin channels require no extra factor.

Source conventions: finite_fermion_state.py (grandcanonical CAR state) and
finite_fermion_response.py:noise_lines/symmetrized_lines (connected zero-gap
terms and opposite-time transpose). The new real-space/full-Fock fixture in
attempts/0023/run_quantum_noise_oracle.py corroborates this construction;
it is not a theorem for general interacting equilibrium or higher moments.
"""
from __future__ import annotations

import numpy as np


_TRANSITION_BLOCK = 4096


def _real_vector(value, name):
    array = np.asarray(value)
    if array.ndim != 1 or np.iscomplexobj(array):
        raise ValueError(f"{name} must be a real one-dimensional array")
    array = np.asarray(array, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def ordered_transition_correlation(
    source_energies, target_energies, source_occupations, target_occupations,
    operators, times, weight=1.0,
):
    """Return C(tau), shape (ntimes,nobs,nobs), retaining every transition.

    ``operators`` has shape (nobs,ntarget,nsource). ``weight=1/Ncells``
    implements operators normalized by 1/sqrt(Ncells) when their supplied
    matrix elements omit that normalization. Do not apply it twice. For
    already normalized operators use weight=1. No occupation cutoff, gap
    division, frequency binning, temperature inference, or fitted amplitude
    is applied. Overlapping source/target sectors (including Gamma) retain
    population fluctuations and all non-diagonal degenerate coherences.
    Contractions are bounded transition blocks, not nobs²-by-ntrans tensors.
    """
    es = _real_vector(source_energies, "source_energies")
    et = _real_vector(target_energies, "target_energies")
    fs = _real_vector(source_occupations, "source_occupations")
    ft = _real_vector(target_occupations, "target_occupations")
    tau = _real_vector(times, "times")
    if fs.shape != es.shape or ft.shape != et.shape:
        raise ValueError("occupation and energy shapes must match in each sector")
    if np.any((fs < 0) | (fs > 1)) or np.any((ft < 0) | (ft > 1)):
        raise ValueError("occupations must be raw probabilities in [0,1]")
    scale = np.asarray(weight)
    if (scale.ndim != 0 or np.iscomplexobj(scale)
            or not np.isfinite(scale) or scale < 0):
        raise ValueError("weight must be a finite nonnegative real scalar")
    op = np.asarray(operators, dtype=complex)
    if op.ndim != 3 or op.shape[1:] != (et.size, es.size):
        raise ValueError("operators must have shape (nobs,ntarget,nsource)")
    if not np.all(np.isfinite(op)):
        raise ValueError("operators must be finite")
    nobs = op.shape[0]
    ntrans = et.size * es.size
    result = np.zeros((tau.size, nobs, nobs), dtype=complex)
    flat = op.reshape(nobs, ntrans)
    for start in range(0, ntrans, _TRANSITION_BLOCK):
        stop = min(start + _TRANSITION_BLOCK, ntrans)
        indices = np.arange(start, stop)
        target, source = np.divmod(indices, es.size)
        gaps = et[target] - es[source]
        probabilities = float(scale) * fs[source] * (1.0 - ft[target])
        block = flat[:, start:stop]
        left = block.conj() * probabilities
        for index, time in enumerate(tau):
            result[index] += (left * np.exp(-1j * gaps * time)) @ block.T
    return result


def paired_correlations(ordered_plus, ordered_minus_at_negative_times):
    """Return connected anticommutator/2 and commutator in that order.

    The second input is C_minus(-tau), with minus operators the channelwise
    adjoints of plus operators and the same channel ordering. Stationarity
    gives <delta A_plus[b](0) delta A_minus[a](tau)> = C_minus[b,a](-tau).
    Therefore the reverse ordered product uses a TRANSPOSE, not an adjoint.
    Inputs must already include the same sector sums/normalizations. Output
    ``commutator`` has no extra i, hbar, retarded step, or force-sign factor.
    """
    plus = np.asarray(ordered_plus, dtype=complex)
    minus = np.asarray(ordered_minus_at_negative_times, dtype=complex)
    if (plus.ndim != 3 or plus.shape[1] != plus.shape[2]
            or minus.shape != plus.shape):
        raise ValueError("paired inputs must share shape (ntimes,nobs,nobs)")
    if not np.all(np.isfinite(plus)) or not np.all(np.isfinite(minus)):
        raise ValueError("paired inputs must be finite")
    reverse = minus.swapaxes(1, 2)
    return {"symmetrized": (plus + reverse) * 0.5,
            "commutator": plus - reverse}
