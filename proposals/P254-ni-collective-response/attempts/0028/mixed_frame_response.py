"""Fixed-D coherent q/-q frame derivatives and borrowed low-rank densities.

All frames retain every PW row and their ORIGINAL source column rank.  The
mixed coefficient is d²/(dλ dη), not half of it.  On the physical slice η=λ*,
the bra λ derivative is the dagger of the ket η derivative (and conversely).
These are coherent mean variations, NOT population-inclusive quantum noise.

The caller supplies genuine Hamiltonian actions at the current moving stage,
including the nonstationary reference.  This module neither evolves the mean
for the caller nor constructs/replaces any Hamiltonian or self-consistent
mixed potential.  Arrays in density descriptors are borrowed: keep them alive
and unchanged while consuming descriptors.  No PW-squared density, inverse,
rank compression, frame copy, or Hermiticity repair is performed here.
"""
from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

Action = Callable[[int, np.ndarray], np.ndarray]
Frames = Sequence[np.ndarray]


def _permutations(count, target_plus, target_minus):
    targets = tuple(tuple(p) for p in (target_plus, target_minus))
    for p in targets:
        if (len(p) != count or any(not isinstance(j, (int, np.integer))
                                  or isinstance(j, (bool, np.bool_)) for j in p)
                or sorted(p) != list(range(count))):
            raise ValueError("targets must be explicit endpoint permutations")
    p, m = targets
    if any(m[p[i]] != i or p[m[i]] != i for i in range(count)):
        raise ValueError("q and -q target permutations must be inverses")
    return p, m


def _layout(psi, plus, minus, mixed, target_plus, target_minus):
    count = len(psi)
    if not count or any(len(x) != count for x in (plus, minus, mixed)):
        raise ValueError("complete nonempty mean/plus/minus/mixed tuples required")
    p, m = _permutations(count, target_plus, target_minus)
    for i, frame in enumerate(psi):
        if not isinstance(frame, np.ndarray) or frame.ndim != 2 or min(frame.shape) < 1:
            raise ValueError("mean frames must be full-row two-dimensional arrays")
    for i, frame in enumerate(psi):
        rank = frame.shape[1]
        for name, value, shape in (
                ("plus", plus[i], (psi[p[i]].shape[0], rank)),
                ("minus", minus[i], (psi[m[i]].shape[0], rank)),
                ("mixed", mixed[i], frame.shape)):
            if not isinstance(value, np.ndarray) or value.shape != shape:
                raise ValueError(f"{name}[{i}] must have original-source shape {shape}")
    return p, m


def _action(action, index, vectors, expected):
    result = action(index, vectors)
    if not isinstance(result, np.ndarray) or result.shape != expected:
        raise ValueError(f"Hamiltonian action at endpoint {index} must return {expected}")
    return result


def frame_rhs_item(index, psi, plus, minus, mixed, *, plus_endpoint,
                   minus_endpoint, h0: Action, hplus: Action, hminus: Action,
                   hmixed: Action):
    """Evaluate ONE source's arrays; never open another derivative payload.

    The owner validates the complete geometric permutation/checkpoint before
    dispatch. Already committed other source edges may have been retired;
    this local API deliberately needs no array or stand-in for them.
    Callbacks still select genuine full-row endpoint actions explicitly.
    """
    for value in (index, plus_endpoint, minus_endpoint):
        if (not isinstance(value, (int, np.integer))
                or isinstance(value, (bool, np.bool_)) or value < 0):
            raise ValueError("explicit nonnegative source/target indices required")
    if not isinstance(psi, np.ndarray) or psi.ndim != 2 or min(psi.shape) < 1:
        raise ValueError("the original full-row source frame is required")
    for value in (plus, minus, mixed):
        if (not isinstance(value, np.ndarray) or value.ndim != 2
                or min(value.shape) < 1 or value.shape[1] != psi.shape[1]):
            raise ValueError("each local derivative must retain the source columns")
    if mixed.shape != psi.shape:
        raise ValueError("mixed Gamma rows must match the source frame")
    if not all(callable(h) for h in (h0, hplus, hminus, hmixed)):
        raise TypeError("all four physical Hamiltonian callbacks are required")
    return (
        -1j * (_action(h0, plus_endpoint, plus, plus.shape)
               + _action(hplus, index, psi, plus.shape)),
        -1j * (_action(h0, minus_endpoint, minus, minus.shape)
               + _action(hminus, index, psi, minus.shape)),
        -1j * (_action(h0, index, mixed, mixed.shape)
               + _action(hplus, minus_endpoint, minus, mixed.shape)
               + _action(hminus, plus_endpoint, plus, mixed.shape)
               + _action(hmixed, index, psi, mixed.shape)))


def frame_rhs(psi, plus, minus, mixed, *, target_plus, target_minus,
              h0: Action, hplus: Action, hminus: Action, hmixed: Action):
    """Return (plus_rhs, minus_rhs, mixed_rhs), each an endpoint tuple.

    psi[i],mixed[i] have shape (Npw_i,r_i); plus[i] has
    (Npw_target_plus[i],r_i), minus[i] has (Npw_target_minus[i],r_i).
    h0(i,V),hmixed(i,V) preserve endpoint i; hplus/hminus map i to
    the specified target.  Every callback must accept arbitrary column counts,
    not just the endpoint's native rank.  No callback may mutate its input.
    h0 is the actual reference action, not a stationary energy eigenvalue.
    The reference RHS itself is -1j*h0(i,psi[i]) and remains caller-owned.
    """
    p, m = _layout(psi, plus, minus, mixed, target_plus, target_minus)
    if not all(callable(h) for h in (h0, hplus, hminus, hmixed)):
        raise TypeError("all four physical Hamiltonian callbacks are required")
    rows = tuple(frame_rhs_item(i, psi[i], plus[i], minus[i], mixed[i],
                               plus_endpoint=p[i], minus_endpoint=m[i],
                               h0=h0, hplus=hplus, hminus=hminus, hmixed=hmixed)
                 for i in range(len(psi)))
    return tuple(tuple(row[j] for row in rows) for j in range(3))


@dataclass(frozen=True, slots=True)
class DensityDyad:
    """Borrowed contribution weight * left @ D @ right.conj().T.

    source_index identifies the original occupied-column source of ALL three
    factors and its quadrature weight, even after the columns have moved.
    No multiplication is performed on construction.
    """
    left: np.ndarray
    D: np.ndarray
    right: np.ndarray
    source_index: int
    weight: float


@dataclass(frozen=True, slots=True)
class DensityBlock:
    """Sum of dyads in the (target_endpoint, source_endpoint) PW block."""
    target_endpoint: int
    source_endpoint: int
    dyads: tuple[DensityDyad, ...]


def density_dyads_item(index, psi, plus, minus, mixed, D, *, target_plus,
                       target_minus, source_weights):
    """Borrow one endpoint's blocks, reading only its incident source jets.

    Sequences may be lazy authenticated readers. Only mixed[index],
    plus[index], minus[index], plus[target_minus[index]] and
    minus[target_plus[index]] are requested. Incoming norm dyads retain the
    incoming source's raw rank, D and weight, even when ranks differ.
    """
    count = len(psi)
    p, m = _permutations(count, target_plus, target_minus)
    if (not isinstance(index, (int, np.integer))
            or isinstance(index, (bool, np.bool_)) or not 0 <= index < count):
        raise ValueError("explicit original source index required")
    if any(len(x) != count for x in (plus, minus, mixed, D, source_weights)):
        raise ValueError("complete source-ordered readers required")
    i = int(index)
    sources = set((i, p[i], m[i]))
    frames = {s: psi[s] for s in sources}
    for s, frame in frames.items():
        if not isinstance(frame, np.ndarray) or frame.ndim != 2 or min(frame.shape) < 1:
            raise ValueError("original full-row mean frame required")
    coefficients = {s: D[s] for s in sources}
    weights = {s: source_weights[s] for s in sources}
    for s in sources:
        d, weight = coefficients[s], weights[s]
        rank = frames[s].shape[1]
        if not isinstance(d, np.ndarray) or d.shape != (rank, rank):
            raise ValueError(f"D[{s}] must retain original rank {rank}")
        if not np.isfinite(d).all() or not np.allclose(d, d.conj().T, atol=1e-12, rtol=1e-12):
            raise ValueError("finite Hermitian fixed D required; no repair is made")
        if not np.isscalar(weight) or not np.isreal(weight) or not np.isfinite(weight) or weight <= 0:
            raise ValueError("explicit finite positive real source weights required")
    local_plus, local_minus, local_mixed = plus[i], minus[i], mixed[i]
    incoming_plus = local_plus if m[i] == i else plus[m[i]]
    incoming_minus = local_minus if p[i] == i else minus[p[i]]
    for name, value, shape in (
            ("plus", local_plus, (frames[p[i]].shape[0], frames[i].shape[1])),
            ("minus", local_minus, (frames[m[i]].shape[0], frames[i].shape[1])),
            ("mixed", local_mixed, frames[i].shape),
            ("incoming plus", incoming_plus, (frames[i].shape[0], frames[m[i]].shape[1])),
            ("incoming minus", incoming_minus, (frames[i].shape[0], frames[p[i]].shape[1]))):
        if not isinstance(value, np.ndarray) or value.shape != shape:
            raise ValueError(f"{name} must retain original-source shape {shape}")

    def dyad(left, right, source):
        return DensityDyad(left, coefficients[source], right, source, float(weights[source]))

    return dict(
        mean=DensityBlock(i, i, (dyad(frames[i], frames[i], i),)),
        plus=DensityBlock(p[i], i, (
            dyad(local_plus, frames[i], i),
            dyad(frames[p[i]], incoming_minus, p[i]))),
        minus=DensityBlock(m[i], i, (
            dyad(local_minus, frames[i], i),
            dyad(frames[m[i]], incoming_plus, m[i]))),
        mixed=DensityBlock(i, i, (
            dyad(local_mixed, frames[i], i), dyad(frames[i], local_mixed, i),
            dyad(incoming_plus, incoming_plus, m[i]),
            dyad(incoming_minus, incoming_minus, p[i]))))


def density_dyads(psi, plus, minus, mixed, D, *, target_plus, target_minus,
                  source_weights):
    """Return mean/plus/minus/mixed tuples of copy-free DensityBlock records.

    Physical eta=lambda* swaps bra directions before dagger. Mixed denotes
    d²/(dλ dη), without a half. Every original D and weight occurs once.
    """
    p, m = _layout(psi, plus, minus, mixed, target_plus, target_minus)
    rows = tuple(density_dyads_item(i, psi, plus, minus, mixed, D,
                 target_plus=p, target_minus=m, source_weights=source_weights)
                 for i in range(len(psi)))
    return {name: tuple(row[name] for row in rows)
            for name in ("mean", "plus", "minus", "mixed")}
