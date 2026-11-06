"""Independent dense finite models for the mixed-frame scientific oracle.

This is oracle-only code: full matrices/exponentials are intentional here,
never production PW density implementations.  Models are seeded asymmetric
complex finite Hamiltonians, not material simulations or physical energy data.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag, expm


@dataclass(frozen=True)
class FiniteModel:
    dimensions: tuple
    ranks: tuple
    target_plus: tuple
    target_minus: tuple
    row_slices: tuple
    column_slices: tuple
    frames: tuple
    D: tuple
    weights: tuple
    h0: np.ndarray
    hp: np.ndarray
    hm: np.ndarray
    hx: np.ndarray
    initial_frame: np.ndarray
    weighted_D: np.ndarray


def _slices(sizes):
    edges = np.cumsum((0,) + tuple(sizes))
    return tuple(slice(int(a), int(b)) for a, b in zip(edges[:-1], edges[1:]))


def finite_model(dimensions=(2, 3, 4), *, gamma=False, seed=25428):
    rng = np.random.default_rng(seed)

    def complex_array(shape):
        return (rng.normal(size=shape) + 1j * rng.normal(size=shape)) / 3

    def hermitian(size):
        a = complex_array((size, size))
        return (a + a.conj().T) / 2

    ranks = (1, 2, 3)
    p = (0, 1, 2) if gamma else (1, 2, 0)
    m = tuple(p.index(i) for i in range(3))
    rows, columns = _slices(dimensions), _slices(ranks)
    frames = tuple(complex_array((n, r)) for n, r in zip(dimensions, ranks))
    # Outer products give exactly rank-one complex Hermitian D; no inverse,
    # positivity repair, eigenvector gauge choice, or rank truncation is used.
    vectors = tuple(complex_array((r, 1)) for r in ranks)
    d = tuple(v @ v.conj().T for v in vectors)
    weights = (0.17, 0.31, 0.52)
    n = sum(dimensions)
    h0 = block_diag(*(hermitian(size) for size in dimensions))
    hx = block_diag(*(hermitian(size) for size in dimensions))
    hp = np.zeros((n, n), complex)
    for i in range(3):
        hp[rows[p[i]], rows[i]] = complex_array((dimensions[p[i]], dimensions[i]))
    initial = block_diag(*frames)
    return FiniteModel(tuple(dimensions), ranks, p, m, rows, columns, frames,
                       d, weights, h0, hp, hp.conj().T, hx, initial,
                       block_diag(*(w * a for w, a in zip(weights, d))))


def physical_actions(model):
    """Actual matrix actions of the declared finite-model Hamiltonians."""
    rows, p, m = model.row_slices, model.target_plus, model.target_minus

    def h0(i, vectors):
        return model.h0[rows[i], rows[i]] @ vectors

    def hplus(i, vectors):
        return model.hp[rows[p[i]], rows[i]] @ vectors

    def hminus(i, vectors):
        return model.hm[rows[m[i]], rows[i]] @ vectors

    def hmixed(i, vectors):
        return model.hx[rows[i], rows[i]] @ vectors

    return dict(h0=h0, hplus=hplus, hminus=hminus, hmixed=hmixed)


def embed_frames(model, stages):
    """Independently embed source-column frames in a single dense frame."""
    result = []
    for blocks, targets in zip(stages, (range(3), model.target_plus,
                                       model.target_minus, range(3))):
        full = np.zeros_like(model.initial_frame)
        for source, target in enumerate(targets):
            full[model.row_slices[target], model.column_slices[source]] = blocks[source]
        result.append(full)
    return tuple(result)


def extract_frames(model, full_stages):
    return tuple(tuple(full[model.row_slices[target], model.column_slices[source]]
                       for source, target in enumerate(targets))
                 for full, targets in zip(full_stages, (range(3), model.target_plus,
                                                       model.target_minus, range(3))))


def exact_frame_jets(model, time):
    """Exponential derivative via a full block-triangular Fréchet lift.

    No endpoint product rules or production density formulas are used.  The
    four-variable lift is obtained by multiplying the polynomial H(λ,η) by a
    polynomial ket and retaining coefficients 1,λ,η,λη.
    """
    n = model.h0.shape[0]
    zero = np.zeros_like(model.h0)
    generator = -1j * np.block([
        [model.h0, zero, zero, zero],
        [model.hp, model.h0, zero, zero],
        [model.hm, zero, model.h0, zero],
        [model.hx, model.hm, model.hp, model.h0]])
    initial = np.concatenate((model.initial_frame,) + (np.zeros_like(model.initial_frame),) * 3)
    solution = expm(time * generator) @ initial
    return tuple(solution[k * n:(k + 1) * n] for k in range(4))


def evolve_production(model, time, frame_rhs):
    """Integrate the supplied production RHS and actual moving reference."""
    zeros = tuple(np.zeros_like(x) for x in model.frames)
    plus = tuple(np.zeros((model.dimensions[t], model.ranks[i]), complex)
                 for i, t in enumerate(model.target_plus))
    minus = tuple(np.zeros((model.dimensions[t], model.ranks[i]), complex)
                  for i, t in enumerate(model.target_minus))
    template = (model.frames, plus, minus, zeros)
    shapes = tuple(tuple(a.shape for a in stage) for stage in template)
    sizes = tuple(tuple(a.size for a in stage) for stage in template)
    actions = physical_actions(model)

    def pack(stages):
        return np.concatenate([a.ravel() for stage in stages for a in stage])

    def unpack(flat):
        cursor, stages = 0, []
        for shape_stage, size_stage in zip(shapes, sizes):
            blocks = []
            for shape, size in zip(shape_stage, size_stage):
                blocks.append(flat[cursor:cursor + size].reshape(shape))
                cursor += size
            stages.append(tuple(blocks))
        return tuple(stages)

    def rhs(t, flat):
        stages = unpack(flat)
        reference = tuple(-1j * actions['h0'](i, a) for i, a in enumerate(stages[0]))
        response = frame_rhs(*stages, target_plus=model.target_plus,
                             target_minus=model.target_minus, **actions)
        return pack((reference,) + response)

    solution = solve_ivp(rhs, (0, time), pack(template), method='DOP853',
                         atol=2e-12, rtol=2e-12, t_eval=[time])
    if not solution.success:
        raise RuntimeError(solution.message)
    return unpack(solution.y[:, -1]), {'nfev': solution.nfev, 'message': solution.message}


def dense_density_ode(model, time):
    """Independent full Liouville product differentiation, not frame dyads."""
    rho = model.initial_frame @ model.weighted_D @ model.initial_frame.conj().T
    n = rho.shape[0]

    def commutator(a, b):
        return a @ b - b @ a

    def rhs(t, flat):
        r0, rp, rm, rx = flat.reshape(4, n, n)
        return (-1j * np.stack((
            commutator(model.h0, r0),
            commutator(model.h0, rp) + commutator(model.hp, r0),
            commutator(model.h0, rm) + commutator(model.hm, r0),
            commutator(model.h0, rx) + commutator(model.hp, rm)
            + commutator(model.hm, rp) + commutator(model.hx, r0)))).ravel()

    initial = np.stack((rho, np.zeros_like(rho), np.zeros_like(rho), np.zeros_like(rho)))
    solution = solve_ivp(rhs, (0, time), initial.ravel(), method='DOP853',
                         atol=2e-12, rtol=2e-12, t_eval=[time])
    if not solution.success:
        raise RuntimeError(solution.message)
    return tuple(solution.y[:, -1].reshape(4, n, n))


def analytic_density(model, time, lam, eta):
    """Holomorphic continuation U rho0 U^-1; NOT U rho0 U† off slice."""
    h = model.h0 + lam * model.hp + eta * model.hm + lam * eta * model.hx
    rho0 = model.initial_frame @ model.weighted_D @ model.initial_frame.conj().T
    return expm(-1j * time * h) @ rho0 @ expm(1j * time * h)


def materialize_descriptors(model, descriptors):
    """Oracle-only dense contraction, deliberately outside production."""
    stages = []
    n = model.h0.shape[0]
    for name in ('mean', 'plus', 'minus', 'mixed'):
        full = np.zeros((n, n), complex)
        for block in descriptors[name]:
            target, source = model.row_slices[block.target_endpoint], model.row_slices[block.source_endpoint]
            for d in block.dyads:
                full[target, source] += d.weight * d.left @ d.D @ d.right.conj().T
        stages.append(full)
    return tuple(stages)
