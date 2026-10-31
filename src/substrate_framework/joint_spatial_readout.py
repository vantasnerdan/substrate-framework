"""Full-grid postprocessing only: no evolution, symmetry completion or repair."""
import numpy as np

from .joint_response_observables import spatial_pair
from .spatial_readout_metrics import periodic_distances, selected_weight

ORIGINS = ('electronic', 'ionic', 'field', 'total')


def complete_grid(q_fractional, mesh):
    """Require exactly the original Cartesian character group, modulo integers."""
    shape = np.asarray(mesh)
    if (shape.shape != (3,) or not np.issubdtype(shape.dtype, np.integer)
            or np.any(shape < 1)):
        raise ValueError('three positive original integer mesh dimensions required')
    n = int(np.prod(shape))
    q = np.asarray(q_fractional, float)
    if q.shape != (n, 3) or not np.isfinite(q).all():
        raise ValueError('complete original q grid required, no selected-q reconstruction')
    scaled = np.mod(q, 1)*shape
    nearest = np.rint(scaled)
    if not np.allclose(scaled, nearest, rtol=0, atol=1e-11):
        raise ValueError('q is not an original mesh character')
    indices = nearest.astype(int) % shape
    if len({tuple(row) for row in indices}) != n:
        raise ValueError('duplicate/missing original q character')
    return q, np.indices(tuple(shape)).reshape(3, -1).T


def coherent_local(q, values, cells):
    """Inverse orthonormal transform: kick phase and 1/sqrtN already in values.

    local_momentum_kick seeds p_q=A exp(-iq.R0)/sqrtN. A unit response
    therefore reconstructs A delta_(R,R0), not A/sqrtN or sqrtN*A.
    """
    z = np.asarray(values, complex)
    if z.shape != (len(q),) or not np.isfinite(z).all():
        raise ValueError('finite coherent force per q required')
    return np.exp(2j*np.pi*np.asarray(cells) @ np.asarray(q).T) @ z/np.sqrt(len(q))


def reconstruct(q, coherent, ab, ba, *, mesh, cell_bohr, kick_cell, reference_cell,
                probe_position_bohr, seed_position_bohr):
    """AB/BA origins share one probe; coherent extent is relative to its seed."""
    q, cells = complete_grid(q, mesh)
    n = len(q)
    a, b = np.asarray(ab, complex), np.asarray(ba, complex)
    if a.shape != (n, 4) or b.shape != (n, 4) or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('all initial origins and totals required for both orderings')
    for cell in (kick_cell, reference_cell):
        raw = np.asarray(cell)
        if raw.shape != (3,) or not np.array_equal(raw, np.rint(raw)):
            raise ValueError('integer original cell indices required')
    local = coherent_local(q, coherent, cells)
    pair_a = np.asarray([spatial_pair(q, a, cell, reference_cell, ncells=n) for cell in cells])
    pair_b = np.asarray([spatial_pair(q, b, cell, reference_cell, ncells=n) for cell in cells])
    pairs = dict(AB=pair_a, BA=pair_b, symmetrized=(pair_a+pair_b)/2,
                 commutator=pair_a-pair_b)
    offset = np.asarray(probe_position_bohr, float)-np.asarray(seed_position_bohr, float)
    force_distances = periodic_distances(cells, kick_cell, mesh, cell_bohr, offset_bohr=offset)
    correlation_distances = periodic_distances(cells, reference_cell, mesh, cell_bohr)
    indices = np.rint(np.mod(q, 1)*np.asarray(mesh)).astype(int) % mesh
    lookup = {tuple(index): i for i, index in enumerate(indices)}
    minus = [lookup[tuple(-index % mesh)] for index in indices]
    coherent_q = np.asarray(coherent, complex)
    residuals = dict(
        coherent_imaginary_max=float(np.max(np.abs(local.imag))),
        coherent_imaginary_l2=float(np.linalg.norm(local.imag)),
        coherent_q_reality_max=float(np.max(np.abs(coherent_q-coherent_q[minus].conj()))),
        origin_sum_q_max={name: float(np.max(np.abs(x[:, :3].sum(axis=1)-x[:, 3])))
                          for name, x in (('AB', a), ('BA', b))},
        AB_minus_conjugate_BA_max_by_origin=np.max(np.abs(pair_a-pair_b.conj()), axis=0),
        symmetrized_imaginary_max_by_origin=np.max(np.abs(pairs['symmetrized'].imag), axis=0),
        commutator_real_max_by_origin=np.max(np.abs(pairs['commutator'].real), axis=0),
        origin_sum_spatial_max={name: float(np.max(np.abs(x[:, :3].sum(axis=1)-x[:, 3])))
                                for name, x in pairs.items()})
    return dict(cells=cells, mesh=list(mesh), ncells=n, origins=list(ORIGINS),
        reference_cell=list(reference_cell), kick_cell=list(kick_cell),
        coherent_distance_probe_minus_seed_bohr=offset,
        coherent_force_ry_bohr=local,
        coherent_components='Producer retains total coherent readout only; no electronic/ionic/field decomposition is fabricated.',
        connected_force_pairs_ry2_bohr2=pairs, reality_residuals=residuals,
        selected_observable_weights=dict(
            coherent_total=selected_weight(local, force_distances, value_units='Ry/bohr'),
            connected={name: {origin: selected_weight(values[:, i], correlation_distances,
                                                      value_units='(Ry/bohr)^2')
                              for i, origin in enumerate(ORIGINS)}
                       for name, values in pairs.items()}))
