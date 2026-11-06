"""Physical density, ionic, field and bare-H0 moments of borrowed blocks.

This is a contraction layer, not an evolution or material-result claim.
The caller owns the actual directed coupling, endpoint H0 and field vertex.
No baseline contraction, self-consistent potential or energy residual is
constructed here. DensityDyad's original source D and weight are passed once.
"""
import numpy as np
from scipy.linalg.blas import zgemm

from substrate_framework.material_physical_operators import PhysicalFieldVertices, _require_physical
from substrate_framework.qe_bloch_finiteq_ionic import BlochFiniteQIonic
from mixed_dyad_grid_moments import grid_moments
from mixed_dyad_projector_moments import field_projector_moments
from mixed_moment_terms import moment_terms


def _geometry(block, coupling, field_vertex, diagonal):
    if not isinstance(diagonal, (bool, np.bool_)):
        raise TypeError('diagonal must be an explicit geometry boolean')
    if not isinstance(coupling, BlochFiniteQIonic):
        raise TypeError('actual directed BlochFiniteQIonic required')
    _require_physical(coupling.source_factors)
    _require_physical(coupling.target_factors)
    same_endpoint = block.target_endpoint == block.source_endpoint
    gamma = np.all(np.asarray(coupling.edge.q_fractional) == 0)
    if diagonal != same_endpoint or diagonal != gamma:
        raise ValueError('block endpoints, directed q and diagonal geometry disagree')
    if not block.dyads:
        raise ValueError('at least one borrowed density dyad required')
    if diagonal:
        if not isinstance(field_vertex, PhysicalFieldVertices):
            raise TypeError('diagonal Gamma requires the actual physical field vertex')
        coupling._check_factors(field_vertex.hamiltonian.factors, coupling.edge.source_basis)
        coupling._check_factors(field_vertex.hamiltonian.factors, coupling.edge.target_basis)
    elif field_vertex is not None:
        raise ValueError('off-diagonal nonGamma has no global field trace')


def _one_body(dyad, bare_action, block_size):
    """weight Tr(H0 left D right.H), using bounded bare-H0 action columns.

    H0 is the actual zero-field action, never a stationary eigenvalue. This
    direct ordering does not rely on Hermiticity of a complex tangent drive.
    No full weighted PW frame or PW-square density is needed.
    """
    left, right, coefficient = dyad.left, dyad.right, dyad.D
    value = 0j
    for start in range(0, left.shape[1], block_size):
        ls = slice(start, start + block_size)
        hleft = bare_action(left[:, ls])
        if not isinstance(hleft, np.ndarray) or hleft.shape != left[:, ls].shape:
            raise ValueError('bare H0 must act on arbitrary full-row column blocks')
        if not np.isfinite(hleft).all():
            raise ValueError('bare H0 returned nonfinite values')
        for other in range(0, right.shape[1], block_size):
            rs = slice(other, other + block_size)
            overlap = zgemm(1., right[:, rs], hleft, trans_a=2)
            value += np.einsum('ab,ba->', coefficient[ls, rs], overlap)
    return dyad.weight * value


def response_blocks(groups, *, coupling, field_vertex, diagonal):
    """Contract independent dyad groups sharing ONE actual directed coupling.

    Original coefficient/weight ownership and all ordered terms are retained.
    Request-local grid/projector factors are shared, never trajectories or
    conjugate responses. NonGamma directions use separate calls. No H0 action
    or one-body energy is evaluated by this response-only API.
    """
    if not groups:
        raise ValueError('at least one original density block required')
    for block in groups.values():
        _geometry(block, coupling, field_vertex, diagonal)
    terms = moment_terms(coupling, groups)
    grid, kinetic = grid_moments(coupling, terms,
                                 field_momenta=field_vertex.momenta if diagonal else None)
    field = field_projector_moments(field_vertex, terms, kinetic) if diagonal else {}
    return {name: dict(density=grid[name]['density'], gradient=grid[name]['gradient'],
                       field_gradient_electronic=field[name] if diagonal else np.zeros(3, complex))
            for name in groups}
