"""Physical density, ionic, field and bare-H0 moments of borrowed blocks.

This is a contraction layer, not an evolution or material-result claim.
The caller owns the actual directed coupling, endpoint H0 and field vertex.
No baseline contraction, self-consistent potential or energy residual is
constructed here. DensityDyad's original source D and weight are passed once.
"""
import numpy as np
from scipy.linalg.blas import zgemm

from substrate_framework.material_physical_operators import _require_physical
from substrate_framework.qe_bloch_finiteq_ionic import BlochFiniteQIonic


def _geometry(block, coupling, field_moments, bare_action, diagonal, energy):
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
        if not callable(field_moments) or (energy and not callable(bare_action)):
            raise TypeError('diagonal Gamma requires physical field moments and energy requires bare H0')
    elif field_moments is not None:
        raise ValueError('off-diagonal nonGamma has no global field trace callback')


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


def _contract(block, *, coupling, field_moments, bare_action, diagonal, energy):
    """One physical dyad loop; bare one-body work is energy-readout only."""
    _geometry(block, coupling, field_moments, bare_action, diagonal, energy)
    density = np.zeros((coupling.bridge.source_fft.nnr, 2), complex)
    gradient = np.zeros(coupling.ndof, complex)
    field = np.zeros(3, complex)
    one_body = 0j
    for dyad in block.dyads:
        n, g = coupling.moments(dyad.left, dyad.right,
                                coefficient=dyad.D, weight=dyad.weight)
        density += n
        gradient += g
        if diagonal:
            ell = np.asarray(field_moments(dyad.left, dyad.right,
                                          dyad.D, weight=dyad.weight))
            if ell.shape != (3,) or not np.isfinite(ell).all():
                raise ValueError('physical field moments must return finite complex Cartesian (3,)')
            field += ell
            if energy:
                one_body += _one_body(dyad, bare_action, coupling.bridge.block_size)
    result = dict(density=density, gradient=gradient,
                  field_gradient_electronic=field)
    if energy:
        result['one_body'] = one_body
    return result


def response_block(block, *, coupling, field_moments, diagonal):
    """Density/positive ionic gradient/L traces, with NO bare H0 actions.

    NonGamma global L traces are geometrically absent. No one_body entry is
    returned: this API does not evaluate or claim a one-body energy.
    """
    return _contract(block, coupling=coupling, field_moments=field_moments,
                     bare_action=None, diagonal=diagonal, energy=False)


def contract_block(block, *, coupling, field_moments, bare_action, diagonal):
    """Return genuine density/gradient/field_gradient_electronic/one_body.

    ``coupling`` is the actual directed physical_finiteq, including q=0 for
    diagonal blocks. ``field_moments`` binds physical_field_dyad_moments to
    matching PhysicalFieldVertices; ``bare_action`` is actual H0.apply.
    NonGamma global H0/L traces vanish by endpoint geometry. All complex
    defects remain; original source D and quadrature weight occur once.
    """
    return _contract(block, coupling=coupling, field_moments=field_moments,
                     bare_action=bare_action, diagonal=diagonal, energy=True)
