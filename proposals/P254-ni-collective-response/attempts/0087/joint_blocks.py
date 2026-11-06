"""A41's same Sigma monomial, restricted to IF, FI and FF coordinate blocks.

Canonical A50 vertices already contain ionic sqrt(N) and field normalization.
The rectangular contraction below therefore has n=1, exactly as A56 common_pair.
No ionic-ionic result is evaluated or subtracted; Pi retains the full 114 bubble.
"""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/n) for n in (
    '0062', '0063', '0060', '0056', '0055', '0054', '0053', '0052',
    '0051', '0050', '0048', '0046', '0045', '0044', '0043', '0042',
    '0041', '0030', '0029', '0028', '0027', '0025', '0024')]
import numpy as np
from electron_ion_functional import EdgePair, evaluate_components, _TERMS
from joint_gamma_adapter import vertices

# Owning primitive is evaluate_components' explicit einsum, not a separate
# private rectangular API. Reuse its literal factor/expression and operands.
_SIGMA_FACTOR, _SIGMA_EXPRESSION = _TERMS[0][1:3]


def sigma_rectangular(left, electron, right, boson):
    """Same i g_i G g_j D_ij monomial with independent rectangular i,j."""
    left, electron, right, boson = map(np.asarray, (left, electron, right, boson))
    if (left.ndim != 3 or right.ndim != 3 or electron.ndim != 2 or
            electron.shape != (left.shape[2], right.shape[1]) or
            boson.shape != (left.shape[0], right.shape[0]) or
            left.shape[1] != right.shape[2] or
            not all(np.isfinite(x).all() for x in (left, electron, right, boson))):
        raise ValueError('finite independent rectangular Sigma operands required')
    return _SIGMA_FACTOR*np.einsum(_SIGMA_EXPRESSION, left, electron, right,
                                  boson, optimize=True)


def sigma_correction(left, electron, right, boson, ionic_directions):
    """IF + FI + FF; preserve D_IF versus D_FI and raw reverse band order."""
    d = ionic_directions
    if not 0 < d < len(left) or len(left) != len(right):
        raise ValueError('nonempty ionic and field blocks required')
    ionic, field = slice(None, d), slice(d, None)
    return (sigma_rectangular(left[ionic], electron, right[field], boson[ionic, field])
            + sigma_rectangular(left[field], electron, right[ionic], boson[field, ionic])
            + sigma_rectangular(left[field], electron, right[field], boson[field, field]))


def joint_pair(c, reference, a, b, canonical_vertices=None):
    """Once-global full Pi plus stable joint-minus-ionic Sigma, all128 slots."""
    forward, reverse = vertices(c) if canonical_vertices is None else canonical_vertices
    dg, dl = reference.ordered(0, a, b)
    gg, gl = reference.electron(a, b)
    bg, bl = reference.electron(b, a)
    result = {'sigma_'+order: np.zeros((c.ns, c.m, c.m), np.complex128)
              for order in ('greater', 'lesser', 'retarded')}
    result.update({'pi_'+order: np.zeros((c.d+c.f, c.d+c.f), np.complex128)
                   for order in ('greater', 'lesser', 'retarded')})
    for source in range(c.ns):
        target = int(c.target[0, source])
        left = forward[source]
        right = reverse[source].transpose(0, 2, 1)
        for order, electron, boson in (('greater', gg[source], dg),
                                        ('lesser', gl[source], dl)):
            result['sigma_'+order][target] += sigma_correction(
                left, electron, right, boson, c.d)
        pair = EdgePair(left, right, left, right, gg[source], gl[source],
                        gg[target], gl[target], bg[source], bl[source], dg, dl)
        for name, value in evaluate_components(
                pair, n=1, t=reference.times[a], s=reference.times[b],
                components=('pi',)).items():
            result[name] += value
    result['sigma_retarded'] = (result['sigma_greater']-result['sigma_lesser']
                               if reference.times[a] >= reference.times[b]
                               else np.zeros_like(result['sigma_greater']))
    return result
