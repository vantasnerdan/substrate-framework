"""One electron-ion skeleton diagram, its tangent, and its bilinear transpose.

hbar=1. H_int=sum(q,k,spin,i) u_qi c†_(k+q) g_i(q,k)c_k/sqrt(N).
g is the POSITIVE Hamiltonian gradient, not the force. G^>=-i<c c†>,
G^<=i<c† c>, D^>=-i<du_q du_-q>, D^<=-i<du_-q du_q>.
All matrices are in original orthonormal Cartesian/band coordinates. No spin
multiplicity or quadrature weight is implicit. Vertices below are UNSCALED.

An edge contributes Sigma at target=k+q and Pi at q. Sum every q exactly once
for Sigma, every explicit spin/k exactly once for Pi. Pi uses the opposite
fermionic ordered component at REVERSED times. Do not also add a mirrored edge
contribution to Sigma: the full q/source enumeration already includes it.
The same formulas apply to any contour pair when the forward/reverse contour
components are supplied; this module performs no KMS/thermal reconstruction.
For an initially factorized Gaussian IVP there is no correlated imaginary-leg
injection. The lesser collision injection is retained, not inferred from R.

Source conventions: NESSi Holstein_impurity_impl.cpp:24-50 (without its local
one-orbital spin2 assumption), Schueler2020 Eqs12,15,19; A29
reduced_electronic_dynamics.py:147-174,181-199. Phi's one closed electron loop
has the usual half symmetry factor; cutting an electron/boson line gives the
+i g G g D / -i tr(g G g G) pair implemented together here. Hartree/mean force
belongs to the existing mean scalar; D is CONNECTED and contains no mean dyad.
"""
from dataclasses import dataclass, fields

import numpy as np


@dataclass(frozen=True)
class EdgePair:
    """Independent legs, including response variations not obeying conjugacy.

    left_t,left_s: [i,target,source]. right_t,right_s: [i,source,target].
    gs_*: G_source(t,s); gt_*: G_target(t,s); back_*: G_source(s,t).
    d_*: D_q(t,s) [i,j]. Stored ionic_reverse[i,a,b] must be TRANSPOSED
    on band axes to obtain right[i,b,a], never conjugated or manufactured.
    """
    left_t: np.ndarray
    right_t: np.ndarray
    left_s: np.ndarray
    right_s: np.ndarray
    gs_greater: np.ndarray
    gs_lesser: np.ndarray
    gt_greater: np.ndarray
    gt_lesser: np.ndarray
    back_greater: np.ndarray
    back_lesser: np.ndarray
    d_greater: np.ndarray
    d_lesser: np.ndarray


@dataclass(frozen=True)
class Kernels:
    sigma_greater: np.ndarray
    sigma_lesser: np.ndarray
    pi_greater: np.ndarray
    pi_lesser: np.ndarray
    sigma_retarded: np.ndarray
    pi_retarded: np.ndarray


# These four monomials are the ENTIRE functional. Differentiation below acts on
# their actual operands (including all four independently supplied vertex legs).
_TERMS = (
    ('sigma_greater', 1j, 'iab,bc,jcd,ij->ad',
     ('left_t', 'gs_greater', 'right_s', 'd_greater')),
    ('sigma_lesser', 1j, 'iab,bc,jcd,ij->ad',
     ('left_t', 'gs_lesser', 'right_s', 'd_lesser')),
    ('pi_greater', -1j, 'iba,ad,jdc,cb->ij',
     ('right_t', 'gt_greater', 'left_s', 'back_lesser')),
    ('pi_lesser', -1j, 'iba,ad,jdc,cb->ij',
     ('right_t', 'gt_lesser', 'left_s', 'back_greater')),
)


def _validate(pair, n):
    if type(n) is not int or n < 1:
        raise ValueError('original positive integer cell count N required')
    shape = np.shape(pair.left_t)
    if len(shape) != 3 or min(shape) < 1:
        raise ValueError('left vertex must have [Cartesian mode,target,source] axes')
    d, a, b = shape
    expected = dict(left_t=(d,a,b), left_s=(d,a,b), right_t=(d,b,a),
                    right_s=(d,b,a), gs_greater=(b,b), gs_lesser=(b,b),
                    gt_greater=(a,a), gt_lesser=(a,a), back_greater=(b,b),
                    back_lesser=(b,b), d_greater=(d,d), d_lesser=(d,d))
    for name, wanted in expected.items():
        value = np.asarray(getattr(pair, name))
        if value.shape != wanted or value.dtype.kind not in 'fci' or not np.isfinite(value).all():
            raise ValueError('invalid finite numeric edge operand: '+name)


def _causal(t, s):
    if not np.isfinite(t) or not np.isfinite(s):
        raise ValueError('finite real times required')
    return float(t >= s)


def _pack(ordered, t, s):
    # Native modified R at equal time; theta(0)=1, no half-valued jump.
    theta = _causal(t, s)
    return Kernels(**ordered,
        sigma_retarded=theta*(ordered['sigma_greater']-ordered['sigma_lesser']),
        pi_retarded=theta*(ordered['pi_greater']-ordered['pi_lesser']))


def evaluate_components(pair, *, n, t, s, components=('sigma', 'pi')):
    """Compute only requested diagrams; absent components are not fake zeros."""
    _validate(pair, n)
    theta = _causal(t, s)
    if components not in (('sigma',), ('pi',), ('sigma', 'pi')):
        raise ValueError('request sigma, pi, or both in canonical order')
    ordered = {name: factor/n*np.einsum(expression,
        *(getattr(pair, key) for key in keys), optimize=True)
        for name, factor, expression, keys in _TERMS
        if name.split('_', 1)[0] in components}
    for component in components:
        ordered[component+'_retarded'] = theta*(
            ordered[component+'_greater']-ordered[component+'_lesser'])
    return ordered


def evaluate(pair, *, n, t, s):
    """Both ordered selfenergies and causal retarded pair from one interaction.

    Sigma units Ry^2; Pi units Ry^2/bohr^2 for g in Ry/bohr and D in bohr^2.
    D Dyson convention D=D0+D0*Pi*D; no oscillator mass/frequency prefactor.
    """
    return Kernels(**evaluate_components(pair, n=n, t=t, s=s))


def jvp(pair, direction, *, n, t, s):
    """Exact complex-linear Frechet derivative, with N and times held fixed."""
    _validate(pair, n)
    _validate(direction, n)
    ordered = {}
    for name, factor, expression, keys in _TERMS:
        operands = [getattr(pair, key) for key in keys]
        total = None
        for k, key in enumerate(keys):
            args = operands.copy()
            args[k] = getattr(direction, key)
            term = factor/n*np.einsum(expression, *args, optimize=True)
            total = term if total is None else total+term
        ordered[name] = total
    return _pack(ordered, t, s)


def transpose(pair, dual, *, n, t, s):
    """True J.T dual under sum(z*x), NOT vdot; includes retarded readouts.

    Independent reversed-time legs remain independent in this local API. A
    history assembler must scatter-add their cotangents into the corresponding
    reversed-time history entries, and sum every edge sharing a variable.
    """
    _validate(pair, n)
    theta = _causal(t, s)
    weights = dict(sigma_greater=dual.sigma_greater+theta*dual.sigma_retarded,
        sigma_lesser=dual.sigma_lesser-theta*dual.sigma_retarded,
        pi_greater=dual.pi_greater+theta*dual.pi_retarded,
        pi_lesser=dual.pi_lesser-theta*dual.pi_retarded)
    result = {f.name: np.zeros_like(getattr(pair, f.name), dtype=complex) for f in fields(pair)}
    for name, factor, expression, keys in _TERMS:
        inputs, output = expression.split('->')
        indices = inputs.split(',')
        for k, key in enumerate(keys):
            others = [j for j in range(len(keys)) if j != k]
            contraction = ','.join([output]+[indices[j] for j in others])+'->'+indices[k]
            args = [weights[name]]+[getattr(pair, keys[j]) for j in others]
            result[key] += factor/n*np.einsum(contraction, *args, optimize=True)
    return EdgePair(**result)


def replace_represented_polarization(new, represented):
    """Dyson insertion replacing, never stacking, an ALREADY included Pi.

    Both arguments have identical contour components/units. A frozen clamped
    contact Hessian <d2H/du2> is not a bubble and must NOT be subtracted here.
    For an explicit-electron tangent driver, its represented electronic bubble
    must be removed before adding this functional (or this difference used).
    This API does not invent an unprovided DFT exchange-correlation counterterm.
    """
    return Kernels(**{f.name: getattr(new, f.name)-getattr(represented, f.name)
                      for f in fields(new)})
