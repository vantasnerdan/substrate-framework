"""Connected cubic force polynomial and its exact bilinear differential.

V=Phi_ijk u_i u_j u_k/6; force=-Phi_ijk u_j u_k/2.
Inputs are RAW ordered W=<delta u delta u>, not -iW Green functions.
Both vertices and both ordered legs are independent; never infer a dagger.
Outputs greater/lesser=-i<delta F delta F>; retarded=theta*(greater-lesser).
Cartesian bohr coordinates, Ry/bohr^3 vertices; no mass/sector/N factors here:
canonical 1/sqrt(N) is ALREADY in the material vertices.
"""
import numpy as np


def _terms():
    terms = [('mean_left', -.5, 'iab,ab->i', ('left', 'equal_left')),
             ('mean_left', -.5, 'iab,a,b->i', ('left', 'mu_left', 'mu_left')),
             ('mean_right', -.5, 'jcd,cd->j', ('right', 'equal_right')),
             ('mean_right', -.5, 'jcd,c,d->j', ('right', 'mu_right', 'mu_right'))]
    for output, leg in [('greater', 'wg'), ('lesser', 'wl')]:
        for equation in ('iab,jcd,ac,bd->ij', 'iab,jcd,ad,bc->ij'):
            terms.append((output, -.25j, equation, ('left', 'right', leg, leg)))
        for a, c, w in [('a','c','bd'), ('a','d','bc'), ('b','c','ad'), ('b','d','ac')]:
            terms.append((output, -.25j, f'iab,jcd,{a},{c},{w}->ij',
                          ('left', 'right', 'mu_left', 'mu_right', leg)))
    return terms


TERMS = _terms()


def evaluate(data, theta, terms=TERMS):
    """Return both means and connected ordered/retarded kernels.

    left[I,A,B], right[J,C,D]; wg/wl[A,C] (same shape), equal_left[A,B],
    equal_right[C,D], mu_left[A], mu_right[C]. theta is causal 0 or 1,
    with theta(0)=1 matching the existing A41 convention.
    """
    out = {}
    for key, factor, equation, names in terms:
        value = factor*np.einsum(equation, *(data[n] for n in names), optimize=True)
        out[key] = out.get(key, 0) + value
    out['retarded'] = theta*(out['greater']-out['lesser'])
    return out


def jvp(data, tangent, theta, terms=TERMS):
    """Analytic Frechet derivative in every independent input (complex linear)."""
    out = {k: np.zeros_like(v) for k, v in evaluate(data, theta, terms).items()}
    for key, factor, equation, names in terms:
        for slot, name in enumerate(names):
            if name not in tangent:
                continue
            operands = [data[n] for n in names]
            operands[slot] = tangent[name]
            out[key] += factor*np.einsum(equation, *operands, optimize=True)
    out['retarded'] = theta*(out['greater']-out['lesser'])
    return out


def transpose(data, cotangent, theta, terms=TERMS):
    """True nonconjugating transpose: sum c*dK = sum (J^T c)*dx.

    This is NOT a Hilbert adjoint; all complex cotangents stay unconjugated.
    Repeated appearances of the same leg/mean are accumulated independently.
    """
    cot = dict(cotangent)
    if 'retarded' in cot:
        cot['greater'] = cot.get('greater', 0) + theta*cot['retarded']
        cot['lesser'] = cot.get('lesser', 0) - theta*cot['retarded']
    out = {name: np.zeros_like(value, dtype=complex) for name, value in data.items()}
    for key, factor, equation, names in terms:
        if key not in cot:
            continue
        inputs, result = equation.split('->')
        labels = inputs.split(',')
        for slot, name in enumerate(names):
            others = [i for i in range(len(names)) if i != slot]
            dual = ','.join([result]+[labels[i] for i in others])+'->'+labels[slot]
            operands = [cot[key]]+[data[names[i]] for i in others]
            out[name] += factor*np.einsum(dual, *operands, optimize=True)
    return out


# Routing the two Wick pairings into the full ordered internal-q sum gives
# 1/2, not 1/4 or 1. This reduction requires a permutation-symmetric physical
# Phi and an UNTRUNCATED internal-q sum, not symmetric/asymmetric W matrices.
# The same polynomial engine supplies JVP and transpose with terms=ROUTED_TERMS.
ROUTED_TERMS = [
    ('mean_left', -.5, 'iab,ab->i', ('left', 'equal_left')),
    ('mean_right', -.5, 'jcd,cd->j', ('right', 'equal_right')),
    ('greater', -.5j, 'iab,jcd,ac,bd->ij', ('left','right','wg_first','wg_second')),
    ('lesser', -.5j, 'iab,jcd,ac,bd->ij', ('left','right','wl_first','wl_second')),
]
