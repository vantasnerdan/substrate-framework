"""One cubic fermion polynomial owns direct/exchange, JVP and bilinear transpose.

Sigma_ij(t,s)=1/2 A_ia;bc(t) A_de;jf(s) G_bd(t,s) G_ce(t,s) G_fa(s,t).
A=V-V(exchanged incoming). Greater uses >,>,< reverse; lesser uses <,<,>.
The four expanded terms retain DD+XX (direct) and -DX-XD (exchange).
Factor legs at both times are independent operands. No conjugation occurs in
this polynomial or its dual. Material factors alone own physical conjugation.
"""
from dataclasses import dataclass
import numpy as np

SIGNATURES = (('rib','rac','sdj','sef'), ('rib','rac','sdf','sej'),
              ('ric','rab','sdj','sef'), ('ric','rab','sdf','sej'))
SIGNS = (.5,-.5,-.5,.5)
_PATHS = {}


@dataclass(frozen=True)
class Term:
    signature: tuple
    factors: tuple
    coefficient: float
    channel: str


def terms_from_vertices(left, right):
    terms = []
    for index, (a,b) in enumerate(((0,0),(0,1),(1,0),(1,1))):
        if left[a] is not None and right[b] is not None:
            terms.append(Term(SIGNATURES[index], (*left[a],*right[b]), SIGNS[index],
                              'direct' if a == b else 'exchange'))
    return tuple(terms)


def contract(labels, operands, output):
    equation = ','.join(labels)+'->'+output
    key = (equation,tuple(x.shape for x in operands))
    if key not in _PATHS:
        _PATHS[key] = np.einsum_path(equation,*operands,optimize=('greedy',2_000_000))[0]
    return np.einsum(equation,*operands,optimize=_PATHS[key])


def ordered(terms, forward1, forward2, reverse):
    shape = (forward1.shape[0],forward1.shape[0])
    result = {name:np.zeros(shape,complex) for name in ('direct','exchange')}
    for term in terms:
        value = contract((*term.signature,'bd','ce','fa'),
                         (*term.factors,forward1,forward2,reverse),'ij')
        result[term.channel] += term.coefficient*value
    result['total'] = result['direct']+result['exchange']
    return result


def leg_keys(legs, order):
    forward = (order+'1',order+'2') if order+'1' in legs else (order,order)
    reverse = 'reverse_lesser' if order == 'greater' else 'reverse_greater'
    return (*forward,reverse)


def evaluate(terms, legs, *, retarded=True):
    """legs={greater,lesser,reverse_greater,reverse_lesser}, no inferred reverse."""
    greater = ordered(terms,*(legs[k] for k in leg_keys(legs,'greater')))
    lesser = ordered(terms,*(legs[k] for k in leg_keys(legs,'lesser')))
    return dict(greater=greater,lesser=lesser,
                retarded={k:(greater[k]-lesser[k]) if retarded else np.zeros_like(v)
                          for k,v in greater.items()})


def jvp(terms, legs, direction, *, factor_direction=None, retarded=True):
    """Directions include independent time-vertex factors and all four G legs.

    factor_direction has one four-factor tuple per expanded term. If physical
    factors are shared by terms/calls, caller gathers directions at their source.
    """
    if factor_direction is None:
        factor_direction = [(None,)*4 for _ in terms]
    result = {}
    for order in ('greater','lesser'):
        f,h,r = leg_keys(legs,order)
        channels = {name:np.zeros_like(legs[f]) for name in ('direct','exchange')}
        for term, df in zip(terms,factor_direction,strict=True):
            operands = (*term.factors,legs[f],legs[h],legs[r])
            changes = (*df,direction.get(f),direction.get(h),direction.get(r))
            labels = (*term.signature,'bd','ce','fa')
            for index, delta in enumerate(changes):
                if delta is not None:
                    changed = list(operands); changed[index] = delta
                    channels[term.channel] += term.coefficient*contract(labels,changed,'ij')
        channels['total'] = channels['direct']+channels['exchange']
        result[order] = channels
    result['retarded'] = {k:(result['greater'][k]-result['lesser'][k]) if retarded else np.zeros_like(v)
                          for k,v in result['greater'].items()}
    return result


def transpose(terms, legs, dual, *, retarded=True):
    """True sum(z*dx) transpose; never vdot, conjugation or reversed-leg merging.

    dual[order][channel] may address direct/exchange/total. Return all G cotangents
    plus per-term four independent factor cotangents. Shared physical operands
    require scatter-add by caller; each occurrence is returned rather than lost.
    """
    gl = {k:np.zeros_like(v) for k,v in legs.items()}
    gf = [[np.zeros_like(x) for x in term.factors] for term in terms]
    for order, sign in (('greater',1),('lesser',-1)):
        f,h,r = leg_keys(legs,order)
        for ti, term in enumerate(terms):
            z = np.zeros_like(legs[f])
            for channel in (term.channel,'total'):
                z += dual.get(order,{}).get(channel,0)
                if retarded:
                    z += sign*dual.get('retarded',{}).get(channel,0)
            labels = (*term.signature,'bd','ce','fa')
            operands = (*term.factors,legs[f],legs[h],legs[r])
            for index in range(7):
                keep = [j for j in range(7) if j != index]
                value = term.coefficient*contract([labels[j] for j in keep]+['ij'],
                    [operands[j] for j in keep]+[z],labels[index])
                if index < 4:
                    gf[ti][index] += value
                else:
                    gl[(f,h,r)[index-4]] += value
    return dict(legs=gl,factors=gf)
