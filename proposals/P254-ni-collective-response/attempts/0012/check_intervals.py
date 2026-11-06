"""Expose each pure arithmetic block with exact and independent optimizers."""
from itertools import product
import numpy as np
import sympy as sp
from scipy.optimize import minimize
from substrate_framework.material_force_sensitivity import (
    printed_interval, difference_box, box_norm_squared, vector_difference_squared,
    max_force_difference, printing_decision, scalar_difference_interval,
)
from substrate_framework.material_primary_boundary import force_rounding_upper


def rejects(call):
    try:
        call()
    except ValueError:
        return True
    return False


def run(ledger, context):
    q = sp.Rational
    ledger.check('decimal lexeme retains sign/exponent/trailing precision',
                 printed_interval('-1.230D-4') == (q('-0.0001230'), q('0.00000005')))
    ledger.check('half quantum changes with trailing zero, not decimal value',
                 printed_interval('1.20')[0] == printed_interval('1.2')[0] and
                 printed_interval('1.20')[1]*10 == printed_interval('1.2')[1])
    for token in ('NaN', 'inf', '1.0 junk', '', 1.0):
        ledger.check('malformed decimal rejected '+repr(token), rejects(lambda: printed_interval(token)))
    a, b = ['1.00','-2.0','0.000'], ['0.00','1.0','0.000']
    box = difference_box(a,b,q(2))
    ledger.check('difference includes BOTH independent printed radii and sign',
                 box == ((q('1.98'),q('2.02')),(q('-6.2'),q('-5.8')),(q('-.002'),q('.002'))))
    ledger.check('component omission and negative conversion rejected',
                 rejects(lambda: difference_box(a[:2],b,2)) and rejects(lambda: difference_box(a,b,-2)))
    ledger.check('full signed subtraction exposes add/absolute/component mutations',
                 vector_difference_squared(a,b,q(2)) == 40 and
                 vector_difference_squared(a,['0.00','-1.0','0.000'],q(2)) == 8)
    witnesses = [box, ((q(-2),q(1)),(q(3),q(4)),(q(-6),q(-5))),
                 ((q(-1),q(1)),)*3, ((q(0),q(0)),)*3]
    rows=[]
    for i,bounds in enumerate(witnesses):
        lo, hi = box_norm_squared(bounds)
        corners = [sum(x*x for x in point) for point in product(*bounds)]
        projected = [min(max(q(0),lower),upper) for lower,upper in bounds]
        ledger.check(f'box {i} farthest corner and projected origin exact extrema',
                     hi == max(corners) and lo == sum(x*x for x in projected))
        fbounds = np.asarray(bounds,float)
        fit=minimize(lambda x: float(np.dot(x,x)), fbounds.mean(axis=1),
                     jac=lambda x:2*x, bounds=fbounds,method='L-BFGS-B',
                     options={'ftol':1e-15,'gtol':1e-12,'maxiter':100})
        ledger.check(f'box {i} independent SciPy convex minimum succeeds',
                     fit.success and abs(fit.fun-float(lo)) < 2e-12)
        rows.append({'box':str(bounds),'lower_squared':str(lo),'upper_squared':str(hi),
                     'scipy_success':bool(fit.success),'scipy_fun':float(fit.fun)})
    ledger.check('lower-bound mutation nearest CORNER fails when box crosses origin',
                 min(sum(x*x for x in p) for p in product(*witnesses[2])) != box_norm_squared(witnesses[2])[0])
    ledger.check('upper-bound mutation center-only fails at exposed corner',
                 sum(((lo+hi)/2)**2 for lo,hi in box) < box_norm_squared(box)[1])
    gate=q('.003'); zero=[['0.000000000']*3]
    for token,decision in [('0.002999990','guaranteed_below'),('0.003000010','guaranteed_above'),
                            ('0.003000000','print_inconclusive')]:
        result=max_force_difference([[token,'0.000000000','0.000000000']],zero,1)
        ledger.check('near-gate exact witness '+token,
                     printing_decision(result['lower_squared'],result['upper_squared'],gate)==decision)
    skew=[['.00220000','.00220000','.00000000']]
    result=max_force_difference(skew,[['.00000000']*3],1)
    ledger.check('Euclidean witness defeats componentwise weaker gate',
                 result['lower_squared'] > gate**2 and max(map(q,skew[0])) < gate)
    two=[['.000000','.000000','.000000'],['.004000','.000000','.000000']]
    result=max_force_difference(two,[['.000000']*3]*2,1)
    ledger.check('last atom omission changes failing maximum to passing',
                 result['lower_squared']>gate**2 and max_force_difference(two[:1],[['.000000']*3],1)['upper_squared']<gate**2)
    ledger.check('wrong atom order changes observable',
                 max_force_difference(two,two,1)['center_squared']==0 and
                 max_force_difference(two,two[::-1],1)['center_squared']==q('.004')**2)
    ledger.check('scalar intervals expose sign and retained printing uncertainty',
                 scalar_difference_interval('-2.00','-1.00',2)==(q('-2.02'),q('-1.98')))
    # Retain the accepted old triangle helper, supplying the full difference radius.
    center=np.array([[(float(l)+float(h))/2 for l,h in box]])
    radius=max(float((h-l)/2) for l,h in box)
    triangle=float(force_rounding_upper(center,radius)[0])
    ledger.check('correct retained triangle encloses exact box upper norm',
                 triangle >= float(sp.sqrt(box_norm_squared(box)[1])))
    tiny=difference_box(['0.00000000']*3,['0.00000000']*3,1)
    correct=float(force_rounding_upper(np.zeros((1,3)),1e-8)[0])
    wrong=float(force_rounding_upper(np.zeros((1,3)),5e-9)[0])
    ledger.check('old radius mistake exposed: one-output half quantum misses difference corner',
                 wrong < float(sp.sqrt(box_norm_squared(tiny)[1])) <= correct*(1+1e-15))
    return {'optimizer_crosschecks':rows,'triangle_example_upper':triangle,
            'rounding_model':'closed +/- half last printed decimal quantum; independent boxes',
            'scope':'Exact extrema of Cartesian printing boxes, not a native numerical error bound.'}
