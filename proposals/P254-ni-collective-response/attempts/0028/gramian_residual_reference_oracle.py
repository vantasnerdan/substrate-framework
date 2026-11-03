#!/usr/bin/env python3
"""Repair the weak oscillatory quadrature reference with an exact two-level identity."""
import argparse
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
sys.path[:0]=[str(ROOT/'src'),str(HERE)]
import numpy as np
import sympy as sp
from gramian_residual_reference import two_level_residual_integral
from gramian_krylov_bounds import gramian_envelope
from substrate_framework.streamed_material_state import file_binding


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    d=sp.symbols('d',positive=True)
    w,t=sp.symbols('w t',real=True)
    h=sp.Matrix([[0,d],[d,w]])
    shifted=h-w*sp.eye(2)/2
    frequency=sp.sqrt(w*w/4+d*d)
    assert (shifted**2-frequency**2*sp.eye(2)).applyfunc(sp.simplify)==sp.zeros(2)
    flow=sp.exp(-sp.I*w*t/2)*(sp.cos(frequency*t)*sp.eye(2)-sp.I*sp.sin(frequency*t)/frequency*shifted)
    assert (flow.diff(t)+sp.I*h*flow).applyfunc(sp.simplify)==sp.zeros(2)
    assert flow.subs(t,0)==sp.eye(2)
    u,r=sp.symbols('u r',real=True)
    assert sp.calculus.util.function_range(sp.sin(u),u,sp.Interval(0,sp.pi))==sp.Interval(0,1)
    assert sp.integrate(sp.sin(u),(u,0,sp.pi))==2
    assert sp.simplify(sp.integrate(sp.sin(u),(u,0,r))-(1-sp.cos(r)))==0
    direct=two_level_residual_integral(1e-8,40.,30.)
    h=np.array([[0,1e-8],[1e-8,40]],complex)
    residual=np.array([[0,1]],complex)
    bound=gramian_envelope(h,residual,1.,0.,30.,ambient_growth=0.)
    assert direct>0 and bound['absolute_error_envelope_estimate']>=direct
    result=dict(schema='P254-exact-weak-residual-reference-v1',status='passed',
        exact_matrix_ODE_and_initial_condition='zero',exact_period_primitive='zero',
        independent_integral=direct,gramian_bound=bound,
        source_bindings=[file_binding(p) for p in (Path(__file__),HERE/'gramian_residual_reference.py',HERE/'gramian_krylov_bounds.py')],
        scope='Exact encoded two-level source/period identity and high-precision scalar evaluation; replaces warning-bearing quadrature, not a Native claim')
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(result))


if __name__=='__main__':
    main()
