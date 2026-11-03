#!/usr/bin/env python3
"""Exact block identity, independent quadrature and full-row error witnesses."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'),str(HERE)]
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.linalg import expm
from gramian_krylov_bounds import gramian_envelope
from gramian_column_mean import bounded_source_columns
from gramian_residual_reference import two_level_residual_integral
from joint_column_mean import ColumnLimits
from substrate_framework.column_krylov_bounds import column_envelope
from substrate_framework.streamed_material_state import file_binding


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    t,k = sp.symbols('t k',real=True)
    a = sp.Matrix([[0,k],[0,0]])
    g = sp.Matrix([[2,1+sp.I],[1-sp.I,2]])
    block = (-a.H).row_join(g).col_join(sp.zeros(2).row_join(a))
    assert block**4 == sp.zeros(4)
    flow = sum((t**j*block**j/sp.factorial(j) for j in range(4)),sp.zeros(4))
    gramian = flow[2:4,2:4].H*flow[0:2,2:4]
    s = sp.symbols('s',real=True)
    independent = ((sp.eye(2)+s*a.H)*g*(sp.eye(2)+s*a)).applyfunc(lambda x:sp.integrate(x,(s,0,t)))
    assert (gramian-independent).applyfunc(sp.simplify) == sp.zeros(2)
    results = []
    projected = np.array([[.2+.03j,.4+.2j],[-.1+.1j,-.7-.04j]])
    residual = np.array([[.1+.3j,-.2j],[.4j,.3-.2j],[-.2,.1+.1j]])
    coefficient,initial,alpha = 1.3,.007,.2
    for time in (.37,-.41,3.):
        tau=abs(time); b=(-1j if time>=0 else 1j)*projected
        source=np.array([coefficient,0],complex)
        integral=quad(lambda u:np.exp(alpha*(tau-u))*np.linalg.norm(residual@expm(u*b)@source),
                      0,tau,epsabs=1e-12,epsrel=1e-12)[0]+np.exp(alpha*tau)*initial
        bound=gramian_envelope(projected,residual,coefficient,initial,time,ambient_growth=alpha)
        assert bound['absolute_error_envelope_estimate']+1e-12 >= integral
        results.append(dict(time=time,independent_Duhamel_integral=integral,bound=bound))
    rng=np.random.default_rng(225393)
    h=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5))
    q=rng.normal(size=(5,2))+1j*rng.normal(size=(5,2))
    small=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
    delta=(rng.normal(size=5)+1j*rng.normal(size=5))*.001
    x0=q[:,0]+delta
    alpha=float(np.linalg.norm((h-h.conj().T)/(2j),2))
    actual_errors=[]
    for time in (.08,-.09):
        actual=np.linalg.norm(expm(-1j*time*h)@x0-q@expm(-1j*time*small)[:,0])
        bound=gramian_envelope(small,h@q-q@small,1.,np.linalg.norm(delta),time,ambient_growth=alpha)
        assert actual <= bound['absolute_error_envelope_estimate']+1e-12
        actual_errors.append(dict(time=time,error=float(actual),bound=bound))
    h=np.array([[0,1e-8,0],[1e-8,40,1],[0,1,80]],complex)
    q=np.eye(3,dtype=complex)[:,:2]; small=q.conj().T@h@q; residual=h@q-q@small
    old=column_envelope(small,residual,1.,0.,30.,ambient_growth=0.)
    new=gramian_envelope(small,residual,1.,0.,30.,ambient_growth=0.,digits=50)
    high=gramian_envelope(small,residual,1.,0.,30.,ambient_growth=0.,digits=80)
    direct=two_level_residual_integral(1e-8,40.,30.)
    assert new['absolute_error_envelope_estimate'] >= direct-1e-13
    assert old['absolute_error_envelope_estimate'] > 1 and new['absolute_error_envelope_estimate'] < 2e-8
    assert abs(new['absolute_error_envelope_estimate']-high['absolute_error_envelope_estimate']) < 1e-20
    frame=np.eye(3,dtype=complex)[:,:1]
    evolved,receipt=bounded_source_columns(lambda x:h@x,frame,30.,
        limits=ColumnLimits(2,1024**2,2e-8),norm_bound=float(np.linalg.norm(h,2)),growth_bound=0.,identity='near-invariant-full-row')
    error=float(np.linalg.norm(evolved-expm(-1j*30*h)@frame))
    assert error <= receipt['envelope']['absolute_error_envelope_estimate']+1e-13
    constant=gramian_envelope(np.zeros((1,1)),np.array([[1j*1e-8]]),1.,0.,30.,ambient_growth=0.)
    assert constant['absolute_error_envelope_estimate'] >= 30e-8*(1-1e-14)
    growth=gramian_envelope(np.zeros((1,1)),np.zeros((1,1)),1.,1.,2.,ambient_growth=.3)
    assert growth['absolute_error_envelope_estimate'] >= np.exp(.6)*(1-1e-14)
    result=dict(schema='P254-source-aware-Gramian-oracle-v1',status='passed',
        exact_noncommuting_complex_block_identity='zero',quadrature=results,
        nonorthogonal_nonHermitian_full_row_errors=actual_errors,
        near_invariant=dict(old_source_blind_remainder=old,new=new,high_precision=high,
            independent_integral=direct,controller_error=error,controller=receipt),
        conjugation_and_sqrt_time_witness=constant,ambient_growth_witness=growth,
        source_bindings=[file_binding(path) for path in (Path(__file__),HERE/'gramian_krylov_bounds.py',
            HERE/'gramian_column_mean.py',ROOT/'src/substrate_framework/column_frame_krylov.py')],
        scope='Exact encoded block identity and audited finite matrix/error evidence; no Native or interval/physical-accuracy claim')
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(status=result['status'],controller_error=error,
        old_bound=old['absolute_error_envelope_estimate'],new_bound=new['absolute_error_envelope_estimate'])))


if __name__=='__main__':
    main()
