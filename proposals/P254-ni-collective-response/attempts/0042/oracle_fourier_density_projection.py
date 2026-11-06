"""Changed-representation oracle: literal grid sums, actual FFT and characters.

No old density/Fock/trajectory suite is replayed. This checks only the exact
Fourier-index projection and its independently applied opposite leg.
"""
from pathlib import Path
import argparse
import json
import sys
from types import SimpleNamespace
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[3]/'src'))
import numpy as np
import sympy as sp
from substrate_framework.qe_material_grid import PhysicalRhoFFT
from density_projection import BorrowedBathFrame,FourierDensityProjector


def character_identity():
    z=sp.symbols('z')
    for n in (1,3,4,5,7):
        polynomial=sp.cyclotomic_poly(n,z)
        for difference in range(n):
            expression=sum(z**((difference*r)%n) for r in range(n))
            expression-=n if difference==0 else 0
            if sp.rem(expression,polynomial,z)!=0:
                raise AssertionError('finite Fourier character identity failed')
    return True


def coordinates(shape):
    axes=np.meshgrid(*(np.arange(n)/n for n in shape),indexing='ij')
    return np.stack([axis.ravel(order='F') for axis in axes])


def literal(source,target,bs,bt,g,fold,n):
    x=coordinates(source.physical)
    ws=np.exp(2j*np.pi*x.T@source.miller)@bs
    wt=np.exp(2j*np.pi*x.T@target.miller)@bt
    multiplier=np.exp(2j*np.pi*(fold-g)@x)/np.sqrt(n)/np.sqrt(source.omega)
    forward=ws.conj().T@(multiplier[:,None]*wt)/source.nnr
    adjoint=wt.conj().T@(multiplier.conj()[:,None]*ws)/source.nnr
    return forward,adjoint


def run():
    rng=np.random.default_rng(4255)
    random=lambda shape:(rng.normal(size=shape)+1j*rng.normal(size=shape))*.2
    shape=(5,3,1);volume=7.3;n=9
    sm=np.array([[0,1,-2,6],[0,1,0,-1],[0,0,0,0]])
    tm=np.array([[0,-1,2,4,1],[0,1,-1,0,-1],[0,0,0,0,0]])
    source=PhysicalRhoFFT.create(shape,sm,np.eye(3),volume)
    target=PhysicalRhoFFT.create(shape,tm,np.eye(3),volume)
    ps,pt=random((4,2)),random((5,3))
    si,ti=random((2,2)),random((3,3))
    es,et=random((4,2)),random((5,1))
    sf,tf=BorrowedBathFrame(ps,si,es),BorrowedBathFrame(pt,ti,et)
    bs,bt=np.concatenate((ps@si,es),axis=1),np.concatenate((pt@ti,et),axis=1)
    fold=np.array([1,-1,0])
    bridge=SimpleNamespace(source_fft=source,target_fft=target,
                           edge=SimpleNamespace(reciprocal_fold=fold))
    projector=FourierDensityProjector(bridge,sf,tf,ncells=n,block_size=128)
    narrow=FourierDensityProjector(bridge,sf,tf,ncells=n,block_size=1)
    defects={};mutations={}
    for index,g in enumerate((np.array([0,0,0]),np.array([1,0,0]),
                              np.array([-1,1,0]),np.array([6,-2,0]))):
        actual,adjoint=projector.project(g)
        expected,expected_adjoint=literal(source,target,bs,bt,g,fold,n)
        small,small_adjoint=narrow.project(g)
        x=coordinates(shape)
        multiplier=np.exp(2j*np.pi*(fold-g)@x)/np.sqrt(n)/np.sqrt(volume)
        fft_reference=bs.conj().T@source.to_g(
            multiplier[:,None]*target.from_g(bt))
        fft_adjoint=bt.conj().T@target.to_g(
            multiplier.conj()[:,None]*source.from_g(bs))
        defects['literal_forward_'+str(index)]=float(np.max(abs(actual-expected)))
        defects['literal_adjoint_'+str(index)]=float(np.max(abs(adjoint-expected_adjoint)))
        defects['actual_FFT_forward_'+str(index)]=float(np.max(abs(actual-fft_reference)))
        defects['actual_FFT_adjoint_'+str(index)]=float(np.max(abs(adjoint-fft_adjoint)))
        defects['column_partition_'+str(index)]=float(max(np.max(abs(actual-small)),
                                                         np.max(abs(adjoint-small_adjoint))))
    g=np.array([-1,1,0]);actual,adjoint=projector.project(g)
    wrong,_=literal(source,target,bs,bt,-g,fold,n)
    mutations['wrong_profile_sign']=float(np.linalg.norm(actual-wrong))
    wrong,_=literal(source,target,bs,bt,g,np.zeros(3,int),n)
    mutations['omitted_reciprocal_fold']=float(np.linalg.norm(actual-wrong))
    mutations['missing_volume']=float(np.linalg.norm(actual*np.sqrt(volume)-actual))
    mutations['missing_rootN']=float(np.linalg.norm(actual*np.sqrt(n)-actual))
    mutations['transpose_instead_of_dagger']=float(np.linalg.norm(adjoint-actual.T))
    # Real q=0 physical interface and an empty extra subspace, not a fake action.
    gamma=SimpleNamespace(factors=SimpleNamespace(fft=source))
    empty=BorrowedBathFrame(ps,si,np.empty((4,0),complex))
    gp=FourierDensityProjector(gamma,empty,empty,ncells=n,block_size=128)
    ga,gd=gp.project(np.array([1,0,0]))
    ge,gde=literal(source,source,ps@si,ps@si,np.array([1,0,0]),np.zeros(3,int),n)
    defects['Gamma_empty_extra_forward']=float(np.max(abs(ga-ge)))
    defects['Gamma_empty_extra_adjoint']=float(np.max(abs(gd-gde)))
    predicates={name:bool(value<2e-11) for name,value in defects.items()}
    predicates.update({name+'_exposed':bool(value>1e-9) for name,value in mutations.items()})
    predicates['symbolic_character_identity']=character_identity()
    return dict(schema='P254-A42-exact-Fourier-projection-oracle-v1',
        defects=defects,mutation_distances=mutations,predicates=predicates,
        passed=bool(all(predicates.values())),
        scope='same finite physical FFT grid/profile/fold/full-PW frame; no new continuum accuracy or material-rate claim',
        unchanged='No state, physical feature, cutoff, N, volume or backward-filter approximation changed')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=run()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(output=str(args.output),passed=result['passed'])))
    raise SystemExit(0 if result['passed'] else 1)
