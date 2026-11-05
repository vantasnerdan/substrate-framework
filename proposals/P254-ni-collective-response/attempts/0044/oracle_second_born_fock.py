"""Explicit CAR/Fock connected force and weak-coupling injection oracle."""
import itertools
import numpy as np
from scipy.linalg import expm
from second_born_functional import terms_from_vertices, evaluate


def annihilators(n):
    out = []
    for mode in range(n):
        c = np.zeros((2**n,2**n),complex)
        for ket in range(2**n):
            if ket & (1<<mode):
                c[ket^(1<<mode),ket] = (-1)**((ket&((1<<mode)-1)).bit_count())
        out.append(c)
    return out


def run():
    n = 4
    c = annihilators(n)
    dag = [x.conj().T for x in c]
    rng = np.random.default_rng(4403)
    # Explicit spin support: first two up, last two down; no implicit spin2.
    spin = np.array([0,0,1,1])
    x = rng.normal(size=(3,n,n))+1j*rng.normal(size=(3,n,n))
    x = (x+x.conj().swapaxes(-1,-2))/2
    x *= (spin[:,None] == spin[None,:])[None]
    weight = np.array([.23,.37,.19])
    v = np.einsum('r,rik,rjl->ijkl',weight,x,x)
    a = v-v.swapaxes(2,3)
    h = np.zeros_like(c[0])
    for i,j,k,l in itertools.product(range(n),repeat=4):
        h += .5*v[i,j,k,l]*dag[i]@dag[j]@c[l]@c[k]
    occ = np.array([.17,.38,.62,.79])
    probability = np.array([np.prod([occ[i] if ket&(1<<i) else 1-occ[i] for i in range(n)]) for ket in range(2**n)])
    state = np.diag(probability)
    trace = lambda op:np.trace(state@op)
    hf = np.einsum('iaka,a->ik',a,occ)
    e = trace(h)
    mean = sum((hf[i,k]*dag[i]@c[k] for i,k in itertools.product(range(n),repeat=2)),np.zeros_like(h))
    residual = h-mean+e*np.eye(2**n)
    force = [ci@residual-residual@ci for ci in c]
    normal = []
    for i in range(n):
        cubic = sum((.5*a[i,j,k,l]*dag[j]@c[l]@c[k]
                     for j,k,l in itertools.product(range(n),repeat=3)),np.zeros_like(h))
        normal.append(cubic-sum((hf[i,k]*c[k] for k in range(n)),np.zeros_like(h)))
    assert abs(trace(residual))<2e-12
    assert all(np.allclose(f,q,atol=2e-12) for f,q in zip(force,normal,strict=True))
    assert max(abs(trace(f@dj+dj@f)) for f in force for dj in dag)<2e-12
    energy = np.array([-.43,.13,.37,.71])
    h0 = sum((energy[i]*dag[i]@c[i] for i in range(n)),np.zeros_like(h))
    t,s = .41,.09
    ut,us = expm(1j*t*h0),expm(1j*s*h0)
    ft = [ut@f@ut.conj().T for f in force]
    fsdag = [us@f.conj().T@us.conj().T for f in force]
    fockg = np.array([[-1j*trace(f@g) for g in fsdag] for f in ft])
    fockl = np.array([[1j*trace(g@f) for g in fsdag] for f in ft])
    phase = np.exp(-1j*(t-s)*energy)
    gg,gl = np.diag(-1j*(1-occ)*phase),np.diag(1j*occ*phase)
    bg,bl = np.diag(-1j*(1-occ)/phase),np.diag(1j*occ/phase)
    factors = (weight[:,None,None]*x,x)
    terms = terms_from_vertices((factors,factors),(factors,factors))
    legs = dict(greater=gg,lesser=gl,reverse_greater=bg,reverse_lesser=bl)
    answer = evaluate(terms,legs)
    assert np.allclose(answer['greater']['total'],fockg,atol=2e-11,rtol=2e-11)
    assert np.allclose(answer['lesser']['total'],fockl,atol=2e-11,rtol=2e-11)
    assert np.linalg.norm(answer['greater']['exchange'])>1e-8
    # Independent weak-coupling interacting Fock injection, no Dyson solver.
    lam,dt = 1e-3,.2
    ci = c[0]
    comm = ci@residual-residual@ci
    double = comm@residual-residual@comm
    target = .5j*dt**2*trace(double@dag[0])
    zero = -1j*trace(ci@dag[0])
    values = []
    for coupling in (lam,-lam):
        u = expm(-1j*dt*coupling*residual)
        values.append(-1j*trace(u.conj().T@ci@u@dag[0]))
    injection = (values[0]+values[1]-2*zero)/(2*lam**2)
    assert abs(injection-target)<1e-6*max(1,abs(target))
    assert np.linalg.norm(.5*answer['greater']['total']-fockg)>1e-8
    assert np.linalg.norm(2*answer['greater']['total']-fockg)>1e-8
    assert np.linalg.norm(answer['greater']['direct']-fockg)>1e-8
    fresh_mean_force = [ci@h-h@ci for ci in c]
    assert max(np.linalg.norm(f-q) for f,q in zip(fresh_mean_force,force,strict=True))>1e-8
    opposite = fockg[:2,2:]
    assert np.linalg.norm(opposite)<2e-12
    return dict(schema='p254-a44-fock-oracle-v1',force_defect=float(max(np.linalg.norm(f-q) for f,q in zip(force,normal,strict=True))),
        greater_defect=float(np.linalg.norm(answer['greater']['total']-fockg)),
        weak_coupling_injection_defect=float(abs(injection-target)),
        mutations_exposed=['missing-antisymmetric-half','spin2','direct-only','fresh-HF-mean'],
        scope='finite explicit-spin Fock connected normal-ordered residual, independently of material applicability')
