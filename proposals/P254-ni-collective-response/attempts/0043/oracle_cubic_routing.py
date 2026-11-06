"""Literal finite-cell Fock forces versus full-q routed cubic functional.

Independent finite-coordinate SymPy scalar establishes 1/sqrt(N) in Phi and
both Wick pairings' routed 1/2 factor. No production Ewald formula is reused.
"""
from itertools import product
import json
import numpy as np
import sympy as sp
from cubic_force_functional import evaluate,jvp,transpose,ROUTED_TERMS


def main():
    n = 3; cutoff = 3; size = cutoff**n; eye = np.eye(cutoff)
    a = np.diag(np.sqrt(np.arange(1,cutoff)),1)
    ops = []
    for site in range(n):
        factors = [a if j == site else eye for j in range(n)]
        ops.append(np.kron(np.kron(*factors[:2]),factors[2]))
    ket = np.eye(size,dtype=complex)[:,0]
    x = [(o+o.conj().T)/np.sqrt(6) for o in ops]
    p = [-1j*np.sqrt(6)/2*(o-o.conj().T) for o in ops]
    xt = [x[r]+.21*p[r]+.13*x[(r+1)%n]+.13*x[(r-1)%n] for r in range(n)]
    xs = [x[r]-.17*p[r]+.09*p[(r+1)%n]+.09*p[(r-1)%n] for r in range(n)]
    expect = lambda o: np.vdot(ket,o@ket)
    fields = sp.symbols('x0:3')
    scalar = sum((fields[r]+sp.Rational(3,10)*fields[(r+1)%n])**3/6 for r in range(n))
    phi = np.array([float(sp.diff(scalar,fields[r],fields[s],fields[t]))
                    for r,s,t in product(range(n),repeat=3)]).reshape(n,n,n)
    transform = np.exp(2j*np.pi*np.outer(np.arange(n),np.arange(n))/n)/np.sqrt(n)
    # Direct finite-coordinate tensor transform, not presumed N prefactors.
    vertex = np.einsum('rst,rq,sk,tl->qkl',phi,transform,transform,transform)
    forces = lambda coordinates: [-sum(phi[r,s,t]*coordinates[s]@coordinates[t]/2
                                      for s in range(n) for t in range(n)) for r in range(n)]
    ft,fs = forces(xt),forces(xs)
    xqt = [sum(transform[r,q].conjugate()*xt[r] for r in range(n)) for q in range(n)]
    xqs = [sum(transform[r,q].conjugate()*xs[r] for r in range(n)) for q in range(n)]
    errors = {}; wrong = []
    for q in range(n):
        direct_t = sum(transform[r,q].conjugate()*ft[r] for r in range(n))
        direct_s = sum(transform[r,(-q)%n].conjugate()*fs[r] for r in range(n))
        direct_mean = expect(direct_t)
        direct_g = -1j*(expect(direct_t@direct_s)-direct_mean*expect(direct_s))
        direct_l = -1j*(expect(direct_s@direct_t)-direct_mean*expect(direct_s))
        accumulated = dict(mean_left=0j,greater=0j,lesser=0j,retarded=0j)
        for k in range(n):
            r = (q-k)%n
            data = dict(left=np.array([[[vertex[-q%n,k,r]]]]),
                right=np.array([[[vertex[q,-k%n,-r%n]]]]),
                equal_left=np.array([[expect(xqt[k]@xqt[-k%n])]]),
                equal_right=np.array([[expect(xqs[-k%n]@xqs[k])]]),
                wg_first=np.array([[expect(xqt[k]@xqs[-k%n])]]),
                wg_second=np.array([[expect(xqt[r]@xqs[-r%n])]]),
                wl_first=np.array([[expect(xqs[-k%n]@xqt[k])]]),
                wl_second=np.array([[expect(xqs[-r%n]@xqt[r])]]))
            result = evaluate(data,1,ROUTED_TERMS)
            for name in accumulated:
                if name != 'mean_left' or q == 0:
                    accumulated[name] += result[name].item()
            errors[f'transpose_{q}_{k}'] = duality(data)
        for name,reference in [('mean_left',direct_mean),('greater',direct_g),
                               ('lesser',direct_l),('retarded',direct_g-direct_l)]:
            errors[f'fock_q{q}_{name}'] = float(abs(accumulated[name]-reference))
        wrong.append(abs(accumulated['greater']*2-direct_g))
    # Onsite scalar sum x_R^3/6 has literal Phi(q,k,l)=delta/sqrt(N).
    onsite = np.einsum('rq,rk,rl->qkl',transform,transform,transform)
    for q,k,l in product(range(n),repeat=3):
        expected = 1/np.sqrt(n) if (q+k+l)%n == 0 else 0.
        errors[f'canonical_{q}{k}{l}'] = float(abs(onsite[q,k,l]-expected))
    mutations = dict(statistical_factor=float(max(wrong)),missing_N=float(abs(1-1/np.sqrt(n))))
    if max(errors.values()) > 3e-10 or min(mutations.values()) < 1e-8:
        raise AssertionError(dict(errors=errors,mutations=mutations))
    print(json.dumps(dict(schema='p254-a43-cubic-routing-oracle-v1',errors=errors,mutations=mutations),indent=2))


def duality(data):
    rng = np.random.default_rng(4354)
    tangent = {k:rng.normal(size=v.shape)+1j*rng.normal(size=v.shape) for k,v in data.items()}
    value = evaluate(data,1,ROUTED_TERMS)
    cot = {k:rng.normal(size=v.shape)+1j*rng.normal(size=v.shape) for k,v in value.items()}
    dx = jvp(data,tangent,1,ROUTED_TERMS)
    jt = transpose(data,cot,1,ROUTED_TERMS)
    return float(abs(sum(np.sum(cot[k]*dx[k]) for k in cot)-sum(np.sum(jt[k]*tangent[k]) for k in data)))


if __name__ == '__main__':
    main()
