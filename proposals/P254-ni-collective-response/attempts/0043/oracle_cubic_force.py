"""Executable independent oscillator/Fock and symbolic cubic-force oracle.

Written only. A fourth moment visits at most occupation 2 from vacuum; cutoff
3 is sufficient. We deliberately use cutoff 5 and compare 3, not inherit A41's
cutoff. Polynomial linear canonical evolution preserves this degree bound.
"""
import json
import numpy as np
import sympy as sp
from cubic_force_functional import evaluate, jvp, transpose


def fock_inputs(cutoff):
    a = np.diag(np.sqrt(np.arange(1, cutoff)), 1)
    ops = [np.kron(a, np.eye(cutoff)), np.kron(np.eye(cutoff), a)]
    ket = np.eye(cutoff**2, dtype=complex)[:, 0]
    x = [(o+o.conj().T)/np.sqrt(2*m*w) for o,m,w in zip(ops,[2.,3.],[3.,5.])]
    p = [-1j*np.sqrt(m*w/2)*(o-o.conj().T) for o,m,w in zip(ops,[2.,3.],[3.,5.])]
    p = [p[0]+.2*x[0]+.7*x[1], p[1]+.7*x[0]-.1*x[1]]
    xt = [x[0]+.31*p[0]+.13*x[1], x[1]-.19*p[1]+.23*x[0]]
    xs = [x[0]-.17*p[0]+.21*x[1], x[1]+.29*p[1]-.11*x[0]]
    expectation = lambda o: np.vdot(ket, o@ket)
    ordered = lambda l,r: np.array([[expectation(a@b) for b in r] for a in l])
    left = np.array([[[.4+.2j,.3-.1j],[.3-.1j,-.2j]],
                     [[-.2,.1+.4j],[.1+.4j,.7]]])
    right = np.array([[[.3-.1j,.2+.3j],[.2+.3j,-.6]],
                      [[.8j,-.4+.1j],[-.4+.1j,.2]]])
    ml = np.array([.17+.03j,-.21]); mr = np.array([-.13,.08-.07j])
    data = dict(left=left, right=right, mu_left=ml, mu_right=mr,
        equal_left=ordered(xt,xt), equal_right=ordered(xs,xs),
        wg=ordered(xt,xs), wl=ordered(xs,xt).T)
    def force(vertex, position, mean):
        shifted = [v+m*np.eye(cutoff**2) for v,m in zip(position,mean)]
        return [-sum(vertex[i,a,b]*shifted[a]@shifted[b]/2
                     for a in range(2) for b in range(2)) for i in range(2)]
    fl, fr = force(left,xt,ml), force(right,xs,mr)
    means = [np.array([expectation(o) for o in f]) for f in (fl,fr)]
    greater = -1j*(ordered(fl,fr)-np.outer(*means))
    lesser = -1j*(ordered(fr,fl).T-np.outer(*means))
    exact = dict(mean_left=means[0], mean_right=means[1], greater=greater,
                 lesser=lesser, retarded=greater-lesser)
    return data, exact


def main():
    errors, mutations = {}, {}
    data, direct = fock_inputs(5)
    low, low_direct = fock_inputs(3)
    actual = evaluate(data, 1)
    for name in actual:
        errors['fock_'+name] = float(np.max(np.abs(actual[name]-direct[name])))
        errors['cutoff_'+name] = float(np.max(np.abs(direct[name]-low_direct[name])))
    errors['causality'] = float(np.max(np.abs(evaluate(data,0)['retarded'])))
    # Independent symbolic differentiation of a noncommuting-leg polynomial:
    # the production differential is not used to construct these derivatives.
    z = sp.symbols('z0:8')
    poly = -sp.I*sp.Rational(1,4)*z[0]*z[1]*(z[2]*z[3]+z[4]*z[5])
    errors['symbolic_factor'] = float(abs(complex(sp.diff(poly,z[0],z[1],z[2],z[3])+sp.I/4)))
    rng = np.random.default_rng(4301)
    tangent = {k: rng.normal(size=v.shape)+1j*rng.normal(size=v.shape) for k,v in data.items()}
    cot = {k: rng.normal(size=v.shape)+1j*rng.normal(size=v.shape) for k,v in actual.items()}
    dv = jvp(data,tangent,1); jt = transpose(data,cot,1)
    lhs = sum(np.sum(cot[k]*dv[k]) for k in cot)
    rhs = sum(np.sum(jt[k]*tangent[k]) for k in data)
    errors['bilinear_transpose'] = float(abs(lhs-rhs))
    # Entrywise SymPy reference for ALL independent legs/vertices, including
    # asymmetric complex perturbations; finite differences are not the oracle.
    symbolic_jvp = symbolic_directional(data,tangent)
    for name in dv:
        errors['symbolic_jvp_'+name] = float(np.max(np.abs(dv[name]-symbolic_jvp[name])))
    mutations['statistical_factor'] = float(np.linalg.norm(2*actual['greater']-direct['greater']))
    disconnected = -1j*np.outer(direct['mean_left'],direct['mean_right'])
    mutations['missing_mean_subtraction'] = float(np.linalg.norm(disconnected))
    wrong = dict(data, wl=data['wg'])
    mutations['ordered_leg'] = float(np.linalg.norm(evaluate(wrong,1)['lesser']-direct['lesser']))
    mutations['hilbert_adjoint'] = float(abs(lhs-sum(np.sum(jt[k].conj()*tangent[k]) for k in data)))
    wrong = dict(data, right=data['right'].conj())
    mutations['conjugated_vertex'] = float(np.linalg.norm(evaluate(wrong,1)['greater']-direct['greater']))
    if max(errors.values()) > 3e-10 or min(mutations.values()) < 1e-8:
        raise AssertionError(dict(errors=errors, mutations=mutations))
    print(json.dumps(dict(schema='p254-a43-cubic-force-oracle-v1', errors=errors,
                          mutations=mutations, cutoff_degree_bound=4),indent=2))


def symbolic_directional(data, tangent):
    e = sp.Symbol('e', real=True)
    def value(name, *index):
        return complex(data[name][index])+e*complex(tangent[name][index])
    out = {name: np.zeros_like(v) for name,v in evaluate(data,1).items()}
    for side, vertex, equal, mu in [('left','left','equal_left','mu_left'),
                                   ('right','right','equal_right','mu_right')]:
        for i in range(2):
            expression = -sum(value(vertex,i,a,b)*(value(equal,a,b)+value(mu,a)*value(mu,b))/2
                              for a in range(2) for b in range(2))
            out['mean_'+side][i] = complex(sp.diff(expression,e).subs(e,0))
    for name,leg in [('greater','wg'),('lesser','wl')]:
        for i in range(2):
            for j in range(2):
                expression = 0
                for a in range(2):
                    for b in range(2):
                        for c in range(2):
                            for d in range(2):
                                w = lambda a,c: value(leg,a,c)
                                ml = lambda a: value('mu_left',a)
                                mr = lambda c: value('mu_right',c)
                                wick = (w(a,c)*w(b,d)+w(a,d)*w(b,c)+ml(a)*mr(c)*w(b,d)
                                        +ml(a)*mr(d)*w(b,c)+ml(b)*mr(c)*w(a,d)+ml(b)*mr(d)*w(a,c))
                                expression += -sp.I*value('left',i,a,b)*value('right',j,c,d)*wick/4
                out[name][i,j] = complex(sp.diff(expression,e).subs(e,0))
    out['retarded'] = out['greater']-out['lesser']
    return out


if __name__ == '__main__':
    main()
