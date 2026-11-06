"""Exposing Fock-operator / exact-SymPy checks of the NEW interaction block.

Execute only by Root after handoff. Fock residual operators eta=[c,H_ei]
produce the lambda^2 electron selfenergy coefficient; direct connected force
products produce Pi. A full Fock double commutator independently checks the
short-time population injection coefficient. No fitted lifetime is tested.
Two explicit spin copies have distinct coherent Slater states/vertices; two
Cartesian oscillator coordinates use their own masses/preparation widths.
Finite oscillator truncation is exact for the tested vacuum two-point products.
SymPy explicitly sums indices and differentiates independent scalar variables,
without importing the production contraction table or derivative machinery.
"""
import argparse
from dataclasses import fields, replace
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm, block_diag
import sympy as sp

from electron_ion_global_functional import (
    GlobalPair, edge_operands, evaluate_global, jvp_global, transpose_global)
from electron_ion_functional import EdgePair, Kernels, evaluate, jvp, transpose


def annihilators(n):
    operators = []
    for a in range(n):
        c = np.zeros((2**n,2**n),complex)
        for bits in range(2**n):
            if bits & (1 << a):
                c[bits ^ (1 << a),bits] = (-1)**((bits & ((1 << a)-1)).bit_count())
        operators.append(c)
    return operators


def expect(ket, operator):
    return np.vdot(ket,operator@ket)


def fock_case(spin, t, s, n):
    cs = annihilators(4)
    hs = np.array([[.2+spin*.13,.31+.17j],[.31-.17j,-.4]])
    ht = np.array([[-.1,.19-.23j],[.19+.23j,.6+spin*.07]])
    h = block_diag(hs,ht)
    hf = sum(h[a,b]*cs[a].conj().T@cs[b] for a in range(4) for b in range(4))
    state = np.eye(16,dtype=complex)[:,0]
    vs = np.array([1,1j*(2+spin)],complex); vs /= np.linalg.norm(vs)
    vt = np.array([2-spin,1-.7j],complex); vt /= np.linalg.norm(vt)
    state = sum(vs[a]*cs[a].conj().T for a in range(2))@state
    state = sum(vt[a]*cs[a+2].conj().T for a in range(2))@state
    def evolved(time):
        u = expm(-1j*time*hf)
        return [u.conj().T@c@u for c in cs]
    ct, c_s = evolved(t), evolved(s)
    def green(left,right,offset):
        greater = np.array([[-1j*expect(state,left[a+offset]@right[b+offset].conj().T)
                             for b in range(2)] for a in range(2)])
        lesser = np.array([[1j*expect(state,right[b+offset].conj().T@left[a+offset])
                            for b in range(2)] for a in range(2)])
        return greater,lesser
    # The two exported vertices are read independently by the production API;
    # this tiny Hermitian Hamiltonian supplies its physically opposite legs.
    left = np.array([[[.4+.1j,.2-.3j],[-.1+.2j,.7]],
                     [[-.2j,.5+.1j],[.3-.4j,-.2]]])*(1+.3*spin)
    right = left.conj().transpose(0,2,1)
    cutoff = 3
    b = np.diag(np.sqrt(np.arange(1,cutoff)),1)
    ib = np.eye(cutoff)
    bs = [np.kron(b,ib),np.kron(ib,b)]
    bosket = np.eye(cutoff**2,dtype=complex)[:,0]
    mass, preparation, frequency = [2.,3.],[3.,5.],[1.,2.]
    xs = [(a+a.conj().T)/np.sqrt(2*m*w) for a,m,w in zip(bs,mass,preparation)]
    ps = [-1j*np.sqrt(m*w/2)*(a-a.conj().T) for a,m,w in zip(bs,mass,preparation)]
    # Symmetric canonical momentum shear gives a nonstationary correlated
    # Gaussian and nonsymmetric two-time mode matrices, without a thermal ansatz.
    shear = np.array([[.2,.7],[.7,-.1]])
    ps = [ps[i]+sum(shear[i,j]*xs[j] for j in range(2)) for i in range(2)]
    mix = np.array([[.8,.6],[-.6,.8]])
    def positions(time):
        raw = [np.cos(w*time)*x+np.sin(w*time)*p/(m*w)
               for x,p,m,w in zip(xs,ps,mass,frequency)]
        return [sum(mix[i,j]*raw[j] for j in range(2)) for i in range(2)]
    xt, x_s = positions(t), positions(s)
    dg = np.array([[-1j*expect(bosket,xt[i]@x_s[j]) for j in range(2)] for i in range(2)])
    dl = np.array([[-1j*expect(bosket,x_s[j]@xt[i]) for j in range(2)] for i in range(2)])
    gg,gl = green(ct,c_s,0); tg,tl = green(ct,c_s,2)
    bg,bl = green(c_s,ct,0)
    pair = EdgePair(left,right,left,right,gg,gl,tg,tl,bg,bl,dg,dl)
    joint = np.kron(bosket,state)
    def interaction(c, x):
        return sum(np.kron(x[i],sum(left[i,a,b]*c[a+2].conj().T@c[b]
                    +right[i,b,a]*c[b].conj().T@c[a+2]
                    for a in range(2) for b in range(2))) for i in range(2))/np.sqrt(n)
    def residual(c, x):
        v = interaction(c,x)
        annihilation = [np.kron(np.eye(cutoff**2),c[a+2]) for a in range(2)]
        return [a@v-v@a for a in annihilation]
    et, es = residual(ct,xt), residual(c_s,x_s)
    sg = np.array([[-1j*np.vdot(et[a].conj().T@joint,es[b].conj().T@joint)
                    for b in range(2)] for a in range(2)])
    sl = np.array([[1j*np.vdot(es[b]@joint,et[a]@joint)
                    for b in range(2)] for a in range(2)])
    def forces(c):
        return [-sum(right[i,b,a]*c[b].conj().T@c[a+2]
                     for a in range(2) for b in range(2))/np.sqrt(n) for i in range(2)]
    ft, fs = forces(ct),forces(c_s)
    pg = np.array([[-1j*(expect(state,ft[i]@fs[j].conj().T)
                         -expect(state,ft[i])*expect(state,fs[j].conj().T))
                    for j in range(2)] for i in range(2)])
    pl = np.array([[-1j*(expect(state,fs[j].conj().T@ft[i])
                         -expect(state,fs[j].conj().T)*expect(state,ft[i]))
                    for j in range(2)] for i in range(2)])
    reference = Kernels(sg,sl,pg,pl,float(t>=s)*(sg-sl),float(t>=s)*(pg-pl))
    # Full fermion-boson Fock weak-coupling population coefficient at time0.
    v = interaction(cs,positions(0.))
    number = np.kron(np.eye(cutoff**2),sum(c.conj().T@c for c in cs[2:]))
    comm = v@number-number@v
    acceleration = -expect(joint,v@comm-comm@v)
    return pair,reference,acceleration


def error(a,b):
    return float(np.linalg.norm(np.asarray(a)-np.asarray(b)))


def explicit_symbolic(pair, direction, dual, n, t, s):
    """Independent entrywise polynomial, directional derivative and transpose."""
    variables, values, directions, arrays = [],{}, {},{}
    for field in fields(pair):
        value = getattr(pair,field.name)
        symbols = np.asarray(sp.symbols(field.name+':'+str(value.size)),object).reshape(value.shape)
        arrays[field.name] = symbols
        for symbol,z,dz in zip(symbols.flat,value.flat,getattr(direction,field.name).flat):
            variables.append(symbol)
            values[symbol] = sp.Rational(str(z.real))+sp.I*sp.Rational(str(z.imag))
            directions[symbol] = sp.Rational(str(dz.real))+sp.I*sp.Rational(str(dz.imag))
    p = arrays; out = {}
    for component in ('greater','lesser'):
        out['sigma_'+component] = sp.Matrix(2,2,lambda a,d:
            sp.I/n*sum(p['left_t'][i,a,b]*p['gs_'+component][b,c]*p['right_s'][j,c,d]
                      *p['d_'+component][i,j] for i in range(2) for j in range(2)
                      for b in range(2) for c in range(2)))
        back = 'lesser' if component == 'greater' else 'greater'
        out['pi_'+component] = sp.Matrix(2,2,lambda i,j:
            -sp.I/n*sum(p['right_t'][i,b,a]*p['gt_'+component][a,d]*p['left_s'][j,d,c]
                       *p['back_'+back][c,b] for a in range(2) for b in range(2)
                       for c in range(2) for d in range(2)))
    theta = int(t>=s)
    out['sigma_retarded'] = theta*(out['sigma_greater']-out['sigma_lesser'])
    out['pi_retarded'] = theta*(out['pi_greater']-out['pi_lesser'])
    eps = sp.Symbol('eps',real=True)
    shift = {x: values[x]+eps*directions[x] for x in variables}
    tangent, numeric = {},{}
    objective = sp.S.Zero
    for name,matrix in out.items():
        numeric[name] = np.array(matrix.applyfunc(lambda z: complex(z.xreplace(values))),complex)
        tangent[name] = np.array(matrix.applyfunc(lambda z:
            complex(sp.diff(z.xreplace(shift),eps).subs(eps,0))),complex)
        weights = getattr(dual,name)
        objective += sum(matrix[a,b]*complex(weights[a,b]) for a in range(2) for b in range(2))
    gradient = {x: complex(sp.diff(objective,x).xreplace(values)) for x in variables}
    cotangent = EdgePair(**{name: np.array([gradient[x] for x in symbols.flat],complex).reshape(symbols.shape)
                            for name,symbols in arrays.items()})
    return Kernels(**numeric),Kernels(**tangent),cotangent


def bilinear(a,b):
    return sum(np.sum(getattr(a,f.name)*getattr(b,f.name)) for f in fields(a))


def run():
    n,t,s = 7,.37,.11
    pairs, refs, kernels, accelerations = [],[],[],[]
    for spin in range(2):
        pair,reference,acceleration = fock_case(spin,t,s,n)
        pairs.append(pair); refs.append(reference); accelerations.append(acceleration)
        kernels.append(evaluate(pair,n=n,t=t,s=s))
    defects = {f'spin{spin}_{f.name}': error(getattr(kernels[spin],f.name),getattr(refs[spin],f.name))
               for spin in range(2) for f in fields(Kernels)}
    pg = sum(k.pi_greater for k in kernels)
    defects['explicit_spin_sum'] = error(pg,sum(r.pi_greater for r in refs))
    zero_pair,_,acceleration = fock_case(0,0.,0.,n)
    zero_kernel = evaluate(zero_pair,n=n,t=0.,s=0.)
    injection = 2*np.trace(zero_kernel.sigma_lesser@zero_pair.gt_greater
                           -zero_kernel.sigma_greater@zero_pair.gt_lesser)
    defects['full_Fock_double_commutator_injection'] = error(injection,acceleration)
    reversed_pair,reversed_reference,_ = fock_case(0,s,t,n)
    reversed_kernel = evaluate(reversed_pair,n=n,t=s,s=t)
    defects['causal_earlier_time'] = error(reversed_kernel.pi_retarded,reversed_reference.pi_retarded)
    # Independent arbitrary complex operands break physical Hermiticity relations,
    # exposing a false conjugating dual that could pass on symmetric states.
    rng = np.random.default_rng(4107)
    def arbitrary(value):
        return (rng.integers(-3,4,value.shape)+1j*rng.integers(-3,4,value.shape))/5
    base = EdgePair(**{f.name: arbitrary(getattr(pairs[0],f.name)) for f in fields(EdgePair)})
    direction = EdgePair(**{f.name: arbitrary(getattr(base,f.name)) for f in fields(EdgePair)})
    sample = evaluate(base,n=n,t=t,s=s)
    dual = Kernels(**{f.name: arbitrary(getattr(sample,f.name)) for f in fields(Kernels)})
    exact,exact_j,exact_t = explicit_symbolic(base,direction,dual,n,t,s)
    actual_j = jvp(base,direction,n=n,t=t,s=s)
    actual_t = transpose(base,dual,n=n,t=t,s=s)
    defects['entrywise_functional'] = max(error(getattr(sample,f.name),getattr(exact,f.name)) for f in fields(Kernels))
    defects['SymPy_JVP'] = max(error(getattr(actual_j,f.name),getattr(exact_j,f.name)) for f in fields(Kernels))
    defects['SymPy_true_transpose'] = max(error(getattr(actual_t,f.name),getattr(exact_t,f.name)) for f in fields(EdgePair))
    defects['bilinear_identity'] = float(abs(bilinear(dual,actual_j)-bilinear(actual_t,direction)))
    # Unrestricted dense spatial/spin/k coherences: independently differentiated
    # entrywise polynomial, with all shared occurrences scatter-added by J.T.
    global_pair = GlobalPair(base.left_t,base.left_s,base.gs_greater,base.gs_lesser,
        base.back_greater,base.back_lesser,base.d_greater,base.d_lesser)
    global_direction = GlobalPair(direction.left_t,direction.left_s,direction.gs_greater,direction.gs_lesser,
        direction.back_greater,direction.back_lesser,direction.d_greater,direction.d_lesser)
    global_exact,global_exact_j,_ = explicit_symbolic(
        edge_operands(global_pair),edge_operands(global_direction),dual,1,t,s)
    global_actual = evaluate_global(global_pair,t=t,s=s)
    global_actual_j = jvp_global(global_pair,global_direction,t=t,s=s)
    global_actual_t = transpose_global(global_pair,dual,t=t,s=s)
    defects['global_entrywise_functional'] = max(error(getattr(global_actual,f.name),getattr(global_exact,f.name)) for f in fields(Kernels))
    defects['global_SymPy_JVP'] = max(error(getattr(global_actual_j,f.name),getattr(global_exact_j,f.name)) for f in fields(Kernels))
    defects['global_true_transpose'] = float(abs(bilinear(dual,global_exact_j)-bilinear(global_actual_t,global_direction)))
    pair,reference = pairs[0],refs[0]
    changed = lambda p: evaluate(p,n=n,t=t,s=s)
    mutations = dict(
        inconsistent_vertex_sign=error(changed(replace(pair,right_s=-pair.right_s)).sigma_greater,reference.sigma_greater),
        wrong_force_loop_sign=error(-reference.pi_greater,reference.pi_greater),
        wrong_vertex_conjugation=error(changed(replace(pair,right_s=pair.right_s.conj())).sigma_greater,reference.sigma_greater),
        wrong_reverse_ordered_leg=error(changed(replace(pair,back_lesser=pair.back_greater)).pi_greater,reference.pi_greater),
        wrong_mode_leg=error(changed(replace(pair,d_greater=pair.d_greater.T)).sigma_greater,reference.sigma_greater),
        dimensionless_oscillator=error(changed(replace(pair,d_greater=6*pair.d_greater)).sigma_greater,reference.sigma_greater),
        disconnected_mean_dyad=error(changed(replace(pair,d_greater=pair.d_greater-1j*np.outer([.3,.7],[.3,.7]))).sigma_greater,reference.sigma_greater),
        wrong_orbital_leg=error(changed(replace(pair,gs_greater=pair.gt_greater)).sigma_greater,reference.sigma_greater),
        wrong_spin2=error(2*kernels[0].pi_greater,sum(r.pi_greater for r in refs)),
        missing_N=error(evaluate(pair,n=1,t=t,s=s).sigma_greater,reference.sigma_greater),
        missing_lesser_injection=error(2*np.trace(-zero_kernel.sigma_greater@zero_pair.gt_lesser),acceleration),
        refilled_diagonal_state=error(changed(replace(pair,gs_lesser=np.diag(np.diag(pair.gs_lesser)))).sigma_lesser,reference.sigma_lesser),
        wrong_causal_order=error(kernels[0].pi_retarded,reversed_reference.pi_retarded),
        conjugating_dual=float(abs(bilinear(dual,actual_j)-sum(np.vdot(getattr(actual_t,f.name),getattr(direction,f.name)) for f in fields(EdgePair)))))
    predicates = {name: bool(value < 3e-10) for name,value in defects.items()}
    predicates.update({'exposes_'+name: bool(value > 1e-8) for name,value in mutations.items()})
    return dict(schema='p254-electron-ion-independent-oracle-v1',defects=defects,
        mutation_distances=mutations,predicates=predicates,passed=bool(all(predicates.values())),
        scope='lambda^2 composite/Fock force coefficients + double-commutator injection + symbolic polynomial derivatives, not full KB/material accuracy',
        sign_limit='simultaneous g->-g is invisible at second order; inconsistent single-leg/loop signs are exposed; mean-force sign is source-bound not inferred here')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(output=str(args.output),passed=result['passed'])))
    if not result['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
