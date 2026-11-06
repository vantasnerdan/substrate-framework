#!/usr/bin/env python3
"""Independent symbolic, literal spherical/radial, finite-operator mutation oracles.

Run explicitly; no material/DFT execution. Output is evidence only if assertions
complete. The angular integral is independent of production solid harmonics
and Bessel identities; it differentiates exp(-i r p.n) under finite quadrature.
"""
from pathlib import Path
import json
import sys
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from substrate_framework.qe_material_field_derivatives import (
    UPFRadial, HermiteRadial, centered_derivatives, product_actions, solid_harmonics)


def symbolic():
    x, y, z, t, u = sp.symbols('x y z t u', real=True)
    e = sp.Matrix([sp.Rational(2, 3), -sp.Rational(1, 5), sp.Rational(3, 7)])
    f = sp.Matrix([-sp.Rational(1, 2), sp.Rational(4, 5), sp.Rational(1, 3)])
    p = sp.Matrix([x, y, z])
    q = p+t*e+u*f
    radial = sp.Function('F')
    polys = [1, z, -x, -y, (2*z*z-x*x-y*y)/2, -x*z, -y*z, x*x-y*y, x*y]
    for poly in polys:
        S = sp.sympify(poly)
        Se = sum(sp.diff(S, p[i])*e[i] for i in range(3))
        Sf = sum(sp.diff(S, p[i])*f[i] for i in range(3))
        Sef = sum(sp.diff(Se, p[i])*f[i] for i in range(3))
        s = sp.symbols('s', real=True)
        F, Fs, Fss = radial(s), sp.diff(radial(s), s), sp.diff(radial(s), s, 2)
        formula = Sef*F+Se*Fs*2*p.dot(f)+Sf*Fs*2*p.dot(e)+S*(Fss*4*p.dot(e)*p.dot(f)+Fs*2*e.dot(f))
        actual = sp.diff(S.subs(dict(zip(p, q)), simultaneous=True)*radial(q.dot(q)), t, u).subs({t:0, u:0})
        assert sp.simplify(actual-formula.subs(s, p.dot(p))) == 0
    # Literal forward cubic obstruction: even smooth q^4 acquires a nonzero
    # radial origin slope; adjoining intervals have unequal knot slopes.
    q = sp.symbols('q', real=True)
    first = sp.interpolate([(j, j**4) for j in range(4)], q)
    following = sp.interpolate([(j, j**4) for j in range(1, 5)], q)
    assert sp.diff(first,q).subs(q,0) == 6
    assert sp.diff(following-first,q).subs(q,1) == 8
    h = sp.symbols('h', positive=True)
    origin_contact = 2*(first.subs(q,h)-first.subs(q,0))/h**2
    assert sp.limit(h*origin_contact,h,0,dir='+') == 12
    # General complex projector legs, independent formal matrix product rule.
    B = sp.Function('B')(t); D = sp.Symbol('D', commutative=False)
    A = sp.Function('A')(t)
    assert sp.expand(sp.diff(B*D*A, t, 2)-sp.diff(B,t,2)*D*A-2*sp.diff(B,t)*D*sp.diff(A,t)-B*D*sp.diff(A,t,2)) == 0
    return {'solid_radial_mixed_identities': 9, 'product_rule_identity': True}


def fixture():
    r = np.array([0., .2, .55, .9, 1.3])
    beta = np.column_stack([np.exp(-r*r)*(1+.2*r), r*np.exp(-r), r*r*np.exp(-r), r*(1+.3*r)*np.exp(-r)])
    ls = np.array([0, 1, 2, 1])
    radial = np.repeat(np.arange(4), 2*ls+1)
    ell = ls[radial]
    lm = np.concatenate([np.arange(l*l, (l+1)**2) for l in ls])
    upf = SimpleNamespace(r=r, rab=np.array([.1,.16,.22,.27,.31]), beta_mesh=5,
                         beta_r=beta, angular_l=ls, nbeta=4, channel_radial=radial,
                         channel_l=ell, channel_lm=lm)
    return SimpleNamespace(upf=upf, omega_bohr3=17.)


def literal_integral(tables, p, e, f):
    # Literal angular functions on the unit sphere, not production polynomial API.
    mu, w = leggauss(32)
    phi = np.arange(64)*2*np.pi/64
    mu, phi = np.meshgrid(mu, phi, indexing='ij')
    nx, ny, nz = (np.sqrt(1-mu*mu)*np.cos(phi), np.sqrt(1-mu*mu)*np.sin(phi), mu)
    n = np.column_stack([nx.ravel(), ny.ravel(), nz.ravel()])
    c0, c1, c2 = 1/np.sqrt(4*np.pi), np.sqrt(3/(4*np.pi)), np.sqrt(5/(4*np.pi))
    Y = np.column_stack([np.full(n.shape[0], c0), c1*n[:,2], -c1*n[:,0], -c1*n[:,1],
                        c2*(3*n[:,2]**2-1)/2, -c2*np.sqrt(3)*n[:,0]*n[:,2],
                        -c2*np.sqrt(3)*n[:,1]*n[:,2], c2*np.sqrt(3)*(n[:,0]**2-n[:,1]**2)/2,
                        c2*np.sqrt(3)*n[:,0]*n[:,1]])
    aw = np.repeat(w, 64)*2*np.pi/64
    upf = tables.upf
    rw = upf.rab*np.array([1.,4.,2.,4.,1.])/3  # literal original odd Simpson
    out = [np.zeros((len(p),len(upf.channel_radial)),complex) for _ in range(4)]
    for j, r in enumerate(upf.r):
        phase = np.exp(-1j*r*(p@n.T))
        legs = [np.ones(len(n)), -1j*r*(n@e), -1j*r*(n@f), -r*r*(n@e)*(n@f)]
        for k, leg in enumerate(legs):
            angular = (phase*(aw*leg))@Y[:,upf.channel_lm]
            out[k] += angular*(rw[j]*r*upf.beta_r[j,upf.channel_radial]/np.sqrt(tables.omega_bohr3))
    return out


def run():
    result = symbolic()
    tables = fixture(); exact = UPFRadial(tables)
    p = np.array([[0.,0.,0.],[.2,-.3,.4],[1.1,.7,-.6],[-.31,0.,.5]])
    e = np.array([.7,-.2,.4]); f = np.array([-.3,.8,.1])
    actual = centered_derivatives(exact,p,e,f)
    literal = literal_integral(tables,p,e,f)
    errors = [float(np.max(np.abs(a-b))) for a,b in zip(actual,literal)]
    assert max(errors) < 2e-12, errors
    # Non-diagonal radial-channel mixing remains in D; use complex input columns.
    rng = np.random.default_rng(22528)
    size = actual[0].shape[1]
    d = np.diag(np.linspace(.4,1.2,size))
    for a,b in [(1,9),(2,10),(3,11)]: d[a,b]=d[b,a]=.23
    v = rng.normal(size=(len(p),2))+1j*rng.normal(size=(len(p),2))
    b, be, bf, bef = actual
    lin, con = product_actions(*actual,d,v)
    dense_l = be@d@b.conj().T+b@d@be.conj().T
    dense_c = bef@d@b.conj().T+be@d@bf.conj().T+bf@d@be.conj().T+b@d@bef.conj().T
    assert np.linalg.norm(lin-dense_l@v)<1e-12
    assert np.linalg.norm(con-dense_c@v)<1e-12
    mutations = {
        'omit_right_current_leg': float(np.linalg.norm(b@d@be.conj().T@v)),
        'omit_mixed_cross_terms': float(np.linalg.norm((be@d@bf.conj().T+bf@d@be.conj().T)@v)),
        'reverse_electron_current_sign': float(np.linalg.norm(lin)),
        'zero_origin_first_derivative': float(np.linalg.norm(literal[1][0])),
        'zero_origin_second_derivative': float(np.linalg.norm(literal[3][0])),
        'discard_offdiagonal_D': float(np.linalg.norm(con-product_actions(*actual,np.diag(np.diag(d)),v)[1])),
    }
    S, Se, Sf, Sef = solid_harmonics(p,e,f)
    rr, rs, rss = exact.evaluate(np.sum(p*p,axis=1))
    angular_leg = Se[:,tables.upf.channel_lm]*rr[:,tables.upf.channel_radial]
    mutations['omit_angular_current_leg'] = float(np.linalg.norm(angular_leg))
    omitted_mixed_radial_angular = (
        Se[:,tables.upf.channel_lm]*rs[:,tables.upf.channel_radial]*2*(p@f)[:,None]
        + Sf[:,tables.upf.channel_lm]*rs[:,tables.upf.channel_radial]*2*(p@e)[:,None])
    mutations['omit_radial_angular_contact_cross_terms'] = float(np.linalg.norm(omitted_mixed_radial_angular))
    assert min(mutations.values()) > 1e-5, mutations
    finite = []
    def operator(shift):
        B = centered_derivatives(exact,p+shift,e,f)[0]
        return np.diag(np.sum((p+shift)**2,axis=1))+B@d@B.conj().T
    reference_l = (np.diag(2*(p@e))+dense_l)@v
    reference_c = (2*(e@f)*np.eye(len(p))+dense_c)@v
    for h in [.04,.02,.01]:
        dl = ((operator(h*e)-operator(-h*e))/(2*h))@v
        dc = ((operator(h*(e+f))-operator(h*(e-f))-operator(h*(-e+f))+operator(-h*(e+f)))/(4*h*h))@v
        finite.append([h,float(np.linalg.norm(dl-reference_l)),float(np.linalg.norm(dc-reference_c))])
    assert finite[-1][1] < finite[0][1]/10 and finite[-1][2] < finite[0][2]/10, finite
    representations=[]
    for spacing in [.8,.4,.2]:
        approximation=HermiteRadial(exact,2.,spacing)
        data=centered_derivatives(approximation,p,e,f)
        representations.append([spacing]+[float(np.linalg.norm(a-b)) for a,b in zip(data,actual)])
        # C2 at an internal knot: stored pieces agree through two derivatives.
        knot=approximation.polynomial.x[1]
        for poly in approximation.derivatives:
            assert np.max(np.abs(poly(knot-1e-10)-poly(knot+1e-10)))<1e-8
    assert representations[-1][-1] < representations[0][-1]/10, representations
    result.update(literal_radial_spherical_errors=errors, mutation_exposures=mutations,
                  finite_operator_refinement=finite,C2_representation_refinement=representations,
                  scope='Synthetic finite operator only; material execution is separate.')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    run()
