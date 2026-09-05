"""Actual common-V source and whole-law physical-current forcing rows."""

import importlib.util
from pathlib import Path

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger

source = Path(__file__).resolve().parents[1]/"0188"/"verify_stationary.py"
spec = importlib.util.spec_from_file_location("frozen_0188_law",source)
frozen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(frozen)


def main():
    checks = CheckLedger("P251-0196-current")
    a,b = sp.symbols("A B",positive=True)
    kap = sp.Matrix(sp.symbols("kx ky kz",real=True))
    disp = sp.Matrix(sp.symbols("Dx Dy Dz",real=True))
    psi = ef.add(ef.scale(ef.trig(1),b),ef.scale(ef.trig(2),a))
    alpha = ef.add(ef.scale(ef.trig(1),b),ef.scale(ef.trig(2),-a))
    u = (psi,ef.scale(ef.trig(2,kind="sin"),a),ef.scale(ef.trig(1,kind="sin"),-b))

    def clean(field):
        return ef.add({wave:sp.cancel(value) for wave,value in field.items()})

    def add(*vectors):
        return tuple(clean(ef.add(*(vec[i] for vec in vectors))) for i in range(3))

    def scale(vector,value):
        return tuple(ef.scale(component,value) for component in vector)

    def law(polynomial):
        return frozen.transverse_law_average(polynomial,kap,disp)

    def pairing(left,right):
        return sp.cancel(ef.mul(left,right).get(ef.ZERO,0))

    def average_field(field,coefficient=1):
        return clean({wave:law(coefficient*value) for wave,value in field.items()})

    translation = tuple(ef.add(*(ef.scale(ef.derivative(u[i],j),-disp[j]) for j in range(3))) for i in range(3))
    kv = tuple({ef.ZERO:kap[i]} for i in range(3))
    q = scale(ef.cross(kv,translation),-1)
    au = ef.add(*(ef.scale(u[i],kap[i]) for i in range(3)))
    ud = ef.add(*(ef.scale(u[i],disp[i]) for i in range(3)))
    au_d = tuple(ef.scale(au,disp[i]) for i in range(3))
    b_source = scale(ef.leray(add(au_d,q)),-1)
    k_t = ef.add(*(ef.scale(translation[i],kap[i]) for i in range(3)))
    projection_derivative = tuple(ef.add({wave:sp.I*wave[i]*value/sum(v*v for v in wave) for wave,value in k_t.items()}) for i in range(3))
    l1 = add(projection_derivative,scale(ef.leray(au_d),-1))
    checks.check("common-V source follows the full Leray derivative for the same curl lift", not any(add(l1,scale(q,-1),scale(b_source,-1))))
    curl_b = ef.curl(b_source)[0]
    beta = clean({wave:-value/sum(v*v for v in wave) for wave,value in curl_b.items()})
    hbeta = clean({wave:sum(v*v for v in wave)*value for wave,value in beta.items()})
    r_b = clean(ef.add(b_source[0],ef.scale(hbeta,-1)))
    d = kap[1]*disp[1]-kap[2]*disp[2]
    e = kap[1]*disp[1]+kap[2]*disp[2]
    beta_d = average_field(beta,d)
    beta_e = average_field(beta,e)
    print("actual aggregate beta_d=",beta_d)
    print("actual aggregate beta_e=",beta_e)
    print("actual axial-plus-vorticity source r_b=",r_b)
    checks.check("the physical velocity source remains mean-zero solenoidal", not ef.divergence(b_source) and all(component.get(ef.ZERO,0) == 0 for component in b_source))
    checks.check("the actual planar source is completely in the first Fourier shell", all(sum(v*v for v in wave) == 1 for wave in beta))
    checks.check("the known whole-law strain normalization is retained", law(d*d) == sp.Rational(1,5))
    cx = ef.add(ef.scale(au,disp[0]),ef.scale(ud,kap[0]))
    basis = (ef.trig(1),ef.trig(2),ef.trig(1,kind="sin"),ef.trig(2,kind="sin"))
    c_coefficients = [2*pairing(cx,field) for field in basis]
    r_coefficients = [2*pairing(r_b,field) for field in basis]
    covariance = sp.Matrix(4,4,lambda i,j:law(c_coefficients[i]*r_coefficients[j]))
    print("complete passive-current forcing covariance on cosY,cosZ,sinY,sinZ=",covariance)
    reconstructions = []
    for original,coefficients in ((cx,c_coefficients),(r_b,r_coefficients)):
        reconstructed = ef.add(*(ef.scale(field,coef) for field,coef in zip(basis,coefficients,strict=True)))
        reconstructions.append(not clean(ef.add(original,ef.scale(reconstructed,-1))))
    checks.check("the current covariance is derived from all actual field components", all(reconstructions))
    checks.check("the aggregate planar V forcing is the actual alpha fifth", not clean(ef.add(beta_d,ef.scale(alpha,-sp.Rational(1,5)))) and not beta_e)
    expected_covariance = -sp.Matrix([[b*b,a*b,0,0],[a*b,a*a,0,0],[0,0,b*b,0],[0,0,0,a*a]])/5
    checks.check("the complete passive response has the derived rank-three physical forcing", covariance == expected_covariance)
    current = tuple(ef.add(ef.scale(au,disp[i]),ef.scale(ud,kap[i])) for i in range(3))

    def minus(vector):
        return scale(add(ef.leray(vector),scale(ef.curl(vector),-1)),sp.Rational(1,2))

    material_rate_shell = add(q,scale(au_d,-1),(ef.scale(alpha,-d/2),{},{}))
    preparation = sp.Symbol("tau",real=True)
    returned = add(scale(minus(material_rate_shell),-1),scale(minus(current),preparation))
    curl_returned = ef.curl(returned)[0]
    beta_returned = clean({wave:-value/sum(v*v for v in wave) for wave,value in curl_returned.items()})
    aggregate_returned = average_field(beta_returned,d)
    aggregate_e_returned = average_field(beta_returned,e)
    print("actual0188 correlated D first-shell d-weighted planar row=",aggregate_returned)
    print("actual0188 correlated D first-shell e-weighted planar row=",aggregate_e_returned)
    psi_jordan_moment = pairing(psi,ef.add(aggregate_returned,ef.scale(alpha,-sp.Rational(1,5))))
    print("actual aggregate necessary psi Jordan moment=",sp.factor(psi_jordan_moment))
    checks.check("all correlated D returns are actual stationary helicity fields", not any(add(ef.curl(returned),returned)))
    checks.check("the whole-law Jordan moment is an exact nonzero preparation defect", sp.factor(psi_jordan_moment-(a*a-b*b)*(preparation+1)/20) == 0)
    projected_current = ef.leray(current)
    transpose_du = tuple(ef.add(*(ef.mul(ef.derivative(u[j],i),projected_current[j]) for j in range(3))) for i in range(3))
    adjoint = ef.leray(add(ef.transport(u,projected_current),scale(transpose_du,-1)))
    adjoint_norm = sp.factor(law(ef.inner(adjoint,adjoint)))
    free_v_initial_slope = sp.factor(law(ef.inner(current,add(q,b_source))))
    print("physical common-V initial slope without microscopic return=",free_v_initial_slope)
    print("exact whole-law squared norm of the adjoint current derivative=",adjoint_norm)
    checks.check("the common-V first derivative retains the actual universal bare slope", sp.factor(free_v_initial_slope+2*(a*a+b*b)/15) == 0)
    checks.check("the adjoint current derivative genuinely vanishes at either one-wave endpoint", all(not clean({wave:sp.simplify(value.subs(a,0)) for wave,value in component.items()}) and not clean({wave:sp.simplify(value.subs(b,0)) for wave,value in component.items()}) for component in adjoint))
    passive_channel = kap[1]*disp[0]+kap[0]*disp[1]
    material_final_shell = add(material_rate_shell,returned)
    phase_density = average_field(ef.add(material_final_shell[0],current[0]),passive_channel)
    mean_slip_density = average_field(ef.add(ef.scale(au,disp[0]),ef.scale(ud,-kap[0])),passive_channel)
    print("actual configuration-return DV phase density=",phase_density)
    print("actual configuration-return V dot J density=",mean_slip_density)
    checks.check("the passive correlated channel has the actual fifth normalization", law(passive_channel**2) == sp.Rational(1,5))
    checks.check("the actual mean-velocity correction contributes no averaged VV energy", not mean_slip_density)
    checks.check("range D fields do not silently enter the passive phase matching row", law(d*passive_channel) == 0)
    checks.check("the complete phase row is the actual sine moment with its tenth normalization", not clean(ef.add(phase_density,ef.scale(ef.trig(2,kind="sin"),-a/10))))
    raw_momentum_density = average_field(ef.scale(ud,kap[0]),passive_channel)
    internal_phase_density = average_field(ef.add(material_final_shell[0],ef.scale(au,disp[0])),passive_channel)
    checks.check("the raw material momentum and internal phase rows remain separately accounted", not clean(ef.add(raw_momentum_density,ef.scale(ef.trig(2,kind="sin"),-a/10))) and not internal_phase_density)
    return int(checks.finish())


if __name__ == "__main__":
    raise SystemExit(main())
