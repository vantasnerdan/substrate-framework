"""Exact triangular-cell Biot--Savart/action algebra; no lattice numerics."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0139")
    theta=s.symbols("theta",real=True)
    angles=[theta+j*s.pi/3 for j in range(6)]
    ledger.check("sixfold second moment", s.trigsimp(sum(s.sin(t)**2 for t in angles)/6-s.Rational(1,2))==0)
    ledger.check("sixfold fourth mixed moment", s.trigsimp(s.expand_trig(sum(s.sin(t)**2*s.cos(t)**2 for t in angles)/6)-s.Rational(1,8))==0)
    ledger.check("sixfold odd shear cross", s.trigsimp(s.expand_trig(sum(s.sin(t)**3*s.cos(t) for t in angles)/6))==0)
    q,reg,area,rho,circulation,k=s.symbols("q reg area rho circulation k",positive=True)
    f=s.exp(-reg*q)/q
    coeff=q**2*s.diff(f,q,2)/4+q*s.diff(f,q)/2
    ledger.check("full regularized shear derivative", s.simplify(coeff-reg**2*q*s.exp(-reg*q)/4)==0)
    radial=s.integrate(coeff,(q,0,s.oo))
    ledger.check("radial boundary integral", radial==s.Rational(1,4))
    reciprocal_density=area/(2*s.pi)**2
    lattice_limit=reciprocal_density*s.pi*radial
    ledger.check("reciprocal normalization", s.simplify(lattice_limit-area/(16*s.pi))==0)
    gamma=circulation/area
    mu=s.simplify(rho*gamma**2*lattice_limit)
    ledger.check("derived physical shear", s.simplify(mu-rho*circulation**2/(16*s.pi*area))==0)
    rr=s.symbols("rr",nonnegative=True)
    poisson_term=area/(4*s.pi*reg)*s.exp(-rr**2/(4*reg))
    poisson_shear=-reg**2*s.diff(poisson_term,reg)/4
    ledger.check("exact Poisson error term", s.simplify(poisson_shear-area/(16*s.pi)*(1-rr**2/(4*reg))*s.exp(-rr**2/(4*reg)))==0)
    cap_h=s.Function("h")(q)
    general_f=cap_h/q
    ledger.check("radial profile-independent total derivative", s.simplify(q**2*s.diff(general_f,q,2)/4+q*s.diff(general_f,q)/2-s.diff(q**2*s.diff(general_f,q),q)/4)==0)

    longitudinal,transverse=s.symbols("longitudinal transverse",real=True)
    j=s.Matrix([[0,-1],[1,0]])
    wave=s.Matrix([k,0])
    delta_omega=-s.I*gamma*k*longitudinal
    measured=-s.I*j*wave*delta_omega/k**2
    ledger.check("actual inverse curl transverse velocity", s.simplify(measured-s.Matrix([0,-gamma*longitudinal]))==s.zeros(2,1))
    ledger.check("actual incompressible mean", s.simplify(wave.dot(measured))==0)
    macro_energy=s.simplify(rho*measured.dot(measured)/2)
    ledger.check("longitudinal pole is fluid kinetic energy", s.simplify(macro_energy-rho*gamma**2*longitudinal**2/2)==0)
    hll=rho*gamma**2-mu*k**2
    htt=mu*k**2
    hlt=s.symbols("hlt",real=True)
    rate=s.symbols("rate",real=True)
    lagrangian=-rho*gamma*longitudinal*rate-(hll*longitudinal**2+2*hlt*longitudinal*transverse+htt*transverse**2)/2
    eliminated=s.solve(s.diff(lagrangian,longitudinal),longitudinal)[0]
    reduced=s.simplify(lagrangian.subs(longitudinal,eliminated))
    expected=(rho*gamma)**2/(2*hll)*(rate+hlt*transverse/(rho*gamma))**2-htt*transverse**2/2
    ledger.check("exact momentum Schur action", s.simplify(reduced-expected)==0)
    measured_rate=s.simplify(-gamma*eliminated)
    ledger.check("physical mean observation after elimination", s.simplify(measured_rate-rho*gamma**2/hll*(rate+hlt*transverse/(rho*gamma)))==0)
    mass=s.diff(reduced,rate,2)
    ledger.check("physical mass is ambient density in long-wave limit", s.simplify(s.limit(mass,k,0)-rho)==0)
    ledger.check("full retained mass second jet", s.simplify(s.diff(mass,k,2).subs(k,0)/2-mu/gamma**2)==0)
    omega_squared=s.simplify(hll*htt/(rho*gamma)**2)
    ledger.check("same-action physical acoustic coefficient", s.simplify(s.diff(omega_squared,k,2).subs(k,0)/2-mu/rho)==0)
    ledger.check("nonzero-k positive matrix sufficient interval", s.simplify(hll.subs(k**2,rho*gamma**2/(2*mu))-rho*gamma**2/2)==0)
    ledger.check("wrong mean-velocity sign is exposed", s.simplify(measured[1]-gamma*longitudinal)!=0)
    translation,phase,amplitude=s.symbols("translation phase amplitude",real=True)
    fourier_velocity=amplitude*s.exp(-s.I*phase*translation)
    tangent=s.diff(fourier_velocity,translation).subs(translation,0)
    curvature=s.diff(fourier_velocity,translation,2).subs(translation,0)
    first_square=s.expand(s.conjugate(tangent)*tangent)
    full_translation_hessian=s.simplify(first_square+s.re(amplitude*curvature))
    ledger.check("actual Fourier translation full energy Hessian cancels", full_translation_hessian==0)
    ledger.check("omitting second variation creates spurious translation energy", s.simplify(first_square-amplitude**2*phase**2)==0 and first_square!=0)
    raise SystemExit(ledger.finish())


if __name__ == "__main__":
    main()
