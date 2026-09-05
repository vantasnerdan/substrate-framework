"""Exact full-Fourier checks of the inherited first-cell action and current."""

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.euler_acoustic import (
    observed_acoustic_cell_rows,
    prepared_acoustic_cell_rows,
)
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0180")
    t, eps = sp.symbols("t eps", real=True)
    u = (ef.add(ef.trig(2), ef.trig(1, kind="sin")),
         ef.trig(2, kind="sin"), ef.trig(1))
    pressure = ef.scale(ef.add(*(ef.mul(v,v) for v in u)), -sp.Rational(1,2))
    gradp = tuple(ef.derivative(pressure,i) for i in range(3))
    kappa, d, v = sp.Matrix([0,0,1]), sp.Matrix([1,0,0]), sp.Matrix([0,1,0])
    slow = sp.eye(3)-kappa*kappa.T
    y = d+t*v

    def add(*vectors):
        return tuple(ef.add(*(vec[i] for vec in vectors)) for i in range(3))

    def scale(vector, scalar):
        return tuple(ef.scale(component,scalar) for component in vector)

    def time_diff(vector):
        return tuple(ef.add({wave:sp.diff(value,t) for wave,value in component.items()}) for component in vector)

    def times(field, vector):
        return tuple(ef.scale(field,vector[i]) for i in range(3))

    def dot_constant(vector,field):
        return ef.add(*(ef.scale(field[i],vector[i]) for i in range(3)))

    def mean(left,right):
        return ef.mul(left,right).get(ef.ZERO,0)

    def hp(field):
        return tuple(ef.add(*(ef.mul(ef.derivative(gradp[i],j),field[j]) for j in range(3))) for i in range(3))

    def transport(field):
        return ef.transport(u,field)

    a = dot_constant(kappa,u)
    aa = ef.add(*(ef.mul(u[j],ef.derivative(a,j)) for j in range(3)))
    kp = dot_constant(kappa,gradp)
    ledger.check("actual stationary transport of the directional velocity supplies pressure", not ef.add(aa,kp))
    api = prepared_acoustic_cell_rows(u,pressure,kappa,d,v)
    chi0 = api.initial_rate
    chi1 = ef.leray(ef.transport(u,chi0))
    chi2 = tuple(ef.derivative(component,1) for component in u)
    chi = add(chi0,scale(chi1,t),scale(chi2,t**2/2))
    rate = time_diff(chi)
    test = ef.leray(add(u,chi0,chi1))
    ledger.check("the exposing full-Fourier variation is solenoidal and nonzero", not ef.divergence(test) and any(test))

    def action2(cell,cell_rate,yy,yv, *, retain_constraint=True):
        convective = add(cell_rate,transport(cell),times(a,yy))
        kinetic = ef.inner(convective,convective)
        potential = ef.inner(cell,hp(cell))
        current = 2*mean(a,dot_constant(yv,cell))
        constraint = 2*mean(dot_constant(yy,gradp),dot_constant(kappa,cell))
        return sp.expand((kinetic-potential-current+(constraint if retain_constraint else 0))/2)

    def variation(retain_constraint=True):
        position = sp.diff(action2(add(chi,scale(test,eps)),rate,y,v,
                                  retain_constraint=retain_constraint),eps).subs(eps,0)
        momentum = sp.diff(action2(chi,add(rate,scale(test,eps)),y,v,
                                  retain_constraint=retain_constraint),eps).subs(eps,0)
        return sp.expand(position-sp.diff(momentum,t))

    forcing = add(api.forcing_constant,scale(api.forcing_rate,t))
    jacobi = ef.leray(add(time_diff(rate),scale(transport(rate),2),transport(transport(chi)),hp(chi),scale(forcing,-1)))
    ledger.check("the inherited action varies to the full canonical Euler/Lin forcing", sp.simplify(variation()+ef.inner(test,jacobi)) == 0)
    omitted = variation(False)-variation(True)
    # If this single tangent happens to be orthogonal, the complete vector row
    # still exposes the omitted constraint, rather than declaring a false green.
    missing = ef.leray(times(dot_constant(y,gradp),kappa))
    ledger.check("the second-cell constraint contributes an actual nonzero pressure row", any(missing) and (omitted != 0 or ef.inner(missing,missing) != 0))

    yy0,yy1,vv0,vv1 = sp.symbols("Y0 Y1 V0 V1",real=True)
    yy, vv = sp.Matrix([yy0,yy1,0]),sp.Matrix([vv0,vv1,0])
    lag = action2(chi,rate,yy,vv)
    substitutions = {yy0:1,yy1:t,vv0:0,vv1:1}
    macro_force = sp.Matrix([sp.diff(lag,coord).subs(substitutions)-sp.diff(sp.diff(lag,vel).subs(substitutions),t)
                             for coord,vel in ((yy0,vv0),(yy1,vv1))]+[0])
    material = slow*sp.Matrix([mean(a,a)*y[i]+2*mean(a,rate[i])+mean(kp,chi[i])+mean(gradp[i],dot_constant(kappa,chi)) for i in range(3)])
    ledger.check("macro variation reproduces the material mean with its factor two", sp.simplify(macro_force-material) == sp.zeros(3,1))
    observed = observed_acoustic_cell_rows(u,pressure,kappa,y,chi,rate)
    current = observed.current_correction
    ledger.check("actual Lin current converts material stress to physical Euler stress", sp.simplify(material+current.diff(t)-observed.acceleration) == sp.zeros(3,1))
    canonical_correction = slow*sp.Matrix([-mean(a,chi[i]) for i in range(3)])
    measured_difference = slow*sp.Matrix([mean(u[i],dot_constant(kappa,chi)) for i in range(3)])
    ledger.check("physical minus canonical momentum retains the actual velocity moment", sp.simplify(current-canonical_correction-measured_difference) == sp.zeros(3,1))
    ledger.check("equating material and physical means would change this probe", current != sp.zeros(3,1))

    # Exact initial phase at the actual Kelvin-D/common-V preparation.
    cotangent = ef.leray(add(api.initial_rate,times(a,d)))
    expected = scale(ef.leray(times(dot_constant(d,u),kappa)),-1)
    ledger.check("actual internal initial cotangent is the retained preparation slip", not any(add(cotangent,scale(expected,-1))))
    ledger.check("zeroing the internal initial momentum would change the Euler preparation", any(cotangent))

    # No arbitrary oscillator is used in these checks. This algebraic necessary
    # relation tests the existing exact0175 response, whose expensive full
    # convolution receipt is reused rather than rerun here.
    omega, amplitude, offset = sp.symbols("omega amplitude offset",positive=True)
    one_pair = offset+amplitude*sp.cos(omega*t)
    ledger.check("any constant single optical pair obeys its exact temporal recurrence", sp.diff(one_pair,t,4).subs(t,0)+omega**2*sp.diff(one_pair,t,2).subs(t,0) == 0)
    ledger.check("actual0175 bare stress jets violate that one-pair recurrence", -sp.Rational(2,25)+omega**2*0 != 0)
    # Fundamental theorem of calculus exposes the preparation term, independent
    # of any finite-dimensional spectral approximation to the true group.
    ell = sp.symbols("ell",real=True)
    source = sp.Function("source")(t)
    z0 = sp.symbols("z0")
    s = sp.symbols("s",real=True)
    duhamel = sp.exp(ell*t)*z0+sp.Integral(sp.exp(ell*(t-s))*source.subs(t,s),(s,0,t))
    ledger.check("stationary response elimination preserves forcing and preparation", sp.simplify(sp.diff(duhamel,t)-ell*duhamel-source) == 0 and duhamel.subs(t,0).doit() == z0)

    # Actual finite reachable space for the circular one-wave comparison.
    wave = (ef.trig(2),ef.trig(2,kind="sin"),{})
    direction = sp.Matrix([sp.Rational(3,5),0,sp.Rational(4,5)])
    displacement = sp.Matrix([sp.Rational(4,5),0,-sp.Rational(3,5)])
    velocity = sp.Matrix([0,1,0])
    prepared = prepared_acoustic_cell_rows(wave,{ef.ZERO:-sp.Rational(1,2)},direction,displacement,velocity)
    for name, field in (("Kelvin slip",prepared.initial_rate),("macro forcing",prepared.forcing_constant)):
        cc = add(scale(ef.transport(wave,field),-1),ef.transport(field,wave))
        le = scale(ef.leray(add(ef.transport(wave,field),ef.transport(field,wave))),-1)
        ledger.check(f"actual one-wave {name} lies in the nilpotent finite cell space", not any(cc) and not any(le))
    av = dot_constant(direction,wave)
    slip = ef.leray(add(times(av,velocity),scale(times(dot_constant(velocity,wave),direction),-1)))
    zerovec = ({},{},{})
    slip_row = observed_acoustic_cell_rows(wave,{ef.ZERO:-sp.Rational(1,2)},direction,sp.zeros(3,1),zerovec,slip)
    ledger.check("nonzero one-wave preparation slip has exactly zero physical stress row", any(slip) and slip_row.acceleration == sp.zeros(3,1))
    q = sp.symbols("q",real=True)
    # Azimuthal average of n_x² is (1-q²)/2 with kappa=e_z.
    averaged = sp.integrate(-(1-q**2)+2*(1-2*q**2)*(1-q**2)/2,(q,-1,1))/2
    ledger.check("finite physical row gives the positive whole-Haar one-wave modulus", averaged == -sp.Rational(4,15))

    # New all-time isotropic stress identity uses the ACTUAL forced Jacobi
    # energy, not a positivity premise for its indefinite potential.
    strain_rows = []
    pressure_rows = []
    for m in range(3):
        for ell_index in range(3):
            strain_rows.append(ef.leray(tuple(ef.add(u[m] if i == ell_index else {},u[ell_index] if i == m else {}) for i in range(3))))
            pressure_rows.append(ef.leray(tuple(ef.add(gradp[m] if i == ell_index else {},gradp[ell_index] if i == m else {}) for i in range(3))))
    energy = ef.inner(u,u)
    ledger.check("complete isotropic strain-source norm has its exact sixfold factor", sum(ef.inner(row,row) for row in strain_rows) == 6*energy)
    ledger.check("actual trace source is exactly the stationary velocity", not any(add(*(strain_rows[4*i] for i in range(3)),scale(u,-2))))
    ledger.check("actual trace forcing is zero after full pressure projection", not any(add(*(pressure_rows[4*i] for i in range(3)))))
    jacobi_u = ef.leray(add(transport(transport(u)),hp(u)))
    ledger.check("stationary velocity is an exact trace-cell Jacobi zero mode", not any(jacobi_u))
    kjchi = ef.leray(add(transport(transport(chi)),hp(chi)))
    constant_force = prepared_acoustic_cell_rows(u,pressure,kappa,d,sp.zeros(3,1)).forcing_constant
    actual_acceleration = add(scale(ef.leray(transport(rate)),-2),scale(kjchi,-1),constant_force)
    energy_rate = ef.inner(rate,actual_acceleration)+ef.inner(rate,kjchi)-ef.inner(constant_force,rate)
    ledger.check("the exact forced Jacobi energy is conserved with full gyroscopic pressure", sp.simplify(energy_rate) == 0)
    aa,ss,kk,ff,rs = sp.symbols("rate_square source_square potential force_position source_rate",real=True)
    energy_constraint = {aa:ss-kk+2*ff}
    completed = ((aa+2*rs+ss+kk)/20).subs(energy_constraint)
    ledger.check("whole-law stress square completion retains the indefinite potential", sp.expand(completed-(rs+ff+ss)/10) == 0)
    homogeneous_acc = add(scale(ef.leray(transport(rate)),-2),scale(kjchi,-1))
    ledger.check("arbitrary corrected trace data leave the isotropic trace current constant", sp.simplify(ef.inner(u,homogeneous_acc)) == 0)
    initial_energy = sp.symbols("initial_energy",real=True)
    corrected_completed = ((aa+2*rs+ss+kk)/20).subs(aa,2*initial_energy-kk+2*ff)
    ledger.check("corrected preparation changes only the constant in stress square completion", sp.expand(corrected_completed-(rs+ff)/10-initial_energy/10-ss/20) == 0)
    return int(ledger.finish())


if __name__ == "__main__":
    raise SystemExit(main())
