"""Exact complete correlated-current reduction on the actual two-wave Euler field."""

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.homogenization import sphere_fourth_moment_isotropic
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0184")
    y,z = sp.symbols("Y Z",real=True)
    kap = sp.Matrix(sp.symbols("kx ky kz",real=True))
    disp = sp.Matrix(sp.symbols("Dx Dy Dz",real=True))
    psi, alpha = sp.cos(y)+sp.cos(z), sp.cos(y)-sp.cos(z)
    u = sp.Matrix([psi,sp.sin(z),-sp.sin(y)])
    c = kap.dot(u)*disp+disp.dot(u)*kap
    planar_row = sp.diff(c[1],z)-sp.diff(c[2],y)
    hc = -sp.diff(c[0],y,2)-sp.diff(c[0],z,2)
    d = kap[1]*disp[1]-kap[2]*disp[2]
    e = kap[1]*disp[1]+kap[2]*disp[2]
    full = sp.expand_trig(sp.expand(hc+planar_row+e*psi+d*alpha))
    ledger.check("complete physical current reduces using only actual transversality", sp.simplify(full-2*kap.dot(disp)*psi) == 0)
    ledger.check("the omitted axial-vorticity term would leave real mixed strain rows", sp.simplify(planar_row+e*psi+d*alpha-2*kap.dot(disp)*psi) != 0)
    ledger.check("the axial current is on the genuine first scalar Fourier shell", sp.simplify(hc-c[0]) == 0)
    adv_psi = sp.sin(z)*sp.diff(psi,y)-sp.sin(y)*sp.diff(psi,z)
    adv_alpha = sp.sin(z)*sp.diff(alpha,y)-sp.sin(y)*sp.diff(alpha,z)
    ledger.check("psi is an actual streamline first integral", adv_psi == 0)
    ledger.check("alpha supplies the nonzero actual forcing harmonic", sp.simplify(adv_alpha+2*sp.sin(y)*sp.sin(z)) == 0)

    fourth = sphere_fourth_moment_isotropic()
    def pair_moment(i,j,a,b):
        # Conditional E[D_j D_b | kappa]=(delta_jb-kappa_j kappa_b)/2.
        return sp.KroneckerDelta(j,b)*sp.KroneckerDelta(i,a)/6-fourth[i,a,j,b]/2

    def average(polynomial):
        result = 0
        for powers,coefficient in sp.Poly(sp.expand(polynomial),*kap,*disp).terms():
            ii = [index for index in range(3) for _ in range(powers[index])]
            jj = [index for index in range(3) for _ in range(powers[index+3])]
            if len(ii) != 2 or len(jj) != 2:
                raise ValueError("the current contraction must have exactly two input factors each")
            result += coefficient*pair_moment(ii[0],jj[0],ii[1],jj[1])
        return sp.simplify(result)

    sigma = average(d**2)
    ledger.check("arbitrary correlated planar forcing retains its nonzero whole-law weight", sigma == sp.Rational(1,5))
    ledger.check("forced axial return cancels only after the complete whole law", average(d*c[0]) == 0)
    ledger.check("the planar and axial cross-sector cancellation is not separability", average(d*c[1]) == sp.sin(z)/5 and average(d*c[2]) == sp.sin(y)/5)

    base = (ef.add(ef.trig(1),ef.trig(2)),ef.trig(2,kind="sin"),ef.scale(ef.trig(1,kind="sin"),-1))
    aa = ef.add(ef.trig(1),ef.scale(ef.trig(2),-1))
    wave_y = (ef.trig(1),{},ef.scale(ef.trig(1,kind="sin"),-1))
    wave_z = (ef.trig(2),ef.trig(2,kind="sin"),{})
    force0 = ef.leray(ef.cross(wave_z,wave_y))
    ledger.check("the forcing is the actual complete solenoidal Euler cross field", not ef.divergence(force0) and force0 == ef.cross(wave_z,wave_y))

    def add(*vectors):
        return tuple(ef.add(*(vec[i] for vec in vectors)) for i in range(3))

    def scale(vector,value):
        return tuple(ef.scale(component,value) for component in vector)

    def transport(field):
        return ef.transport(base,(field,{},{}))[0]

    def pairing(left,right):
        return sp.simplify(ef.mul(left,right).get(ef.ZERO,0))

    def multiplier(field,symbol):
        return ef.add({wave:symbol(sp.Integer(sum(q*q for q in wave)))*value for wave,value in field.items() if wave != ef.ZERO})

    translation = tuple(ef.add(*(ef.scale(ef.derivative(base[i],axis),-disp[axis]) for axis in range(3))) for i in range(3))
    directional_translation = tuple(ef.add(*(ef.scale(ef.derivative(translation[i],axis),kap[axis]) for axis in range(3))) for i in range(3))
    full_forcing = ef.leray(ef.cross(base,directional_translation))
    ledger.check("full translation-leading Bloch forcing reduces to the same d sector", not any(add(full_forcing,scale(force0,-d))))

    phi = ef.add(ef.mul(ef.trig(1,kind="sin"),ef.trig(2,kind="sin")),ef.trig(1,2),ef.scale(aa,3))
    b = ef.add(ef.trig(2,2),ef.trig(1))
    velocity = (b,ef.scale(ef.derivative(phi,2),-1),ef.derivative(phi,1))
    eta = ef.curl(velocity)[0]
    rr = ef.add(b,eta)
    velocity_t = add(scale(ef.leray(add(ef.transport(base,velocity),ef.transport(velocity,base))),-1),scale(force0,sigma))
    rt = ef.add(velocity_t[0],ef.curl(velocity_t)[0])
    ledger.check("the full forced axial-plus-vorticity current has its retained source", not ef.add(rt,transport(rr),ef.scale(transport(aa),sigma/2)))
    ledger.check("declaring the forced return passive would change the exact equation", bool(ef.add(rt,transport(rr))))

    def qq(field):
        return multiplier(field,lambda norm:sp.sqrt((norm-1)/norm))

    def generator(field):
        return ef.scale(qq(transport(qq(field))),-1)

    bphi = multiplier(phi,lambda norm:norm-1)
    phit = multiplier(ef.add(ef.scale(transport(bphi),-1),ef.scale(transport(aa),sigma)),lambda norm:1/norm)
    ledger.check("the planar vorticity equation is the full forced Euler curl", not ef.add(ef.curl(velocity_t)[0],multiplier(phit,lambda norm:norm)))
    zz = multiplier(phi,lambda norm:sp.sqrt(norm*(norm-1)))
    zt = multiplier(phit,lambda norm:sp.sqrt(norm*(norm-1)))
    ff = ef.scale(qq(transport(aa)),sigma)
    ledger.check("the averaged physical field obeys the actual forced unitary range equation", not ef.add(zt,ef.scale(generator(zz),-1),ef.scale(ff,-1)))
    current_t = -pairing(aa,phit)
    energy_t = pairing(zz,zt)
    ledger.check("the complete correlated current has the exact positive range-energy factor", sp.simplify(current_t+energy_t/sigma) == 0 and current_t != 0)
    ledger.check("using the one-orientation normalization would fail this whole-law current", sp.simplify(current_t+energy_t) != 0)
    ledger.check("the nonzero averaged forcing retains the exact natural norm", pairing(ff,ff) == sigma**2/2)
    ledger.check("first-shell physical current is retained outside the range coordinate", not multiplier(aa,lambda norm:sp.sqrt(norm*(norm-1))) and pairing(aa,aa) == 1)
    other = ef.add(ef.trig(2,2,kind="sin"),ef.mul(ef.trig(1),ef.trig(2)))
    ledger.check("actual range generator is skew on separate full-Fourier probes", sp.simplify(pairing(generator(zz),other)+pairing(zz,generator(other))) == 0)
    psi_fourier = ef.add(ef.trig(1),ef.trig(2))
    ledger.check("arbitrary corrected psi current is actually conserved", pairing(psi_fourier,phit) == 0)
    ledger.check("streamline parity obstruction retains its nonzero whole-law coefficient", pairing(aa,ef.scale(aa,-sigma)) == -sigma)
    reversed_phit = multiplier(ef.add(transport(bphi),ef.scale(transport(aa),sigma)),lambda norm:1/norm)
    reversed_zt = multiplier(reversed_phit,lambda norm:sp.sqrt(norm*(norm-1)))
    ledger.check("whole time reversal preserves the current-energy sign", sp.simplify(pairing(aa,reversed_phit)+pairing(zz,reversed_zt)/sigma) == 0)
    print(f"whole-law weight={sigma}; current_t={current_t}; range_energy_t={energy_t}")
    return int(ledger.finish())


if __name__ == "__main__":
    raise SystemExit(main())
