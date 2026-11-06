"""Rational Machin enclosures for the source conversion's exact pi dependence."""
import sympy as sp


def arctangent_bounds(denominator,terms):
    """Alternating Taylor remainder for atan(1/denominator), denominator > 1."""
    if denominator<=1 or terms<1:
        raise ValueError('Convergent positive arctangent series required')
    partial=sum(sp.Rational((-1)**n,(2*n+1)*denominator**(2*n+1)) for n in range(terms))
    remainder=sp.Rational(1,(2*terms+1)*denominator**(2*terms+1))
    return (partial-remainder,partial) if terms%2 else (partial,partial+remainder)


def machin_bounds(terms):
    a,b=arctangent_bounds(5,terms);c,d=arctangent_bounds(239,terms)
    return 16*a-4*d,16*b-4*c


def run(ledger,bohr,ry):
    # (5+i)^4/(239+i) has equal positive real/imaginary parts. Since
    # 0 < 4 atan(1/5)-atan(1/239) < 1 < pi/2, its argument is pi/4.
    z=sp.expand((5+sp.I)**4*(239-sp.I))
    ledger.check('Machin identity exact complex tangent and principal-angle enclosure',
                 sp.re(z)==sp.im(z) and sp.re(z)>0 and
                 4*(sp.Rational(1,5)-sp.Rational(1,3*5**3))-sp.Rational(1,239)>0 and
                 4*sp.Rational(1,5)<1)
    lo,hi=machin_bounds(32);coarse_lo,coarse_hi=machin_bounds(8)
    ledger.check('alternating-series rational pi bounds nested and positive',0<coarse_lo<lo<hi<coarse_hi)
    ledger.check('omitted alternating remainder exposed by next partial sum',
                 arctangent_bounds(5,2)[0]<arctangent_bounds(5,3)[0] and
                 arctangent_bounds(5,4)[1]<arctangent_bounds(5,2)[1])
    a_coeff=sp.simplify(bohr*sp.pi**2);ry_coeff=sp.simplify(ry/sp.pi**2)
    factor_coeff=sp.simplify((ry/bohr)/sp.pi**4)
    ledger.check('source Bohr/Ry/force pi powers reduce to exact positive rational coefficients',
                 all(x.is_Rational and x>0 for x in (a_coeff,ry_coeff,factor_coeff)))
    force_bounds=(factor_coeff*lo**4,factor_coeff*hi**4)
    ledger.check('certified conversion relative enclosure narrower than 1e-40',
                 (force_bounds[1]-force_bounds[0])/force_bounds[0]<sp.Rational(1,10**40))
    return {'force_bounds':force_bounds,'pi_bounds':(lo,hi),'force_coefficient':factor_coeff,
            'proof':'Machin identity with principal angle, alternating atan series with signed first-omitted-term remainder; 32 terms each.',
            'scope':'Rational enclosure of the fixed exact source convention, not uncertainty in physical constants.'}
