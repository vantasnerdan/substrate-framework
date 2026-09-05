"""Full finite-Fourier Bloch Kelvin operator and its actual first pressure jet."""

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0223-first-pressure")
    wave = sp.Symbol("k",real=True,nonzero=True)
    coefficient = sp.Rational(1,100)
    psi = ef.add(ef.trig(2),ef.scale(ef.trig(1),coefficient))
    background = (psi,ef.trig(2,kind="sin"),ef.scale(ef.trig(1,kind="sin"),-coefficient))
    vorticity = ef.curl(background)

    def clean(field):
        return {q:value for q,c in field.items() if (value := sp.cancel(c)) != 0}

    def clean_vector(vector):
        return tuple(clean(field) for field in vector)

    def vector_add(left,right,factor=1):
        return clean_vector(tuple(ef.add(a,ef.scale(b,factor)) for a,b in zip(left,right,strict=True)))

    def inverse(field):
        return {q:c/(q[1]**2+q[2]**2) for q,c in field.items() if q[1]**2+q[2]**2 != 0}

    def lift(field):
        return {(wave,q[1],q[2]):c for q,c in field.items()}

    def coefficient_jet(field,order):
        return clean({(0,q[1],q[2]):sp.diff(c,wave,order).limit(wave,0)/sp.factorial(order) for q,c in field.items()})

    def transport(field):
        return ef.transport(background,(field,{},{}))[0]

    def mean_free(field):
        return {q:c for q,c in field.items() if q != ef.ZERO}

    streams = [ef.trig(1,2),ef.mul(ef.trig(1,2),ef.trig(2,2)),ef.trig(1)]
    for index,stream in enumerate(streams):
        inverse_stream = inverse(stream)
        generator = (lift(ef.scale(stream,-1)),
            lift(ef.add(ef.scale(ef.derivative(stream,2),-1),ef.scale(ef.derivative(inverse_stream,1),-sp.I*wave))),
            lift(ef.add(ef.derivative(stream,1),ef.scale(ef.derivative(inverse_stream,2),-sp.I*wave))))
        ledger.check(f"the actual axial Bloch generator is exactly solenoidal {index}", not clean(ef.divergence(generator)))
        velocity = clean_vector(ef.leray(ef.cross(generator,vorticity)))
        lin_rate = vector_add(vector_add(velocity,ef.transport(background,generator),-1),ef.transport(generator,background))
        canonical_rate = clean_vector(ef.material_kelvin_operator(background,generator))
        ledger.check(f"full pressure and the canonical Kelvin generator give the SAME Lin rate {index}", not any(vector_add(lin_rate,canonical_rate,-1)))
        ledger.check(f"the induced actual velocity retains three-dimensional divergence {index}", not clean(ef.divergence(velocity)))
        zero_velocity = tuple(coefficient_jet(field,0) for field in velocity)
        ledger.check(f"the exact zero-wave field is the positive passive supplier {index}", zero_velocity == (transport(stream),{},{}))
        h_rate1 = coefficient_jet(lin_rate[0],1)
        expected_h_rate1 = ef.scale(mean_free(ef.mul(psi,stream)),2*sp.I)
        ledger.check(f"the actual axial generator develops the non-Doppler first jet {index}", not clean(ef.add(h_rate1,ef.scale(expected_h_rate1,-1))))
        normal_rate1 = ({},coefficient_jet(lin_rate[1],1),coefficient_jet(lin_rate[2],1))
        s_rate1 = ef.scale(inverse(ef.curl(normal_rate1)[0]),-1)
        f = ef.scale(inverse_stream,-1)
        f_minus_s = ef.add(f,ef.scale(stream,-1))
        weighted_gradient = ({},ef.mul(psi,ef.derivative(f_minus_s,1)),ef.mul(psi,ef.derivative(f_minus_s,2)))
        symmetric_gradient = ({},*tuple(ef.add(*(ef.mul(ef.add(ef.derivative(background[i],j),ef.derivative(background[j],i)),ef.derivative(f,j)) for j in (1,2))) for i in (1,2)))
        pressure_source = ef.add(ef.divergence(weighted_gradient),ef.curl(symmetric_gradient)[0])
        expected_s_rate1 = ef.scale(inverse(pressure_source),-sp.I)
        ledger.check(f"curling the FULL Lin operator gives the stated scalar pressure correction {index}", not clean(ef.add(s_rate1,ef.scale(expected_s_rate1,-1))))
        first_complement_force = clean(ef.add(h_rate1,s_rate1))
        ledger.check(f"the passive s+h=0 relation is not a finite-wave invariant {index}", bool(first_complement_force))
        # Independently reconstruct the pressure-free planar complement of v1.
        velocity1 = tuple(coefficient_jet(field,1) for field in velocity)
        potential = ( {},ef.scale(ef.derivative(inverse(transport(stream)),1),sp.I),
            ef.scale(ef.derivative(inverse(transport(stream)),2),sp.I))
        complement = vector_add(velocity1,potential,-1)
        complement = ({},complement[1],complement[2])
        expected_complement = ef.scale(ef.derivative(inverse_stream,2),-1),ef.derivative(inverse_stream,1)
        expected_complement = ef.leray(({},ef.scale(ef.mul(psi,expected_complement[0]),-sp.I),ef.scale(ef.mul(psi,expected_complement[1]),-sp.I)))
        ledger.check(f"the defining full pressure gives the actual horizontal complement {index}", not any(vector_add(complement,expected_complement,-1)))
        ledger.check(f"the complement is genuinely nonzero and planar solenoidal {index}", any(complement) and not clean(ef.divergence(complement)))
        # The exact full physical mean is exposed, not identified with a tag.
        ledger.check(f"normal reflection protects the actual transverse mean {index}", all(not field.get((wave,0,0),0) for field in velocity[1:]))
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
