"""Independent full-pressure Euler check of the physical mean bending row."""

import sympy as s

from substrate_framework import euler_fourier as f
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0144-axial-Euler-pressure")
    eps = s.Symbol("epsilon", positive=True)
    psi = f.mul(f.trig(0), f.trig(1))
    base = (f.scale(f.derivative(psi, 1), -1), f.derivative(psi, 0), {})
    omega = f.curl(base)
    carrier = (0, 0, eps)
    generator = ({carrier: s.Integer(1)}, {}, {})
    kelvin_velocity = f.leray(f.cross(generator, omega))

    def linear_euler(velocity):
        first = f.transport(base, velocity)
        second = f.transport(velocity, base)
        return f.leray(tuple(f.scale(f.add(first[j], second[j]), -1) for j in range(3)))

    def mean(velocity):
        return s.Matrix([s.factor(component.get(carrier, 0)) for component in velocity])

    initial_mean = mean(kelvin_velocity)
    initial_acceleration = mean(linear_euler(kelvin_velocity))
    covariance_column = s.Matrix([f.mul(base[j], base[0]).get(f.ZERO, 0) for j in range(3)])
    second_jet = initial_acceleration.applyfunc(lambda value: s.limit(value/eps**2, eps, 0))
    ledger.check("actual three-dimensional Kelvin preparation has zero initial physical mean velocity",
                 initial_mean == s.zeros(3, 1))
    ledger.check("complete nonlocal pressure gives the positive axial covariance stiffness",
                 second_jet == -covariance_column and covariance_column[0] > 0)
    common_velocity = generator
    first_time = linear_euler(common_velocity)
    second_time = linear_euler(first_time)
    ledger.check("actual common initial velocity retains the full physical mean",
                 mean(common_velocity) == s.Matrix([1, 0, 0]))
    ledger.check("common-V phase has no immediate acoustic second time derivative",
                 mean(second_time).applyfunc(s.simplify) == s.zeros(3, 1))
    ledger.check("the two preparations are physically distinct, not a renamed momentum",
                 mean(common_velocity) != initial_mean
                 and initial_acceleration != s.zeros(3, 1))
    print("Actual full-pressure initial mean acceleration:", initial_acceleration.T)
    print("Its second spatial coefficient:", second_jet.T)
    print("Direct velocity covariance:", covariance_column.T)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
