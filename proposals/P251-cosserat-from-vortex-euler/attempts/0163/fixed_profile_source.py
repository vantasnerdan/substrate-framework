"""Exact source/projection for the registered fixed-profile corrector route."""

import sympy as s

from substrate_framework import euler_fourier as f


def calculate(power):
    half = s.Rational(1, 2)
    psi = f.add(f.trig(0), f.trig(1),
                {(1, 1, 0): half, (-1, -1, 0): half})
    axial = {f.ZERO: s.Integer(1)}
    for _ in range(power):
        axial = f.mul(axial, psi)
    axial_prime = {f.ZERO: s.Integer(power)}
    for _ in range(power-1):
        axial_prime = f.mul(axial_prime, psi)

    def dx(field):
        return f.add(f.derivative(field, 0), f.scale(f.derivative(field, 1), -half))

    def dy(field):
        return f.scale(f.derivative(field, 1), s.sqrt(3)/2)

    def green(field):
        return {wave: value/s.Integer(wave[0]**2+wave[1]**2-wave[0]*wave[1])
                for wave, value in field.items() if wave != f.ZERO}

    grad_psi = (dx(psi), dy(psi))
    grad_axial_x = (dx(green(dx(axial))), dy(green(dx(axial))))
    speed_square = f.add(*(f.mul(entry, entry) for entry in grad_psi))
    # lambda=Psi=1: the exact general-profile F_x, not a fitted drift.
    source = f.add(*(f.mul(grad_psi[j], grad_axial_x[j]) for j in range(2)),
                   f.mul(f.add(axial, f.scale(f.mul(psi, axial_prime), -1)), dx(psi)),
                   f.scale(f.mul(axial_prime, dx(speed_square)), -half))

    def integral_exponential(n):
        return s.pi if n == 0 else (s.Integer(-1)**n-1)/(s.I*n)

    def triangle_mode(n, m):
        if m:
            return s.Integer(-1)**m*(integral_exponential(n)
                                     -integral_exponential(n-m))/(s.I*m)
        if not n:
            return s.pi**2/2
        return s.pi*s.Integer(-1)**n/(s.I*n)+(s.Integer(-1)**n-1)/s.Integer(n*n)

    upper = s.simplify(sum(value*triangle_mode(wave[0], wave[1])
                           for wave, value in source.items()))
    raw_x_row = s.simplify(sum(value*s.Integer(-1)**wave[0]/(s.I*wave[0])
                               for wave, value in source.items()
                               if wave[1] == 0 and wave[0] != 0))
    drift = s.simplify(-raw_x_row+upper/(3*s.pi))
    print(f"Axial profile W=psi^{power}, lambda=Psi=1")
    print("Computed first-order horizontal source modes:", len(source))
    print("Uncentered <a F_x>:", raw_x_row)
    print("Upper polygon integral F_x:", upper)
    print("Computed conserved-row axial frame nu:", drift)
    print("Full cell mean W:", axial.get(f.ZERO, 0))
    print("Residual advection mean(W)-nu:", s.simplify(axial.get(f.ZERO, 0)-drift))


def main():
    for power in (1, 3, 4):
        calculate(power)


if __name__ == "__main__":
    main()
