"""Actual transport and both-parity source algebra for0243."""
import sympy as s

from substrate_framework import euler_fourier as ef
from substrate_framework.euler_displacement_preparation import finite_displacement_cell
from substrate_framework.euler_joint import passive_output_weights
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger('P251-0243-hybrid-acoustic-control')
    theta, time = s.symbols('theta t', real=True)
    omega, amp, density, ss, aa, weight = s.symbols('omega G mu s1 A w', positive=True)
    cell = finite_displacement_cell()
    scalar = ef.mul(ef.trig(1, 3), ef.trig(2, 2, kind='sin'))
    axial = (scalar, {}, {})
    transported = ef.transport(cell.background, axial)
    checks.check('Eq3: actual three-component C016 passive sector has zero pressure',
                 not any(ef.transport(axial, cell.background))
                 and not ef.divergence(transported) and ef.leray(transported) == transported)
    for parity in (0, 1):
        initial = amp*(s.sin(theta) if parity == 0 else s.cos(theta))
        velocity = initial.subs(theta, theta-omega*time)
        displacement = time*velocity
        checks.check(f'Eq3 parity{parity}: exact Euler transport and zero-initial Lin history',
                     s.simplify(velocity.diff(time)+omega*velocity.diff(theta)) == 0
                     and s.simplify(displacement.diff(time)+omega*displacement.diff(theta)-velocity) == 0
                     and displacement.subs(time, 0) == 0)
        angular = s.integrate(ss*s.sin(theta)*velocity, (theta, 0, 2*s.pi))/(2*s.pi)
        output = s.simplify(aa*density*angular/5)
        waveform = s.cos(omega*time) if parity == 0 else s.sin(omega*time)
        checks.check(f'Eq6 parity{parity}: physical angular integral yields the prescribed quadrature',
                     s.simplify(output-amp*aa*density*ss*waveform/10) == 0)
        checks.check(f'Eq7 parity{parity}: full source norm remains quadratic in signed amplitude',
                     s.simplify(s.integrate(initial**2,(theta,0,2*s.pi))/(2*s.pi)-amp**2/2) == 0)
        ws = s.symbols('w1:4', positive=True)
        moments = s.Matrix([[w**(2*j+parity) for w in ws] for j in range(3)])
        determinant = s.prod(w**parity for w in ws)*s.prod(
            ws[b]**2-ws[a]**2 for a in range(3) for b in range(a+1,3))
        checks.check(f'Eq8 parity{parity}: determinant is the signed positive-band Vandermonde',
                     s.factor(moments.det()-determinant) == 0)
        # Actual even target has a nonzero constant; the old odd-only family cannot supply it.
        solved = passive_output_weights(moments.subs(dict(zip(ws,[1,2,3]))), [1,-2,3], parity=parity)
        response = sum(w*(s.cos(n*time) if parity == 0 else s.sin(n*time))
                       for w,n in zip(solved.weights,[1,2,3]))
        checks.check(f'Eq8 parity{parity}: source output derivatives reproduce the finite target',
                     all(s.simplify(response.diff(time,2*j+parity).subs(time,0)
                                    -s.factorial(2*j+parity)*c) == 0
                         for j,c in enumerate([1,-2,3])))
    arbitrary_weights = s.symbols('a:3')
    old = sum(w*s.sin(n*time) for w,n in zip(arbitrary_weights,[1,2,3]))
    checks.check('odd-only negative control cannot repair a nonzero even initial acceleration',
                 s.simplify(old.subs(time,0)-1) != 0)
    k = s.Symbol('k', real=True)
    aa1, bb0, xx0 = s.symbols('A1 B0 X0')
    # Difference formula: delta A = O(k^2), delta B=O(k), delta X has no zero jet.
    dA2, dB1 = s.symbols('delta_A2 delta_B1')
    change = -s.I*k*(k**2*dA2)-k**2*(k*dB1)/2
    checks.check('Eq9: zero first tag response makes the full hybrid correction cubic',
                 s.series(change,k,0,3).removeO() == 0)
    checks.check('retaining a nonzero first tag moment exposes a second-order hybrid defect',
                 s.expand(change-s.I*k**2*aa1).coeff(k,2) == -s.I*aa1)
    return checks.finish()


if __name__ == '__main__':
    raise SystemExit(main())
