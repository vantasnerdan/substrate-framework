# Attempt 0031: repaired Rankine branch oracle

Parent objective: P251 / issue #198, unchanged. Active obligation: N2.
Route: derive momentum from Cartesian Euler (0029), eliminate amplitudes,
differentiate the Bessel series rather than transcribing its derivative, then
solve the exact scalar pressure-matching equation.

## Analytic specification

With x=ka, s=omega/Omega-m, D=4-s^2, lambda=x sqrt(D)/|s|,
Jratio=lambda J_m'(lambda)/J_m(lambda), Kratio=x K_m'(x)/K_m(x),
the pressure boundary equation is F=D*Kratio+s^2*Jratio-2*m*s=0.
It follows by eliminating the interior pressure and exterior velocity-potential
amplitudes from radial kinematics and pressure continuity. Kratio is the
potential logarithmic derivative, not the pressure logarithmic derivative:
the latter also differentiates the radially varying exterior Doppler factor.
The tangential jump follows from these equations and
the azimuthal momentum equation; it is a dependent check.

For m=1, Jratio=1-3x^2/4+O(x^4 log x), and differentiating
K1=1/x+x/2*(log(x/2)+gamma-1/2)+... gives
K1'=-1/x^2+(log(x/2)+gamma+1/2)/2+.... The old symbolic
derivative omitted 1/4, manufacturing its disputed subleading constant.
For m=2, Jratio=2-x^2/2+... and Kratio=-2-x^2/2+....
The exact expansion supplies the frequency coefficients; numerical roots
check the remainder and branch identity, not select the expansion.

## Numerical remainder and small-ratio prescriptions

Measure convergence of the exact Bessel roots to the derived small-x
expansions on x=10^-2, 10^-3, 10^-4. Use dimensionless mpmath arithmetic at
40 and 60 digits, signed residuals, and a bracket around the physical
s=-1 branch. The numerical output has numeric-evidence scope.
Error budget: print |F|, precision-ladder frequency discrepancy and observed
asymptotic remainder. There is no mesh, domain truncation or quadrature:
regularity and exterior decay are encoded in exact Bessel functions. Translation
neutrality is checked analytically at x=0. Eigenpair residual and Morse-index
tests do not apply to this scalar dispersion root. No tiny difference is
evaluated in machine float arithmetic. Wrong derivative constant and reversed
Coriolis operator are exact mutation tests, independent of tolerance.

## Continuation

This repairs the N2 mode oracle and its constant. It does not by itself derive
a twist modulus, core inertia or the Euler-to-Cosserat action. Those remain
active alongside the angle-dependent patch interaction in attempt 0030.
