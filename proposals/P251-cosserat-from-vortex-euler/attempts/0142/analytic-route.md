# Analytic route receipt: trapping, physical observation, and a repair

The parent registered0142 and validated263/12 before opening0140 or
new source bodies. This remains candidate construction, not a finished
optical theorem. All calculations below are exact or analytic asymptotic
ones; no radial box or floating-point spectrum is being selected.

Use U,lambda>0, Omega=lambda U/2, negative axial carrier k=-p, p>0,
delta=sqrt(lambda/p), and m>=2. The Lundquist column has

    O(r)=U J1(lambda r)/r,
    Z(r)=lambda U J0(lambda r), W(r)=U J0(lambda r).

Its exact radial transfer equations are0140(1), with no exterior
vorticity cutoff. In particular s=omega-mO-kW and the full radial,
azimuthal and axial pressure equations are retained.

## 1. Counterpropagating axial-Doppler trap

For p/lambda large, the Doppler function mO-pW has its global minimum
at r=0. Indeed 1-J0(x)>0 for x>0, and the continuous ratio
[1/2-J1(x)/x]/[1-J0(x)] is bounded on (0,infinity), with a finite
limit1/4 at zero and bounded limit behavior at infinity. The p term
dominates that fixed m ratio. A frequency below mOmega-pU therefore
has no critical radius s=0 anywhere.

Write omega=kU+mOmega+sigma, sigma=-2Omega+Omega delta c. The
axis expansion gives

    s=-2Omega+Omega delta[c-R²/sqrt(2)]+higher terms,
    r=ell R, ell=(2/(lambda p³))^(1/4).

Keeping the full incompressibility and pressure elimination, the leading
radial oscillator is

    [-partial_R²-R^-1 partial_R+(m-1)²/R²+R²]f=sqrt(2)c f.

Its regular decaying radial states are

    f_n=R^(m-1)exp(-R²/2)L_n^(m-1)(R²),
    c_n=sqrt(2)(2n+m), n=0,1,... .

Thus the leading intrinsic branch is
sigma=-2Omega+sqrt(2)Omega(2n+m)delta+O(delta²). Its intrinsic
sigma² second carrier derivative is negative. That observation is not
silently replaced by a favorable Floquet winding.

The pressure formulation's apparent near-axis small denominator
m²+k²r² is best handled in Cartesian form, not by discarding m²/r².
In an axial Fourier fiber, the complete pressure inverse is
(p²-Delta_perp)^-1. Under r=ell R its dimensionless expansion parameter
is 1/(p²ell²)=delta/sqrt(2); the exact resolvent supplies a bounded
remainder on Schwartz functions at every finite expansion order.
This gives a route to globally smooth Gaussian/Laguerre quasimodes,
without an artificial radial wall or an assumed isolated pole.

## 2. Why actual material-angle curvature differs from sigma curvature

For n=0 and a fixed physical tag weight chi(r) proportional to r^(2l)
on the packet scale, the angle row has radial weight
r^(2m+2l-1)exp(-g r²/2), g=ell^-2. Its transported axial phase is
exp[-i a p r² t], a=U lambda²/4. Put d=m+l and cD=sqrt(2)Omega delta.
The exact leading normalized moment is

    A(t,p)/A(0,p)=(1+i cD t)^(-d).

Therefore the physical angle phase and amplitude connection are

    gamma_tag=sigma-d cD/(1+cD²t²),
    connection=-d cD²t/(1+cD²t²).

In particular gamma_tag=-2Omega-l cD+O(delta²) in the fixed optical
window. For l>0 its squared-frequency carrier curvature is positive.
The tag is an actual fixed nonnegative material fraction, not a selected
logarithm branch. Its full time-dependent one-form pullback must retain
both displayed connections and its mass, as in0140(8).

## 3. A named ground-packet momentum/jet mismatch

The leading ground-packet spin row has the form

    T(p)=const*Omega delta integral chi b(r)r^(m+1)exp(-g r²/2)dr.

Its action matching target beta/A0 scales as g^l. Since g'/g=3/(2p),
matching both T and its carrier derivative forces the signed spin-row
mean of g r² to equal -2/3-2l. But the angle-row mean is 2(m+l).
Matching the transported spin to the same canonical momentum through
first order in delta would require those means to agree. In this
monomial ground-state class that forces l=-m/2-1/6, which has the wrong
leading physical curvature. Extra copies of the same ground-state
radial rows cannot fix their exact linear dependence.

This refutes that particular simultaneous ground-packet jet/time matching
route, not trapped optical packets or the positive physical-angle action
itself. It motivates the following append-only candidate expansion.

## 4. Failure-derived radial excited packets

For n>=1, the pressure/spin row carries the Laguerre polynomial

    P_n(x)=L_n^(m-1)(x)-2 partial_x L_n^(m-1)(x), x=g r².

Its carrier derivatives act by -1/2+(3/2)x partial_x on
exp(-x/2)P_n(x). The time-Doppler row is multiplication by x. Register
the n=1 and n=2 states before evaluating their rank. The candidate is
to match three carrier jets of the spin row and three jets of its first
time-Doppler moment, with fixed radial marker controls plus the nonzero
reference-shape row. Six independent polynomial-exponential rows would
permit this; a dependence must instead be compared with the exact
canonical targets. This rank calculation is the next analytic step,
followed by the full Euler residual and fixed-tag error bounds.
