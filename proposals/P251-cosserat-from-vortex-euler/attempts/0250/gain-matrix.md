# Thin-ring physical gain matrix and vector adjoint criterion

This continuation evaluates the leading physical moment rows of
`isovortical-channel.md` on the exact 0211 axisymmetric annulus. It preserves
the whole-space Leray evolution and separates the finite local gain question
from arbitrary vector-history density.

## 1. Three explicit helical source branches

Use the 0217 thin-ring orientation

    r=R+s cos(theta),       z=-s sin(theta),
    u=V(s)e_theta+W(s)e_varphi,       V=s Omega_p,          (1)

and the complex transverse laboratory component `Y_+=Y_x+iY_y`. Fix the
toroidal number `n=-1`; it is nonzero, so the cohomological construction in
`isovortical-channel.md` applies after removing finitely many action
resonances. Consider the three poloidal harmonics

    g_m=exp[i(m theta-varphi)],       m=0,+1,-1.             (2)

Their leading laboratory frequencies are

    nu_0=-W/R,
    nu_+= V/s-W/R,
    nu_-=-V/s-W/R.                                          (3)

On any compact core annulus away from `s=0`, each exact action-angle
frequency is a smooth perturbation of (3). Since `W'=-lambda s Omega_p` and
the 0211 rotation number has strict twist, smaller open action intervals give
nonconstant bands separated from zero and from the finite cohomological
resonances. Every branch has both real phase quadratures.

For the fixed radial tag ladder, divide out the common nonzero resonant
factor `rho epsilon N_j c_j A J/(2 lambda)` and harmless nonzero phase factors.
The complex centroid numerator uses `F_X=x+i y=r exp(i varphi)`. The physical
0217 covariance tilt uses

    F_theta=(x+i y)z=r z exp(i varphi),
    theta_+=delta M_Ftheta/(i D_tag),       D_tag>0.         (4)

Direct angular integration of (1)--(2) gives the three centroid/tilt rows,
up to one common positive factor `pi`,

    X:      (2R,              s,                s),
    theta:  (0,        -R s/D_tag,        R s/D_tag).       (5)

Changing the fixed real phase reverses the displayed signs but not their
rank. In particular `m=0` is centroid-dominant and angle-null, while the two
first poloidal harmonics give opposite nonzero transverse covariance tilts.

The cohomological solution has
`beta=N f/(lambda nu_m)+O(1)`. In the same orientation the transverse
angular-displacement kernel is

    (x cross u)_+
      =exp(i varphi)[-W z+iV(r cos(theta)-z sin(theta))].    (6)

Its three angular coefficients give the scaled leading `G_+` row

    G: (-g_0,-g_+,-g_-),                                   (7)

where

    g_0=2Vs/nu_0,
    g_+=(VR+Ws)/nu_+,
    g_-=(VR-Ws)/nu_-.                                      (8)

The signs in (5)--(8) use one common column phase. The literal current is not
silently identified with `G_t`. For each exact finite carrier,

    S_+=-i nu_m G_+
          +2rho integral chi (xi cross u)_+ dx
          + lower full-pressure sidebands,                 (9)

which is the exact identity `G_t-S=-2rho integral chi xi cross u` rewritten
for the leading clock. Since `xi cross u=a f/lambda`, the displayed correction
is `O(c_j)` while the resonant `G` term is `O(N_j c_j)`. Thus
`S_+/G_+ ->-i nu_m` after the fixed-tag normalization. Exact compatibility
`S_+=-i nu_m G_+` is earned only by a same-frequency source combination that
annihilates the explicit correction in (9); it is not assumed here.

## 2. The finite local matrix is nondegenerate

The three leading direction columns for `(X_+,theta_+,G_+)` are therefore

    M_s=[[2R,             s,             s],
         [ 0,    -R s/D_tag,     R s/D_tag],
         [-g_0,          -g_+,          -g_-]].             (10)

Substitution of (3) into (8) yields

    det M_s=
      4 R^2 V s^3 [R^2 V^2+2R^2 W^2-W^2 s^2]
      /[D_tag W(R^2 V^2-W^2 s^2)].                         (11)

The omitted common phase and `pi` factors are nonzero. On the positive 0211
core, `V,W,D_tag>0`; choosing one sufficiently large but finite ring with
`RV>Ws` makes (11) strictly positive. Smooth dependence of the exact
action-angle chart, exact finite-`R` moments and the cohomological inverse
then preserves this rank on smaller action neighborhoods. This identifies
three actual source harmonics/action branches with independent finite local
centroid, transverse angle and angular-displacement/current directions.

Equation (11) is not a pointwise same-frequency matrix: its three columns
have the distinct bands (3). Whole-law rotations can supply laboratory
vector components, but rotate each column's angle and `G` together and hence
do not change this distinction.

## 3. Exact same-frequency angle/G inverse from positive fractions

The accepted C-CST-017 input includes the already reviewed two fixed positive
tag fractions. They apply without changing the 0211 field or source. Let
`chi_sigma=f_sigma chi`, with fixed `0<f_1!=f_2<1`, and retain the continuous
ambient complement `1-f_sigma chi`. For any actual displacement history,

    theta_sigma=theta,
    G_sigma=f_sigma G,
    S_sigma=f_sigma S.                                    (12)

These identities are exact at finite carrier and with full pressure: the
reference covariance and its variation both scale by `f_sigma`, whereas the
unnormalized material currents scale once.

Choose one of the `m=+1` or `m=-1` branches on a smaller band where its exact
finite-carrier angle gain and `G` gain are smooth and nonzero. Denote their
ratio by `j_N(omega)=Gamma_G/Gamma_theta`. For two signed physical source
profiles `F_1,F_2`, averaged with equal positive species probabilities, the
actual leading spectral rows are

    A=(F_1+F_2)/2,
    B=j_N(omega)(f_1 F_1+f_2 F_2)/2.                      (13)

The determinant is `j_N(omega)(f_2-f_1)/4`, uniformly separated from zero
on a smaller compact band. Consequently

    F_1=2[f_2 A-B/j_N]/(f_2-f_1),
    F_2=2[B/j_N-f_1 A]/(f_2-f_1)                          (14)

sets the angle and `G` spectral profiles independently at the same frequency.
All derivatives of the actual `j_N`, its inverse and the full finite source
norm enter the 0250 diagonal. The probabilities, material fractions and
ambient density remain positive; only the linear initial amplitudes are
signed. Thus the three-harmonic determinant (11) is a useful independent
check and a one-tag fallback, but angle/`G` separation itself is closed by
(12)--(14).

The spin row is preserved rather than promoted to an exact oscillator
identity. Write the exact same-branch ratio as

    Gamma_S/Gamma_G=-i omega+delta_N(omega).                (15)

Equation (9) and the all-order pressure expansion give
`delta_N=O(1/N)` with every fixed frequency/time derivative on the selected
band: the `G` gain is `O(N c_N)`, while the material correction is `O(c_N)`.
The inverse (14) therefore realizes `S=G_t+o(1)` after normalization. This is
the compatibility needed when the target spin is the derivative of its
chosen current primitive. Literal finite-`N` equality requires an additional
same-frequency column annihilating the correction in (9); it is not claimed.

The 0257 isotropic reconstruction acts after these per-realization physical
rows. It gives a condition-one full-vector angle map and uniformly rotates
the `G/S` costs, but does not supply the inverse (14); the two positive
fractions are what separate those outputs.

## 4. Vector adjoint criterion for any enlarged output family

Let `Gamma_m(omega)` denote the exact finite-`R`, full-pressure gain vector
for branch `m`, including every retained physical row, parametrized by its
exact action-frequency inverse on an open band `J_m`. For a compactly
supported vector distribution `mu` on the requested time window, write its
componentwise entire Fourier transform as `mu_hat(z)`. The exact annihilator
criterion for the combined source family is

    Gamma_m(omega)^T mu_hat(omega)=0
       for every omega in J_m and m=0,+1,-1.                (16)

The vector-valued output span is dense precisely when (16) implies
`mu_hat identically 0`. A full-rank gain matrix at every common frequency is
one sufficient proof, but it is not necessary. In the present construction
the bands need not overlap, so (11) by itself does not prove (16). Narrowing
the bands and treating `Gamma_m` as constant is also insufficient without an
inverse-cost estimate: the finite exponential coefficient sum can grow
faster than the direction variation shrinks.

A sufficient analytic continuation repair is explicit. If each rescaled
`Gamma_m(z)` has a common connected meromorphic continuation, clear its
denominators in (16). The three resulting entire row identities then hold for
all complex `z`. If their determinant is not identically zero—(11) supplies
the thin-ring nonzero comparison—then `mu_hat=0`. The remaining exact task is
to construct that common continuation for any output columns not covered by
the exact inverse (14), or prove (16) directly by a vector adjoint uniqueness
argument. Finite Taylor
rank or continuity of (11) does not replace this step.

## 5. Ambient hybrid row

The centroid row in (5) is an actual material moment, but it is not the full
0241 point-to-hybrid coefficient. For an exact source `xi_m,v_m`, define the
additional row without approximation by

    Gamma_H,m(omega)=coefficient of K tensor K in
      X_pt,tt-rho^-1 partial_t Delta J[xi_m,v_m],            (17)

where `Delta J` includes the complete first momentum, shape rate, second
moment and moving-centroid phase of 0241. The continuous ambient and the
whole-space pressure in `v_m=P(xi_m cross omega)` are part of (13). Every
term is a finite action-angle Fourier integral and can be appended as a
fourth component of `Gamma_m` in (16).

One diagnostic piece of (17) can be evaluated for the `m=0,n=-1` branch.
With `k=I'(s)e_s`, its leading velocity component is

    v_+=exp(iNI)[-W sin(theta)-iV]/I'.                     (18)

For a macroscopic wave vector parallel to the ring axis, the
Eulerian-velocity contribution to the second moment in 0241 contains

    B_(+,zz)=integral rho chi v_+ z^2 dx
      =-i C_N integral A(s)chi_res(s)
                    [V(s)/I'(s)] R s^3 ds+O(a/R)+O(1/N),   (19)

where `C_N>0` contains the explicit nonzero angular factor and the conjugate
fixed-tag coefficient. The `W sin^3(theta)` integral vanishes, while the
`V sin^2(theta)` integral has one sign. Choose the radial envelope where
`V/I'` and the resonant tag coefficient have fixed sign; then (19) is nonzero
on one sufficiently large fixed ring and for every sufficiently large finite
carrier. The whole-space pressure correction is the displayed `O(1/N)`
relative WKB row, not omitted.

For this axial contraction, the source has `delta X_z=0` by its toroidal
harmonic, so the moving-centroid phase in 0241 does not cancel (19).  But the
literal material variation of `B` also contains the moved-domain and
background-velocity terms in equation (8) of `axisymmetric-hybrid.md`.
Equation (19) alone therefore does not prove the full coefficient; calling it
the complete material row would be the velocity-only mutation exposed by the
new symbolic check.

The quadratic form behind (19) is not special to axial `K`. Before the final
theta integration, the `V` contribution is proportional to

    pi r^3 |K_perp|^2+2pi r z^2 |K_parallel|^2,            (20)

and is positive for nonzero `K`; the odd `W` contribution integrates to zero
in the circular thin-ring model.  This identifies the sign of one Eulerian
contribution, not the complete covariant hybrid tensor.  The full calculation
is instead carried out for the axisymmetric `n=0,m=+/-1` pair in
`axisymmetric-hybrid.md`.  There the Lagrangian variation retains `D_t xi`,
moved positions, centroid phase and the untagged point row.  Its
same-frequency `G`-cancelled column has nonzero parallel and transverse hybrid
gains of one sign, and its frequency-dependent whole-law inverse supplies the
full acoustic vector.  An untagged ambient branch would be useful only after
its own actual `Gamma_H` is shown nonzero; disjoint support alone does not
create that gain.

## 6. Strongest scope

`route_verdict: full-vector covariance angle and same-frequency angle/G
separation established here; axisymmetric-hybrid.md establishes the complete
acoustic hybrid band and its reflected common-band calculation closes the
coupled fixed-ring physical observation density`

`evidence_scope: explicit 0211 thin-ring harmonics m=0,+1,-1, nonzero leading
determinant, exact finite-N spin correction, and necessary-and-sufficient
vector annihilator criterion with full pressure retained; exact two-positive-
fraction angle/G inverse, asymptotic full-spin compatibility, and an explicit
nonzero Eulerian contribution whose velocity-only limitation is exposed`

`axisymmetric-hybrid.md` supplies the minimum repair: reflected O(3) parity
cancels the forbidden even cross blocks, retains the explicit odd curl/current
rows and gives an overlapping-band finite-`K` Schur inverse.  The completed
0257 transaction supplies the longitudinal covariance reconstruction with
uniform rotated costs. Full KKS/Jacobi normalizer columns remain separate and
are appended only after that physical gain family is fixed. None of these
local results changes the geometry scope or closes the parent objective.
