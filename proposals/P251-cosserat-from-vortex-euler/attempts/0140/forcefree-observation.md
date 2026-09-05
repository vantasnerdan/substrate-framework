# Actual force-free pole and one transported material-angle action

## 1. Fixed geometry, carrier and order of choices

Use0137's smooth compact-vorticity ordinary column, with axial vorticity
Z(r), exactly 2O on r<a, and zero outside b. Its swirl is V=r O(r).
First fix a sufficiently small nonzero carrier k*, a compact interval I
around it, the smooth profile, the axial action measure Lz, and a fixed
optical observation window |t|<=T. Lz is held fixed in the Fourier/Bloch
normalization; it is not replaced by 2pi/k when differentiating k.
All constants below may depend on these fixed choices, including the
small but nonzero tag-moment determinant.

Now use the actual globally stationary generalized force-free column
of0136:

```
D(r)=2 integral_0^r Z(s)V(s)ds,
W(r)=sqrt(U0²-D(r)), h=W-U0,
u=V e_theta+W e_z,
curl u=(0,-W',Z)=(Z/W)u.
```

This is smooth, has compact vorticity, and has no wall or vortex sheet.
Its factor Z/W is not constant; no constant-lambda API or global Runge
theorem is imported. An axial Galilean observation frame gives h in
place of W, with `||h||_C^j<=C_j/U0` for each fixed j. The field itself
is the stationary force-free laboratory field. Put delta=U0^-1.

## 2. The actual force-free Euler pole

Let nu(k) be the eigenfrequency in the axial comoving frame. The exact
Euler transfer system from0137 is

```
s(r,k)=nu(k)-m O(r)-k h(r),
f'=-[1/r+2m O/(r s)]f+(m²/r²+k²)p/s²,
p'=(s²-2O Z)f+2m O p/(r s).                 (1)
```

Here p is pressure divided by density, f is actual radial material
displacement. Axial shear has not been discarded: its terms cancel in
(1) after full incompressibility is used. Reconstructed fields are

```
v_r=-i s f, v_theta=m p/(r s)-Zf, v_z=k p/s-h'f,
eta=(f,iB,iC),
B=[m p/(r s)-2O f]/s, C=k p/s².             (2)
```

On the fixed nonresonant neighborhood selected in0137, the coefficient
and parameter-derivative changes in (1) are O(delta). The origin is a
regular singular endpoint. In the rigid core h is analytic in r²;
the regular solution has p=p0 r^m+O(r^(m+2)) and
`f=m p0 r^(m-1)/[sigma(2O+sigma)]+O(r^(m+1))`, where
`sigma=nu-mO` uses h(0)=0. Factoring these powers gives the usual
regular Frobenius solution, with parameter derivatives continuous in
delta. Integrate to the annulus and match the full decaying exterior
potential exactly as in0137, now with exterior constant h_infinity.

The real Evans determinant therefore differs in C³(nu,k) by O(delta)
from the smooth ordinary-column determinant. Its simple root persists:

```
||nu_delta-nu_ordinary||_C²(I)<=C delta.      (3)
```

This is an actual smooth Euler eigenvelocity, with compact perturbed
vorticity and exponentially decaying exterior velocity in its fixed
nonzero Fourier fiber. It is not a generic isolated-eigenvalue assertion
against arbitrary exterior-vorticity perturbations. The exact radial
matching construction supplies the pole on the specified invariant
compact-vorticity class. Its laboratory frequency is k U0+nu_delta.

The same transfer bounds hold for p,f,B,C and their required k/r
derivatives on every fixed compact core interval. No steep-profile
vorticity-generator operator-norm limit is used. Since the profile is
fixed before U0, the derivative constants are finite.

## 3. Full force-free KKS and positive mode action

The azimuthal component of base vorticity now contributes. For the two
real Fourier columns in (2), direct cross multiplication gives

```
Re eta cross Im eta=(0,-f C,f B),
beta(k)=2pi rho Lz integral_0^b f(ZB+W'C)r dr.       (4)
```

The W'C term is part of the same Euler KKS action, not an optional
small correction. It is O(delta), and all other terms converge to the
positive ordinary-column beta in C²(I). Thus beta remains positive.
The physical SO(2) generator and axial Galilean generator give the
rotating-core/comoving quadratic Hamiltonian
`H''=-beta sigma I2`, where `sigma=nu_delta-mO<0`.
The squared intrinsic frequency and its first k² derivative and second
k derivative retain their positive ordinary-column curvature by (3).
Generally a nonzero first k derivative remains; no standing-carrier
closure is inferred.

## 4. One actual transported tag, without a constant radial Doppler

Choose m=2 for an ordinary covariance orientation, or any m>=2 for
the central m-lobed shape angle. Let the actual material fraction at
t=0 be

```
w0=mu chi(r) chi_z(z)[1+eps b_tag(r)cos(m theta)],
|eps b_tag|<1, 0<mu<=1/2,
```

with chi supported inside r<a, even chi_z, and both between zero and
one. The signed b_tag describes radially opposed lobe orientations;
the material fraction is everywhere nonnegative. The same tag is used
for all times and all k in I. Define

```
Z0=integral chi_z, Zc(k)=integral chi_z cos(kz),
Q=integral chi b_tag r^(m+1)dr !=0,
F(r,k)=chi r^m[p'+m p/r]/[s(2O+s)],
A0(k)=integral F, A1(k)=integral s F,
T_tag(k)=integral chi b_tag r p/s dr.          (5)
```

The tag and I can be chosen with Zc and A0 nonzero. Its base central
shape moment is `rho mu pi eps Z0 Q exp(i m O t)`: axial shear does
not alter this reference transverse moment. Actual material axial
transport is `z(t)=z0+W(r)t`. Using (2), both its angle and complete
axial angular momentum for the prepared sine-angle phase are exactly

```
theta(t,k)=Zc/(eps Z0 Q) integral F(r,k)sin(s(r,k)t)dr,
S(t,k)=rho mu m pi eps Zc
       integral chi b_tag r p/s cos(s(r,k)t)dr.     (6)
```

Spin includes the position variation as well as velocity. Its pointwise
coefficient is `rho r(2O eta_r+v_theta)=rho m p/s`; omitting the first
term would change the result. Axial shear enters these genuine material
observations through both s and the profiles. The other mode phase
gives the corresponding cosine angle and sine spin rows. Consequently
the results below apply to the entire real two-dimensional mode, not
only to one trajectory.

## 5. Exact time-dependent pullback to this physical angle

Factor out the constant reference sigma and define

```
A(t,k)=integral F(r,k)exp[-i k h(r)t]dr,
a_complex=A/A0, delta_phase=arg(a_complex),
c(t,k)=c0(k)|a_complex|,
c0=Zc A0/(eps Z0 Q),
gamma=sigma+partial_t delta_phase,
connection=partial_t log|c|.                 (7)
```

c0 may be signed; it never vanishes. At t=0, a_complex=1. For large
enough U0 it stays away from zero throughout I times [-T,T], so the
phase and logarithm are well defined continuously. In particular
`gamma(0)=A1/A0`, not necessarily sigma, and `connection(0)=0`.

Let y be the real canonical mode coordinates, with `ydot=sigma J y`.
Rotate them by delta_phase, write z=R_delta_phase y, and use the actual
material coordinate `theta=c z2`. Rotating the full one-form and
eliminating z1 gives the **exact scalar physical-angle action**

```
L_theta=M(t,k)/2 [(theta_dot-connection theta)²-gamma² theta²],
M=-beta/[gamma c²]>0,
Pi=M(theta_dot-connection theta).            (8)
```

The positivity follows from gamma<0, retained uniformly at large U0.
Both temporal connection terms are necessary; (8) is not obtained by
putting a measured phase into an unaltered oscillator. For unit prepared
sine phase, with psi=sigma t+delta_phase,

```
theta=c sin psi,   Pi=-beta cos psi/c.        (9)
```

At t=0 this gives `Pi(0)=-beta/c0` and
`theta_dot(0)=c0 A1/A0`. These facts fix the matching target below.

## 6. Four fixed radial controls match the physical momentum jet

Import the ordinary-column rank proof and its first 8/8 receipt from
0138, `radial-jet-rank.md`, `radial_rank_verify.py` and
`radial-rank-first-run.txt`. The four ordinary radial rows are
Q and the first three k jets of `integral chi b_tag r p/s` at k*.
Their independence is proved there from the Bessel recurrence. It is
not reverified here.

In this force-free field s depends on r; dividing by s is NOT an exact
scalar triangular row transformation. Instead fix the four bump
supports furnished by0138, and use the C² profile/Doppler estimates from
(1)--(3). Their actual moment matrix is M0+Delta, with
`||Delta||<=C delta`. Choose U0 such that
`||M0^-1 Delta||<1`. This proves invertibility for the actual field.

The exact physical-action matching target is

```
R(k)=-beta(k) Z0/[rho mu m pi Zc(k)² A0(k)].  (10)
```

Solve the four finite moment equations

```
Q=Q_selected !=0,
partial_k^j T_tag(k*)=Q partial_k^j R(k*), j=0,1,2.
```

The target comes from (6) and the exact canonical momentum (9), so
there is no fitted frequency or introduced inertia. Scale all four
bump coefficients and Q_selected by one common positive factor if
needed to enforce |eps b_tag|<1. The target ratios remain unchanged.
Freeze that resulting SINGLE tag. Then

```
partial_k^j [S(0,k)-Pi(0,k)]_(k*)=0, j=0,1,2.      (11)
```

For the other phase both initial momentum rows vanish. At t=0 the
spin/rate inertia equals the exact action mass because connection=0.
The ordinary-column limit of (10) agrees with0135/0138 after converting
the p/s convention. The four-bump gap can be small (the explicit
ordinary coefficient minor scales as k*^9); k* is fixed before U0.

Append-only correction to0135: its final normalized-pressure series
has denominator 4(m+1), not 2(m+1). The latter is the logarithmic
derivative's coefficient. The corrected pressure recurrence and its
nonzero rank determinant are already checked independently in0138.
This changes no dispersion coefficient or O(k²) conditioning verdict.

## 7. Finite optical window and two carrier derivatives

All estimates here are at fixed geometry, carrier interval, tag and
finite T. Since `s-sigma=-k h`,

```
max_(j<=2) ||partial_k^j(s-sigma)||_infinity <= C delta.
```

Differentiating the actual exponential inside (6), rather than
replacing it by a constant, gives bounds by
`C delta (|t|+t²+|t|³)` for its first two k derivatives relative to
exp(i sigma t). Coefficient derivatives are bounded by the radial
transfer estimates. The fixed positive lower bounds for |Q|, |A0|,
|Zc| and |sigma| then give, for the two-dimensional unit mode amplitudes,

```
max_(j<=2) sup_|t|<=T |partial_k^j(S-Pi)(t,k*)| <= C_T delta.       (12)
```

The proof splits each row into its t=0 coefficient times cos(sigma t)
or sin(sigma t), plus the explicit phase difference. Equation(11)
cancels the t=0 coefficient mismatch and both derivatives. In (9),
c-c0 and delta_phase have the same controlled bounds. No retuning of
the tag with k or t occurs.

Likewise (7)--(8), including the needed time derivatives, imply a C²
carrier/C^j time bound of order delta on every fixed j for the exact
physical-action coefficients relative to the frozen intrinsic scalar
action `M0=-beta/(sigma c0²)`. In particular the connection, the change
of gamma from sigma, and the change of M from M0 are O_T(delta).
The measured equation has those computed connection terms, not an
unqualified autonomous frequency. Its second-carrier spatial action
jet retains the positive intrinsic curvature whenever U0 is chosen
so this coefficient error is below the already fixed positive margin.

These are coefficient/jet error bounds. An O(delta) absolute signal
error is not o(K²) for arbitrarily small macro sideband K at fixed U0.
Consumers either retain the actual zero/first/second coefficients and
their error intervals, or declare a joint limit delta=o(K²). Neither
limit order nor its physical interpretation is silently exchanged.

## 8. Carrier-coherence normalization is chosen before fixing the tag

This theorem uses the single-carrier target (10). More generally,
replace R by eta_weight R for a prescribed positive constant
eta_weight. The same four fixed controls then give (11)--(12) with
S-eta_weight Pi. No extra rank proof or redefinition of the measured
angle is needed. This is a different tag choice, frozen once before
the sideband expansion, not a retuning of the original tag per carrier.

The parent's0138 same-fluid coherent standing pair has carrier weights
1/2: quadratic action scales by 1/4 per carrier while linear measured
spin scales by 1/2. Its matching choice is consequently eta_weight=1/2,
not this theorem's single-carrier eta_weight=1. Time-reversal population
averaging does not remove that coherent-amplitude normalization. The
parent's full standing action determines which target is used.

## 9. Scientific result and next parent construction

Route verdict: **established as stated** for the actual stationary
smooth generalized-force-free column mode, its exact positive
time-dependent physical-angle action, and ONE nonnegative material tag
whose physical spin/action mismatch and first two carrier derivatives
are controlled over the fixed optical window. At t=0 those three
momentum mismatches vanish exactly. The actual radial pressure, full
KKS and every material phase remain in this statement.

This is a controlled stationary-column construction, not an arbitrary
global EPS knot or an autonomous two-field homogenized theorem. The
global ambient extension of the local closed torus remains the parent's
separate geometry problem. Coherent counterpropagating carrier and
translation assembly consumes the explicitly bounded physical jets
rather than assuming equality at all wavelengths or all times.
