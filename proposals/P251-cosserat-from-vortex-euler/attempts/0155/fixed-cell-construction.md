# Finite-time calibrated optical construction in the actual fixed periodic cell

This is the calibrated, odd-lobed candidate of0155. The observation is
the actual strain-tensor-calibrated material moment(4) in
`physical-observer-repair.md`, not an absolute Cosserat director.
The following construction is on the physical periodic Euler cell.
A Bloch continuation and a whole-space finite packet have different
normalizations, identified explicitly in section6.

## 1. Full pressure operator and the executed second-order repair

Normalize time by Omega and use scaled transverse coordinates R=X/ell,
delta=sqrt(T/p), T=trG, ell^4=T/p³. After removing only the axial
translation pW(0), the exact linear/quadratic core pressure equation is

    (p²-Delta_G)P=2Omega curl_X V+2ip grad_X W dot V.

It follows by taking the divergence of the TWO transverse equations
and substituting v_s=i div_X V/p into the longitudinal equation.
In particular no longitudinal pressure or Leray component is omitted.
With derivatives in R, the transverse generator through second order is

    L=L0+delta L1+delta² L2+O(delta³),
    L0=partial_theta+J,
    L1 V=i R² V/2-(2/T)G grad curl V,
    L2 V=(2i/T)G grad(R dot V)
                  -(2/T²)G grad Delta_G curl V.             (1)

The remaining cubic transverse jet of the ACTUAL cosine field starts
at delta³. Its axial quartic term starts at delta4. All coefficients
are analytic with fixed d; these powers follow from ell²=O_d(delta³)
and p ell²=delta, not from dropping the strain.

Let D_+=partial_X+i partial_Y, D_-=partial_X-i partial_Y and
DeltaG=G11-G22. The possibly resonant opposite polarization at angular
index l-2 receives at second order the following two contributions:

    from L1 V1: +i DeltaG²/(8T²) D_-² Delta F,
    from L2 V0: -i DeltaG²/(8T²) D_-² Delta F.              (2)

They cancel identically for EVERY smooth radial profile, since
D_-D_+=Delta. The term G grad(R dot V) has no such angular component.
This is why the earlier possible second-order spectral obstruction
does not occur. The explicit rational-metric polynomial calculation
in `verify_continuation.py` exposes the same cancellation. It is not
used as evidence of an isolated spectrum or all-order decoupling.

Choose m=5,n=8. The first pressure eigenvalue is -(2n+m)delta;
the nonresonant V1 is(6) of the observer receipt. Solve the second-
order nonresonant equations by the nonzero integer L0 denominators.
Retain the order-delta² resonant-plus source rather than claim a
localized spectral inverse for it. Its angular index remains l=m-1.
Thus the constructed two real initial data have:

    angle-relevant residual O(delta²) in the plus l sector;
    all other residual sectors O(delta³).                  (3)

These are actual velocity/pressure/Lin residuals. Preparing the data
by the exact Kelvin push v=P(xi cross omega0), with the pressure
corrections to xi included through the same order, changes neither
power in(3). One can implement this preparation by the exact
divergence-free vector potential of the formal Lin displacement;
the remainder is estimated, not equated to zero.

## 2. Why(3) suffices for the actual spin observable

The principal axial material-spin density of a plus-l displacement
has angular index l-1=m-2. The pure-m tag misses it exactly. The same
tag sees the first pressure torque at index m, which is order delta.
Therefore the retained plus-l residual in(3) changes the physical
angle by O(delta²) but changes its observed axial spin only after
another pressure order: O(delta³). Nonresonant residuals already have
that order. The exact material spin here is

    S=rho integral_tag [r cross xi_t+2 xi cross u0]_s dV,

including the moving position and actual centered reference parcel.
It is never replaced by a KKS coordinate before its moment solve.
The pressure torque formula(9) in the observer receipt proves this
angular selection directly. The actual nonlinear core changes the
selection only at order delta³, already covered by the remainder.

The distinction is important: a raw velocity residual of order delta²
would ordinarily be too large relative to the order-delta spin. Here
its explicitly computed angular type is what gains the extra order.
An unspecified small residual would not license the conclusion.

## 3. Fixed physical tag and the eight-row solve

Use a smooth nonnegative label density with an m-lobed modulation,
radial cutoff chi proportional to R² over the fixed scaled mode region,
and an axial label cutoff of width c/p_star. Here p_star is a fixed
large INTEGER and c>0 is fixed. The initial angular/axial modulation
is cos(m theta-p_star s), with a bounded signed radial amplitude b.
The complete label density, not b alone, is nonnegative. Smooth radial
cutoffs and finite bump supports are chosen before p_star.

Let Q0 be its calibrated reference moment. In action-normalized units
the fixed-cell spin occupies an axial fraction O(1/p_star), and starts
at order delta. Consequently its matching reference moment is of size
Q0=O(delta/p_star)=O_d(delta³), not O(delta). This normalization keeps
the full cell KKS; it does not silently truncate the action to the tag.

Besides the reference moment impose one vanishing higher reference
moment, proportional to integral b R^(m+2) chi dR². It cancels the
entire first cubic-core correction to the reference moment, including
its angular harmonics. The remaining reference error is O(delta6);
division by Q0 gives O(delta³). Because m=5 is odd, the unmodulated
reference odd moment vanishes by inversion symmetry. The exact centered
moment and its centroid correction are used in the finite-dimensional
map; their first possible reference correction is also covered by this
order. There is no assertion that every displacement/current moment
of the tag is zero.

The leading spin functional is e^(-x/2)P_n(x),
P_n=L_n^(m-1)-2(L_n^(m-1))', against the positive radial reference
measure. Fit its value and first slow-time coefficient, together with
their first two relative carrier jets. These are the six rows

    P, DP, D²P, xP, D(xP), D²(xP),
    D=c0+(3/2)x partial_x,

where D acts also on the Gaussian. The scalar c0 is the actual
dimensional prefactor. Changing c0 only performs a triangular row
operation, so it does not change the rank. Together with reference
rows1,x, the exact zero-point Wronskian for n=8,m=5 is

    471442749716799919349866575/65536 !=0.                  (4)

Thus eight sufficiently small disjoint radial bumps give an invertible
matrix, with finite constants fixed before p_star. The exact normalized
map includes the second-order pressure spin, the Kelvin corrections,
the true KKS, the finite axial filter, and the calibrated angle row.
After dividing the spin rows by their explicit nonzero order-delta
factor, this map is an O_d(delta) perturbation of(4). The ordinary
finite-dimensional implicit-function theorem provides b. Its common
amplitude can be reduced together with Q0 so the full label density
remains nonnegative. No frequency or inertia is fitted: the data in
this moment solve are the actual action and observation integrals.

The target may be eta=1 or a separately declared standing-pair eta=1/2.
Unlike the stationary m2 class, the odd reference moment is independently
tunable. Matching the first slow-time coefficient leaves a relative
O(delta²) error over any fixed Omega T_time. Matching the two carrier
jets is necessary; an initial-value match alone is not this theorem.

## 4. Actual finite-time Euler control with the nonlocal pressure

The exact field is x-independent and periodic. Its axial p fiber has
the pressure inverse (p²-Delta_perp)^(-1) on the transverse torus, not
a pinned core wall. Its periodic Green kernel is a sum of modified-
Helmholtz kernels. Away from the source each derivative is bounded by
a polynomial in p and distance times exp(-p distance). Retain all
images and the complete outside velocity when computing KKS and spin.

A uniform local estimate follows without using local C^k closeness as
a substitute for global dynamics. In the invariant elliptic core put
w=exp(sigma sqrt(E)/ell), smoothly flattened as a function of E before
the separatrix. Then u_perp dot grad w=0 exactly and
|grad log w|<=C/ell. Since p ell tends to infinity, Schur's kernel
bound for the conjugated pressure inverse is uniform:

    exp(|log w(x)-log w(y)|) exp(-p dist(x,y))
       <=exp[-(p-C/ell)dist(x,y)].                         (5)

The usual pressure-energy integration gives a weighted fiber estimate
with exponent C_d,T independent of p; the additional commutator cost
is bounded by C/(p ell). Polynomial R moments are controlled with
nested exponential weights. Relative carrier derivatives p partial_p
add the centered transport factor p[W-W(0)]=O_d(delta R²) in this region.
Outside it the fixed positive weighted distance supplies an exponential
tail. Two such derivatives therefore satisfy the same fixed-time
estimate using finitely many nested weights. This proves a local/full-
pressure estimate, rather than multiplying a local error by a growing
uncontrolled global norm. The global background and all its derivatives
are fixed bounded periodic functions.

Apply this estimate to the explicitly typed residual(3), then use the
exact Lin transport to follow the material labels. The angular selection
of section2 and the eight-row solve give uniformly on fixed Omega T_time
and through the two relative carrier jets:

    relative angle/action error O_d(delta²),
    S=eta Pi_theta + O_d(delta²) on the same action scale.  (6)

The constants depend on the fixed d, cutoffs, time, eight bumps and
Gaussian tail margin. Choose those first, then p_star. An exact global
Euler solution, not just the truncated polynomial equations, supplies
the history in(6). The result is a controlled finite-time statement;
it is not an invariant finite-dimensional nonlinear Euler ansatz.

## 5. Positive physical action and its natural carrier scale

Pull the same two-dimensional Kelvin phase form into the actual
calibrated angle chart, retaining the time-dependent amplitude
connection. Its KKS sign is negative; the phase rate is

    gamma=2Omega+(5/21)Omega delta+O_d(delta²)>0,
    M=-beta/[gamma |c_theta|²]>0.

The scalar action is

    L=M/2[(theta_dot-(partial_t log|c_theta|)theta)²
                                                    -gamma² theta²].

All parameter derivatives act on c_theta and M as well. The actual
observation-clock curvature satisfies

    p² partial_p² gamma²=(5/7)Omega² delta+O_d(delta²)>0.    (7)

The error bound is compared to Omega² delta, not to an order-one
frequency. The observation connection is part of the action; equation
(7) is not by itself a statement that an autonomous scalar PDE has
that gap. This is the same calibrated, physically measured shape
throughout, with literal spin supplied by(6).

## 6. Precisely what the fixed-cell license does and does not supply

At integer p_star the field and the two prepared perturbations are
smooth periodic finite-action data on the physical three-torus. KKS
is integrated over its full volume. Continuous neighboring carrier
values are genuine Bloch fibers; their quadratic cell averages and
their local material observations have the two jets used above. They
are not individually periodic fields on the original torus.

One may form a Bloch wave packet to obtain finite action on the covering
space, but its envelope, total normalization and observations must then
be retained. This file does not quietly equate that construction with
a single cell. Nor does it identify this contractible transverse tube's
marked shape with a distant knotted tube or an absolute director.
The fixed periodic field is globally O(d) from0151's one-wave field;
its independently constructed mean-response sector remains an input at
its own time and wave-number scope. A joint continuum or a uniform
infinite-time acoustic/optical theorem does not follow from(1)--(7).
