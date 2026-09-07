# Two-label forced locking on the charged Cao carrier

## 1. Distinct labels remove the old monodromy inference

Let `chi_c` be the real electromagnetic material tag and let `Theta` be an
independently transported circle phase on a punctured material domain. Put

    q_A=omega dot grad Theta,        q_A=lambda chi_c.     (1)

The forced vorticity and label equations give

    D_t q_A=(curl f) dot grad Theta,
    D_t chi_c=0,

and therefore, without dividing on the zero set,

    chi_c D_t lambda=(curl f) dot grad Theta.             (2)

The P253/0101 closed-line obstruction concerned
`omega dot grad chi=lambda chi`, where contraction along a vorticity line
directly differentiates the same label on the right. Here the two labels are
different, so (1) does not imply `d chi_c/dtau=lambda chi_c`. That monodromy
no-go cannot be imported.

The exact axial continuity current is

    partial_t q_A+div(q_A u-f cross grad Theta)=0.        (3)

It is a Galilean axial current for a true scalar phase. It is not yet a
Lorentz chiral current.

## 2. The steady Cao variable lock is automatic

On the reviewed charged Cao equilibrium,

    W=r^(-1)grad P cross e_theta,
    f=-(g/rho_m)chi_c(grad Phi+H grad P),
    chi_c=chi_c(P).                                       (4)

Direct cylindrical curl calculus gives

    (curl f) dot grad theta
      =(g/rho_m) W dot grad(chi_c' Phi-chi_c H).          (5)

Indeed, the terms with two gradients of `P` vanish and

    curl f=(g/rho_m)grad(chi_c'Phi-chi_c H) cross grad P.

Since
`W dot grad S=r^(-1)e_theta dot(grad S cross grad P)`, (5) follows with the
displayed sign. The modified Grad--Shafranov first integral is

    zeta=h(P)+(g/rho_m)(chi_c'Phi-chi_c H).               (6)

At the equilibrium `Theta=theta`, hence `q_A=zeta`. Wherever `chi_c` is
nonzero,

    lambda=zeta/chi_c                                     (7)

satisfies (2), because `W dot grad chi_c=0` and the `h(P)` term drops out.
Thus the **variable** forced lock is an established equilibrium identity once
the phase domain is supplied. It is not an additional compatibility gate and
does not select a constant ratio.

## 3. Support prevents a global constant lock

The reviewed tag is compactly supported strictly inside a regular band where
the stabilizer proof has `zeta,zeta'` bounded away from zero. The axial density
`q_A=zeta` occupies the full vorticity core. At any core point outside tag
support,

    q_A-lambda_0 chi_c=zeta !=0.                          (8)

Consequently the existing carrier cannot have a global pointwise constant
lock. Dividing by `chi_c` hides this obstruction. The two remaining honest
targets are a band-local lock with its boundary current, or a redesigned tag
and stabilizer covering the full required domain.

Because `q_A` is a pseudoscalar while `chi_c` is a true scalar, nonzero
`lambda_0` is a pseudoscalar carrier/defect datum. It cannot be called an
`O(3)`-scalar coupling.

## 4. First constant-ratio obstruction on a tagged band

Require on a band where `chi_c` is bounded away from zero

    h(P)+(g/rho_m)(chi_c'Phi-chi_c H)=lambda_0 chi_c(P).  (9)

Use the reviewed parity parameter `tau=g^2`:

    chi_c=chi_0+tau chi_2+...,
    h=h_0+tau h_2+...,
    lambda=lambda_0+tau lambda_2+...,
    Phi=g Phi_hat+O(g^3),       H=g H_hat+O(g^3).         (10)

At order one, (9) forces

    chi_0=zeta_0/lambda_0.                               (11)

At order `tau`, all `h_2`, `chi_2` and `lambda_2` terms depend only on the
streamline label. Let `Pi_P` be the period/action weighted mean on one regular
contour and `Q_P=I-Pi_P`. Applying `Q_P` leaves the necessary condition

    Q_P[chi_0' Phi_hat-chi_0 H_hat]=0.                   (12)

This is the first physical constant-lock oracle. An order-`tau` tag correction
cannot change it. If (12) is nonzero on one contour, the fixed-profile route is
refuted at that contour.

The next profile route must vary the uncharged carrier and tag together at
leading order while retaining (11), or change the phase geometry so that
`q_A` itself changes. Varying `chi_2(P)` alone is not a control of (12).

## 5. The scalar Maxwell block already carries angular structure

The reviewed Lorenz elimination gives, after writing
`Phi=g Phi_hat+O(g^3)`,

    L_c Phi_hat=(chi_0/epsilon_EM)
       [1-c^2/c_EM^2-c W_z/c_EM^2],                      (13)

with the frozen field/tag normalization. In the straight core coordinates the
principal cross-section operator is

    L_c=-partial_1^2-a partial_2^2,
    a=1-c^2/c_EM^2>0.                                    (14)

Before taking a straight-column limit, divide the exact reviewed Ampere
primitive by `g`:

    epsilon_EM a c_EM^2 r H_hat
       +epsilon_EM c partial_r Phi_hat=K(P)/r.            (13a)

Thus the exact magnetic-eliminated response is

    chi_0' Phi_hat-chi_0 H_hat
      =chi_0' Phi_hat
       -chi_0 K(P)/(epsilon_EM a c_EM^2 r^2)
       +chi_0 c partial_r Phi_hat/(a c_EM^2 r).           (13b)

This identity fixes the relative electric/magnetic normalization.  In a
curved ring the explicit `r^(-2)` row is itself not constant on a `P`
contour, so it must remain inside `Q_P`; it becomes a streamline-only row only
after the leading straight-column replacement `r=R`.

For a radial test function `phi(s)`, its values on the two coordinate axes
differ by

    (L_c phi)|_(alpha=0)-(L_c phi)|_(alpha=pi/2)
      =(1-a)[phi'(s)/s-phi''(s)].                        (15)

When `c!=0`, a radial solution with radial right-hand side is possible only if
`phi''=phi'/s`, hence `phi=C s^2/2+D` and the radial source is constant.
The Cao charge profile is not generally constant. Moreover the `c W_z` source
contains an odd first harmonic. Thus the scalar potential is not generically
a streamline function even at the straight-column rung.

Equation (15) is a useful positive diagnostic, but it does not by itself prove
that the full combination in (12) is nonzero: the toroidal magnetic response
`H_hat` can cancel it. The exact next calculation is therefore the coupled
Lorenz/Ampere Green response followed by the `Q_P` projection, with the affine
Coulomb row and physical band measure retained.

### 5.1 The exact straight-column first-harmonic oracle

Expand the straight-column Maxwell system for small nonzero physical
translation speed `c`, retaining the fixed positive-vorticity column and a
compact radial tag.  This computes the first coupled response exactly, but its
sign must be read in the physical Cao orientation rather than imposed by a
choice of symbol.  Write

    Phi_hat=Phi_0+c v(s)cos(alpha)+O(c^2).                (16)

At `c=0`,

    Phi_0'(s)=-(epsilon_EM s)^(-1)
       integral_0^s t chi_0(t)dt.                        (17)

If `W_z=Omega(s)s cos(alpha)`, differentiating (13) at `c=0` gives

    v''+s^(-1)v'-s^(-2)v
       =chi_0(s)Omega(s)s/(epsilon_EM c_EM^2).            (18)

The regular-at-zero, decaying-dipole solution is exactly

    v(s)=-1/2[s^(-1) integral_0^s t^2 f(t)dt
                    +s integral_s^infinity f(t)dt],
    f(t)=chi_0(t)Omega(t)t/(epsilon_EM c_EM^2).           (19)

The reviewed toroidal Ampere primitive, after writing
`H=g H_hat+O(g^3)` and dividing the exact equation by `g`, has leading
straight-column form

    (1/mu_EM-epsilon_EM c^2)R H_hat
      +epsilon_EM c partial_1 Phi_hat=K(P)/R.             (20)

Since `1/mu_EM=epsilon_EM c_EM^2`, its first `c` derivative at fixed
`P,chi_0,K(P)` yields

    H_1=-(c_EM^2 R)^(-1)partial_1 Phi_0.                  (21)

The first angular coefficient of the lock correction is therefore

    S_1(s)=chi_0'(P(s))v(s)
       +chi_0(s)Phi_0'(s)/(c_EM^2 R).                    (22)

The exact fixed-profile obstruction at this order is therefore

    Q_P[chi_0' Phi_hat-chi_0 H_hat]
       =c S_1(s) cos(alpha)+O(c^2).                        (23)

The two summands in (22) cannot be assigned the same sign from positive
vorticity alone.  With the convention `W_z=P_s cos(alpha)/R`, a positive Cao
vorticity has `P_s<0` near its core maximum.  Thus the sign convention
`Omega>0` used in an earlier scratch argument would reverse the source in
(18): the electric-dipole and Ampere terms can oppose one another.  The exact
criterion is `S_1(s_0)!=0` on one tagged contour, not a universal sign claim.

Equations (16)--(23) establish a one-dimensional Volterra evaluation of that
criterion.  They do not yet evaluate it for the normalized Lane--Emden/Cao
profile and the reviewed annular tag.  That evaluation must retain the Cao
core/major-radius scaling, the tag cutoff, and the physical `P` orientation.
If it is nonzero, the straight-column fixed-profile route is refuted; if it
vanishes, the first curvature response becomes the next test.  Promotion to a
finite Cao member additionally requires uniform coupled Maxwell Green
convergence on a joint subluminal carrier path and control of the band
boundary.

### 5.2 An exact exponential-profile cancellation

The sign failure in Section 5.1 generates a positive representation change.
In physical straight-column radius `s`, the uncharged streamfunction and the
order-zero electric potential obey

    -Delta P=R^2 zeta_0(P),
    -Delta Phi_0=zeta_0(P)/(lambda_0 epsilon_EM).          (24)

With the same radial Green normalization, their radial derivatives satisfy

    Phi_0'=P'/(lambda_0 epsilon_EM R^2).                  (25)

Both two-dimensional potentials carry the corresponding logarithmic affine
tail when their total source is nonzero.  Equation (25) is a statement about
their globally matched derivatives; it does not place `Phi_0` in a decaying
zero-monopole potential space.

Differentiating the first equation in a Cartesian transverse direction gives

    L_1 P'=-R^2 zeta_0'(P)P',
    L_1=d_s^2+s^(-1)d_s-s^(-2).                          (26)

Choose the exact profile law

    zeta_0'(P)=a_0 zeta_0(P),       a_0!=0,               (27)

This law has the explicit regular radial Liouville realization

    P_b(s)=P_0-(2/a_0)log(1+b s^2),
    zeta_b(s)=8b/[a_0 R^2(1+b s^2)^2],       a_0,b>0.     (27a)

Direct radial differentiation gives
`-Delta P_b=R^2 zeta_b` and `d zeta_b/dP_b=a_0 zeta_b`.
Its cross-sectional vorticity is integrable, but `P_b'=O(s^(-1))`; the
kinetic energy per unit length of the infinite column diverges
logarithmically.  Thus (27a) realizes the response kernel exactly while
exposing its localization cost.

and retain `chi_0=zeta_0/lambda_0`. Equations (18) and (26) then give the
regular/matched dipole solution

    v=-P'/(lambda_0 epsilon_EM c_EM^2 R^3 a_0).           (28)

There is no added homogeneous `s` row when (24)--(28) use the same global
Green/matching normalization. Since

    chi_0'=a_0 zeta_0/lambda_0,
    H_1=-P'/(lambda_0 epsilon_EM c_EM^2 R^3),              (29)

the first harmonic cancels exactly:

    S_1=chi_0'v-chi_0H_1=0.                              (30)

Thus the first-speed response operator has a genuine profile kernel in the
global straight-column representation. This prevents the failed sign
argument from becoming a class no-go.

The candidate does not yet meet the Cao carrier domain. Equation (27) gives
`zeta_0=C exp(a_0P)`, which stays nonzero at every finite threshold and cannot
be the compact positive-part free-boundary law. Truncating it changes the
global dipole matching and generally restores a homogeneous `s` row in (28).
The next compact construction is therefore to solve the profile/free-boundary
Fredholm map with that matching coefficient and the circulation, mean-radius,
tag-leaf and stabilizer rows included. A nonzero cokernel there activates a
redesigned full-core tag or phase geometry rather than a global no-go.

### 5.3 Exact annular matching rows

The compact-support mismatch can be typed exactly before attempting the Cao
map. On an annulus where (27) holds, the difference between the actual dipole
and the particular solution (28) solves

    L_1 w=0,             w=A s+B/s.                       (31)

For a source perturbation `delta f`, the Green formula (19) shows that a
perturbation supported strictly inside the annulus produces

    delta A=0,
    delta B=-(1/2) integral t^2 delta f(t)dt,              (32)

whereas one supported strictly outside produces

    delta A=-(1/2) integral delta f(t)dt,
    delta B=0.                                            (33)

The source-to-matching matrix is therefore

    [[0,-1/2],[-1/2,0]],      determinant=-1/4.           (34)

This is an exact two-row surjectivity theorem for independent radial Maxwell
source variations. It makes compact matching constructive rather than merely
naming a Fredholm problem. It is not yet a Cao-profile theorem: admissible
variations must satisfy
`delta chi_0=delta(zeta_0/lambda_0)` and simultaneously move `P`, the free
boundary, circulation, mean radius and Casimir/material distribution. The next achievement is to lift
(34) through that joint linearized map and exhibit enough profile controls to
satisfy its finite rows without changing the declared leaf silently.

## 6. The phase defect is not carrier-local

For `Theta=theta`, the removed set is the entire symmetry axis. It is a
noncompact defect, not a compact puncture carried by the vortex ring. The
Euclidean norm `|dtheta|=1/r` also makes a naive unweighted phase-gradient
energy divergent unless both its axis core and outer behavior are treated.

If Route D replaces the symmetry axis by a compact ring defect, it must
construct a new circle map. The identities `Theta=theta` and `q_A=zeta` then
do not transfer automatically. This topology/energy row is part of the moving
material-phase construction and prevents the equilibrium integral from being
misread as localized particle structure.

## 7. Route verdicts at this rung

- **Route A — established at equilibrium identity scope.** Equations
  (2)--(7) prove the exact variable two-label forced lock wherever the tag is
  nonzero. The moving material phase remains a separate construction.
- **Route B0 — established reduction.** Equations (9)--(12) reduce
  fixed-profile constant locking to one exact zero-mean response.
- **Route B1 — blocked at the normalized straight-column response.**
  Equations (16)--(23) reduce the first harmonic to the explicit Volterra
  coefficient `S_1`.  Its nonvanishing for the actual Lane--Emden/Cao profile
  and tagged band has not yet been proved.
- **Route B2 — blocked at finite-curvature Cao promotion.** If B1 is nonzero,
  uniform coupled
  Maxwell Green convergence, the band boundary, and an actual joint carrier
  path are still required before B1 becomes a finite-Cao refutation.
- **Route C0 — established at global straight-column first-speed scope.**
  The explicit Liouville column (27a) realizes the exponential law and the
  exact cancellation (30). Its noncompact vorticity/free-boundary and
  logarithmically divergent infinite-column energy costs are explicit.
- **Route C1 — established at annular Maxwell-source matching scope.**
  Equations (31)--(34) give independent inner/outer control of the two
  homogeneous dipole coefficients.
- **Route C2 — blocked at a compact joint uncharged-profile/tag Fredholm
  control.** An
  order-`tau` tag perturbation is proved insufficient; the valid leading
  control must preserve (11) while moving the carrier/profile and every leaf
  row.
- **Route D — blocked at the global advected circle phase and moving defect.**
  The equilibrium phase is not yet a nonaxisymmetric material-domain theorem.
- **Route E — blocked at the shared P4 action and Lorentz chiral flavor map.**
  No flavor or P6 inference is made from the axial lock.

The next analytic achievement is to solve the coupled order-`g` Maxwell block
for the fixed Cao column/band and evaluate (12). A nonzero result activates the
joint profile/carrier control; a zero result activates the first curvature or
finite-core response rather than closing the constant-lock route.
