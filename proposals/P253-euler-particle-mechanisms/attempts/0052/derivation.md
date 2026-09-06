# Fixed-domain and full-sector Cao graph-Riesz continuation

## 1. Correct inversion of the resonant index

Use one common core-time normalization and write the two-term high-index
Kelvin law and bending law as

    sigma_n(delta)=A delta/(n+beta_0)
       +B delta/(n+beta_0)^2
       +O(delta/n^3+delta^3/n),                            (1)

    tau(delta)=c delta^2 L+d delta^2+o(delta^2),
    L=log(1/delta),       A,c>0.                           (2)

Put `y=n+beta_0`. Solving (1)--(2), first without discarding either finite
term, gives

    y=A/[delta(cL+d)]+B/A
       +o(1/(delta L^2))+O(1),                             (3)

and hence

    n_*(delta)=A/(c delta L)
      -A d/(c^2 delta L^2)
      +(B/A-beta_0)
      +o(1/(delta L^2))+O(1).                              (4)

The `d` term is parametrically larger than the phase correction. Thus the
correct concise remainder is

    n_*=A/(c delta L)+O(1/(delta L^2)+1),                  (5)

not `O(L^-1)`. Since `n_*^{-1}=Theta(delta L)`, (4) leaves

    g_(n_*)=O(delta^3 L),
    Delta sigma_(n_*)=Theta(delta^3 L^2),
    g_(n_*)/Delta sigma_(n_*)=O(L^-1)                      (6)

unchanged. The finite window and asymptotic bad-set fraction from 0048 survive
with shifted window centers.

In physical units, if

    F(a)=integral_0^a W(r)r dr,       Gamma=2*pi*F(a),
    Omega_a=F(a)/a^2,                                      (7)

then the universal leading bending frequency is

    nu_b=Gamma L l sqrt(l^2-1)/(4*pi*R^2),                 (8)

so in the time unit `Omega_a^-1`,

    tau=nu_b/Omega_a
       ={l sqrt(l^2-1)/2} delta^2 L.                       (9)

Thus `c=l sqrt(l^2-1)/2` in this normalization. Using instead
`Gamma/(4*pi*a^2)` as the clock doubles `c`; it cannot be mixed with (7).

## 2. Why a literal fixed-to-shrinking volume map is wrong

A unit solid torus and the physical Cao core have unequal volumes. Therefore
no volume-preserving diffeomorphism between them exists. The common-domain
construction must separate:

1. physical center/core scaling, whose constant Jacobian is retained; and
2. a normalized shape correction, which is exactly volume preserving.

This distinction is also required for the density and KKS normalization.
The resulting cross-`delta` map is an orientation-preserving operator
conjugacy, not a fluid motion on one coadjoint orbit. Only a displacement of
one already fixed physical carrier is required to have determinant one.

Let

    Psi_delta(r,z)=psi_delta(r,z)-c_delta r^2/2            (10)

be the translating-frame streamfunction, with sign chosen so

    U_delta^r=-r^-1 partial_z Psi_delta,
    U_delta^z= r^-1 partial_r Psi_delta.                   (11)

The meridional weighted area form is

    dA=r dr wedge dz,        i_(U_delta)dA=-dPsi_delta.    (12)

The Cao core is a small simply connected perturbation of the Lane--Emden
disk. For sufficiently small `delta`, its interior level sets are nested
closed streamlines with one elliptic center and a regular boundary streamline.
Define the action and angle by

    I_delta(h)=(2*pi)^-1 integral_{Psi_delta>=h} r dr dz,
    beta_delta=2*pi t/T_delta(I) mod 2*pi,                 (13)

where `t` is travel time along the streamline and `T_delta(I)` its period.
Direct differentiation of enclosed weighted area gives

    dA=dI wedge d beta,
    U_delta=omega_delta(I) partial_beta,
    omega_delta(I)=partial_I Psi_delta.                    (14)

The sign of `beta` is chosen to agree with (12). Normalize
`s^2=I/I_a(delta)` and keep the constant `I_a(delta)` explicit. Equations
(13)--(14) are an exact weighted action--angle chart, not a polar relabeling.
They diagonalize the principal core transport for every `delta`.

## 3. Explicit Hanzawa--Moser and Piola map

Let `K_*=[0,1]_s x S^1_beta x S^1_theta`, with the regular polar condition at
`s=0`. The raw Hanzawa map is

    H_delta(s,beta,theta)
      =(r_delta(s,beta) e_r(theta),z_delta(s,beta)),       (15)

where `(r_delta,z_delta)` is obtained by following (13) from the center and
the boundary value at `s=1` is the Cao free-boundary graph. Extend (15) to a
collar using a fixed cutoff and to the identity relative to the circular
physical reference torus outside that collar.

Its volume pullback is

    dmu_*=s ds d beta d theta,
    H_delta^*(dx)=Jbar_delta dmu_*,
    Jbar_delta=2 I_a(delta).                               (16)

The coordinate determinant relative to `ds d beta d theta` is
`2I_a s`, while the Jacobian relative to the fixed polar volume `dmu_*` is
the constant `Jbar_delta`. The factor is constant in angle because (14) is canonical. If a raw
free-boundary Hanzawa chart is used before (13), write its normalized density
as `j_delta/bar j_delta`. Solve the explicit Moser equation

    div(mu_t Y_t)=bar j_delta-j_delta,
    mu_t=(1-t)j_delta+t bar j_delta,                        (17)

with zero normal trace using the Bogovskii right inverse on the collar, and
take the time-one flow. The two densities in (17) have equal total volume by
the definition of `bar j_delta`. It produces a constant **normalized**
Jacobian after the physical dilation/volume factor has been extracted. It
does not make a fixed unit torus and `K_delta` equal-volume. Thus (13) is the
intrinsic version of the normalized Hanzawa--Moser correction and (17) is its
fixed-domain construction.

For the raw orientation-preserving map, before the normalized correction,
write the full Jacobian ratio as `J_Phi(y)>0`. The uniform geometry target is

    0<c<=J_Phi(y)/Jbar_delta<=C,
    ||D Phi_delta||+||D Phi_delta^-1||<=C_scale(delta),    (17a)

where the known physical dilation is removed from `C_scale` when fixed norms
are compared. The general contravariant Piola map is

    P_Phi eta=J_Phi^-1 D Phi eta o Phi^-1.                 (17b)

Let `Phi_delta` denote the resulting map and `Jbar_delta` the constant physical
scale in (16). On the closed vorticity core, where the action--angle map has
exactly this constant Jacobian, the contravariant Piola map is

    P_delta eta=Jbar_delta^-1
                  D Phi_delta eta o Phi_delta^-1.          (18)

It obeys the exact distributional identities

    div(P_delta eta)=Jbar_delta^-1
       (div_* eta)o Phi_delta^-1,                          (19)

    P_delta[curl_*(xi cross Omega_*)]
      =curl[(Phi_delta)_*xi cross P_delta Omega_*].        (20)

Hence it maps fixed DA tangents modulo their stabilizer onto the physical DA
tangent **spaces** and preserves closed-core support. This cross-carrier map
does not say that the two base vorticities lie on one coadjoint orbit. For a
same-carrier displacement `g` at fixed `delta`, impose `det Dg=1`, use
`g=id` outside the physical collar, and then (20) is the actual coadjoint
tangent covariance. The associated velocity is not a
Piola transform of a reference velocity; it is exactly

    v_delta[eta]=BiotSavart_R3(P_delta eta),               (21)

which retains the collar, exterior harmonic tail, and pressure.

The collar extension of `Phi_delta` need not have constant Jacobian. There
the operator conjugacy uses the general formula (17b), with its full
`J_Phi(y)` and cofactor factors; (18)--(20) are the constant-Jacobian core
restriction. Since every physical vorticity tangent has distributional curl
support in the closed core, the collar enters only through coordinate traces,
the pulled-back metric, and the exact whole-space completion (21). No
constant-Jacobian or volume-preserving assertion is made for the cross-delta
collar map.

## 4. The actual fixed energy space

Pull (21) back to `K_*` and define the exact metric operator

    <eta,G_delta zeta>_*=rho_0 integral_R3
       v_delta[eta] dot conjugate(v_delta[zeta]) dx.       (22)

The thin-tube Green split has one logarithmically large centerline block. Let
`P_c` be the finite physical circulation/displacement projection and put

    S_delta=L^-1/2 P_c+(1-P_c),
    mathcal U_delta=P_delta S_delta.                       (23)

The near-diagonal Biot--Savart integral gives

    S_delta^*G_delta S_delta=G_*+O(L^-1)+O(delta)          (24)

as a coercive form on the energy completion; the far kernel, collar, and
exterior pieces are uniformly bounded. Therefore

    C^-1||eta||_(E,*)<=||mathcal U_delta eta||_(E,delta)
                         <=C||eta||_(E,*)                  (25)

with `C` independent of small `delta`. Without `S_delta`, (25) is false on
the displacement block. This is the explicit energy-domain correction absent
from a raw Hanzawa pullback.

The action--angle transport domain is

    D_*={eta in X_*:
         omega_0(I) partial_beta eta+K_0 eta in X_*},      (26)

where `K_0` is the source-defined nonlocal vorticity/velocity block. Since
`omega_0(I)>=omega_min>0` on the compact core, its `m!=0` restriction is
equivalent to one `beta` derivative. The `m=0` restriction is the energy
completion of the compact Kelvin block.

The exact conjugated generator is

    Ahat_delta=mathcal U_delta^-1 A_delta mathcal U_delta
      =T_delta+K_delta.                                    (27)

Here `T_delta` is the exact block-diagonal local coadjoint
transport/stretching operator. In streamline Fourier sector `m` it has the
triangular form

    T_(delta,m)=-i m omega_delta(I) Identity+S_delta(I),
    S_delta(I)^2=0,                                        (27a)

in the transported vorticity frame. The `m=0` DA reduction is kept in its
canonical displacement variables, where it is the positive Kelvin wave
system rather than an unconstrained shear Jordan block. `K_delta` contains
the exact nonlocal global Biot--Savart/Leray response. Its principal
action--angle order is `-1`; apparent order-zero cross-`m` curvature mixing
in polar coordinates is absorbed into the exact streamline chart (14), while
the local shear remains explicitly in (27a).

Equations (13)--(27) construct a common energy space and the candidate common
graph domain. The **uniform graph equivalence** and a second-order expansion
of (27) require the full coefficient regularity recorded as `HJ2` in section
8; Cao's published elliptic remainder alone does not prove them.

## 5. Exact KKS, energy, and frequency normalization

Use the right-reduced convention

    Omega_KKS(X_xi,X_zeta)
      =-<m,[xi,zeta]>
      =rho_0 integral Omega_delta dot (xi cross zeta) dx,
    i_(X_H)Omega_KKS=dH.                                  (28)

Let `L_delta` be the physical Hessian and
`A_delta=mathbb J_delta L_delta`. If

    A_delta q=i nu q,                                     (29)

then complex-bilinear extension of (28) gives exactly

    L_delta(q,conjugate(q))
       =Omega_KKS(A_delta q,conjugate(q))
       =i nu Omega_KKS(q,conjugate(q)).                    (30)

For the canonical oscillator `Omega=dq wedge dp`,
`H=(p^2+nu^2q^2)/2`, the positive-frequency vector `(1,i nu)` satisfies

    Omega(q_+,conjugate(q_+))=-2i nu,
    L(q_+,conjugate(q_+))=2nu^2,                           (31)

which checks (30). Thus an energy-unit vector
`L(e_E,conjugate(e_E))=1` has

    Omega(e_E,conjugate(e_E))=-i/nu,                      (32)

and the KKS-unit vector with pairing `-i` is

    e_K=sqrt(nu)e_E.                                      (33)

The sign in (32) corrects the positive-`i` scratch convention in 0048.

For an axisymmetric column mode normalized by
`integral Phi|u_n|^2 r dr=1`, the exact velocity equations give

    ||v_n||_cross^2=sigma_n^-2 M_n,
    M_n=integral(Phi+W^2)|u_n|^2 r dr,   1<=M_n<=2.        (34)

After the physical ring integration, the Hessian norm contains
`2*pi*R_delta*a_delta^2*rho_0` times (34), plus the exact
`1+delta s cos beta` correction. Consequently

    e_(n,E)=sigma_n
      [2*pi*rho_0 R_delta a_delta^2 M_n]^-1/2 v_n
      (1+O(delta)),                                       (35)

and (33) gives its KKS scale. Translation modes are normalized from (28)
using their impulse duals, then by the same Hessian convention; no radial
`L2` dual is substituted.

The chain of frequency units is now fixed:

    lambda_n=k_delta^2/sigma_n^2,
    nu_(n,physical)=Omega_a sigma_n,
    eigenvalue(A_delta)=i nu_(n,physical).                 (36)

Equations (7)--(9) define `Omega_a`; all `2*pi` factors and density cancel
only after (28), (34), and (35) have been formed.

## 6. Full-sector complement reduction

Decompose the fixed action--angle space as

    X_*=direct_sum_(m in Z) X_m,
    partial_beta|_(X_m)=i m.                              (37)

Let `P_delta` contain the two real bending coordinates and every real/KKS
partner of an `m=0` Kelvin mode in the corrected four-coupling window; let
`Q_delta=1-P_delta`.

### 6.1 Nonzero streamline modes

For `|z|<omega_min/2`, the triangular local part satisfies

    ||[T_(delta,m)-z]^-1||
      <=C{(|m|omega_min)^-1+(|m|omega_min)^-2}.            (38)

The nonlocal block `K_delta` is order `-1` in the fixed core variables. Thus
`K_delta[-i m omega-z]^-1=O(|m|^-1)`, and a Neumann series gives the complete
large-`|m|` inverse and summation.

Only finitely many nonzero sectors remain. In `m=1`, the ground-state
factorization from 0044 gives a simple translation kernel and a positive
complement gap. For `|m|>=2`, the exact scalar form obeys

    t_m=t_1+(m^2-1) integral |psi|^2/r dr
                 +(m-1)|psi(a)|^2,                        (39)

so it has no zero mode. Small `k_delta` and fixed-domain coefficient
perturbations preserve these finite-sector inverses. Negative modes follow
by the physical real conjugation transported through (18), not by an
unidentified ambient domain. Center, axial-impulse, and stabilizer zeros are
removed by the same finite KKS duals used in `P_delta`.

### 6.2 Exact `m=0` effective block

Let `Q_perp` collect all `m!=0` complement sectors. Their inverse defines the
Grushin/Feshbach operator on `m=0`:

    A_eff(z)=A_00-z
      -A_0perp[Q_perp(Ahat_delta-z)Q_perp]^-1A_perp0.      (40)

In action--angle coordinates `A_00` is order `-1`, with principal eigenvalue
law `A delta/n`. Each off-diagonal curvature block is order `-1`; hence the
Schur term in (40) is order `-2`. Its eigenvalue correction is
`O(delta^2/n^2)`, not an order-zero band. This is precisely why the corrected
spacing and the leading coupling result (6) survive the all-sector reduction.

With the one conjugate pair in the coupling window retained, the remaining
`m=0` resolvent has separation

    r_delta=gamma delta^3 L^2                             (41)

for a fixed `gamma>0` on the good parameter set. The leading translation-tail
Schur sum remains `O(delta^3 L)`, but (40) now proves that this is the Schur
complement of **all** nonzero sectors, rather than a scalar-tail surrogate.

Combining (38)--(41) gives the formal full complement bound

    sup_(z in Gamma_delta)
      ||[Q_delta(Ahat_delta-z)Q_delta]^-1||_(X_* -> D_*)
        <=C/(delta^3 L^2).                                (42)

The collar and exterior velocity enter `K_delta` through (21); fixed-core
support removes passive exterior vorticity but not this global operator.

The limiting estimates (38)--(41) are unconditional. Upgrading (42) to the
exact Cao generator requires the modewise remainder and common graph control
in `HJ2`; an `o(delta^2)` unstructured operator norm is insufficient compared
with (41).

## 7. The exact hypothesis needed for the Cao graph jet

The exact rescaled elliptic equation determines the local cells, and the
action--angle chart determines the common principal transport. What the
published Cao estimates do not supply is the following full Euler statement.

The direct route does evaluate the coefficients hidden by the first printed
Green big-`O`. With

    rho_G=((r-r')^2+(z-z')^2)/(r r'),
    kappa^2=4/(4+rho_G),       kappa_prime^2=1-kappa^2,    (43a)

the exact angular integral in Cao (2.2) is

    G_1=sqrt(r r')/(2*pi)
      {(2/kappa-kappa)K(kappa)-(2/kappa)E(kappa)}.         (43b)

The defining series of the complete elliptic integrals at
`kappa_prime=0` give

    G_1=sqrt(r r')/(4*pi) {
       log(1/rho_G)+2log(8)-4
       +rho_G[3log(1/rho_G)/16+3log(8)/8-1/8]
       +O(rho_G^2 log(1/rho_G))}.                          (43c)

Thus the logarithmic and finite second Green cells are source-defined, not
free constants. Likewise, expanding the exact rescaled Cao equation

    -Delta w+[delta/(1+delta x)]partial_x w
       =(1+delta x)^2(w_+)^p                              (43d)

for `w=U+delta V+delta^2 Z+...` yields

    L_U V=2xU^p-U_x,                                      (43e)

    L_U Z=p(p-1)U^(p-2)V^2/2+2pxU^(p-1)V+x^2U^p
             -V_x+xU_x.                                  (43f)

Equations (43a)--(43f) construct the coefficient recurrence for the profile,
shape, singular Green, and finite Green pieces. What remains is not an absent
Taylor coefficient; it is the uniform mapping estimate obtained after those
pieces are inserted into the moving-boundary Euler generator.

**Hypothesis HJ2.** For one fixed Cao member with sufficiently large finite
`p`, the maps (13)--(23) obey uniform graph equivalence

    C^-1||eta||_(G,*)<=||mathcal U_delta eta||_(G,delta)
                         <=C||eta||_(G,*),                 (43)

and on the common domain (26)

    Ahat_delta=A_col(k_delta)+delta C_1
      +delta^2 L C_(2,log)+delta^2 C_2+R_delta.            (44)

After subtracting its diagonal order-`-1` and `-2` symbols into the exact
renormalized Kelvin frequencies, the remainder satisfies the modewise bounds

    |<e_n,R_delta e_j>_E|
       <=epsilon_delta delta^2/[(1+n)(1+j)],
    |<t,R_delta e_n>_E|
       <=epsilon_delta delta^2/(1+n),
    epsilon_delta ->0,                                    (45)

and the analogous weighted all-`m` bounds. The KKS form and dual projections
are twice differentiable in the same chart.

The exact local Cao equation yields the candidate first/second profile cells,
while the exact Green integral yields the logarithmic and finite kernels.
However, the source proves only an elliptic `O(epsilon^2|log epsilon|)`
remainder on its moving support. It neither proves (43) nor the global
Biot--Savart/Leray matrix estimates (45). No primary result found in the
registered inventory closes that implication. Therefore `HJ2` remains an
explicit unproved construction, not a fact attributed to Cao.

Route D2 does not remove it: form convergence from the printed elliptic bound
controls low fixed modes, but the resonant index tends to infinity as (4), so
norm-resolvent convergence without the modewise rate (45) does not control the
contour (41). The failure is a rate/domain gap, not physical nonexistence.

## 8. Conditional full graph-Riesz theorem

Assume `HJ2`. Equations (38)--(45) make the all-sector remainder smaller than
the local gap (41), after the retained window is enlarged by any pole whose
actual coupling reaches four times its separation. The exact complement
inverse (42) then exists. The physical Feshbach map is

    K_delta(z)=P_delta(Ahat_delta-z)P_delta
      -P_delta Ahat_delta Q_delta
       [Q_delta(Ahat_delta-z)Q_delta]^-1
       Q_delta Ahat_delta P_delta.                         (46)

The same-sign calculation uses (30)--(33); it cannot inherit the incorrect
0048 sign. When the constrained Hessian is positive on the retained Kelvin
and bending coordinates, exact resonance is an elliptic avoided crossing.
The reconstructed Riesz projection

    Pi_delta=(2*pi*i)^-1 integral_(Gamma_delta)
                    (z-Ahat_delta)^-1 dz                  (47)

is finite rank, bounded on `X_*` and `D_*`, KKS nondegenerate, and transported
by (23) to compact DA vorticity with the exact exterior velocity (21).

Thus the full graph-Riesz theorem is **established conditional on HJ2**. The
unconditional achievements are the corrected crossing law, explicit
action--angle/Hanzawa--Moser/Piola energy chart, exact normalization, and
all-sector limiting Grushin reduction. The missing construction is exactly
(43)--(45).

## 9. Nonlinear scope

The 0048 expressions labelled there `(52)--(53)` remain only

    V^(N)(a)=V_0+sum_(r=1)^N a^r V_r,
    ||F(V^(N),c^(N),Omega^(N))||<=C_(N,delta)|a|^(N+1),    (48)

for the finite regularity/order actually earned. Equation (48) is an
approximate residual expansion, not an exact branch. Nothing in the
conditional Riesz theorem changes that boundary. Exact nonlinear existence
still requires the separately identified critical-layer range/transparency
and convergence construction.

## 10. Route verdicts

- Correct resonance inversion and physical frequency clock: **established**
  by (1)--(9).
- Fixed-domain chart and energy equivalence: **established** by the exact
  action--angle/Hanzawa--Moser/Piola construction (13)--(25). Uniform graph
  equivalence is isolated inside `HJ2`, not silently included.
- KKS/energy/generator normalization: **established** by (28)--(36), including
  the canonical negative-`i` sign.
- All-sector complement: **established for the source-defined limiting and
  polyhomogeneous operator**, with its exact Grushin structure (37)--(42);
  transfer to the exact Cao generator is conditional on `HJ2`.
- Direct moving-boundary/global-Green graph jet: **blocked by the named
  modewise common-domain construction (43)--(45)**. Cao's elliptic remainder
  is not enough at `n_*(delta)`.
- Full physical graph-Riesz theorem: **established conditional on HJ2**, not
  unconditionally.
- The fixed-order nonlinear calculation remains **approximate only** as in
  (48); no exact branch is claimed.

This attempt does not refute the rotating-wave carrier, infer a global
impossibility, or license stability, quantization, an electron, or a neutrino.
