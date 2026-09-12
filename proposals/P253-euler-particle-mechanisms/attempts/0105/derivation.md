# Derivation: normalized power-law response and the true leaf tangent

## 1. Frozen response convention

On the straight circular column, let `s` be transverse radius, let `P=P(s)`
be positive for `0<=s<a`, strictly decreasing for `0<s<a`, and vanish at the
free edge. Use the reviewed physical orientation

    -Delta_2 P=R^2 zeta(P),       Omega(s)=P_s/(R s)<0.       (1)

For the full-core power-law comparator set

    zeta(P)=C P^p,       chi(P)=zeta(P)/lambda_0,             (2)

with `C,p,P(0),epsilon_EM,c_EM,R>0`, `p>1`, and nonzero `lambda_0`.
This is a redesigned-tag comparator. It is not the existing strict-band tag
and it is not assumed to satisfy the 0077 stabilizer construction.

The first-speed Maxwell equations and global regular/decaying dipole Green
normalization reviewed in 0102 give

    Phi_0'(s)=-(epsilon_EM s)^(-1) integral_0^s t chi(t)dt,

    L_1 v=f,
    f=chi Omega s/(epsilon_EM c_EM^2),
    L_1=partial_s^2+s^(-1)partial_s-s^(-2),                 (3)

and

    v=-1/2[s^(-1) integral_0^s t^2 f(t)dt
               +s integral_s^a f(t)dt],
    S_1=chi_P v+chi Phi_0'/(c_EM^2 R).                       (4)

All formulas below use this one matching normalization.

## 2. Exact Lane--Emden full-core response

Define

    M_j(s)=integral_0^s t P(t)^j dt.                         (5)

Equation (1)--(3) makes the dipole source a total derivative:

    f(s)=C P(s)^p P_s(s)
          /(lambda_0 epsilon_EM c_EM^2 R).                  (6)

Since `P(a)=0`,

    integral_s^a f(t)dt
      =-C P(s)^(p+1)/[(p+1)lambda_0 epsilon_EM c_EM^2 R],    (7)

while integration by parts gives

    integral_0^s t^2 f(t)dt
      =C[s^2 P(s)^(p+1)-2M_(p+1)(s)]
        /[(p+1)lambda_0 epsilon_EM c_EM^2 R].               (8)

The two local `s P^(p+1)` rows cancel in (4), leaving the exact positive
Volterra field (with sign carried by `lambda_0`)

    v(s)=C M_(p+1)(s)
          /[(p+1)lambda_0 epsilon_EM c_EM^2 R s].           (9)

The electrostatic row is

    Phi_0'(s)=-C M_p(s)/(lambda_0 epsilon_EM s).             (10)

Substitution into (4) yields the full response.  Factoring the positive
profile power isolates the sign function

    D(s)=p M_(p+1)(s)/(p+1)-P(s)M_p(s),                    (11)

    S_1(s)=C^2/[lambda_0^2 epsilon_EM c_EM^2 R s]
      P(s)^(p-1) D(s).                                     (12)

This is a physical normalized response, not a copied sign.

At a smooth center,

    M_p=P(0)^p s^2/2+O(s^4),
    M_(p+1)=P(0)^(p+1)s^2/2+O(s^4),                        (13)

and therefore

    D(s)=-P(0)^(p+1)s^2/[2(p+1)]+O(s^4).                  (14)

Consequently,

    S_1(s)=-C^2 P(0)^(2p)s
      /[2(p+1)lambda_0^2 epsilon_EM c_EM^2 R]+O(s^3).      (15)

It is strictly negative on a sufficiently small punctured center annulus.
Near the edge, `M_p` and `M_(p+1)` approach positive constants and

    D(s) -> p M_(p+1)(a)/(p+1)>0.                          (16)

Thus `S_1>0` sufficiently near the edge and `S_1->0` there when `p>1`.
Continuity gives at least one zero.  For the actual Cao integer-power
Lane--Emden interior, the radial ODE and regular center data make `P`, the
moments, `D`, and `S_1` analytic on the connected open core.  Since (14) is
strictly negative for small positive `s`, the identity theorem gives the
stronger conclusion

    S_1 cannot vanish on any open full-core tagged band.     (17)

Hence the full-core power-law fixed-profile tag does not inherit the
Liouville cancellation. This route is refuted as an exact constant-lock
profile at first speed. The sign change also proves that no universal one-sign
response theorem should replace the exact formula.

For the existing strict-band tag, derivatives of the cutoff enter `f` and its
total moment in (7). Equations (7)--(12) do not evaluate those transition
rows. The actual strict-band normalized response remains open.

### 2.1 Classification of every regular full-core cancellation profile

The same calculation classifies the whole regular full-core constitutive
class. Let `chi=h(P)=zeta(P)/lambda_0`, assume `h(P)>0` throughout the core,
`P_s!=0` for `s>0`, `h` is smooth, and the matched dipole satisfies `v=O(s)`
at the center. The radial equations are

    -(P_ss+P_s/s)=R^2 lambda_0 h(P),
    -(Phi_0''+Phi_0'/s)=h(P)/epsilon_EM.                    (C1)

Regularity at the center removes the integration constant, so

    Phi_0'=P_s/(R^2 lambda_0 epsilon_EM).                   (C2)

If `S_1=0`, then `h_P` cannot vanish for `s>0`: the second term of (4) is
nonzero there. Define `G=h/h_P`. Equation (4) gives

    v=-P_s G/(R^3 lambda_0 epsilon_EM c_EM^2).              (C3)

Differentiating the first equation in (C1) gives

    L_1 P_s=-R^2 lambda_0 h_P P_s.                         (C4)

Substitute (C3) into `L_1v=hP_s/(R epsilon_EM c_EM^2)` and use
`G h_P=h`. The Lane--Emden terms cancel exactly and leave

    G_PP P_s^2+G_P(3P_ss+P_s/s)=0.                         (C5)

Equivalently,

    partial_s[s P_s^3 G_P]=0,
    G_P(P(s))=C_*/[s P_s(s)^3].                            (C6)

There is no unexamined degenerate-center case under these full-core
hypotheses: (C1) and `lambda_0 h(P(0))!=0` give
`P_s=-k s+O(s^3)` with nonzero `k`. If `C_*!=0`, (C6) gives
`G_P=O(s^-4)`, hence `G=O(s^-2)` and (C3) gives `v=O(s^-1)`, contradicting
the regular dipole condition. Therefore `C_*=0`, `G` is constant, and

    h_P=a h,        h(P)=A exp(aP).                         (C7)

The constant-source case `h_P=0` cannot evade the argument because it makes
the nonzero electric term in (4) uncancellable. Thus (C7) is the only
nontrivial regular full-core cancellation family. A compact Cao core has
finite `P(a)=0` and requires its constitutive vorticity/tag source to vanish
at the free boundary, whereas `A exp(aP(a))=A!=0`. Consequently

    no smooth positive compact full-core fixed-ratio tag cancels S_1.       (C8)

The noncompact Liouville exponential of 0102 is therefore the unique regular
constitutive kernel at this order, rather than a representative of a compact
family. This classification does not apply through a strict tag cutoff:
there `chi` is not the same full-core constitutive source as
`zeta/lambda_0`, and its transition moments remain Route A2.

## 3. Fixed-material-leaf mean-flux obstruction

### 3.1 Straight circular column

Let `xi=(xi_s,xi_alpha)` be a smooth compactly supported area-preserving
cross-sectional displacement. In polar coordinates

    partial_s(s xi_s)+partial_alpha xi_alpha=0.             (18)

Angular averaging gives

    partial_s[s <xi_s>_alpha]=0.                            (19)

Regularity at `s=0` or compact support fixes the constant to zero, hence

    <xi_s>_alpha=0.                                         (20)

For a radial base `zeta=F(P(s))`, its coadjoint tangent is

    delta zeta=-xi dot grad zeta=-zeta_s xi_s,
    <delta zeta>_alpha=0.                                   (21)

If the tag is a fixed function of the vorticity on the declared leaf, then
`delta chi=chi_zeta delta zeta` and its angular mean vanishes as well.
Therefore a nonzero axisymmetric radial source perturbation is not the tangent
of this compact fixed material leaf. The two radial inner/outer controls that
produced the unconstrained `-1/4` matrix in 0102 do not lift through this
first-order coadjoint map.

### 3.2 Finite torus

In the physical meridional half-plane the axisymmetric divergence row is

    partial_r(r xi_r)+partial_z(r xi_z)=0.                  (22)

For a regular closed `P` contour `C_I`, the divergence theorem gives

    integral_(C_I) r xi dot n ds=0.                         (23)

Since `delta zeta=-F'(P)xi dot grad P`,

    integral_(C_I) r delta zeta/|grad P| ds
      =-F'(P) integral_(C_I) r xi dot n ds=0.               (24)

This is the exact `r dr dz` version of (21). Boundary displacement,
circulation, impulse, center and tag-distribution rows remain part of the
finite-Cao steady map; none can turn a forbidden radial mean into the two
formal source controls without changing the leaf.

Thus the proposed radial `A/B` lift is refuted by a named material-flux
mechanism. This does not decide the full zero-mean band operator.

## 4. A positive raw zero-mean material tangent

The same divergence equation shows what survives. Let `T(s,alpha)` be smooth,
compactly supported inside a regular band, have zero angular mean at every
`s`, and suppose `chi_s` is bounded away from zero there. Set

    xi_s=-T/chi_s.                                          (25)

Then the right side of

    partial_alpha xi_alpha=-partial_s(s xi_s)               (26)

has zero angular mean, so it has a smooth periodic primitive. Choosing its
mean smoothly makes `xi` compactly supported, and

    delta chi=-xi_s chi_s=T.                                (27)

When `chi=zeta/lambda_0`, the same displacement gives
`delta zeta=lambda_0 T`.

The same construction is exact on a finite-Cao regular foliation.  Take
coordinates `(I=P,alpha)` with physical meridional measure

    r dr dz=J(I,alpha)dI dalpha.                            (28)

Then

    partial_I(J xi^I)+partial_alpha(J xi^alpha)=0,
    xi^I=-T/chi'(I).                                        (29)

Integration around a contour shows that a smooth periodic `xi^alpha` exists
precisely when

    integral J(I,alpha)T(I,alpha)dalpha=0                  (30)

for every `I` in the support.  Indeed, `chi'` depends only on `I`, and compact
support in `I` fixes the otherwise constant weighted normal flux to zero.
Conversely, under (30), integrate `-partial_I(J xi^I)` in `alpha` and choose
the periodic mean smoothly; this constructs compactly supported `xi^alpha`.
Thus the raw material-label tangent maps onto smooth weighted-zero-mean
contour functions on every regular finite-Cao band. This is a genuine positive
alternative to the forbidden radial mean.

The construction has one normal derivative of loss.  On a fixed compact
regular band, if `T`, `J`, `J^-1`, and `(chi')^-1` have the required bounded
`C^(k+1,alpha)` norms, then `xi^I=-T/chi'` has that regularity, while the
right side `-partial_I(J xi^I)` in the periodic primitive has only
`C^(k,alpha)` regularity.  Division by `J` therefore gives the tame estimate

    ||xi||_(C^(k,alpha))
      <=C(J,J^-1,(chi')^-1)||T||_(C^(k+1,alpha)).           (30a)

The corresponding fixed-band Sobolev map is `H^(k+1)->H^k`.  Smooth input
still gives a smooth compact displacement, but no same-order bounded inverse
is asserted.

Equations (25)--(30) do not solve the linearized steady Euler--Maxwell map.
The induced streamfunction, pressure, Maxwell field, free boundary and finite
rows can change or obstruct the response. Consequently they establish raw
material-tangent surjectivity, not surjectivity of `L_resp` or `L_full`.

## 5. Exact versus approximate solvability

Let the frozen spaces and operators be

    L_resp:X_leaf->Y_band,
    L_AB:X_leaf->R^2,
    L_full=(L_resp,L_AB),
    T_0=(S_Cao,m_0).                                        (31)

Exact cancellation is `T_0 in Ran L_full`. For every
`ell in ker L_full^*`, `<ell,T_0>=0` is necessary. It is sufficient only after
a closed-range/Fredholm estimate or bounded right inverse is proved. If the
range is nonclosed, these pairings characterize only membership in
`closure(Ran L_full)`; controls with diverging norm do not construct an exact
steady carrier.

No such full steady/free-boundary estimate is obtained here. The exact result
is instead a sharper decomposition:

- the full-core power-law fixed profile has a nonzero, sign-changing response;
- the only regular full-core cancellation law is exponential and cannot
  vanish at a finite compact-core edge;
- the radial source controls do not lie in the compact fixed leaf;
- the raw zero-mean material tangent is locally onto; and
- the passage from that raw tangent to the steady operator remains the next
  Fredholm/free-boundary construction.

## 6. Route verdicts and continuation

- **A0 — regular full-core cancellation classification: established.**
  Equations (C1)--(C7) prove that the only nontrivial regular constitutive
  cancellation law is exponential. Equation (12) independently gives the
  exact power-law response, with negative center coefficient (15), positive
  edge sign (16), and open-band nonidentity (17).
- **A1 — compact full-core fixed-ratio cancellation: refuted.** The only
  regular kernel is exponential and cannot vanish at the finite free edge;
  the power-law formula is an exposing special case. This remains a
  redesigned-tag route until a charged carrier with such a tag is constructed.
- **A2 — existing strict-band normalized response: blocked.** Its cutoff
  transition moments must be evaluated in the reviewed physical pairing.
- **B0 — radial annular source lift on the fixed material leaf: refuted.**
  Equations (18)--(24) force its radial mean to vanish.
- **B1 — raw zero-mean material tangent: established.** Equations (25)--(30)
  give an explicit divergence-free local realization on straight circular
  and finite-Cao regular bands, with the one-normal-derivative tame estimate
  (30a).
- **B2 — full steady `L_full` solvability: blocked.** The common-domain
  steady/free-boundary derivative, closed range and target pairing remain to
  be constructed.
- **C — redesigned tag/defect carrier: blocked.** It must re-enter the charged
  existence and stabilizer map; the formal full-core comparator is not a
  carrier.
- **D — curvature and second-order nonaxisymmetric continuation: blocked at
  its activated construction.** Use the exact zero-mean tangent as the input
  to the finite-Cao steady map; if the full residual class vanishes, compute
  the first curvature response, then the next harmonic at a zero.

The attempt advances P6 source selection but establishes no P2 persistence,
P4 action, Lorentz chirality, flavor dynamics, particle or parent claim.
