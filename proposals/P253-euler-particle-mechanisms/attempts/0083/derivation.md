# Fixed-mean Cao scale path, IVT coverage, and finite charged window

## 1. Exact conventions and the author-scoped input

Write `x=(r,z)` and

    dnu=r dr dz,
    dm_epsilon=zeta_epsilon dnu,
    kappa=integral dm_epsilon,
    R=kappa^(-1) integral r dm_epsilon.                    (1)

Thus the exact mean-radius row is

    B_R(zeta)=integral (r-R)zeta dnu=0.                    (2)

The vorticity law and threshold are

    zeta=epsilon^(-2)P_+^p,
    P=G_1 zeta-c r^2/2-mu,       p>=6.                    (3)

P253/0080 is used at author scope for one precise input: on the even compact
source space, the derivative of `(3), kappa, B_R` with respect to
`(zeta,mu,c)` is an index-zero isomorphism with a uniform rescaled inverse for
sufficiently thin cores.  The mean-radius border is distinct from 0080's
circulation--impulse border.  Nothing below replaces that input by Cao's
auxiliary equation (3.36), a uniqueness statement, or the full Euler
generator.

## 2. The coefficient is fixed by the PDE and physical mass

Let `U` be the positive radial Lane--Emden profile

    -Delta U=U_+^p,
    {U>0}=B_1,
    Lambda_p=integral_(B_1) U^p dy.                       (4)

Choose a provisional physical scale `s_0=C_s epsilon` and amplitude `H_0`,
and set

    r=R+s_0 y_1=R(1+delta_0 y_1),
    z=s_0 y_2,
    P(r,z)=H_0 w(y),
    delta_0=s_0/R.                                        (5)

For Cao's operator

    L=-r^(-1) div_(r,z)(r^(-1) grad_(r,z)),               (6)

the exact interior equation becomes

    -q^(-1) div_y(q^(-1) grad_y w)
       =R^2 C_s^2 H_0^(p-1) w_+^p,
    q=1+delta_0 y_1.                                      (7)

The leading coefficient is one precisely when

    R^2 C_s^2 H_0^(p-1)=1.                                (8)

The circulation row is exactly

    kappa=H_0^p R C_s^2 integral q w_+^p dy.              (9)

At `delta_0=0`, `w=U`; hence (8)--(9) require

    H_0^p R C_s^2 Lambda_p=kappa.                         (10)

Eliminating `H_0` gives the unique positive coefficient

    C_s(kappa,R)
      =(Lambda_p/kappa)^((p-1)/2) R^(-(p+1)/2),
    H_0=(R C_s)^(-2/(p-1)).                               (11)

This is also the coefficient obtained by solving Cao (3.35),

    Lambda_p x_epsilon,1^(-(p+1)/(p-1))
      (epsilon/s_Cao)^(2/(p-1))
      =kappa+O(epsilon^2 |log epsilon|),                  (12)

after the exact mean row is used to replace the core center by `R` at leading
order.  Equation (11) is derived from (7)--(10), not fitted from the
Kelvin--Hicks answer.

The normalized finite rows are equally explicit:

    Lambda_p^(-1) integral q w_+^p dy=1,
    integral y_1 q w_+^p dy=0.                            (13)

The second row is the exact pullback of (2).  Axial evenness supplies the
other translation row.

## 3. Route A derivative ledger and its exact proof boundary

The positive set of `U` has nonzero normal derivative at `|y|=1`.  A natural
Route-A chart uses its normal graph plus a collar Hanzawa map and the unknown

    Z=(w,hat_mu,hat_c,boundary_graph)                      (14)

together with (13).  The exact local Green expansion contains
`delta_0`, `delta_0 log(delta_0)`, and higher polylogarithmic rows; formally
their logarithmic derivative `epsilon partial_epsilon` has the same vanishing
size.  The limiting bordered derivative is the 0080 Lane--Emden border.

This is not yet a complete weighted-`C^1` theorem.  To infer

    ||Z_epsilon-Z_0||
      +||epsilon partial_epsilon Z_epsilon||
        <=C epsilon |log epsilon|                        (15)

one must differentiate the *entire* pulled-back whole-space Green map,
positive-part interface, Hanzawa collar, center definition, and far/axis rows
and prove a uniform nonlinear inverse radius.  The value-level source
remainders do not license that differentiation.  Route A is therefore kept
open at precisely this map-derivative ledger; equation (15) is a target, not
an input to Route B.

## 4. Route B from Cao's exact auxiliary pair

Let `x_epsilon=(x_epsilon,1,0)` and `s_Cao,epsilon` be Cao's refined center and
scale.  For a speed written as `c_epsilon=W_epsilon log(1/epsilon)`, define
`(r_epsilon^*,s_epsilon^*)` by the exact two-equation system (Cao (3.36)):

    W_epsilon r^* log(1/epsilon)
      -kappa/(4 pi)[log(8r^*/s^*)+(p-1)/4]=0,
    Lambda_p (r^*)^(-(p+1)/(p-1))
      (epsilon/s^*)^(2/(p-1))=kappa.                     (16)

The second equation is solved exactly by

    s_epsilon^*=C_s(kappa,r_epsilon^*)epsilon.           (17)

Cao Proposition 3.13 gives

    |x_epsilon,1-r_epsilon^*|<=C epsilon^2,
    |s_Cao,epsilon-s_epsilon^*|
       <=C epsilon^3 |log epsilon|.                       (18)

Although Cao states the proposition for fixed `W`, its proof uses no
derivative of `W`: (3.31), the local Pohozaev identity, and Appendix-B Green
estimates are pointwise in `epsilon`.  On the fixed-mean bordered path,
(3.34) gives

    W_epsilon=c_epsilon/log(1/epsilon)
      ->kappa/(4 pi R).                                   (19)

Thus `W_epsilon` lies in one compact positive interval.  Every coefficient in
the estimates is uniform when `W` and the center range over those compact
sets; re-running the proof with `W=W_epsilon` gives (18) with one constant.
This compact-`W` uniformization is derived here and is not attributed to the
literal fixed-`W` statement of Proposition 3.13.

The missing mean-versus-center conversion is obtained from the same refined
density decomposition.  Write `xi=r-x_epsilon,1`.  The leading density is
radial in `(xi,z)`; the first correction `F+phi_epsilon^o` in Cao (3.29)--
(3.33) is odd in `xi` and is `O(epsilon)` relative to the core amplitude; the
remaining potential error is `O(epsilon^2|log epsilon|)`.  Hence

    integral xi zeta_epsilon r dr dz
       =integral xi^2 zeta_lead dr dz
          +x_epsilon,1 integral xi zeta_odd dr dz
          +O(kappa epsilon^3|log epsilon|)
       =O(kappa epsilon^2).                               (20)

The exact mean row (2) says the left side is `kappa(R-x_epsilon,1)`.
Therefore

    |R-x_epsilon,1|<=C epsilon^2,
    |R-r_epsilon^*|<=C epsilon^2.                         (21)

This is stronger than the convex-hull `O(epsilon)` bound and does not identify
the mean, maximum, and auxiliary radius.

Define the physical geometric core scale by the meridional area radius

    s_epsilon=sqrt(|Omega_epsilon|/pi).                   (22)

Cao Lemma A.2 writes the actual boundary as a normal graph over the
`s_Cao,epsilon` circle and proves only a uniform `O(epsilon)` graph at the
scope used here.  Without an explicit Fourier calculation of the constant
mode of `F+phi_epsilon^o` and the parameter row, that lemma does not prove
cancellation of the full first-order angular mean.  The safe area consequence
is therefore

    |s_epsilon-s_Cao,epsilon|
       <=C epsilon^2.                                     (23)

Combining (17)--(18), (21), and (23), and differentiating only the explicit
smooth coefficient `C_s(kappa,r)` with respect to `r`, gives the value theorem

    s_epsilon=C_s(kappa,R)epsilon+O(epsilon^2),
    delta(epsilon)=D_delta epsilon+O(epsilon^2),
    D_delta=C_s(kappa,R)/R.                               (24)

No derivative of either remainder in (18) or (24) is taken.

Start with one thin fixed-`(kappa,R)` seed furnished by the 0080 bordered
construction.  Its uniform rescaled inverse and the uniform `C^2`
positive-part map give a locally unique continuous solution chart on a fixed
radius in *relative* epsilon.  At the edge of that chart, the same uniform
inverse applies to the endpoint solution; local uniqueness identifies the
new chart with the old one on their nonempty overlap.  Iterating this
seed-plus-overlapping-relative-`epsilon` continuation in either direction
constructs one thin connected component whose epsilon projection contains a
punctured interval.  Thus connectedness is a continuation conclusion using
the seed, uniform chart radius, overlap, and local uniqueness—not a formal
consequence of having unrelated local charts.  This continuity plus (24),
rather than an unproved sign for `delta'`, is the Route-B coverage theorem.

## 5. Exact massive-fiber bracket coverage

To avoid collision with the Lane--Emden exponent `p`, write the rational ray
as `Q/P` with coprime positive integers `P,Q`.  Let the reviewed column
mismatch have opposite signs at

    k_N^-=k_*-h_N,       k_N^+=k_*+h_N,                  (25)

or use fixed opposite-sign endpoints when no shrinking bracket is needed.
Put

    n_1=NP,       n_2=NQ,
    delta_N^plusminus=k_N^plusminus/(NP).                (26)

Let `A_cov` dominate the uniform constant in (24), and define buffered
external endpoints

    epsilon_N^-
      =[delta_N^- -A_cov/N^2]/D_delta,
    epsilon_N^+
      =[delta_N^+ +A_cov/N^2]/D_delta.                   (27)

For all sufficiently large `N`, (24) makes the actual geometric scales at
these parameters lie below `delta_N^-` and above `delta_N^+`, respectively.
The connected fixed-`(kappa,R)` path therefore covers every geometric delta
between the target endpoints.  No inverse or monotonicity of `delta(epsilon)`
is used.

Let `rho_N` be a monotone envelope of the reviewed uniform massive graph/Riesz
error over the two fixed compact `k` neighborhoods.  To include the path
modulus, set

    E_N=max(rho_N,1/N),
    h_N=sqrt(sup_(M>=N) E_M).                             (28)

Then `h_N->0`, `rho_N=o(h_N)`, and the geometric-scale error, after
multiplication by `n_1=NP`, obeys

    NP |delta(epsilon)-D_delta epsilon|
       <=C/N=o(h_N).                                     (29)

Uniform reviewed graph/Riesz convergence gives, at every carrier in the
buffered interval, operator and simple-frequency error at most `rho_N`; no
derivative of that error is used.  Thus both the scale error (29) and the
spectral error are `o(h_N)`, smaller than the column endpoint signs.
Continuous labeling of the two simple Riesz eigenvalues follows from the
common graph domain and the continuous bordered path.  The frequency mismatch
therefore has opposite signs at (27), so IVT supplies an exact Cao crossing
between them.

The 0079/0082 review owns the final crossing claim.  The new 0083 supplier is:
the missing same-family connected carrier path exists, keeps `(kappa,R)` and
the limiting column fixed, exactly covers the massive geometric brackets, and
has path error `o(h_N)`.

## 6. Full response as a constrained dual functional

This section is independent of the scale proof.  On one regular core cell let
`h=[xi,Omega]` in the adopted DA sign convention, with `xi` compact,
divergence free, tangent to the cell, and in toroidal character
`ell=n_1-n_2`.  On the fixed translating slice

    V_h e_2=-[B h,e_2]-[B e_2,h].                        (30)

Let `a_1` be a smooth covector representative of the physical functional
`q -> Omega_KKS(e_1^#,q)` on this cell.  Direct integration by parts, with
compact support and the decaying self-adjoint Hodge solve, gives

    Omega_KKS(e_1^#,V_h e_2)=integral h dot K_12,         (31)

where, writing `v_2=B e_2`,

    K_12=v_2 dot grad a_1+(grad v_2)^T a_1
       -B^*[(grad e_2)^T a_1+e_2 dot grad a_1].           (32)

No boundary term at the vorticity edge occurs in this calculation because the
test displacement is compactly supported inside the regular cell, while
`B h` is still the global decaying Hodge velocity.  Extending the test to the
edge would require the separate coefficient-interface trace and is not used.
Substituting `h=[xi,Omega]` and integrating once more yields

    Omega_KKS(e_1^#,V_[xi,Omega] e_2)
       =integral xi dot F_12,
    F_12=(grad Omega)^T K_12+Omega dot grad K_12,          (33)

up to the single global sign fixed by choosing `h=+[xi,Omega]` rather than
`-[xi,Omega]`.  Finite physical rows add their explicit finite-dimensional
dual covectors before projection.

Let `P_DA` be the constrained Helmholtz projection on the
`(Delta m=1,ell)` slice.  Vanishing for every admissible compact `xi` implies

    P_DA F_12=0

in the quotient by gradients and finite-row annihilators.  Therefore

    curl(P_DA F_12)(x_*)!=0                               (34)

at one regular-cell point is a sufficient exposing predicate.  A ball on
which its selected component has fixed sign, followed by the local
divergence/right-row correction, gives a compact DA seed and a lower bound in
the dual norm.  Replacing `(e_1,e_2)` by the two diagonal pairs gives the same
construction for `F_11-F_22`.

Equations (31)--(34) are a genuine local calculation, not the unproved global
`ad/ad^*` formula removed from 0079.  However the actual column/Cao mode
functions and absolute KKS covector `a_1` have not been normalized, so (34)
is not evaluated here.  Nonvanishing, its high-`ell` scaling, two autonomous
histories, and the long gate-time return kernels remain open and enter
`N_response` below.

## 7. Finite subluminal charged hierarchy

For the carrier endpoint with the smallest geometric scale in (24), Cao's
local Pohozaev identity (3.34), transferred through (20), gives

    c_N
      =alpha[log N+C_ray]+r_N,
    alpha=kappa/(4 pi R),
    C_ray=log(8P/k^-)+(p-1)/4,
    r_N=o(1).                                             (35)

Here `k^-` is the smallest retained positive column endpoint.  Fix
`N_Pois` and `C_rem` so that `|r_N|<=C_rem` for `N>=N_Pois`.  For a fixed
Maxwell speed and strict margin `m_EM>0`, a sufficient explicit ceiling is

    N_EM^suf
      =exp((c_EM-m_EM-C_rem)/alpha-C_ray).                (36)

Every integer `N>=N_Pois` with `N<N_EM^suf` is subluminal.  Conversely the
matching lower Pohozaev bound gives a finite necessary ceiling, so no
`N->infinity` charged transfer exists at fixed `c_EM`.

Define

    N_0=max(N_graph,N_IVT,N_response,N_Pois),
    N_min=floor(N_0)+1.                                   (37)

The exact finite-window test is

    N_min<N_EM^suf,                                      (38)

or equivalently the transparent sufficient foundation-scale inequality

    c_EM>m_EM+C_rem+alpha[log N_min+C_ray].               (39)

The path construction proves `N_graph,N_IVT<infinity`.  It does not prove
`N_response<infinity`, because absolute KKS normalization and a nonzero
quantitative response remain open.  Hence (38) is now an explicit finite
criterion but is not established for the charged analyzer.  A later
`O(g^2)` speed correction must fit inside `m_EM`; charged existence and
persistence are not inferred here.

## 8. Route verdicts and continuation

**Route A (differentiated fixed-domain map): blocked at a named proof
construction.**  The coefficient (11) and exact logarithmic row structure are
known, but the full differentiated Green/interface/Hanzawa/far-field ledger
needed for (15) has not been proved.  No `C^1` scale theorem is claimed.

**Route B (derivative-free coverage): established at author scope** under the
0080 bordered continuity input.  Cao's exact auxiliary system, Proposition
3.13, the derived mean-versus-center row (20)--(21), and the geometric scale
(24) give coverage without differentiating a source remainder.  Equations
(27)--(29) provide the exact massive bracket and `o(h_N)` error.

**Route C (finite charged hierarchy): blocked at a named construction.**
Equations (35)--(39) quantify the finite window and prove why the asymptotic
sequence is not a charged theorem.  The missing object is a physically
normalized nonzero mixer/diagonal response with two-sided gate bounds, hence
a finite `N_response`, followed by one integer satisfying (38).

If (38) is empty after that response is constructed, the two active
failure-derived routes are: (i) fixed `ell` with `J=Theta(delta^(-1))`, which
needs a joint high-`J`/thin-ring Riesz error below `J^(-2)` and physical KKS
normalization; and (ii) a genuinely slower compact Euler carrier family with
the same coercive massive-column structure.  Neither is claimed here.

The strongest positive 0083 bridge is the actual fixed-`(kappa,R)` connected
Cao path and exact geometric-bracket coverage.  It closes the carrier-path
dependency in the 0079 rational-ray construction at author scope.  Physical
KKS coupling, analyzer control, nonlinear persistence, P2/P4, electron,
neutrino, quantum, and relativity claims remain active.
