# Finite-`delta` enlarged KKS/Feshbach calculation

## 1. Physical DA realization and sign

Let `Omega_delta` be one fixed polynomial Cao ring, and let

    eta_xi=curl(xi cross Omega_delta),  div xi=0,            (1)

with compactly supported smooth `xi`. Complete these tangents after the
circulation, impulse, center, and stabilizer reductions in

    ||eta||_E^2=||BiotSavart eta||_2^2+||eta||_H^-1^2,
    ||eta||_G^2=||eta||_E^2+||A_delta eta||_E^2.             (2)

Convergence in (2) implies distributional convergence. Since every (1) is
supported in the closed vorticity core, the closure has no passive exterior
vorticity. Its velocity and pressure remain whole-space fields.

Use the right-reduced convention of 0005,

    X_xi=-ad*_xi m,
    Omega_KKS(X_xi,X_zeta)=-<m,[xi,zeta]>
      =rho_0 integral Omega_delta dot (xi cross zeta) dx,   (3)

so `i_(X_H) Omega_KKS=dH`. The energy--impulse Hessian `L_delta`
and Poisson map satisfy `A_delta=J_delta L_delta` on the quotient. This fixes
the sign and density of every finite block below; an Euclidean projection is
not used.

## 2. Exact physical normalization of the `m=0` Kelvin modes

For the limiting compact column, put

    Phi=2 Omega W,  partial_r^*=partial_r+1/r,
    H_k=-partial_r partial_r^*+k^2.                         (4)

The decaying exterior field is part of the domain. For `k!=0`, the physical
axisymmetric eigen-equation is exactly

    H_k u_n=lambda_n(k) Phi u_n,
    lambda_n=k^2/sigma_n^2,                                (5)

with `u_n(0)=0` and exterior `K_1(|k|r)` matching. Normalize

    integral Phi |u_n|^2 r dr=1.                           (6)

The Gallay--Smets first-order equations give, for `s=i sigma_n`,

    u_z=i partial_r^*u_n/k,
    u_theta=i W u_n/sigma_n.                               (7)

Thus the exact kinetic norm (apart from the common azimuthal/axial Fourier
volume factor) is

    ||v_n||_2^2=sigma_n^-2 M_n,
    M_n=integral (Phi+W^2)|u_n|^2 r dr.                    (8)

For the decreasing compact profile,

    W(r)<=2 Omega(r),

because `2 Omega` is the weighted disk average of `W`. Equations (6)--(8)
therefore imply

    1<=M_n<=2.                                             (9)

The energy-unit mode is

    e_n^E=|sigma_n| v_n/sqrt(M_n).                         (10)

The positive wave equation behind (5) is

    H_k partial_t^2 u+k^2 Phi u=0.                         (11)

Consequently every nonzero axisymmetric pair has positive constrained-energy
signature. If `e_n^+` denotes its positive-frequency energy-unit vector, the
Hamiltonian identity gives

    Omega_KKS(e_n^+,conjugate(e_n^+))=i/sigma_n            (12)

up to the fixed real-to-complex factor used for both translation and Kelvin
coordinates. Hence `sqrt(|sigma_n|)e_n^E` is the corresponding unit KKS
coordinate. Equation (12), rather than radial `L2` normalization, is used in
the enlarged block.

## 3. Discreteness, Weyl scale, and spacing

Multiplying (5) by `r` gives the separated generalized Sturm--Liouville
problem

    -(r u')'+(r^-1+k^2 r)u=lambda r Phi u.                 (13)

The weight is positive on `(0,a)` and vanishes like `(a-r)^p` at the compact
edge. The Liouville length

    L_Phi=integral_0^a sqrt(Phi(r)) dr                     (14)

is finite and nonzero. Dirichlet--Robin bracketing away from the two endpoints
and then letting the collars shrink gives

    sqrt(lambda_n(k))=pi n/L_Phi+O(1),                     (15)

uniformly for `|k|a<=k_0`. The singular `r^-1` term at the axis and the
vanishing edge weight change the `O(1)` phase, not the Weyl coefficient.
In particular, for constants independent of large `n`,

    c n^2<=lambda_n(k)<=C n^2,
    sigma_n(k)=|k|L_Phi/(pi n)+O(|k|/n^2),                 (16)

and

    sigma_n-sigma_(n+1)=|k|L_Phi/(pi n^2)+O(|k|/n^3).      (17)

This is a discrete internal compact-core accumulation at zero. It is not the
passive exterior essential endpoint excluded by (1)--(2).

## 4. The `W/2` coefficient and the correction to the scratch estimate

The exact first toroidal scalar block found in 0044 contains

    C[psi_0 exp(i alpha)]|_(m=0)=W/2.                      (18)

Let

    a_n=integral (W/2) conjugate(u_n) r dr.                (19)

Since `Phi=2 Omega W`, weighted Cauchy--Schwarz and (6) give the explicit
bound

    |a_n|^2
      <=integral W/(8 Omega) r dr=:C_W^2<infinity.         (20)

After physical energy normalization, (10) yields

    |d_n^E|<=C_W |sigma_n|/sqrt(M_n)
             <=C |k|/n.                                   (21)

The preliminary `O(k n^-2)` claim is **not** used. It attempted to write
`W/2=Phi q`, `q=1/(4 Omega)`, and transfer `H_k` to `q`; but `q(0)!=0`, while
the radial `m=0` form contains the `r^-1` axis term. Thus `q` is not in that
form domain and the transfer is not licensed without an axis corrector. The
unconditional physical estimate is (21). It is already square summable and
is sufficient for the finite-window construction.

The remaining pieces of the full `C_1` forcing are the terms in 0044
`P_1N_0+P_0N_1`. On the translation mode they are fixed compact-core
functions plus whole-space harmonic pressure fields. The cylindrical
Biot--Savart estimate and the same energy normalization give

    |c_n^E|<=C |k|/n                                      (22)

provided their scalar source divided by `sqrt(Phi)` is in `L2(r dr)`.
For `W=U_+^p`, `p>=2`, every undifferentiated `W` term satisfies this directly;
the apparently marginal `W'` terms are first integrated by parts in the
physical KKS/Hessian form, where the translation field and `W'` vanish to the
required order at the edge. The boundary term is zero for `p>=2`. This proves
(22) for the complete source-defined first block, while (20) is the explicit
constant for its noncancelling `W/2` piece.

In KKS-unit coordinates the Hessian coefficient is obtained from (12):

    h_n=sqrt(tau_l sigma_n) b_n^E,                          (23)

where `b_n^E` is the energy-metric off-diagonal coefficient. This conversion
is important: KKS normalization changes individual matrix entries, but the
Schur product and the spectral window are invariant.

## 5. Actual second Cao core cell and the logarithmic split

Let

    q=s_epsilon/x_epsilon=delta

be the core-to-center ratio in the exact rescaling (3.11). Equation (3.12) of
the Cao source becomes, with `x=y_1` and `h=1+q x`,

    -Delta w+[q/(1+q x)] partial_x w
       =(1+q x)^2 (w_+)^p.                                (24)

This equation, rather than the `O(epsilon^2|log epsilon|)` estimate printed in
the source, fixes the higher local jets. Write

    w=U+q V+q^2 Z+o(q^2).                                 (25)

With

    L_U=-Delta-p U_+^(p-1),                               (26)

the first and second cells are exactly

    L_U V=2x U^p-partial_x U,                              (27)

    L_U Z={p(p-1)/2}U_+^(p-2)V^2
          +2p x U_+^(p-1)V+x^2U_+^p
          -partial_x V+x partial_x U.                     (28)

Here `V` is odd in `x` and even in the other core coordinate; `Z` is even in
both. Translation is the entire bounded kernel of (26). Centering removes the
odd kernel, and the even equation (28) is uniquely invertible with the fixed
far logarithmic coefficient. The circulation and center conditions are
differentiated along with (24); equivalently they determine the physical
scale `s_epsilon` and center `x_epsilon` after the unique dimensionless cells
have been found.

Any putative local `q^2 log q Z_log` term obeys

    L_U Z_log=0.                                           (29)

It is even, centered, and has zero far logarithmic coefficient at fixed
circulation, so Cao's limiting nondegeneracy forces

    Z_log=0.                                               (30)

Thus the **profile** second core jet has no logarithmic cell. The actual
whole-space velocity/operator does have a `q^2 log q` term: it comes from the
small-`rho` Green kernel and is the local-induction term whose translation
compression was independently evaluated in 0044. Splitting the exact Green
integral into `|x-y|<r_0` and its regular complement gives

    A_delta=A_col(k_delta)+delta C_1^curv
       +delta^2 log(delta) C_(2,log)+delta^2 C_2+R_delta,
    ||R_delta||_(G->E)=o(delta^2),                         (31)

for an integer `p>=3`. The straight-column operator in (31) already contains
the exact longitudinal wave number `k_delta=l delta`; this prevents the
`i l delta` derivatives from being counted a second time as curvature.
Here `C_(2,log)` is the universal local-induction
operator, while `C_2` is the sum of: the cell (28), the finite part of the
exact Green integral (2.2), the second toroidal connection, and the exact
second Leray expression of 0044. Calderon--Zygmund estimates on the fixed
core chart, the `C^2` Nemytskii map `w -> (w_+)^p`, and the exterior harmonic
trace estimate give the graph remainder in (31). For `p=2` the cell remains a
valid weak second derivative, but this argument earns only a fractional graph
remainder; 0048 therefore selects any fixed Cao member with integer `p>=3`.

This constructs the source-defined second jet as explicit invertible cell
problems. It does not replace it by the universal filament matrix.

## 6. Resonant scale and the finite window

In core time units the positive bending frequency is

    tau_l(delta)=c_l delta^2 L_delta+O(delta^2),
    c_l=l sqrt(l^2-1),   L_delta=log(1/delta).              (32)

Equations (16) and (32) place a crossing at

    n_*(delta)=L_Phi/(pi c_l delta L_delta)+O(L_delta^-1).
                                                                  (33)

The full first curvature map from the translation pair into the energy-unit
`m=0` sector has norm `O(delta |k|)=O(delta^2)` and coefficients

    g_n(delta)=delta c_n^E=O(delta^2/n).                   (34)

At (33),

    |g_n|=O(delta^3 L_delta),
    sigma_n-sigma_(n+1)=Theta(delta^3 L_delta^2),          (35)

so

    |g_n|/(sigma_n-sigma_(n+1))=O(1/L_delta).              (36)

The second core jet does not spoil this comparison. After Liouville
conjugation, its diagonal `m=0` part is an order-zero coefficient operator on
the finite interval. The two-term WKB diagonal matrix element has the form
`beta_infinity+beta_1/n+O(n^-2)`; its constant part is absorbed into the exact
renormalized `sigma_n(delta)`, and successive finite-part shifts are
`O(delta^2 n^-2)`. At (33) this is `O(delta^4L_delta^2)`, one factor `delta`
below (35). Off-diagonal entries obey the corresponding Fourier-coefficient
decay supplied by the `p>=3` cell regularity. Thus (35)--(36) apply to the
source-defined renormalized levels, not only the unperturbed column list.

This is the decisive finite-`delta` result: for sufficiently small `delta`,
the frozen four-coupling window contains at most one positive-frequency
`m=0` mode and its real/Hamiltonian conjugate. Therefore `P_delta` has real
symplectic rank at most four: one bending oscillator and one Kelvin
oscillator. Its complexification contains the associated two positive- and
two negative-frequency eigenvectors. Its identity may change with `delta`,
and it does not converge to one fixed rank-two projection as `delta->0`.

Away from that window, the scalar tail in the translation Schur complement
obeys

    sum_n |g_n|^2/|sigma_n-tau_l|
       <=C delta^3 L_delta,                                (37)

where the sum is split at `n_*` and the one window mode is removed. This is
`O(delta)` relative to the bending scale `delta^2L_delta`.

## 7. Quantitative nonresonance measure

On a dyadic shell `Delta_j=(d/2,d]`, `d=2^-j`, differentiate the leading
crossing equation. At a crossing,

    |partial_delta(sigma_n-tau_l)|
       >=c delta L_delta                                  (38)

for all sufficiently small `delta`; the `k` derivative of `lambda_n(k)` and
the finite part of (32) are lower order. By (34), the bad interval belonging
to one crossing has length at most

    C delta^2.                                             (39)

More explicitly, its spectral half-width is `4|g_n|=O(delta^3L_delta)`
at (33), and division by (38) gives `O(delta^2)`. The number of crossing
indices on the shell is `O((delta L_delta)^-1)`. Hence

    |B_j|<=C d/L_d,       |B_j|/|Delta_j|<=C/L_d ->0.      (40)

The union of good sets therefore has asymptotic relative measure one at the
thin-ring endpoint. On a good `delta`, the enlarged window is empty; on a bad
one, retaining its unique conjugate Kelvin pair produces the finite enlarged
block. This is a measure alternative, not deletion of exact resonances.

## 8. Enlarged Hamiltonian block and exact resonances

Let `z_0` denote a positive bending coordinate and `z_n` the selected
axisymmetric coordinate, each normalized by (3) and (12). Up to the real
conjugate block, the quadratic Hamiltonian is

    H_2=tau_l |z_0|^2+sigma_n |z_n|^2
          +2 Re(h_n z_0 conjugate(z_n))+R_2.               (41)

Both diagonal signs are positive: the first follows from the length--impulse
Hessian for `l>=2`, and the second from (11). Equations (31), (34), and (36)
give

    |h_n|^2<tau_l sigma_n                                 (42)

for sufficiently small `delta`, including at exact resonance. Thus (41) is
positive definite. Its two frequencies are

    nu_+/-=(tau_l+sigma_n)/2
       +/-sqrt(((tau_l-sigma_n)/2)^2+|h_n|^2),             (43)

and remain real and nonzero. The Euler generator has the simple imaginary
pairs `+/- i nu_+`, `+/- i nu_-`, except at the codimension-one equality
where the exact finite `C_2` entry must be used to split an accidental double
root. Opposite Krein signs would replace the plus sign under the square root
by a minus sign; (11) rules out that mechanism for this `m=0` channel.

Define `P_delta` using the KKS duals of the translation pair and every mode in
the four-coupling window. On the complement, (17), (31), and (37) give the
weighted graph resolvent after the diagonal `m=0` finite correction is
absorbed into `sigma_n(delta)`. Analytic Fredholm reconstruction then yields

    Pi_delta=(2 pi i)^-1 integral_Gamma
                   (z-A_delta)^-1 dz                      (44)

as a finite-rank projection bounded on both `X_E` and `D_G(A_delta)`. Its
range consists of compact-core DA vorticity and the exact global exterior
velocity. The finite part of its matrix is defined by (28), (31), and the
Feshbach formula, not by (2).

The alternative infinite-block route is also controlled: (34) is in weighted
`ell^2_1`, and (37) is precisely the one-step tame Schur estimate. Repeating
the block elimination loses one Kelvin weight but gains one factor `delta`.
Indeed, the first homological generator has coefficients

    x_n=g_n/(sigma_n-tau_l),
    sum_(n notin I_delta)|x_n|^2<=C delta/L_delta.          (45)

For `n<n_*`, each summand is `O(delta^2)`; for `n>n_*`, it is
`O((n L_delta)^-2)`. Thus the transformation is
`O(sqrt(delta/L_delta))` in the physical energy sequence norm. It converges
on the scale relevant to (32) unless a finite resonant pair is retained.
Route B is not needed to establish the finite-window spectral family, but it
supplies the same complement when parameters are followed through successive
window exchanges.

## 9. Nonlinear branch bridge and its remaining construction

For every good `delta`, and for every resonant `delta` where (43) is simple,
(44) supplies the exact finite-dimensional kernel/range splitting for the
relative-equilibrium map

    F(V,c,Omega)=P_L[(V dot grad)V]
                  -c partial_z V-Omega[R,V].               (46)

The solid-torus displacement chart preserves compact support and the DA leaf;
the impulse/center slice removes translation and the physical rotation
quotient leaves one real amplitude. The KKS crossing is nonzero because the
mode in (43) has nonzero positive signature.

The elementary tame estimate is not the obstruction. For `s_0>5/2`, the
whole-space Leray projector has order zero and the Sobolev product estimate
gives

    ||D^2F(V)[h_1,h_2]||_(s-1)
      <=C_s(||h_1||_s||h_2||_s0
             +||h_1||_s0||h_2||_s)                       (47)

for the fixed `delta` chart. Smooth the divergence-free displacement, apply a
collar Bogovskii correction, and subtract the finite impulse/center
components. With the fixed spectral projector, replace this smoother by

    S_N^Pi=Pi_delta+(1-Pi_delta)S_N(1-Pi_delta).            (48)

Finite-rank commutators have the required tame bounds, and exterior velocity
is reconstructed after smoothing by whole-space Biot--Savart. For a chosen
integer `p` larger than the finite Sobolev index used in the iteration, the
free-boundary coefficient `(w_+)^p` has the required finite differentiability.
Thus (47)--(48) construct the local tame chart and compatible smoother.

The actual obstruction appears in the inverse at the nonlinear harmonics.
In straight-core action--angle variables, a cross-sectional mode `m` in the
`j`th generated ring harmonic has transport divisor

    D_(m,j)(r)=m Omega(r)-j l Omega_pat,                    (49)

where `Omega_pat=tau_l/l+O(delta^2)` is the pattern speed. Since
`Omega(r)` continuously covers `[Omega(a),Omega(0)]` and
`Omega_pat=Theta(delta^2 L_delta)`, every fixed nonzero `delta` admits
integers `j` satisfying

    m Omega(a)/(l Omega_pat)
       <=j<=m Omega(0)/(l Omega_pat).                       (50)

The interval in (50) has length `Theta((delta^2L_delta)^-1)` and therefore
contains integers. At those orders (49) has a genuine internal critical
radius. It belongs to the compact DA core and cannot be removed using the
excluded passive exterior spectrum.

This does not obstruct the finite-order construction allowed by the carrier's
regularity. Fix `N<=p-2` (or first choose one Cao exponent `p>=N+2` and then
keep that carrier fixed), and choose `delta_N` so that

    N l |Omega_pat|<Omega(a)/2.                            (51)

Then every `m!=0`, `|j|<=N` transport divisor is bounded away from zero;
the `m=0` equations use the enlarged Riesz/Feshbach inverse above. Induction
on the quadratic equation (46) constructs source-defined coefficients

    V(a)=V_0+sum_(r=1)^N a^r V_r,
    c(a)=c_0+sum_(r=1)^N a^r c_r,
    Omega(a)=Omega_pat+sum_(r=1)^N a^r omega_r,             (52)

on the same DA leaf, with circulation, impulse, center, and phase conditions,
and

    ||F(V(a),c(a),Omega(a))||_(s-N)
       <=C_(N,delta)|a|^(N+1).                             (53)

Thus 0048 earns a relative-equilibrium expansion through every preassigned
finite order permitted by the one fixed carrier, and an arbitrary
preassigned order can be earned by choosing the polynomial exponent once
before constructing that carrier. This is not an all-orders statement for a
single finite `p`. For fixed `delta`, the first unavoidable critical harmonic is

    j_crit=Theta((delta^2 L_delta)^-1),                    (54)

so its forcing is beyond all orders in the thin-core asymptotic hierarchy.

An exact branch now requires a **critical-layer transparency theorem**: at
every radius where (49) vanishes, the recursively generated numerator must be
divisible by `D_(m,j)` in the physical graph domain, or an additional
coadjoint/Casimir profile variable must cancel its trace. The two executable
continuations are (i) prove this divisibility from the vorticity transport and
Bernoulli equations in streamline coordinates, retaining the fixed Cao leaf,
or (ii) enlarge the nonlinear unknown by the permitted streamline-vorticity
function and solve the resulting infinite range conditions by a tame
Lyapunov--Schmidt scheme. A generic Nash--Moser inverse cannot cross (50)
without one of these mechanisms.

Thus 0048 closes the finite-`delta` enlarged spectral mechanism and its exact
linear bridge, including true resonances, and constructs the nonlinear branch
through every preassigned carrier-regularity-compatible order under (51). It
does **not** yet establish the
exact nonlinear Cao rotating branch: the next executable dependency is the
critical-layer transparency/range construction following (54), followed by
convergence and physical period/KKS-action reconstruction. No quantum or
particle inference follows.

## 10. Route verdicts

- Route A is **established** at the linear level: the coupling decay, finite
  window, asymptotically full-measure alternative, enlarged positive
  Hamiltonian block, graph-domain Riesz reconstruction, and its fixed-`delta`
  quantifiers are explicit.
- Route B is **established as a tame linear continuation**: the infinite tail
  Schur sum is one weight smoothing and smaller than the bending scale by
  `O(delta)`. It is not claimed as a completed nonlinear KAM theorem.
- Route C is **established for the physical axisymmetric resonance**: equal
  positive Krein signatures turn a crossing into an elliptic avoided crossing,
  not a Hamiltonian--Hopf quartet.
- The nonlinear solid-torus bridge is **blocked by the named construction
  following (54)**. Equations (47)--(53) establish the tame chart and every
  carrier-regularity-compatible fixed-order branch; only the fixed-`delta`
  critical-layer range condition and convergence remain. Formal CR
  consequences are not called an exact branch.

The parent P253/LP2/P4 objective remains active. The 0044 fixed-rank
refutation and the limited scope of the universal filament matrix are
unchanged.
