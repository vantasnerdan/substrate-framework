# P253/0048: finite-window KKS/Feshbach continuation for the Cao ring

Status: **README-only preregistration awaiting central registration and schema
activation.** Until `0048/activation-schema.exit` exists and contains exactly
`0`, this README is the sole 0048 artifact. Root owns `proposal.yaml`, the
central attempt registry, activation receipts, commits, review, and promotion.
No derivation body, source download, verifier, numerical run, API, or claim is
authorized before activation.

## Failure-derived objective and inherited boundary

Attempt 0044 established the compact-core dynamically accessible (DA) support
closure, the simple finite-core `m=1` translation response, the nonzero
crossing

    <psi_0,T'(0)psi_0>_(r dr)=-2F(a)<0,                    (1)

the corrected exterior streamfunction matching, and the exact first toroidal
Leray blocks. It also found an actual `m=1 -> m=0` forcing with leading scalar
piece `W/2`. The physical `m=0` Kelvin spectrum accumulates at zero. This
refutes **only** a uniform fixed-rank, rank-two Riesz reduction of the curved
ring. It does not refute a finite-`delta` enlarged Riesz family, the Cao
carrier, or a nonlinear rotating branch.

For one fixed azimuthal ring number `l>=2` and every sufficiently small but
nonzero core ratio `delta=a_epsilon/R`, construct either:

1. a finite-dimensional physical KKS/Feshbach block containing the two real
   `m=1` translations and every `m=0` mode in a quantitatively declared
   coupling window, with a graph-domain Riesz projection and a simple signed
   bending pair; or
2. a tame Hamiltonian infinite-block normal form on the same DA domain when
   no useful nonresonant finite-window family exists.

Then prove the exact analytic bridge from that spectral family to one
compact-vorticity, finite-energy, nonaxisymmetric solid-torus relative
equilibrium of Euler. A finite matrix without domain reconstruction or a
formal bifurcation implication does not meet the objective.

The universal thin-filament compression

    Gamma log(R/a)/(4*pi*R^2)
      [[0,l^2-1],[-l^2,0]]                                 (2)

is retained only as a leading asymptotic comparator. It is not the full
finite-core matrix, does not contain the `m=0` response, and supplies no Cao
Riesz theorem by itself.

## Fixed carrier, physical leaf, and graph domains

Fix the polynomial-profile, no-swirl Cao ring used in 0040 and 0044, with
`p>=2`, circulation `Gamma>0`, major radius `R`, core radius `a_epsilon`, and

    delta=a_epsilon/R,             k_delta=l*delta,
    omega_core=Gamma/(4*pi*a_epsilon^2).                   (3)

Let `Omega_delta` and `U_delta` be the exact compact vorticity and translating-
frame velocity. Start with the integrable DA tangents

    eta_xi=curl(xi cross Omega_delta),
    div xi=0,   xi in C_c^infinity(R^3),                   (4)

modulo the stabilizer kernel `eta_xi=0`, and impose the circulation, axial
impulse, and centering slices used in 0040. The vorticity in (4) remains
distributionally supported in the closed core; its velocity
`v_eta=BiotSavart(eta)` and exterior harmonic field remain global.

Freeze the physical KKS form, with mass density `rho_0`, as

    Omega_KKS(eta_xi,eta_zeta)
      =rho_0 integral Omega_delta dot (xi cross zeta) dx,   (5)

with its sign checked against `i_X Omega_KKS=dH`. Let `L_delta` be the exact
constrained energy--impulse Hessian and

    A_delta=J_delta L_delta                                (6)

the full whole-space Leray/Hodge Euler generator. Complete (4), after the
finite symmetry quotient, in the energy space

    ||eta||_E^2=||v_eta||_L2^2+||eta||_H^-1^2              (7)

and use the closed graph domain

    D_G(A_delta)={eta in X_E : A_delta eta in X_E},
    ||eta||_G^2=||eta||_E^2+||A_delta eta||_E^2.            (8)

If (5) is weak on this completion, construct its strong symplectic chart on
the spectral subspace and state the precise continuous dual pairing used for
the complement. An arbitrary `L2`-vorticity or full-enstrophy exterior space
is not substituted for (7)--(8). All operator expansions retain the global
Biot--Savart field, the exterior Dirichlet-to-Neumann map, pressure/Leray
terms, the core boundary, and the distributional support closure.

## Translation pair and the finite coupling window

Let `t_c,t_s` be the real `m=1` translation pair, normalized by their physical
KKS/impulse duals, and let `P_tr` be the corresponding biorthogonal
projection. In the limiting `m=0` sector use the exact compact positive
operator found in 0044,

    B_k=Phi^(1/2) H_k^(-1) Phi^(1/2),
    H_k=-partial_r partial_r^*+k^2,                         (9)

with physical boundary conditions and eigenpairs
`B_k phi_n=mu_n(k)phi_n`. Its Kelvin frequencies are

    sigma_n^+/- (k)=+/- |k| sqrt(mu_n(k)),
    mu_n(k)>0,       mu_n(k) -> 0.                         (10)

Reconstruct each vorticity--velocity eigenvector `e_n^+/-` in (7)--(8) and
normalize it in the exact KKS/Krein metric, not merely in radial `L2`.

Derive from the full first curvature operator `C_1`, including its Leray and
base-profile pieces, the graph-norm couplings

    c_n^+/- (delta)=Omega_KKS(e_n^-/+,C_1 t_+/-),           (11)

whose scalar source contains `W/2`. Define the physical bending scale
`tau_l(delta)` by the actual translation compression, not by declaring (2)
exact, and define the earned interaction size

    e_delta=max{||Q_0 delta C_1 P_tr||_(G->E),
                ||P_tr delta C_1 Q_0||_(G->E),
                ||delta^2 C_2||_(G->E), ||R_delta||_(G->E)}. (12)

The preregistered finite window is

    I_delta={ (n,s) : s in {+,-},
              |sigma_n^s(k_delta)-tau_l(delta)|
                 <=4 e_delta |c_n^s(delta)|_* },           (13)

where `|c|_*:=min(1,|c|/c_ref)` and `c_ref` is the norm of the full
translation-to-`m=0` block fixed by (11). If exact normalization shows that
the coupling enters quadratically, replace the right side only by the derived
Schur scale `4 e_delta^2 |c_n|_*^2`; record the derivation before selecting the
window. The factor four is fixed now. Prove that `I_delta` is finite for every
fixed `delta>0`, contains every pole whose omission would invalidate the
contour bound, and is invariant under the real/Hamiltonian conjugations.

Set

    P_delta=P_tr + sum_((n,s) in I_delta) P_(n,s),
    Q_delta=1-P_delta.                                     (14)

The enlarged projection must be KKS-biorthogonal and bounded simultaneously
on `X_E` and `D_G(A_delta)`.

## High-mode coupling decay and tail control

The first analytic rung is the exact decay of (11). Put the complete radial
forcing, after pressure reconstruction and KKS normalization, in the form
`g_delta`; its leading scalar term is `W/2`, but no other `C_1` contribution
is discarded. Derive the largest integer or fractional `q=q(p)>0` for which

    g_delta in D(B_k^(-q))

uniformly on the declared small-`delta` interval, including all compatibility
conditions at the Lane--Emden free boundary. Repeated self-adjoint transfer
must then give the earned bound

    |c_n^s(delta)| <= C_q(delta) mu_n(k_delta)^q            (15)

or its exactly equivalent Sturm--Liouville eigenvalue form. Derive the Weyl
law and KKS normalization before translating (15) into a power of `n`.
Determine the maximal decay actually allowed by `W=U_+^p`; do not assume
smooth compact support across the free boundary. Establish an explicit tail
estimate, for example

    sum_(n>N,s) |c_n^s|^2/dist(z,sigma_n^s)^2
       <= Tail_q(N,delta),       Tail_q ->0,                (16)

uniformly on the proposed contour. Failure of a guessed exponent activates a
weighted/fractional-domain version of (15), not a claim that coupling does not
decay.

## The actual Cao second core jet

Derive, from the exact Cao free-boundary equation and not from the filament
model, the second core/operator jet in the graph norm (8):

    U_delta=U_0+delta U_1
      +delta^2 log(delta) U_(2,log)+delta^2 U_2+r_delta,
    ||r_delta||_(profile graph)=o(delta^2),                 (17)

and consequently

    A_delta=A_0+delta C_1
      +delta^2 log(delta) C_(2,log)+delta^2 C_2+R_delta,
    ||R_delta||_(G->E)=o(delta^2).                          (18)

The cell equations must include the moving core boundary, circulation and
impulse constraints, center gauge, toroidal metric/connection, exact exterior
Hodge field, and both derivatives of the field-dependent Leray projector.
State the regularity earned from the chosen finite `p` and use a fractional
graph remainder if a classical second boundary derivative is unavailable.
The 0044 formulas for `P_1,P_2` are inputs to be rederived in this physical
operator convention; they are not a substitute for the Cao profile cells in
(17).

## Route A: quantitative nonresonance and enlarged finite block

Prove differentiable or analytic dependence of `sigma_n(k_delta)`,
`tau_l(delta)`, and (11) on `delta` at the exact regularity earned above. For
the dyadic interval `Delta_j=(2^(-j-1),2^(-j)]`, define the bad set

    B_j={delta in Delta_j :
         exists (n,s) notin I_delta with
         |sigma_n^s(k_delta)-tau_l(delta)|
            <4 e_delta |c_n^s(delta)|_*}.                  (19)

Derive a transversality bound for each crossing and use (15) to prove either

    |B_j| <= C 2^(-j) epsilon_j,       epsilon_j ->0,       (20)

or a summable full-measure estimate of equal strength. On the good set prove
the graph-domain complement resolvent estimate on a contour of radius fixed
by the right side of (19), form the exact Feshbach map

    K_delta(z)=P_delta(A_delta-z)P_delta
      -P_delta A_delta Q_delta
       [Q_delta(A_delta-z)Q_delta]^(-1)
       Q_delta A_delta P_delta,                             (21)

and evaluate all entries through the second jet (18). Exact resonances are
retained inside `P_delta`, not excluded by fiat. Classify their finite
Hamiltonian/Krein blocks and prove whether a simple imaginary bending pair
survives.

## Route B: tame infinite-block normal form

If crossings are too dense for (20), retain all `m=0` modes. On the weighted
sequence space dictated by (7)--(8) and (15), write the coupled generator as

    A_diag + S_delta + N_delta,                             (22)

where `A_diag` contains the exact Kelvin frequencies, `S_delta` is the
smoothing translation--axisymmetric coupling, and `N_delta` contains the
controlled second jet and nonlinear-coordinate commutators. Construct a
KKS-symplectic block transformation by a convergent normal form, or a tame
Nash--Moser inverse with its exact derivative loss. The target estimate is a
uniform right inverse on the symmetry-reduced complement with a stated
weighted graph loss, not a boxed eigenvalue. A finite resonant cluster may be
split off at each step, but the tail cannot be discarded without (16).

## Route C: exact-resonance Hamiltonian block

When (19) selects one or several true resonances, calculate the enlarged
finite KKS and constrained-energy matrices directly. Determine Krein signs,
Jordan structure, and the splitting produced by (18). A hyperbolic quartet is
a route-scoped obstruction to the desired Lyapunov-center branch but is also
a concrete Euler instability candidate; an elliptic signed pair activates
the nonlinear construction below. Neither outcome is inferred from the
universal matrix (2).

## Riesz family and nonlinear solid-torus bridge

For whichever spectral route succeeds, construct the physical graph-domain
Riesz family

    Pi_delta=(2*pi*i)^(-1) integral_(Gamma_delta)
                (z-A_delta)^(-1) dz,                       (23)

prove its rank, real/Hamiltonian symmetries, KKS nondegeneracy, graph
boundedness, and reconstruction to compact-core vorticity plus global
finite-energy velocity. Track all quantifiers: a full-measure or Cantor
sequence of `delta` values is acceptable only if it accumulates at zero and
each selected Cao carrier is exact.

Then use the volume-preserving solid-torus displacement chart, the fixed
circulation/impulse/center slice, and physical axial-rotation quotient to
solve the exact relative-equilibrium equation

    F(V,c,Omega)=P_L[(V dot grad)V]
                 -c partial_z V-Omega[R,V]=0.              (24)

The spectral result must supply the kernel, range, and KKS crossing; (18),
(21) or (22) must supply the tame complement inverse and nonlinear remainder.
Execute Crandall--Rabinowitz/Lyapunov--center when the inverse is bounded, or
an explicitly loss-balanced Nash--Moser scheme otherwise. The target branch
has compact deformed toroidal vorticity, whole-space finite kinetic energy,
physical `SO(2)` period and `C_l` isotropy, and nonzero KKS/angular-impulse
pullback. Prove those properties from the reconstructed field. This attempt
does not claim nonlinear stability, quantization, half-integer spin, an FR
sign, an electron, or a neutrino.

## Sources, verification, and route verdicts

Primary source inventory begins with Cao--Lai--Qin--Zhan--Zou,
arXiv:2206.10165v2, for the exact compact-core carrier and its earned
regularity; Gallay--Smets, arXiv:1805.05064v3, only for source-scoped column
equations; and the 0044 source receipts for the exterior matching and
thin-filament comparator. Full papers remain in `/tmp/primary-source-cache`.
Every receipt records URL, version, SHA-256, and exact theorem/equation
location. Cao's elliptic uniqueness operator is not used as the full Euler
generator, and ambient passive exterior-vorticity spectrum is not imported
into the DA space.

The strongest analytic oracle is agreement of three independently exposed
objects: the high-mode bound (15)--(16), the graph-norm Cao jet (17)--(18), and
the poles/zeros of the reconstructed resolvent and Feshbach map (21)--(23).
Replayable exact symbolic checks may verify geometry, signs, selection rules,
and Schur identities. **No production numerics are authorized.** If an
irreducible small divisor or small splitting survives the analytic ladder,
the small-ratio-numerics skill is read before a separate numerical design is
frozen.

Each of Routes A--C and the nonlinear bridge receives exactly one verdict:
`established`, `refuted_with_mechanism`, or
`blocked_by_named_construction`. A failure of Route A activates Route B; a
finite exact resonance activates Route C. No route verdict closes P253, LP2,
or P4. Success here is one actual same-carrier classical rotating branch and
its classical KKS action, not a quantum or particle identification.
