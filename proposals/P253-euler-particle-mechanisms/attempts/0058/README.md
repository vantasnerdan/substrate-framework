# P253/0058: nonlinear same-leaf rotating Cao branch

Status: **README-only preregistration awaiting central registration and schema
activation.** Until `0058/activation-schema.exit` exists and contains exactly
`0`, this README is the sole 0058 artifact. Root owns `proposal.yaml`, the
central attempt registry, activation receipts, commits, review, and promotion.
No derivation body, verifier, source download, numerical run, API, or claim is
authorized before activation.

## Author-input boundary and positive target

This attempt uses the following only as **author-stage inputs pending their
independent reviews**:

1. P253/0052: the corrected finite-core bending/Kelvin crossing, exact KKS and
   physical-frequency normalization, enlarged finite spectral window, and
   full-sector Feshbach cluster conditional there on `HJ2`; and
2. P253/0054: the `p>=6` common dynamically accessible (DA) graph domain,
   moving-boundary whole-space Green/Leray jet, modewise order-minus-one and
   order-minus-two bounds, and contour-relative graph-Riesz transfer.

Registration of 0058 does not adjudicate either input. Any correction from
their reviews is inherited before a 0058 theorem may be promoted.

The positive target is one exact nonaxisymmetric, compact-vorticity,
finite-energy Euler branch on the same physical coadjoint leaf as an actual
centered `p>=6` Cao vortex ring. In the frame translating with the ring, the
branch is a rotating relative equilibrium and therefore a relative-periodic
Euler solution in the laboratory frame. It must have:

- an exact SO(2) period and real field;
- nonzero physical KKS and angular-momentum pullback;
- convergence of the nonlinear range construction through and beyond the
  first internal critical harmonic;
- whole-space pressure and exterior velocity reconstruction; and
- a separately proved restoring or orbital-stability statement on a genuine
  same-leaf perturbation neighborhood.

No particle, spin-statistics, quantum, electron, or neutrino inference is part
of this attempt.

## Fixed carrier, cluster, and symmetries

Fix one sufficiently thin Cao member with integer `p>=6`, circulation
`Gamma>0`, density `rho_0>0`, major/core radii `R_delta,a_delta`, and

    delta=a_delta/R_delta,       L=log(1/delta).           (1)

Work in the centered frame translating at the exact carrier speed `c_0`.
Let `A_delta` be the exact graph-domain linearized Euler generator on the
physical DA quotient and let `Pi_delta` be the author-stage 0052/0054 Riesz
projection. Select one simple positive-frequency azimuthal pair

    A_delta q=i nu_delta q,
    q=q_c+i q_s,       ell>=2,                            (2)

from the finite bending/Kelvin cluster, with positive constrained Hessian and
the exact convention

    Omega_KKS(q,conjugate(q))=-i/nu_delta                (3)

for an energy-unit complex mode. The real basis and its KKS duals are fixed
before nonlinear reduction. With

    [R,q_c]=ell q_s,
    X_J=-[R,omega],                                      (4)

the physical rotation sign is

    X_J(q_c)=-ell q_s,
    Omega_pat,0=nu_delta/ell,
    J(a)=J(0)+(ell sigma/2)a^2+O(a^3),                   (5)
    sigma=Omega_KKS(q_c,q_s)!=0.

Equation (5) preserves the supervisor-corrected 0040 sign. Neither `nu_delta`
nor `sigma` may be replaced by a filament normalization.

Real fields use the conjugate `+/-ell` pair. Let

    S:(r,theta,z)->(r,-theta,-z),
    mathscr R u=-S_*u.                                   (6)

The centered translating Cao velocity is fixed by `mathscr R`, and Euler's
quadratic vector field obeys

    X(mathscr R u)=-mathscr R X(u).                      (7)

The SO(2) action commutes with Euler, while the reflection in (6) reverses the
rotation parameter. The cosine slice is `mathscr R`-fixed and the sine mode is
its phase partner. The derivation must verify (6)--(7) componentwise in the
actual translating frame; generic Hamiltonian reversibility is not imported.

## Exact same-leaf nonlinear map

Let `g=Exp(xi)` be a volume-preserving `H^s` diffeomorphism equal to the
identity outside the physical collar, modulo the stabilizer and with the
impulse, center, and phase rows fixed. Define

    omega_xi=g_* omega_0,
    v_xi=BiotSavart_R3(omega_xi).                         (8)

This is the nonlinear DA chart. An arbitrary compact vorticity perturbation or
an independently changed streamline profile is not substituted for (8).

The rotating-wave equation has equivalent vorticity and velocity forms

    F_omega(xi,c,Omega)
      =curl[(v_xi-c e_z-Omega R) cross omega_xi]=0,       (9)

    F_v(xi,c,Omega)
      =P_L[(v_xi dot grad)v_xi]
         -c partial_z v_xi-Omega[R,v_xi]=0.              (10)

The exact whole-space Leray operator in (10) supplies pressure. Prove the
equivalence of (9) and (10), including the harmonic/circulation row and decay
at infinity. The support in (8) may move but remains a smooth solid torus; the
velocity is never required to be compact.

Use the phase/amplitude decomposition

    xi=a q_c+w,       Pi_delta w=0,
    <q_s^*,w>=0,                                        (11)

and allow exact corrections `c(a),Omega(a)`. The range and bifurcation
equations are

    Q_delta F_omega(aq_c+w,c,Omega)=0,                   (12)
    Pi_delta F_omega(aq_c+w,c,Omega)=0.                  (13)

Their derivative on the complement is the physical rotating operator

    L_rot=Q_delta(A_delta-Omega_pat,0 X_J)Q_delta.        (14)

The translation and phase duals determine `c'(a)` and `Omega'(a)`; no fitted
frequency is allowed. Equation (13) must compute the first nonzero nonlinear
frequency/action coefficient and distinguish the full SO(2) group period from
the shorter period allowed by the `Z_ell` isotropy:

    T_SO2(a)=2 pi/|Omega(a)|,
    T_min(a)=2 pi/[ell |Omega(a)|].                       (15)

Both close the appropriate real velocity and vorticity fields. In the
laboratory frame the corresponding axial displacements are
`c(a)T_SO2(a)` and `c(a)T_min(a)` and must be retained in the relative-periodic
statement.

## Critical layers and exact range divisibility

In the 0054 action--angle variables `(I,beta,theta)`, the `(m,j)` harmonic of
the complement range equation has principal divisor

    D_(m,j)(I;a)=m omega_delta(I)-j ell Omega(a).          (16)

For the Cao core, `omega_delta(I)` covers a positive interval. The first
unavoidable generated critical harmonic satisfies

    j_crit=Theta((delta^2 L)^-1).                         (17)

Smallness of the `j`th coefficient does not invert a zero of (16). At each
simple critical action `I_*`, the exact numerator `N_(m,j)` must satisfy

    N_(m,j)(I_*)=0,
    N_(m,j)(I)=D_(m,j)(I) G_(m,j)(I)                     (18)

in the physical graph space. The divided-difference inverse is then

    G_(m,j)(I)
      =[N_(m,j)(I)-N_(m,j)(I_*)]/D_(m,j)(I),             (19)

with its Hardy/tame estimate derived from the actual lower bound on
`|partial_I D_(m,j)(I_*)|`. The source profile must be used to prove
monotonicity or to classify all degenerate roots. A root of multiplicity `r`
requires the first `r` numerator traces to vanish; it cannot be deleted by a
generic nonresonance condition.

There are two allowed mechanisms for (18):

1. **structural transparency:** the exact vorticity transport, incompressible
   constraint, Bernoulli equation, and DA commutators force the numerator to
   contain `D_(m,j)`; or
2. **same-leaf range variables:** a volume-preserving streamline displacement
   supplies the missing trace while preserving every Casimir and circulation.
   Its induced profile change must be written as a Lie derivative of the
   original vorticity. An arbitrary new function `f(Psi)` would change the
   coadjoint leaf and is excluded.

Freeze the trace operator

    T_(m,j)N=(partial_I^r N(I_*))_(0<=r<mult(I_*))        (20)

and derive its adjoint pairing with the KKS range. Every critical equation is
recorded as `T_(m,j)N=0`; calling it exponentially small is not a solution.

## Route A: nonlinear KKS Lyapunov--Schmidt normal form

Route A solves (12)--(13) directly on the same-leaf chart.

1. Prove that composition, Piola push-forward, whole-space Biot--Savart, and
   the nonlinear Leray map are tame on the fixed 0054 graph scale.
2. For `|j|<j_crit`, invert (14) with the exact finite-window Riesz projector
   and track the deterioration of the divisor down to `Theta(delta^2 L)`.
3. At and above `j_crit`, prove (18) structurally or solve the trace equations
   with the permitted DA range variables.
4. Construct a KKS-symplectic near-identity normal form. Its homological
   inverse must retain the full finite-core frequencies and all continuous
   streamline divisors.
5. Prove tame estimates for the coefficient sequence through the scale
   `j_crit`, not only at each fixed order. A permitted target is

       ||w_j||_(s-mu)
          <=C_delta A_delta^j <j>^tau,                   (21)

   with the derivative loss `mu`, exponent `tau`, and amplitude radius derived
   from (19), followed by an explicit convergent majorant.
6. Sum the series or run a Nash--Moser iteration with the 0054
   projector-compatible smoother. Record the actual radius `a_0(delta)>0`;
   no uniform thin-limit radius is presumed.

The route succeeds only with an exact solution of (9)--(10), not the finite
formal residual already available in 0048.

## Route B: reversible invariant-manifold/Lyapunov-center construction

Route B seeks a real two-dimensional invariant symplectic manifold tangent to
`span{q_c,q_s}`. Parameterize it by

    K(z,conjugate(z))=z q+conjugate(z)conjugate(q)
       +sum_(r+s>=2) K_(r,s)z^r conjugate(z)^s,           (22)

and solve the exact invariance equation

    X_Euler(K)=D K X_red,
    X_red=i partial_(conjugate z) H_red.                  (23)

Require SO(2) equivariance, the reality condition
`K(conjugate z,z)=conjugate(K(z,conjugate z))`, and the reverser (6)--(7).
Pull back the physical KKS form and prove it is nondegenerate with the sign in
(3)--(5).

The coefficient equations in (23) contain the same divisors (16). Route B
must prove the trace identities (18)--(20); finite-dimensional invariant-
manifold terminology does not bypass a critical layer. Success gives a
convergent invariant manifold and a Lyapunov-center family on it. If normal
hyperbolicity fails because the complement is center spectrum, use a
reversible parameterization/KAM scheme with explicit tame bounds rather than a
spectral-gap theorem whose hypotheses are absent.

## Route C: failure-derived steady range reformulation

If both coefficient schemes encounter the same nonzero critical trace, keep
the carrier and rewrite (9) as a coupled streamline-profile/free-boundary
equation. The extra unknown is licensed only when it is generated by a
volume-preserving displacement as in (8). Derive the map from that displacement
to the trace data (20) and test whether it is onto with a tame right inverse.
This is a representation change of the range problem, not permission to move
to another Casimir leaf.

A failure of onto-ness earns a verdict only for this Cao same-leaf branch and
must identify the nonzero adjoint obstruction. It does not refute other
localized Euler carriers or P253's particle objective.

## Pressure, exterior field, and convergence oracle

Every nonlinear iterate is a compact-vorticity DA field from (8). Its velocity
and pressure are reconstructed globally:

    v=B omega,
    p=-Delta_R3^-1 partial_i partial_j(v_i v_j),           (24)

with the translating/rotating affine terms included in the Bernoulli constant.
Prove finite kinetic energy, decay, interface trace regularity, and convergence
of (24) in the same norm used for (9)--(10). A core-only pressure or exterior
Dirichlet condition cannot establish the branch.

The strongest practical oracle must expose at least:

- the physical rotation/KKS sign in (4)--(5);
- a deliberately nondivisible numerator at a simple critical layer;
- loss of the harmonic circulation row;
- a fake inversion based only on `|N_(m,j)|` being small;
- failure of reality or period closure;
- replacement of the global pressure by a core solve; and
- convergence of a majorant or Nash--Moser residual to zero, rather than a
  fixed-order tally.

No production numerics are preregistered. If one irreducible coefficient or
small divisor remains after exact analysis, a separate append-only expansion
must first freeze its scale, error enclosure, and bridge to the continuum
theorem.

## Separate restoring and perturbation-neighborhood obligation

Existence of a rotating branch does not establish persistence. On the exact
branch define the augmented Hamiltonian

    H_aug=H-c(a)P_z-Omega(a)J.                            (25)

Compute its constrained second variation on the full same-leaf tangent,
remove translations and rotation by the physical modulation conditions, and
separate the branch tangent. The positive route proves either:

1. coercivity

       delta^2 H_aug[h,h]>=kappa_delta ||h||_E^2          (26)

   on the modulated complement, followed by orbital stability; or
2. a reduced restoring coefficient from the exact KKS normal form together
   with a controlled invariant-manifold neighborhood and bounded transverse
   leakage.

The perturbation set must be an actual open ball of volume-preserving
same-leaf `H^s` diffeomorphisms (with its induced physical energy distance),
not only the two-mode plane or a future-history-prepared set. State whether
`kappa_delta` and the neighborhood radius are fixed-carrier or uniform as
`delta->0`. If the Hessian has mixed signature, continue with Hamiltonian
normal-form or invariant-manifold control; mixed signature alone is not called
instability.

## Source and dependency inventory

Full papers remain in `/tmp/primary-source-cache`.

- Cao--Lai--Qin--Zhan--Zou, arXiv:2206.10165v2: exact `p>=6` carrier,
  free-boundary regularity, circulation/translation gauges, and constrained
  stability at the source's axisymmetric scope. It does not supply a
  nonaxisymmetric rotating branch.
- Gallay--Smets, arXiv:1805.05064v3: straight-column mode equations and Hodge
  normalization only; its ambient profile/spectrum is not imported.
- P253/0040 and 0048: physical rotation sign, nonlinear map, formal tame
  estimate, and the actual `j_crit` obstruction at their author-record scope.
- P253/0052 and 0054: the exact finite cluster and graph transfer are consumed
  only as pending-review author inputs. A review correction updates the 0058
  dependency before execution or promotion.

No generic Crandall--Rabinowitz, Lyapunov-center, KAM, or invariant-manifold
theorem is cited as if it supplied the Euler range, critical-layer traces,
global pressure, or same-leaf chart.

## Obligation and verdict ledger

| Node | Positive intent | Pass licenses | Does not license | Maximum verdict | Failure scope | Unlocks |
| --- | --- | --- | --- | --- | --- | --- |
| N1 | Exact real KKS cluster, same-leaf nonlinear chart, and rotating-wave map | Nonlinear kernel/range split | Range solvability | Exact analytic | 0052/0054-dependent cluster | N2 |
| N2 | Critical trace divisibility or licensed same-leaf trace right inverse | Tame inverse through all harmonics | Convergence without estimates | Exact analytic | Selected transparency/range mechanism | N3 |
| N3 | Convergent LS, normal-form, or invariant-manifold solution with global pressure | Exact rotating/relative-periodic Euler branch | Stability or particles | Exact analytic | Selected nonlinear construction | N4, P253 LP2 branch existence |
| N4 | Coercive or controlled restoring dynamics on a real same-leaf neighborhood | Classical localized persistent carrier evidence | Quantum, electron, or neutrino identification | Exact analytic, or explicitly norm-bounded dynamical | This branch/neighborhood | P253 LP2 review |

Each executed route receives exactly one route-scoped verdict: established,
refuted with the mechanism, or blocked with its missing construction. Failure
of one nonlinear representation activates the next route and does not inherit
upward as a no-go.

## Frozen success contract

0058 succeeds only when all of the following hold for the same actual Cao
carrier:

1. the reviewed-or-explicitly-pending exact KKS cluster and signs are retained;
2. the real/reversible and SO(2) symmetry reductions are proved in the physical
   translating frame;
3. the nonlinear range equation crosses every critical harmonic by (18)--(20)
   or a proved same-leaf right inverse;
4. the tame construction converges to an exact solution of (9)--(10);
5. the SO(2) period, angular momentum, pressure, exterior velocity, support,
   finite energy, and Casimir preservation are reconstructed; and
6. the separate restoring/stability obligation holds on a genuine real
   same-leaf perturbation neighborhood.

A fixed-order residual, formal normal form, finite Galerkin branch, sampled
spectrum, prequantization, or analogy does not meet this contract. The parent
P253 campaign remains active and no particle or quantum conclusion is licensed.

## Post-execution Sol-High scope correction

The pendulum equations in the executed body are a conditional generic normal
form. An Euler inner-layer claim first requires the full-Hodge
limiting-absorption adjoint polarization and symplectic measure, a
`sqrt(epsilon)` support/derivative/Bogovskii ledger, and an allowed
axisymmetric centralizer field with compatible critical-set topology. The
exact carrier-map, positive-core centralizer, Hodge/DA symbols,
nondivisibility counterexample, reverser, and KKS sign remain unchanged.
