# P253/0044: curved Cao-ring translation Feshbach construction

Status: **README-only preregistration awaiting central registration and schema
activation.** Until `0044/activation-schema.exit` exists and contains exactly
`0`, this README is the sole 0044 artifact. Root owns `proposal.yaml`, central
integration, commits, review, and activation. Particle-foundations owns only
post-activation append-only bodies in `attempts/0044`.

## Failure-derived positive objective

Attempt 0040 left one concrete spectral input for its Cao rotating-wave route.
For one fixed azimuthal number `l>=2` and sufficiently small core ratio
`delta=a_epsilon/R`, construct the actual two-polarization translation-bundle
operator of the full Leray-linearized Cao ring. The attempt must derive:

1. the exact compact-support dynamically accessible domain;
2. the finite-core response of the Lane--Emden column, including its simple
   translation kernel and nonzero parameter pairing;
3. the full toroidal-curvature coupling of cross-sectional Fourier fibers,
   especially `m=1` into `m=0,2`, before choosing a spectral contour;
4. the resulting finite or enlarged Feshbach matrix and an earned complement
   inverse/Riesz projection; and
5. an actual simple imaginary eigenpair, or a route-scoped obstruction after
   all registered representations have been executed.

The target is the spectral/Fredholm input for 0040, not a renamed missing
theorem. Success licenses the next solid-torus CR/Nash--Moser calculation; it
does not itself construct the nonlinear branch or quantum structure.

## Fixed carrier, domain, and scale

Use the fixed polynomial-profile no-swirl Cao ring of 0040 with circulation
`Gamma>0`, major radius `R`, core radius `a_epsilon`, and

    delta=a_epsilon/R,        k_delta=l delta,
    omega_core=Gamma/(4*pi*a_epsilon^2).                    (1)

In local Frenet coordinates `y=(r-R,z)/a_epsilon`, the limiting compact column
has vorticity `W_p(s)=U(s)_+^p`, where `-Delta U=U^p` on the unit disk. The
exact operator is the `e^(i l theta)` sector of the whole-space Leray/Hodge
generator in the translating frame, with all pressure and cylindrical/toroidal
connection terms retained.

The physical tangent space is the energy-plus-graph closure of

    delta Omega=curl(xi cross Omega_0),
    div xi=0,  xi smooth and compactly supported.             (2)

Because `Omega_0` is compactly supported, (2) is distributionally supported in
the same closed core. Prove that the chosen graph closure embeds continuously
in distributions; closedness of fixed distributional support then excludes
passive exterior vorticity. The exterior velocity is retained exactly through
whole-space Biot--Savart/Hodge matching. The full-enstrophy
`Omega(r)->0` endpoint of Gallay--Smets is therefore a representation
comparator, not an obstruction imported into this physical orbit.

## Exact local translation response

For the planar `m=1`, `k=0` core, freeze the source-defined scalar pencil

    T(b) psi=[-d_rr-r^-1 d_r+r^-2
              +W'(r)/(r(Omega(r)-b))] psi=0,                (3)

regular at the origin and matched at the core radius `a` to the exterior Robin
condition

    psi'(a)=-psi(a)/a.                                      (4)

The translation solution is

    psi_0(r)=r Omega(r)>0.                                  (5)

Execute the ground-state transform of `T(0)` in `L2((0,a),r dr)` to prove a
simple zero and a positive complement gap after removing translation. Evaluate
the parameter crossing exactly:

    <psi_0,T'(0)psi_0>_(r dr)
       =integral_0^a r^2 W'(r)dr
       =-2 integral_0^a r W(r)dr=-2F(a) != 0.              (6)

This is the actual finite-core translation response; Cao's elliptic uniqueness
operator is not substituted for (3). For `k!=0`, let

    D(x)=x K_1'(x)/K_1(x),       x=|k|a.                    (7)

The exterior potential has Dirichlet-to-Neumann multiplier `D/a`, but the
scalar unknown in (3) obeys the distinct matched boundary relation

    psi'(a)/psi(a)=(1+x^2)/(a D(x)),
    -1/D(x)=1+x^2[log(x/2)+gamma_E]
               +O(x^4 log(x)^2).                           (8)

Use the exact `T(k,b)` core pencil and the gap following (6) to construct its
scalar matching determinant. Its frozen projected expansion is

    d(k,b)=-2F(a)b+F(a)^2 k^2 log(|k|a)+O(k^2)
      +O((|b|+k^2|log k|)^2),                              (9)

uniformly in the declared small neighborhood. The nonzero derivative in (6)
then gives the actual compact-column Kelvin branch

    b(k)=-(F(a)/2)k^2 log(1/(|k|a))+O(k^2).                (10)

This earns the finite-core inverse and leading longitudinal response. It does
not decide the distinct curved-ring cross-`m` complement.

## Curved block and asymptotic target

Let `T_1,T_2` be the two real translations. Construct their duals from the
physical KKS/impulse pairing and let `P_1` be the resulting biorthogonal
projection, `Q_1=1-P_1`. It is not an arbitrary Euclidean projection. Expand
the exact toroidal operator on (2) as

    A_delta=A_0+delta C_1+delta^2 C_2+R_delta.              (11)

The first curvature coefficients contain `cos(alpha),sin(alpha)` and obey
`C_1:m -> m+/-1`. Thus `P_1 C_1 P_1=0`, while the second-order Schur term runs
through `m=0,2`:

    M_l^(2)(z)=P_1 C_2 P_1
      -P_1 C_1 Q_(0,2)[Q_(0,2)(A_0-z)Q_(0,2)]^-1
         Q_(0,2) C_1 P_1.                                  (12)

Derive `C_1,C_2` from the toroidal metric, base flow, connection, Hodge
projector, and pressure. Audit the accessible `m=0` spectrum before using (8):
internal Kelvin frequencies may accumulate at zero. If the impulse/Casimir
slice cancels every zero or near-zero coupling, prove it and evaluate (8).
Otherwise enlarge `P_1` by the coupled resonant subspace and derive its normal
form. No contour radius or bending gap is asserted before this calculation.

After that gap audit, the target two-polarization matrix is

    M_l(delta)=omega_core*delta^2
      {log(1/delta) M_l^log+M_l^core(p)+M_l^curv(p)}
      +R_l(delta),                                          (13)

The straight-column restriction of `M_l^log` must reproduce (10) at
`k=l delta`; its full radial/binormal entries are to be derived from (11)--(12).
The cutoff-filament comparator

    [[0,l^2-1],[-l^2,0]]                                   (14)

is an exposing limit target, not an imported Cao coefficient. Toroidal
curvature contributes at the same logarithmic order and must produce or
correct the `l^2-1` entry from the physical operator.

and an earned remainder

    ||R_l(delta)||=o(omega_core*delta^2*log(1/delta))        (15)

on the eventual contour. The finite matrices must be evaluated as convergent
Lane--Emden integrals or explicit radial cell problems. They may not be copied
from a hollow filament.

The exact Feshbach object is

    K_l(delta,z)=P(A_delta-z)P
      -P A_delta Q[Q(A_delta-z)Q]^-1Q A_delta P,             (16)

where `P=P_1` only if (12) proves no enlargement is necessary. Both terms,
including field-dependent Leray pressure, must be evaluated through (13).

## Route A: regular physical-core reduction

Prove the support/closure statement for (2), derive the ground-state gap and
(6), and execute the implicit continuation (7)--(10) from the exact `T(k,b)`
pencil. Retain the exterior velocity and distinguish its potential DtN from
the scalar boundary relation (8). Failure scope is an extra `m=1` core zero
mode or a boundary-domain defect, not the carrier.

## Route B: curvature-coupled enlarged Feshbach reduction

Derive the selection rule and all nonzero `m=1 -> m=0,2` blocks in (11)--(12).
Prove their impulse/Casimir cancellations or enlarge the resonant projection.
Only then choose a scaled contour and prove the complement/Riesz bound. The
ambient-enstrophy critical layer is inapplicable here because it adds passive
exterior-vorticity states absent from (2).

## Route C: matched inner/exterior Evans determinant

Match the regular `m=1` core solution to the decaying exterior potential
`K_1(|k|s)e^(i alpha+i k z)`. Use the corrected variable map (7)--(8) and the
core determinant (9), then include the `m=0,2` curvature blocks before applying
Rouche's theorem. This route must reconstruct a physical graph-domain
eigenfunction and control the collar and exterior velocity.

## Route D: Hamiltonian/Krein reduction

Derive the same matrix by constrained second variation and KKS reduction on
the physical leaf. A positive reduced quadratic form plus an actual complement
inverse may establish a signed imaginary pair. Agreement of its finite-core
and curvature entries with Routes A--C is an exposing cross-check; a formal
finite matrix without domain reconstruction does not establish the branch.

## Sources, verification, and verdict boundary

- Cao--Lai--Qin--Zhan--Zou, arXiv:2206.10165v2: the exact compact-core ring and
  Lane--Emden blow-up; its axisymmetric stability and elliptic uniqueness
  operator are not the required 3D spectral theorem.
- Gallay--Smets, arXiv:1805.05064v3: the exact column equations and `m=1`
  translation/Kelvin analysis. Its full ambient enstrophy domain is kept
  distinct from (2).
- Fischer--Schopohl, arXiv:cond-mat/0008215v3: the cutoff-filament logarithmic
  matrix used only as the limit comparator (14), never as a finite-core Cao
  coefficient.

Downloaded bodies remain in `/tmp/primary-source-cache`, with URL, version,
hash, and exact equation/theorem locations recorded in the attempt. Existing
abstract APIs may supply conventions but cannot substitute for (12).

The strongest oracle is exact source-equation reduction plus a replayable
symbolic check of support preservation, ground-state factorization, the
pairing (6), curvature selection rules/signs, finite matching, and any
determinant/Riesz implication. No production numerics are authorized. If one
analytic scalar remainder survives, freeze it and load the applicable
numerics skill before designing a computation.

Each route receives exactly one of `established`, `refuted_with_mechanism`, or
`blocked_by_named_construction`, with `evidence_scope` separate. Success
requires an actual physical graph-domain eigenpair/Riesz projection, not only
(13). This attempt does not construct the nonlinear rotating branch, prove its
stability, select a quantum, infer an FR sign, or identify an electron or
neutrino. P253 and LP2/P4 remain active.
