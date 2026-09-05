# A defined volume-preserving material mean, and what its action retains

## 1. Primary import and the fluctuation variation

[Holm, nlin/0103043](https://arxiv.org/html/nlin/0103043v1), sections 3–4,
provides the factorization and averaged Euler–Poincare structure. Its mean
momentum is a pulled-back one-form, not generally the mean velocity; section
4.2 identifies the pseudomomentum. The theorem does not close fluctuation
statistics. Its displayed mean equations omit delta-Xi terms for the
prescribed-fluctuation description; section 4.2 explicitly separates that
description from self-consistent fluctuation dynamics. The construction here
therefore retains and varies Xi, including the physical angle. It does not
use the prescribed-Xi equations as an angle equation or assume the paper's
mean is volume preserving.

## 2. An actual idempotent averaging map, selected without a target modulus

Use common material labels a and probability/phase label alpha, and let each
g^alpha be an actual volume-preserving map. Assume its arithmetic material
mean F(a)=E_alpha[g^alpha(a)] is a smooth orientation-preserving diffeomorphism
in the considered mean-coordinate chart. This assumption can be checked on
the actual maps; a folded F is not repaired by calling it a mean diffeomorphism.
Work on a common periodic domain, or on R³ with the stated decay and solvable
Poisson correction. This defines the averaging operation, not a periodic
replacement for an EPS microscopic source.

Let rho_F(y)=1/det DF(F^-1(y)) be the density F pushes forward from unit
material volume. Define a physical-space volume correction by

    rho_tau=(1-tau)rho_F+tau,
    Delta psi=rho_F-1,
    v_tau=grad psi/rho_tau,
    dC_tau/dtau=v_tau composed with C_tau,     C_0=id.

The auxiliary tau is NOT physical time. Positivity of rho_tau and the
regularity/decay assumptions give the ordinary smooth flow. Its continuity
identity is exact:

    partial_tau rho_tau+div(rho_tau v_tau)=0,
    (C_1)_*(rho_F dy)=dy.

Set

    g_mean=C_1 composed with F,
    Xi^alpha=g^alpha composed with g_mean^-1.

Both g_mean and every Xi^alpha are volume preserving, and the factorization
g^alpha=Xi^alpha composed with g_mean is exact. This is a specified averaging
procedure: solve one density Poisson problem and one auxiliary transport
flow. No desired inertia, spin, dispersion or Cosserat coefficient occurs
in its definition.

It is a projection: when every realization equals one volume-preserving h,
F=h, rho_F=1, psi=0 and g_mean=h. Averaging that deterministic output again
does not change it. It is also right-equivariant under common volume-
preserving relabeling eta: F changes to F composed with eta, but rho_F and
its physical-space correction do not change; hence g_mean changes to
g_mean composed with eta. The full particle-relabeling symmetry is not
replaced by a label-dependent Poisson metric. Euclidean covariance holds
with the covariantly transformed domain/Poisson data.

This corrected mean need not equal the uncorrected arithmetic mean and need
not satisfy E[Xi]=id as a pointwise vector equality. Those differences and
the correction's variation are retained. Its mean density is exactly rho
and div u_mean=0, unlike assuming those properties for every GLM definition.

## 3. Exact velocity, tag and cotangent formulas

Write u_mean=g_mean,t composed with g_mean^-1. The physical velocity pulled
back to mean position is

    u^alpha composed with Xi^alpha
       =Xi_t^alpha+D Xi^alpha u_mean.

For a material tag with initial value chi0^alpha,
chi_mean^alpha=chi0^alpha composed with g_mean^-1 and
chi_actual^alpha=chi_mean^alpha composed with (Xi^alpha)^-1. The factorization
therefore transports actual tube tags and their collars; it does not freeze
an Eulerian mask while allowing a crossing velocity. The material commutator
identity in 0082 is its infinitesimal form.

At unit reference volume the SAME Euler kinetic action is exactly

    T=rho/2 integral E |Xi_t+D Xi u_mean|² dx.

Pressure, the volume constraints on Xi and g_mean, the defining mean gauge,
and the selected phase/coherence constraints remain in the action. If the
mean is treated as an independent variable, its defining constraint and
multiplier are varied too; prescribing an arbitrary g_mean while ignoring
that constraint is not this averaging procedure.

For independent fluctuation coordinates z, write Xi_t=Xi_,a zdot_a, including
all retained microscopic reaction and physical-angle coordinates. The exact
material kinetic blocks are

    A=rho E[(D Xi)^T D Xi],
    B_a=rho E[(D Xi)^T Xi_,a],
    C_ab=rho E[Xi_,a.Xi_,b].

Thus the mean momentum and fluctuation momenta are

    m=A u_mean+B zdot,
    Pi=B* u_mean+C zdot.

The entire positive Gram block is retained. A legitimate reaction reduction
takes its full constrained Schur complement, not a copied isolated-cell
inertia. In particular m is not set equal to rho u_mean: the difference is
the full pseudomomentum, including the mean-volume correction. A gradient
part has zero circulation but still belongs in the pressure/constraint map.

Varying the actual physical angle changes Xi_,a, D Xi, the pressure constraint
and the mean gauge. Its Euler–Lagrange equation and momentum Pi are kept.
Dropping those variations would remove the very angle dynamics being sought.
The source's mean-only equations are not used to fill in that missing equation.

Each realization retains its material circulation one-form
`g_mean^*[(Xi^alpha)^*(u^alpha flat)]`. Fixing all these initial one-forms
is stronger than fixing only their ensemble average. The averaged Kelvin
momentum is the corresponding pulled-back mean one-form m/rho, with the
full fluctuating/reaction equations and their momentum constraints present.
This construction does not equate a geometric fluctuation metric to an
already reduced Euler coadjoint inertia without performing that reduction.

## 4. Direct calculation of the physical mean-coordinate shift

An exact factorization now exists, so its observation map can be computed
without choosing the answer. The following small-displacement calculation
also checks the sign and factor independently of the parcel argument.

Let raw volume-preserving fluctuation maps be flows exp(epsilon xi^alpha),
with div xi^alpha=0 and E xi^alpha=0. This flow family is specified explicitly;
no assertion that every diffeomorphism has a Lie logarithm is used. Put

    A0=E[(xi.grad)xi]=div E[xi tensor xi],
    P=I-grad Delta^-1 div.

Expanding the actual mean operation of section 2 gives

    F=id+epsilon² A0/2+O(epsilon³),
    C_F=id-epsilon² (I-P)A0/2+O(epsilon³),
    g_mean=id+epsilon² P A0/2+O(epsilon³),
    Xi^alpha=id+epsilon xi^alpha
        +epsilon²[(xi^alpha.grad)xi^alpha-P A0]/2+O(epsilon³).

The normalizing term is calculated by the density Poisson problem, not a
postulated spin-dependent shift. Direct differentiation and inverse-map
composition around a zero mean velocity give

    E u_actual-u_mean
       =epsilon² curl E[xi cross xi_t]/2-epsilon² P A0,t/2
          +O(epsilon³).

Define the actual fluctuation spin and covariance at this order by
S=rho epsilon² E[xi cross xi_t] and
C=epsilon² E[xi tensor xi]. The result is

    u_E-u_mean=curl S/(2rho)-P div(C_t)/2+higher terms.

It is precisely a spin/shape current, but derived from a defined material
mean rather than assigned to a normal-form variable. The sign gives
`u_E=u_mean+curl S/(2rho)` when the relevant shape-rate term vanishes.

## 5. The same shift in the original slow-spatial scope

The observation identity is not restricted to a formal Lie-flow example.
For bounded fluctuation displacements eta=Xi-id, start from the exact
pushforward distribution

    rho u_E(y)=rho E integral [Xi_t+D Xi u_mean](x)
                                      delta(y-Xi(x)) dx.

Expanding this distribution in displacement moments, with its integral
Taylor remainder, gives the full first and second spatial multipoles. The
first fluctuating momentum moment decomposes into a symmetric covariance-
rate term and an antisymmetric fluctuation-spin term. The mean-volume
correction supplies its longitudinal component; the transverse leading
map is the P-projected identity in section 4. Translation/covariance and
higher internal moments are retained when computing second spatial jets.

Thus the parent can use the original slow-affine/spatial-gradient ordering,
with an actual bounded-displacement moment hypothesis, instead of claiming
unrestricted finite-k closure. Neither a neglected third amplitude order
nor a spatial Taylor remainder is declared zero in an exact full-flow claim.

If the FULL varied material fluctuation action independently gives
`S=j Phidot` at its leading local order and supplies the covariance response,
then the actual averaging map yields

    U_E=U_mean+j curl Phi/(2rho)+shape/second-gradient terms.

This realizes the leading physical map behind 0066's coupled normal
displacement. The observable is the trajectory of the explicitly averaged
volume-preserving material map. It is not the whole-fluid Eulerian point
mean, so the cancellation in 0072 does not identify the two coordinates.

The condition S=j Phidot is a physical material spin calculation, not a
license to insert an old coadjoint coefficient under the same symbol j.
For a truly rotating isotropic displacement covariance C=c I, direct
material geometry gives S=2rho c Phidot, while the SAME kinetic norm gives
j=2rho c. That elementary identity explains how the needed map can arise;
an actual divergence-free core/collar fluctuation must retain its returns,
transport and Kelvin reaction before adopting those simple covariance
values. This route supplies the averaging license, not that final microscopic
identification for every previously constructed sector.

## 6. What survives and the next concrete construction

Established: an explicit idempotent, mass-preserving, relabeling-equivariant
material mean; exact factorization and varied-action/pseudomomentum formulas;
and its calculated spin/covariance observation map. A genuine Lagrangian mean
can realize the desired leading canonical displacement without changing the
point-mean physics or omitting its pseudomomentum.

The remaining coupling identification is concrete: evaluate A,B,C, the actual
fluctuation spin and covariance response for a chosen material core/collar
family, impose its individual Kelvin momenta, and compare the resulting
full Schur action and physical angle. The existing 0059 relative-action
coefficient cannot be pasted into S by name. The factorization handles tag
transport automatically once those actual material maps are supplied; the
mean-only LAEP theorem is not substituted for their dynamics. This is a
positive alternative averaging construction, not parent completion or an
exhaustion conclusion.

### Executable material-spin continuation with the new positive material mode

The new material-Jacobi positive stiffness of 0089/0090 is the appropriate
input class for the full kinetic Gram above. It still does not by itself
prove a fluctuation spin law: for xi=q Y, xi cross xi_t is identically zero.
A nonzero registered background/phase fluctuation is needed.

The coordinated 0091 construction supplies a useful exact spin/metric
mechanism. On divergence-free displacement fields let A_e=P(e cross .).
For a divergence-free eta0, Y=A_e eta0 satisfies

    integral (e cross eta0).Y = integral |Y|².

This is orthogonal projection, not an assumption that vector rotation
preserves divergence. A_e is skew-adjoint on that space, so
eta(q)=exp(q A_e)eta0 preserves its displacement norm and its material
spin coefficient equals the computed tangent kinetic norm. This is not
automatically a physical SO(3) action or a stationary fluctuation ensemble;
the actual core-angle and Kelvin registration are retained by 0091.

Spatially translated registered pairs (u0,eta0) provide an explicit way to
preserve the same tracked-core response across phases while averaging a
zero-mean eta0. A naive plus/minus profile pair need not do so. Initial
material maps and their individual circulation one-forms are specified
together; their subsequent covariance evolution is varied, rather than
assumed stationary because the Eulerian background field is stationary.
This is the next concrete material-family computation to feed the actual
mean and pseudomomentum formulas established here.
