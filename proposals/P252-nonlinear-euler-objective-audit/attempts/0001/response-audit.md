# Fresh response/action audit of P251 and the nonlinear objective

## Boundary and verdict

This is an equation-level audit of the P251 response/action/current chain at
base `b6fc902`, release `v0.183.0`.  I formed the conclusions below from the
proof sources before reading the favorable reviews in 0264, 0267 and 0269.
Reviewer labels, review counts and green tallies are not evidence used here.

The strongest supported positive result is genuinely useful but narrower than
a nonlinear Euler-to-Cosserat constitutive derivation:

> Subject to the fixed compact-ring existence hypotheses, P251 constructs, for
> each sufficiently small **nonzero** macroscopic wave number, compact time
> window, finite derivative inventory and requested accuracy, a finite linear
> Euler/Lin initial-data map on one fixed stationary ensemble.  Its selected
> material observations approximate a prescribed linear micropolar solution,
> and the pullbacks of the conserved Euler KKS and Jacobi forms can be made to
> equal the prescribed quadratic canonical forms on that prepared section.
> The exact material momentum and angular-current balances are retained, and
> their difference from a canonical stress representative vanishes in bulk
> divergence/virtual work through `o(|K|^2)`.

That is a prepared finite-window linear response and second-variation theorem.
It is not an autonomous constitutive closure, a nonlinear invariant or slow
manifold, a prediction of the moduli from unrestricted Euler evolution, or a
finite-rotation Cosserat theory.  The exclusions in C-CST-017/018 say this
correctly; their exclusion of the stronger goal is not an in-scope defect.

I found no new contradiction in the final response, form or current identities
at their explicitly prepared scope.  A superficially contradictory pair of
claims about bounded `K` derivatives and a `|K|^-2` inverse is resolved by
tracking which map is normalized: 0250 bounds the source that produces the
quadratic acoustic-acceleration coefficient, whereas 0265 later inverts that
quadratic factor to produce a full order-one acoustic history.  The latter map
is necessarily singular and is used only after fixing nonzero `K`.  The
nonlinear objective remains open because several positive objects have never
been constructed, as listed below.

## Exact checks of the accepted proof chain

### 1. The Kelvin preparations are actual admissible Lin data

On the constant-curl patch, 0260 defines the full cohomological inverse

```text
xi^I     = lambda^-1 D^-1 w^I,
xi^a     = lambda^-1 D^-1 w^a + D^-1(xi^I Omega^{a prime}),
w        = curl(a_b f).
```

The shear term in the second line is load bearing.  Substitution into the
actual bracket gives `lambda [u,xi]=w`; taking divergence gives
`D(div xi)=0`, and the nonresonant Fourier sector (or the explicitly mean-zero
`n=0` extension) gives `div xi=0`.  Since the two compact fields
`xi cross omega` and `a_b f` have equal curl, their difference is a decaying
gradient and the whole-space Leray projections agree.  These steps are in
`0260/vector-lift.md:44-75` and `:148-179`.  They establish field-changing
isovortical initial data rather than a target matrix.

The cellwise Bloch rephasing is also exact.  For

```text
xi_K = exp[-i K.(x-X_c)] xi_0,
```

one has `div_K xi_K=exp[-iK.(x-X_c)] div xi_0=0`; multiplication by the
physical Bloch factor makes the phase constant inside each ring.  This is the
calculation in `0266/cellwise-phase.md:11-39`.  Periodic pressure is not exactly
cell-local: its image and mean terms are smooth high-carrier pairings.  The
claimed smallness is therefore finite-window, finite-inventory and ordered
after fixed nonzero `K`, exactly as the final claim states.

The arbitrary-order pressure construction has a real analytic route rather
than a local-pressure substitution.  The WKB recursion in
`0112/eps-core-floquet.md:141-215` enforces the divergence constraint at every
order, uses an exact curl potential for solenoidality, and leaves a compact
linear-Euler residual of arbitrarily high inverse-carrier order.  The periodic
transfer in `0265/assembly-and-chart.md:96-108` retains the periodic image
kernel and zero mode as smooth tests.  Its maximum verdict is a fixed-time
Lin estimate; it gives no nonlinear or acoustic-time stability.

### 2. The KKS/Jacobi signs and the simultaneous Gram algebra are consistent

For a solenoidal generator supported where `omega=lambda u`, let
`v_xi=P(xi cross omega)`.  Direct triple-product algebra gives

```text
Omega(xi,zeta) = rho integral omega.(xi cross zeta)
                = -rho integral v_xi.zeta.
```

Together with `curl v_eta=lambda[u,eta]`, the Lin generator satisfies
`L eta=v_eta-curl(v_eta)/lambda` and therefore

```text
H(xi,eta)=-Omega(xi,L eta)
 =rho integral v_xi.(v_eta-curl(v_eta)/lambda).
```

Thus the signs in `0262/local-normalizer.md:11-39` match the actual
Hamiltonian/Lin equations.  The high-carrier helicity blocks have opposite
definite `H` signatures and opposite KKS orientation under the stated
constant-curl and fixed-sign-patch assumptions.

The finite common-null construction in `0262/exact-gram.md:9-48` is an exact
identity.  With `V=e+ + e- + p`, `W=(e+-e-)/2`, `Q=q`, the column matrix
`Y=V+(Wh+Qo)/2` has `Y*HY=h` and `Y*Omega Y=o`.  The later affine contraction
correctly retains the baseline/auxiliary cross forms and the actual observation
pullback.  This realizes arbitrary quadratic forms on a chosen prepared
section; it does not show those forms arise as the restriction of an invariant
nonlinear Euler manifold.

### 3. The physical gain exists only with a singular acoustic preparation

The optical determinant is not inferred from nonzero polarization alone.
`0250/two-polarization.md:145-165` retains the material correction
`R_perp` and obtains

```text
det Gamma_opt = theta_u f1(f2-f1)
                (G_u R_perp-R_u G_perp),
```

which has a nonzero leading large-carrier coefficient for two distinct
positive tag fractions.  The acoustic gain likewise uses the complete central
second-momentum tensor.  After the centering repair, the parity-even entries in
`0266/cellwise-phase.md:132-185` yield the whole-law transverse coefficient
`(4 b_perp+2 b_parallel)/15`, with a nonzero fixed-large-radius limit.  The
cellwise phase makes the acoustic-to-optical first row vanish after the whole
orientation/reflection average; it does not say every realization has zero
spin.

Consequently the final leading gain is triangular:

```text
M_inf(K,omega) = [[h_T(K,omega) I_T, B_H(K,omega)],
                  [0,                    B_O(omega)]],
h_T = |K|^2 h_2 + O(|K|^4),     inf |h_2|>0.
```

Its inverse necessarily costs `O(|K|^-2)`.  This is explicitly admitted in
`0265/joint-continuum.md:79-105` and is the correct finite-`K` statement.
If `M_N=M_inf+E_N`, the exact Neumann condition is
`||M_inf^-1 E_N||<1`; hence invertibility requires the raw carrier error to be
smaller than order `|K|^2`, and a final post-inversion `o(|K|^2)` history error
requires raw error `o(|K|^4)`.  The final review states this stronger budget at
`0269/review.md:123-140`, and the arbitrary-order WKB diagonal can meet it after
each nonzero `K` is frozen.  This is the load-bearing order of limits.

### 4. Critical normalization audit: the bounded 0250 map is not the singular 0265 map

`0250/measured-third-jets.md:120-138` assumes a physical gain-inverse source
map that is `C^3` with bounded derivatives on `|K|<=K_0`.  Read alone, that
looks inconsistent with the final `|K|^-2` cost.  The definitions resolve the
issue.  In `0250/axisymmetric-hybrid.md:182-203`, 0250 first divides the
quadratic gain by `|K|^2`, calls the resulting nonzero coefficient `q_N`, and
uses bounded `q_N^-1` to produce the **quadratic acceleration**
`|K|^2 D`.  Its measured-third-jet source is normalized by the material/carrier
gain `N c_N`; it has not yet been multiplied by `|K|^-2` to produce an
order-one displacement history.

The later full-history inverse really is singular.  Already in the scalar
acoustic/optical leading subblock,

```text
M(K) = [[K^2 h_2, K b_1], [0,B_0]],
M(K)^-1_AA = 1/(K^2 h_2),
d M(K)^-1_AA/dK = -2/(K^3 h_2).
```

The exact check is `response-singular-gain-check.py`.  The active execution
returns zero and prints the pole.  Its first draft used a sign-sensitive
symbolic infinity assertion that SymPy could not decide for an unspecified
nonzero `h_2`; replacing it by the exact Laurent identities changed no
scientific predicate.  This check forbids transferring 0250's bounded-jet
language to the final full-history encoder.  P251 does not make that transfer:
`0265/joint-continuum.md:139-169` fixes each nonzero `K`, includes the finite
`|K|^-2` inverse and synthesis costs, and only then chooses the physical
carrier, WKB order and auxiliary carrier.  The mathematical repair is already
present.  The documentation would be clearer if 0250 named its bounded map
explicitly as `q_N^-1` and 0265 named the later multiplication as a separate
operator.

Three other real mistakes were already exposed and repaired in the campaign:

* `0250/schur-correction.md:1-37` correctly rejects the earlier claim that two
  `O(K)` off-diagonal rows are automatically a small perturbation of an
  `O(K^2)` acoustic block.  The cellwise-phase route constructs `c_1=0`.
* The original 0266 all-18-entry tensor omitted the central-moment terms
  `-delta X_j A_il-delta X_l A_ij`; 0267 found this.  The active 0266 source
  includes them, and they change only the reflection-odd chiral sector, not the
  parity-even Schur coefficient (`0267/review.md:120-150,171-208`).
* The first 0270 materialization incorrectly made `j` depend on kinetic energy
  and had a test that copied the error.  The repaired API uses the tag second
  moment (`0270/repair-review.md:1-20`).

These are useful exposing corrections.  Their repaired status should not be
used to erase the failed route evidence.

### 5. Initial form matching and the evolving observation chart

The initial equations in `0265/joint-continuum.md:180-200` use the same actual
state map `E` on both sides:

```text
H_phys=E* H_can E,       Omega_phys=E* Omega_can E.
```

For the evolved source columns `B(t)`, put `Z=B E(t)^-1`.  Since the stationary
Lin flow preserves `H_phys` and `Omega_phys`, and
`H=-Omega(.,L.)`, direct differentiation gives

```text
Omega(Z,Z_t)=-H_Z-Omega_Z E_t E^-1.
```

If the prepared history estimate holds in operator norm across the full
initial-state basis and its inverse remains bounded, then
`E_tE^-1=A_can+o(|K|^2)`, while pullback of the two conserved forms gives
`H_Z=H_can+o(|K|^2)` and `Omega_Z=Omega_can+o(|K|^2)`.  The leading connection
cancels because `Omega_can A_can=-H_can`.  This is a valid implication of the
history estimate plus both form matches.  Initial Gram matching alone would
not establish it, and P251 does not rely on initial matching alone.

More explicitly, write the observed state propagator as
`E(t)=F_can(t)E(0)+R(t)`, where the prepared operator estimates give
`R,R_t=o(|K|^2)`.  Since `F_can,t=A_can F_can`,

```text
E_t E^-1-A_can = (R_t-A_can R)E^-1=o(|K|^2).
```

The bounded inverse is therefore exactly the hypothesis needed to transport
the initial two-form/energy equality into the evolving observation chart; no
unobserved fast motion is discarded in this algebra.

The word “bounded” here applies to the **observed-state** inverse at a fixed
prepared stage.  It must not be confused with the singular source/gain inverse
as `K -> 0`; 0269 makes that distinction explicitly.

### 6. What “literal complete current” establishes

The exact material balances used in `0246/joint-current-bridge.md:86-120` are

```text
J_H,t=div F_full,
S_full,t=div N_full-ax F_full.
```

The superpotential change

```text
S_int=S_full-div Q,       N_int=N_full-Q_t
```

preserves the angular balance identically.  `0235/current-improvement.md:93-145`
also retains its boundary charge, lower endpoint and `q_t` memory and verifies
that the boundary action does not change translational momentum.  The prepared
output control, rather than conservation alone, supplies
`S_int=j Phi_t+o(|K|^2)`.

Subtracting the canonical equations then yields only

```text
P_T div(F_int-F_can)=o(|K|^2),
div(N_int-N_can)-ax(F_int-F_can)=o(|K|^2).
```

This supports bulk periodic virtual-work equivalence.  It does not give
pointwise equality of stresses/couple stresses or arbitrary boundary
tractions.  The nonstationary coefficient `q(t)` and its memory are actual
microscopic current data, not an autonomous material modulus.  C-CST-017/018
state these limitations accurately.

## Why the prepared result is not the nonlinear constitutive law

The key logical direction in P251 is

```text
choose a target linear branch and its coefficients
-> approximate its desired finite-window observation histories by frequencies
   from an actual Euler/Lin response band
-> invert the finite physical gain
-> add nearly unobserved Euler columns to match chosen quadratic forms.
```

This is a sophisticated controllability/section construction.  It shows that
Euler has enough prepared linear directions to imitate the selected reduced
law in the selected observations.  It does not establish the reverse causal
statement needed for constitutive emergence: that one fixed coarse-graining of
an open class of Euler states evolves autonomously according to a law whose
coefficients are forced by the microscopic dynamics.

The distinction is visible in `0265/joint-continuum.md:109-137`: `a`, `j`,
`ell_tag` and one band frequency are measured, but the Cosserat branch
histories and the recipe for `cT,cL` are fixed before the controls synthesize
them.  Different admissible target forms can be realized by the common-null
normalizer.  Therefore agreement does not select that constitutive form or
make its moduli predictive.

The historical record preserves another scope break.  The original issue-198
contract in `0271/historical-route-state.yaml:1-14` asked for effective moduli
from Biot-Savart self-energy, while its N4 representation explicitly excluded
nonlinear finite rotations (`:109-140`).  The terminal closure map acknowledges
that the original `alpha=L_v T/6` route is unsupported and that C-CST-018 uses
the different prepared relation `alpha=j nu^2/4`
(`0271/closure-map.md:13-21`).  The prepared theorem can be true without
closing either the original self-energy-modulus route or the stronger nonlinear
rotational objective.

## Missing positive objects for Federico's full nonlinear objective

These are absent load-bearing constructions, not objections to the accepted
linear claim:

1. **Finite rotation variable.**  `Phi` is the differential of a covariance
   eigenframe, reconstructed linearly after whole-law averaging.  A nonlinear
   theory needs an objective `R(X,t) in SO(3)` (or an explicitly equivalent
   director bundle), a chart-independent angular velocity, and treatment of
   eigenvalue collisions and global frame/gauge ambiguity.
2. **One encoder independent of a desired future history.**  Construct an
   Euler initial-data map from macroscopic initial data
   `(U_0,U_t0,R_0,Omega_0)` that is fixed by the medium and scale.  It may depend
   on the homogenization parameter, but not on a target solution, target
   coefficients chosen after the fact, compact time window, or desired output
   accuracy.
3. **Approximate invariance for nonlinear Euler.**  Show that the unforced
   Euler solution launched by that encoder stays near the encoded macroscopic
   family for a stated open set of amplitudes and a stated time scale.  P251
   has only Lin evolution around one stationary state and allows a different
   preparation at each stage.
4. **A nonlinear action/constitutive functional.**  Pull back the full Euler
   action, including incompressibility and pressure, to the nonlinear
   `(F,R,nabla R)` family and derive the energy and kinetic metric.  Matching a
   quadratic Hessian by invisible positive/negative-energy columns does not
   determine this finite-amplitude functional.
5. **Closure of unresolved modes.**  Prove that the complement produces a
   negligible or controlled memory term.  P251's explicit `q_t` and
   lower-endpoint memory show why this cannot be assumed.  If memory survives,
   the honest target is a nonlinear Cosserat law with memory/internal variables,
   not a local autonomous law.
6. **Nonlinear current and traction statement.**  Derive the coarse momentum
   and angular-current weak limit for arbitrary admissible nonlinear test
   histories.  Identify the stress/couple-stress representative and its
   boundary data.  Divergence agreement on prepared linear modes is not this
   constitutive closure.
7. **Uniform estimates and stability.**  Supply an error estimate over the
   intended macroscopic time scale, with dependence on amplitude, wave number,
   cell scale and source norms.  The existing fixed-window superoscillatory
   controls and `|K|^-2` encoder cost do not provide such a bound.

## Constructive continuation routes

### Route NL-A: nonlinear modulated-ring manifold

Use the compact stationary ring assembly as a microscopic reference, but
replace target-history synthesis by a single volume-preserving nonlinear
ansatz.  Parameterize each cell by a macroscopic placement and
`R(X,t) in SO(3)`, add the unique pressure/solenoidal correctors required by
slow `F=nabla y` and `nabla R`, and compute the full Euler residual before any
coefficient comparison.  The first achievements are:

1. prove the ansatz is a smooth map into `SDiff` and that its physical
   covariance frame equals `R` with a uniform spectral gap;
2. compute the exact pullback of Euler kinetic action and KKS form through
   cubic order in amplitude and second order in cell scale;
3. solve the nonlinear cell corrector equations and derive, rather than assign,
   `W(F,R,nabla R)` and the angular current;
4. prove a modulation estimate for unforced Euler on a stated amplitude/time
   regime.

The decisive failure test is simple: if two Euler states with the same
`(U,R,U_t,Omega)` generate different leading coarse accelerations after the
correctors are minimized, these variables do not close autonomously; append
the missing internal variable or move to NL-B.

### Route NL-B: exact projected dynamics with memory, then a Markovian limit

Start from the exact material momentum/angular balances already earned in
0232/0235.  Fix one nonlinear observation map and derive its exact projected
Euler equation by a Mori-Zwanzig/conditional-expectation or equivalent
Hamiltonian projection.  Keep the orthogonal dynamics, memory kernel and noise
as explicit objects.  Then test two positive subroutes:

* prove a scale-separation/mixing estimate that contracts the memory to a local
  objective Cosserat constitutive term with coefficients fixed by equilibrium
  correlations; or
* retain the measured memory as a nonlinear generalized Cosserat medium and
  prove that its linear short-memory limit reproduces C-CST-018.

This route fits the existing `q(t)` evidence naturally and prevents a boundary
current or hidden fast mode from being renamed a modulus.

### Route NL-C: nonlinear relative equilibria and symplectic reduction

Construct an actual finite-dimensional family of rotating/translating Euler
relative equilibria for the compact cores, with a nondegenerate Legendre map
from conserved impulse/angular momentum to rates.  Reduce the full KKS action
at fixed charges, then couple cells through the global pressure corrector and
homogenize the resulting modulated family.  This route can make finite
orientation and angular momentum intrinsic rather than observation-selected.
It earns the nonlinear objective only after transverse Euler modes are shown to
remain controlled; elliptic particle or WKB Floquet behavior is not sufficient
for nonlinear PDE stability.

All three routes should use the same success statement.  For an open bounded
set of macroscopic initial data, one fixed encoder produces unforced Euler
solutions whose fixed observation map converges to a unique nonlinear
objective rotational continuum, while the action, momentum and angular current
converge in compatible weak norms and every coefficient is fixed before the
solution is known.  That statement would turn P251's prepared linear evidence
into a real constitutive theorem.

`route_verdict: established for the response/action audit; the 0250 bounded
quadratic-coefficient map and the singular 0265 full-history inverse are
consistent once their normalizations are separated; the final C-CST-018
pointwise-nonzero-K diagonal remains supported at its prepared scope; the full
nonlinear rotational-Euler objective is blocked on the seven named positive
constructions above`

`evidence_scope: exact algebra and source-level implication audit of the
C-CST-017/018 response, form, observation-chart and current chain; no verdict
on the separately audited finite-R existence proof and no nonlinear no-go`
