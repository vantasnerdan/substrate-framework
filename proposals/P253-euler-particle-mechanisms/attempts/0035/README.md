# P253/0035: independent review of the 0030 column action and propagation claim

This README preregisters one coherent fixed-claim, non-author review of
root-owned attempt `0030`. The review joins four statements only because
`0030` presents them as one same-background construction: the physical
column/solitary energy--impulse--Casimir variation on its actual accessible
tangents; a positive weighted column energy and exact unitary full-pressure
axisymmetric evolution; exact phase- and group-speed ceilings by the critical
speed `c0`; and a fixed-mode spectral-stability transfer from the
Gallay--Smets class to an explicit smooth compact-vorticity profile through a
proved `C1` closure argument.

The target is the strongest theorem the actual equations and cited sources
support. A concrete failure in one link earns the minimum correction and does
not erase independently valid identities, operator results, or source-transfer
subclaims. Conversely, those surviving subclaims do not inherit completion of
the joined statement.

## Independence, ownership, and activation boundary

`particle-balance-review` did not author or implement `0030` and owns only
append-only review attempt `0035`. Root owns `0030`, central proposal and
registry files, source/API/test changes, memory, validation, and commits. The
reviewer will not edit `0030`, source modules, tests, central records, generated
files, or commits.

Before this freeze, the reviewer opened only `0030/README.md`, SHA-256
`1728735f7d6aa91179f3b264ee7b1b52c1430c186b9df03efdc80fabbb35180e`,
plus filename/tree metadata and Git blob hashes at the pinned commit below.
No `0030` scientific body, equation, result, validation text, implementation,
test body, or captured output has been opened. The exact two primary papers
were inspected only to freeze their applicability hypotheses before the
favored construction was visible.

No `0030` scientific body may be opened until root:

1. registers `0035` centrally with this ownership and frozen scope;
2. records the exact repository-schema invocation and its outputs in `0035`;
3. makes `attempts/0035/activation-schema.exit` exactly zero.

After activation, pin this README and all reviewed artifacts by SHA-256,
perform one substantive review pass, and allow at most one bounded correction
check if a concrete finding requires repair. The accepted base is `v0.183.0`;
the preregistration head and frozen implementation commit are
`1d378f9aebd10a8ff0dd91f26bf15d371f3e243e`.

## Frozen artifact inventory

The following exact Git blobs at `1d378f9` form the maximum review surface.
The SHA-256 values were computed from `git show 1d378f9:<path>` without
displaying their contents.

| Artifact | Git blob | SHA-256 |
|---|---|---|
| `0030/action-construction.md` | `bd01b5f8e8811d5dac0f03b324a083293881afab` | `56ee231804c1ccfbd60b44546a7943f6518731f52bdd7a8ab95f777ce1aa1fe6` |
| `0030/column-propagation.md` | `76b439df61b414afb3c0f8a204648cfe808f44c5` | `d50facb41126c1b0f31609dee04d733dfff720d34df4e83b59c5ca707f69de2e` |
| `0030/column-profile-transfer.md` | `d0b7576603a283f6b8bac04ec66cbcb4a3892c46` | `c4c9acc3ec29e709a0124a402e85ad688055cc945365b1788ceea23a95bd5fba` |
| `0030/result.yaml` | `12e96b0be4bc21bfeb876b9b47908723e08a2a84` | `43b741c6550807a0c1ae97307daa75c6abb7e849d1b21587e9162cc128debac8` |
| `0030/validation.md` | `a75d3c1c6062cfe2076b76d8d208ce9460f675fb` | `c13fb4c529e456fd48210b55d47d4a2380cc990b6ebb99a5ed1cf2682976e305` |
| `0030/source-inventory.md` | `90358ee35f24880bd4113591612e09e66670cb12` | `d47f5faf690c99f4e0703e184f4e015e135479a094a734d04e880ba2aeceff15` |
| `src/substrate_framework/euler_column_wave.py` | `7c10dd4664d89098eee58afd10f26c7a2c08dced` | `b66328cfafe82fb1f56eb2b52911786ddb5686718befc10ad34dfb0dc8c11dc0` |
| `tests/test_euler_column_wave.py` | `f06f281c3851231b746b370cb40d6f6fbeaabacd` | `d106dc0af6b1db2cb9575df7ea4ba5b856583280da66fa0f6358c35b76ca7974` |

The existing focused/full-suite receipts may corroborate only predicates that
the opened tests actually expose. The unchanged oracle need not be rerun unless
the scientific predicate is invalidated. Any current untracked or later module,
including `euler_axisymmetric_action.py`, is outside this frozen inventory.

Accepted `0027/0028` is an unchanged supplier at its pinned scope: it provides
the exact whole-space background and solitary-wave existence/exterior estimates
that `0030` explicitly imports. This review does not reopen that proof and does
not let existence supply action, spectral, unitary, speed, or stability facts.

Attempts `0034`, `0031`, `0032`, and `0033` are excluded in full. In
particular, no conclusion authored in `0031` and no canonical continuation in
`0034` may repair or prejudge `0030`.

## Frozen primary-source domains

The review freezes the exact source versions, hypotheses, and theorem domains
before opening the construction.

### Gallay--Smets 1805.05064v3

T. Gallay and D. Smets, *Spectral stability of inviscid columnar vortices*,
[`arXiv:1805.05064v3`](https://arxiv.org/abs/1805.05064v3), revised
2019-06-18, studies the stationary whole-space column

    u=V(r)e_theta,  omega=W(r)e_z,
    Omega=V/r,  W=r Omega'+2 Omega,

mode by mode in the divergence-free vorticity/enstrophy space
`X_{m,k} subset L2((0,infinity),r dr)^3`. Its stated profile assumptions are:

- `W` is `C1`, positive, `W'(0)=0`, strictly decreasing for every `r>0`,
  and has finite circulation `integral_0^infinity W(r) r dr`;
- `J=Phi/(Omega')^2`, `Phi=2 Omega W`, is positive and `C1`, with
  `J'(r)<0` for every `r>0` and `r J'(r)->0` at infinity.

Under those hypotheses the fixed-mode spectrum of the vorticity generator is
contained in the imaginary axis for each integer `m` and nonzero axial
wavenumber `k`. The `k=0` enstrophy issue and the paper's separately stated
periodic/zero-vertical-mean result must not be silently promoted to unrestricted
whole-space velocity evolution.

### Gallay--Smets 1811.07584v2

T. Gallay and D. Smets, *On the linear stability of vortex columns in the
energy space*, [`arXiv:1811.07584v2`](https://arxiv.org/abs/1811.07584v2),
revised 2019-06-18, treats the full linearized velocity-pressure equations in
the whole-space divergence-free energy space `L2(R3)^3`. Its `H1` strengthens
the profile regularity to `C2` and includes the same positivity and strict
monotonicity, finite circulation, `W'(0)=0`, and the stated decay
`r^3 W'(r)->0`; `H2` is the strict positive/decreasing Richardson-function
condition above.

The theorem gives a strongly continuous group with spectrum on the imaginary
axis and, for every `epsilon>0`, a bound
`||exp(tL)|| <= C_epsilon exp(epsilon |t|)`. It does **not** by itself give an
exact unitary group, a uniform-in-time bound, a polynomial bound, nonlinear
stability, or orbital stability. An exact unitary axisymmetric result in
`0030`, if true, must follow from its own positive conserved form, full
pressure evolution, and operator-domain argument.

### Applicability discontinuity to audit

The `0027`-compatible target described before activation is smooth, constant
near the axis, and compactly supported in vorticity. Such a profile does not
directly satisfy either paper's everywhere-positive, everywhere-strictly-
decreasing hypotheses. Therefore the claimed result can only be a separately
proved closure transfer, not a direct invocation of either theorem.

Local `C1` dependence of a radial ODE on bounded intervals is not alone a
global spectral-closure theorem. The review will require the actual topology
and quantifiers of the approximating profiles, control of the regular-origin
and infinity solutions, and an argument that any off-axis isolated eigenvalue
of the limiting fixed-mode operator would persist for approximants (for
example through an Evans-function/Rouche or uniform-resolvent argument). A
fixed-`(m,k)` vorticity statement does not inherit the all-mode velocity-space
group theorem.

## Frozen calculus and operator checks

### Physical action and accessible tangents

1. Derive the complete physical relative energy--impulse--Casimir functional,
   its first variation at the exact column/solitary field, and the claimed
   second variation. Track `rho_m`, `2*pi*r`, background subtraction, the
   translating-frame sign, all pressure/exact terms, and convergence of every
   whole-space integral.
2. Identify the actual phase space and distinguish arbitrary field variations,
   constrained streamfunction variations, and dynamically accessible
   coadjoint tangents such as `delta omega=curl(eta cross omega)`. Track the
   Leray/Hodge pressure tail of `delta u`; compact `eta` or `delta omega` does
   not generally make `delta u` compact.
3. State the kernel, stabilizer, Casimir leaf, translation/phase directions,
   and boundary/axis conditions. Positivity on a declared tangent subspace is
   retained at that exact strength; it is not positivity on an arbitrary
   finite-excess neighborhood unless the domain and closure are proved.
4. Separate the stationary column functional from any solitary-wave
   evaluation or asymptotic transfer. Finite relative observables and a local
   weak/KKS pairing do not by themselves produce a finite global symmetry
   action or orbital stability.

### Positive weighted energy and full-pressure evolution

5. Derive the `m=0` axisymmetric linearized Euler equations, including the
   elliptic pressure row, divergence constraint, `r=0` regularity, and
   infinity boundary terms, in the exact radial measure and Fourier convention.
6. Verify positivity, nondegeneracy/coercivity, and conservation of the stated
   weighted quadratic form. Identify whether it is equivalent to physical
   `L2` energy and on which mode/domain; a positive formal integral is not yet
   a Hilbert norm on the generator domain.
7. Prove the closed full-pressure generator is skew-adjoint, or give an
   equivalent well-posed two-sided isometric evolution argument. Exact
   conservation plus existence only forward in time does not alone establish
   a unitary group. Audit surjectivity, strong continuity, dense domain, and
   all conjugation weights.
8. Keep the result axisymmetric and in its declared weighted space unless an
   actual theorem proves more. Gallay--Smets subexponential growth is not a
   substitute for exact unitarity.

### Phase and group speed ceiling

9. Pin the dispersion/eigenvalue parameter, temporal branch, sign of `k`, and
   the definitions of phase speed and group speed. Treat `k=0`, equality,
   branch crossings, and differentiability separately.
10. Derive the phase ceiling from the exact variational/Rayleigh identity and
    the critical speed `c0`, with its parameter and normalization. Then derive
    the group ceiling from an actual derivative/Hellmann--Feynman or monotonicity
    identity. A phase-speed bound alone does not imply a group-speed bound.
11. Check that `c0` is the same supplier quantity used by the nonlinear
    solitary branch, rather than a rescaled or wall-domain eigenvalue. State
    whether the ceiling concerns discrete modes, all real `k`, or only the
    regular branch on its existence interval.

### Compact profile and fixed-mode spectral closure

12. Verify the explicit profile is smooth in Cartesian variables, nonnegative,
    constant to the necessary order near `r=0`, compactly supported, finite-
    circulation, and compatible with the accepted `0027` construction. Compute
    its `Omega`, `Phi`, and `J`, including transition and exterior regions.
13. Construct the actual Gallay--Smets approximants and verify their strict
    `H1/H2` conditions, normalization, circulation/tail behavior, and global
    convergence in the exact claimed topology. Mere mollification or addition
    of a small tail need not preserve `J'<0`.
14. For each declared fixed `(m,k)`, prove convergence of the relevant closed
    operators or radial Evans data at both endpoints. Audit critical layers,
    essential spectrum, isolated-eigenvalue multiplicity, and the compact
    profile's changed far field. Establish the contradiction that a limiting
    off-axis eigenvalue would generate a forbidden approximant eigenvalue.
15. Match the transferred conclusion to the correct source space. A `C1`
    fixed-mode enstrophy closure does not establish the `C2`, all-mode,
    velocity-energy group theorem or exact unitary evolution.

## Oracle, sensitivity, and claim boundary

The strongest oracle is the direct physical variation and PDE calculation,
the complete Hilbert-generator proof, the exact speed derivative identity, and
the source-hypothesis plus global spectral-closure argument. The pinned API and
tests are regression evidence only after their implemented predicates are
matched to those theorems.

Representative exposing mutations are: reverse the translating-frame impulse
sign; omit the Leray pressure tail on an accessible tangent; drop the radial
measure or density; remove the full pressure row; preserve a quadratic form
without proving a two-sided group; infer group speed from phase speed; treat a
compact profile as strictly positive; or replace global Evans/resolvent
convergence by local coefficient convergence. Each load-bearing predicate must
be exposed by calculus, an operator/source argument, or a genuinely sensitive
test.

No new soft eigenvalue, stability-window edge, tiny energy difference, or
small force is registered, and no production numerical verifier is designed;
the small-ratio-numerics prescriptions therefore do not bind this frozen
review. If an unanticipated numerical spectral margin becomes load-bearing,
it requires a separate error-budget freeze before use.

This review cannot establish nonlinear or orbital stability, all-time control
of arbitrary perturbations, a full three-dimensional all-mode unitary theorem,
a particle or mechanical-force interpretation, finite quantum spin, statistics,
electron/neutrino completion, or any result in excluded attempts.

## Verdict and completion contract

The joined `0030` route receives exactly one verdict: `established` at the
strongest fully proved scope; `refuted` only with a concrete contradictory
equation or source-domain mechanism; or `blocked` with the exact missing
action, operator, derivative, or closure construction named. The review will
also record each of the four components as established or remaining so that a
single weak link does not recursively weaken the valid mathematics.

After activation, write only `0035/review.md` and `0035/verdicts.yaml`, with
the frozen hashes, exact supported statement, evidence roles, four-axis status,
minimum concrete correction if needed, exclusions, and next construction.
Allow one bounded correction check at most, without reopening unchanged
claims or rerunning an unchanged oracle. Send artifact/verdict/remaining
dependency to `herdr agent prompt w3:p1`.
