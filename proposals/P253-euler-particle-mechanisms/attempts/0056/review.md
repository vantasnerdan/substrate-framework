# P253/0056 independent review of the 0051 Schwinger--Hopf construction

## Frozen transaction and correction provenance

This transaction independently reviews root-owned attempt `0051`. The
reviewer authored or implemented neither the target, correction, API, tests,
nor receipts. The frozen review README has SHA-256
`5548cab7500fed2f1318a1a1c4847c1461a3c8cf6888204a525ea1ab4149a622`.
Central activation exited exactly `0`; its command, stdout, empty stderr, and
exit hashes are respectively
`5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f`,
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`,
and `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.
The accepted authority base is release `v0.183.0`.

The substantive pass found one factor-of-two physical-action defect and
requested one bounded correction. The receipt
`0051/action-moment-map-correction.md`, SHA-256
`993c360d234ef466c485e8261e838b9918861ed0fec1239ac89aa84e1d0ed638`,
pins the following complete map:

| Artifact | Frozen SHA-256 | Corrected SHA-256 |
|---|---|---|
| `0051/README.md` | `321f8429919a3b6c201693ff2220ef00114fcf44958dc910d0a00be6420c448b` | `d8e7b92b388dc1a5f1efdbbf00913de32012b99ff58ef22a3fab0309f1cfdc50` |
| `0051/construction.md` | `8df69d4c88f675af7fae4446656856d1b5f4b5069c544bff64bc3bdf8024d25d` | `841ca0caf6516e5955b2de89630bcf05402cf54473773c70b679308a23cd0103` |
| `0051/result.yaml` | `a889c0fb325bdb08db640dcc35dd5244c8a54be8d2db8fd471db48ce145dbe07` | `f47398cb566bd3d3b4edc6952b6f5682ae4ee8c1f3b4c2161fedc73471948931` |
| `0051/source-audit.md` | `e355dd7046e93595ff1f8062f391bd9f28380e11391ebe240c382408d11db50b` | `82efd88350342fbcd726c0dfe262b35338ccd778b30fb16ad513ce6cd0f22da5` |
| `0051/validation.md` | `3cd41bd877efcdc881be28eca92b1784c10caa53beda8729a181a7fe1be0fae8` | `57e5ff886141c7e2077ceea89ee58d4f6256d1c4bbd9401b569a20d04b0e3139` |
| `euler_schwinger_hopf.py` | `3e3a9981e8be6ca7f21cc6e30087dca2924c624e40c4ac019f08d7dc5b6f9e4e` | `b91bbb6f39750cea2e846574c97d5c1f74ca86e91231b35db164869a25eea3c3` |
| `test_euler_schwinger_hopf.py` | `9868f214c652d8629b27ea7fe524e6fdd344e9d7022cd31bb4137c3fac257ee6` | `185ecb8fe1cdcb99202e913d94ac28abb38d562677862ecf5cc7e5409575f65c` |

The single correction check was limited to the changed normalization and its
direct evidence edges. No second substantive pass was performed.

## Accepted supplier applicability

The actual Euler suppliers support exactly the boundary stated after the
correction:

- `C-CST-008` supplies one positive two-real-dimensional dynamically
  accessible KKS/Hessian tangent block with physical observation rows. It is
  one classical oscillator, not two complex modes, and is not an unrestricted
  invariant Euler subsystem.
- `C-CST-011`, `C-CST-017`, and `C-CST-018` supply actual prepared physical
  histories, inherited action/current forms, full pressure, and same-ensemble
  observations on each fixed compact time window. Their two real optical
  polarizations are a canonical pair for one oscillator. Preparation may
  depend on history and accuracy; these claims do not supply a four-real-
  dimensional autonomous invariant doublet or physical `U(2)` mixer.
- The independently reviewed `0030/0035` column theorem supplies a positive
  full-pressure axisymmetric linear energy space and exact unitary evolution.
  Sine and cosine at one Fourier/radial mode are its two real quadratures, not
  two complex modes. Distinct radial modes are genuine linear spectral
  candidates, but no accepted result turns them into a localized persistent
  full-Euler doublet with protected degeneracy and controllable off-diagonal
  mixing.
- The completed `0046/0050` review supplies one exact positive classical KKS
  oscillator, its affine Heisenberg algebra, and the finite-CCR obstruction.
  It expressly leaves a separately realized compact `su(2)` reduction open.
- The completed `0043/0047` review supplies exact positive-sphere KKS/Chern
  calculus and conditional `CP^1` quantization while leaving action, state,
  and character selection open.
- The completed `0048/0053` review preserves exact Cao cells and finite
  spectral-window sequence scaling, but not the graph-Riesz/full-complement
  construction required for a physical invariant ring doublet. Active `0052`
  is not evidence here.

Thus `0051` neither overlooks an accepted localized doublet nor weakens an
existing one. Its algebra is a sufficiency theorem conditional on a future
supplier.

## 1. Physical energy--action normalization

With

    Omega=B sum_a dq_a wedge dp_a,
    {q_a,p_b}=delta_ab/B,

the positive diagonal oscillator Hamiltonian is

    H_a=(B/2)(p_a^2+nu_a^2 q_a^2).

The symplectic change

    Q_a=sqrt(nu_a)q_a,
    P_a=p_a/sqrt(nu_a),
    z_a=sqrt(B/2)(Q_a+iP_a)

gives

    {z_a,conjugate(z_b)}=-i delta_ab,
    H_a=nu_a |z_a|^2,
    J_a=H_a/nu_a=|z_a|^2.

Consequently `z_a` has square-root-action units and

    J=J_1+J_2=sum_a |z_a|^2

is the physical total action. It generates the common phase by
`{z_a,J}=-i z_a`, with period `2*pi` up to orientation. For unequal positive
frequencies,

    H=sum_a nu_a J_a

and `J` remains conserved, but energy is not proportional to `J`; a generic
unitary mixer preserves action while preserving the free Hamiltonian only if
it commutes with the frequency matrix. Exact degeneracy makes
`H=nu J` and turns free evolution into one common phase.

The original API returned `J/2` as “total action.” The correction now uses
`J` consistently and adds `diagonal_mode_energy` to expose the
unequal-frequency bridge.

**Verdict: established after the bounded factor correction, conditional on
the supplied physical Williamson/KKS modes. KKS scale alone still does not
fix these energy--action coordinates or select a numerical action.**

## 2. Exact Schwinger--Hopf reduction

Define

    S_i=z^dagger sigma_i z/2.

Direct expansion with the preceding Poisson bracket gives

    {S_i,S_j}=epsilon_ijk S_k,
    {J,S_i}=0,
    |S|^2=(J/2)^2.

For `J=J_0>0`, the level is the radius-`sqrt(J_0)` three-sphere. The common
phase action is free, and its quotient is globally `CP^1`, diffeomorphic to
`S^2`. In the inherited positive orientation the reduced KKS sphere has
radius `J_0/2` and area

    integral Omega_red=2*pi*J_0.

The prequantization condition is therefore `J_0/hbar=N` in the integers;
ordinary holomorphic quantization would give spin `N/2` and dimension `N+1`
only after an external action unit, integral class, and quantization rule are
supplied. No such selection follows from the classical quotient.

**Verdict: the classical `C^2 // U(1)=CP^1=S^2` algebra, Poisson brackets,
radius, and area are exactly established after correction. Quantization is a
separate conditional construction.**

## 3. Complex-linear and antilinear split

For any real-linear compressed generator `V` and complex structure `J` with
`J^2=-1`, set

    V_C=(V-JVJ)/2,    V_A=(V+JVJ)/2.

Then `V=V_C+V_A`, and direct multiplication gives

    V_C J=J V_C,    V_A J=-J V_A.

The projections are unique onto the commuting and anticommuting real-linear
subspaces. On an actual compatible Hamiltonian mode block, `V_C` is the
complex-linear, number-preserving component; `V_A` mixes `z` with
`conjugate(z)` and is the squeezing/action-drift component. The latter
interpretation presupposes that `V` is the physical compressed Hamiltonian
generator, not merely an arbitrary matrix.

**Verdict: established exactly at the declared algebraic operator scope. It
does not itself prove that a physical Euler perturbation closes on the
doublet.**

## 4. Physical Pauli axes

A Hermitian `2x2` Hamiltonian gives a mathematical Stokes rotation and
preserves `J`. The API correctly accepts only a supplied unitary matrix. It
does not materialize that matrix as an Euler operation.

No current supplier constructs two independently tunable, number-preserving
Euler interactions whose traceless compressions have noncollinear Pauli
vectors while preserving the same constraints, action level, pressure-coupled
domain, and complement estimate. One avoided-crossing coupling gives at most
one physical axis; free detuning supplies another algebraic axis only when its
same-carrier control and interval are physically available.

**Verdict: blocked by two noncollinear physical number-preserving Euler
controls. The abstract `su(2)` representation is a sufficiency result, not
controllability or an analyzer.**

## 5. Gate-time leakage and action drift

`0051` correctly names but does not claim a full gate estimate. It supplies no
physical projected Euler operator `PAP`, complement propagator, residual
`QAP`, gate duration, or Duhamel bound. In particular, a nominal mixing time
that grows like the inverse coupling can amplify a small instantaneous
leakage or antilinear term. Nor is preservation of the finite matrix norm a
bound on the physical action of the full state when complement leakage is
present.

The required positive result remains an estimate, in one declared physical
energy/graph norm, of phase error, `Q U(T)P` leakage, the integrated `V_A`
effect, and absolute or relative drift of the physical `J`, with unequal
frequencies and full-pressure propagator growth retained over the actual gate
time.

**Verdict: blocked by the missing physical gate generator and finite-time
leakage/action-drift estimate. No all-time or nonlinear conclusion follows.**

## 6. Euler realization and downstream boundary

The route-specific statements are accurate:

- Route A has genuine positive spectral ingredients but no established
  localized persistent invariant full-Euler doublet, protected degeneracy or
  controlled splitting, physical off-diagonal mixing, or complement bound.
- Route B conditionally supplies the classical product `C^2` and Hopf
  reduction for two independent cells. It is a collective relative-amplitude
  sphere for two carriers, not an internal state of one carrier, and no
  same-field coherent coupling or exchange path is constructed.
- Route C supplies one physical oscillator and useful prepared classical
  readout rows on fixed windows. Prepared source-to-history control is not an
  autonomous `U(2)` operation or invariant doublet.

The exact classical construction therefore evades the narrow finite-CCR
obstruction if two positive modes are physically supplied, without asserting
that present Euler evidence supplies them. Continuous Euler scaling selects
neither `J_0` nor `hbar`. The quotient supplies no invariant probability
measure, Born/event law, reset/repreparation, two-copy exchange character,
particle identity, finite-speed band, or Lorentz cone. Elliptic pressure and
Galilean kinematics remain unchanged.

**Verdict: the Euler nontransfer boundary is established. The physical
Schwinger--Hopf analyzer remains blocked by a localized invariant doublet,
two physical axes, and gate-time complement control; this is not a P4 no-go.**

## API and oracle scope

The corrected API has SHA-256
`b91bbb6f39750cea2e846574c97d5c1f74ca86e91231b35db164869a25eea3c3`.
It implements the conditional canonical doublet, physical total action,
unequal-frequency diagonal energy, Stokes vector and Hopf residual, reduced
area, and supplied unitary mixing. The corrected tests have SHA-256
`185ecb8fe1cdcb99202e913d94ac28abb38d562677862ecf5cc7e5409575f65c`.

The focused direct-consumer command is preserved with SHA-256
`d8d75be024e6ecd2c6056ce7c1abcbd82e9fa2c7936b1c10172a25edde9a4f0c`.
Its stdout reports `14 passed in 2.00s`, SHA-256
`ef3f62b56e136ea2c6fe5868700df626f91d87ca9ec74e69f125119686e4eae7`;
stderr is empty and exit is `0`. The symbolic test now checks unequal
frequencies, `H=sum nu_a|z_a|^2`, `{z_a,J}=-i z_a`, the Stokes Casimir with
radius `J/2`, and the area `2*pi*J`. The review reused this corrected receipt
and did not rerun it.

These predicates establish finite algebra and domain behavior only. They do
not prove the Euler supplier, invariant projection, physical controls,
leakage, quantization, detector, or propagation claims.

## Final verdict and next dependency

`0051` is **established after one bounded normalization correction** at scope
`EXACT_CONDITIONAL_SCHWINGER_HOPF_SU2_AND_EULER_SUPPLIER_BOUNDARY`. Two
physically normalized positive canonical modes give an exact classical
Schwinger--Hopf sphere with physical total action `J`, Stokes `su(2)`, and
reduced area `2*pi*J`; the complex-linear/antilinear split gives the exact
number-preserving versus squeezing criterion.

Separately, current Euler suppliers do not provide the required localized
invariant doublet, two noncollinear physical mixing axes, or a gate-time
leakage/action-drift bound. The next construction is to produce those three
objects on one persistent carrier with full pressure and exact physical
normalization. Only after that result should action selection, detector
dynamics, exchange, or a finite-speed/Lorentz mechanism be joined.

The correction check found no stale half-action, doubled-area, or old
prequantization claim in the changed claim-bearing artifacts. The result YAML
parses, terminal-newline and scoped diff checks pass, and no further correction
is needed.
