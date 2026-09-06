# P253/0038 independent joined fixed-carrier review

## Frozen transaction and independence

This transaction independently reviews the two links of one fixed-carrier
claim. Attempt `0032` derives a hyperbolic returned Kelvin polarization
cocycle from the Baldi/Gavrilov source recurrence. Attempt `0039` lifts that
same cocycle to whole-space linearized Euler operator and essential-norm lower
bounds on finite-energy dynamically accessible data. The reviewer authored or
implemented neither unit and opened no target scientific body before central
activation.

The frozen README hash is
`5a7d63ac6b4ad62d9320ffb7ac52f9d07d208c6e062f4523ad11ea6ae4f3d7c1`.
Central registration names `particle-balance-review`; activation contains
exactly `0`, with SHA-256
`9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.
Activation stdout is
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`
and stderr is the empty-file hash
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`.
The accepted base is `v0.183.0`; the observed preregistration head is
`1d378f9aebd10a8ff0dd91f26bf15d371f3e243e`.

The substantive-pass scientific boundary was:

| Unit | Artifact | SHA-256 |
|---|---|---|
| A | `0032/derivation.md` | `a1c252d729716656d758f8449711463a671b41ba17c3eaf4d15351a180d86a73` |
| A | `0032/result.yaml` | `a25166e2ff0009313d3ba36ca8d09ac105b86c5e0a806891ac56a3a650865c8f` |
| A | `0032/validation.md` | `717b03be39c581b7b23b440ed4c9cb5272cfdef76f4c5f7d167e2784d58566a6` |
| A | `0032/verify_trace_metric.py` | `585c66b7f1eb4e71e56f80bba3ea2fa14866d789c80753ecf5840a04d4ef894f` |
| A | `0032/verify_pullback_coefficient.py` | `d2615996405ca5c75c72b8a5498e26dcbf25942729cd0bd6465a2fe313782ab3` |
| A | `0032/verify_c2_recurrence.py` | `4dae3367b62dbaf7641973ac6edc0b56513531deaadaaa3bb9a6d5c65b4d4d2b` |
| A | `0032/verify_c2_source_recurrence.py` | `194fa4b83f08689bd77dae3cc0da6f0a26b12ffe72dd570e716367754a5d95c8` |
| A | `0032/c2-first-success.stdout.txt` | `6fb56ad2baf1f7ebdbf2f67e88a000e80bca753c11c0a7859edfd0949c71e878` |
| B | `0039/access-inventory.md` | `72595f4e2fca7374509da106ecbdce21dcd4ddc26c2cfde6d927413f4649133f` |
| B | `0039/derivation.md` | `0a47bea9732796586208d31f38490eee16e2d72c8a4a315ffee4e7a7bdea4fbb` |
| B | `0039/result.yaml` | `9b2e33dfaaac5bf5d85f9b96f9738069bfbfcbe3ad11b1c2e432cb5055a1ff1b` |
| B | `0039/validation.md` | `b2d932b7470424fc32a433953b15d7d56c3651e047ede5df6081f0aaca9c557f` |
| B | `0039/verify_accessible_semigroup_lift.py` | `375f4c94c6098098baabdceb6beb3495471e20087f1c720082e798d9d5f762ba` |
| B | `0039/first-run.stdout.txt` | `4599bc20a997562ec7f1909b339e2488a3380b32ce848d2fb48aa84bf5c874b9` |

All recorded scientific executions exited zero. No canonical source module or
test is part of the actual `0032`/`0039` inventory. Existing executions were
not rerun merely to reproduce their tallies. Active `0037` was not opened.

## Primary-source applicability

The exact primary PDFs and hashes are:

- Gavrilov, *A steady Euler flow with compact support*,
  `arXiv:1810.08020v1`, SHA-256
  `fcaca85faa77e3876b11d16718037169fc026112dfd3e4248e03e963a0ebc3c9`;
- Baldi, *The dynamics of Gavrilov's compactly supported Euler flows*,
  `arXiv:2302.02982v1`, SHA-256
  `4fe82be0db15a7f117017e505e6682b3438ed14ebcf56225bd7e74bc759136f3`;
  and
- Shvydkoy--Vishik, *On Spectrum of the Linearized 3D Euler Equation*,
  Dynamics of PDE 1 (2004), 49--63, cached journal PDF SHA-256
  `9307efc6096565f74c5a099af25c834af09ba93b3e56ff34dd076ba289551592`.

Baldi Theorem 1.1 and equations (3.52)--(3.68), (4.29)--(4.61)
establish the analytic action-angle map, pressure/action relation, and source
recurrences. In particular,

    K(I)=I+(1065/1024)I^3+O(I^4),
    Omega_2/Omega_1=sqrt(I)(1+(7/4)I+O(I^2)).

They do not state a Kelvin trace theorem. That trace must be obtained by
differentiating Gavrilov's actual velocity, including the cylindrical
orthonormal-frame connection.

Shvydkoy--Vishik states its Euler results on `T^3` and the full
`L^2_div(T^3)` space. Its sentence that free space should not present major
difficulties is not a theorem on `R^3`, and the source does not prove
invariance of this proposal's dynamically accessible closure or descent to
its Euclidean-symmetry quotient. Attempt `0039` correctly uses the published
bicharacteristic-amplitude and weak-packet construction as corroboration and
supplies a direct fixed-time whole-space argument instead of silently
importing the torus essential-spectrum theorem.

## Unit A: source recurrence and physical cocycle

The source recurrence supports the action-angle and resonance data. On a flat
cutoff plateau,

    Omega_1=K'=1+(3195/1024)I^2+O(I^3),
    Omega_2'=1/(2sqrt(I))+O(sqrt(I)),

so derivative resonance `n Omega_1'=m Omega_2'` has

    m/n=(3195/256)I^(3/2)(1+O(I)).

The continuous positive ratio assumes arbitrarily small rational values.
One may therefore choose one sufficiently small exact resonant shell, choose
a flat cutoff plateau around its pressure, and then freeze one compact
carrier. The later repeated-circuit statement does not change that carrier.

The physical reduction has the correct geometry. In Baldi coordinates the
columns of `D Phi` are expressed in the rotating cylindrical orthonormal
frame, `g=D Phi^T D Phi`, and
`k=D Phi^(-T)(n,-m,0)`. The azimuthal covector contribution is relative
`O(epsilon^4)` when `epsilon=sqrt(2I)`, so it does not affect
`C_0,C_1,C_2`. The returned coefficient subtracts both the physical
cylindrical-frame angular velocity and the moving transverse-frame
connection. This repairs the earlier exploratory coefficient that omitted
frame rotation.

The dynamically accessible symbol is also correct. For a compact
divergence-free displacement with principal amplitude `a`,

    delta omega_pr=i[(k dot omega_*)a-(k dot a)omega_*].

On the selected small shell, the normalized leading geometry
`u=(Z,s/sqrt(2),-X)+O(s^2)` gives
`omega_* dot khat=-1/sqrt(2)+O(epsilon)`, hence the scalar is genuinely
nonzero. For every transverse velocity amplitude `A`, choosing
`a=(k cross A)/(k dot omega_*)` and applying the order-minus-one
Biot--Savart symbol returns `A`, with the displayed sign. Thus the expanding
Kelvin eigenpolarization lies in the principal image of true coadjoint
tangents; transversality alone is not being mistaken for accessibility.

The physical-frame source calculation gives

    C0=[[0,0],[-1/sqrt(2),0]],
    tr M=2+11*pi^2*epsilon^2+O(epsilon^3)
        =2+22*pi^2*I+O(I^(3/2)).

The order-`epsilon^2` trace separates into a single-`C2` contribution
`9*pi^2` and a double-`C1` contribution `2*pi^2`. Polarization-area
conservation, together with return of `|k|` and the oriented physical frame
under the axial symmetry identification, gives `det M=1`. Therefore

    Delta=(tr M)^2-4=88*pi^2*I+O(I^(3/2))>0

on all sufficiently small members of the exact resonant sequence. The
returned map is hyperbolic with reciprocal positive real multipliers
`lambda_+>1>lambda_->0`. This is stronger than the old near-Jordan
comparison and does not infer a sign from determinant one alone.

The generic metric classification is exact: an `SL(2,R)` return admits a
positive periodic transported metric precisely in the elliptic and scalar
`+/-I` cases, not in the hyperbolic or nontrivial Jordan cases. The separate
pressure-normal sector has an explicit identity return. Neither statement is
global carrier stability.

**Unit A verdict: established after the bounded source-oracle correction
described below.** The correction changes the pointwise displayed `C2` jet,
but an independent order-aware replay shows that it does not change either
trace contribution, the `22*pi^2` trace coefficient, or the `88*pi^2`
discriminant coefficient.

## Unit B: whole-space accessible Euler lift

The whole-space operator is the correct full-pressure linearization

    Gv=-P_L[(u_* dot grad)v+(v dot grad)u_*]

on `X=L^2_sigma(R^3)`. The skew-adjoint projected-transport realization and
bounded shear were already independently reviewed in `0019/0026`; the
maximal domain is not replaced by an unproved global `H^1` domain. Bounded
perturbation gives a strongly continuous group and the physical energy bound
used in Duhamel. Smooth coadjoint tangents evolve equivariantly, so their
kinetic closure `X_DA` is invariant. Smooth oscillatory tangents belong to the
part domain, while their energy-normalized graph norm is allowed to grow like
the frequency.

The packet scaling is consistent. A displacement of size `N^-1` has
order-one vorticity and order-`N^-1` velocity. Normalizing the actual
Biot--Savart velocity preserves exact accessibility and produces a unit
finite-energy input. The whole-space Leray order-zero principal symbol gives
the full Kelvin pressure term. Its lower-order composition, the separated
Hodge/Biot--Savart tails, tube-cutoff commutators, flat collar, and exterior
are estimated on each fixed interval; compact vorticity is never identified
with compact velocity. The source-independent content needed by the theorem
is a relative error tending to zero for fixed `T,delta`; the displayed
`N^-1 delta^-3 C_T` estimate supplies it.

For each integer circuit count `j`, `t_j=jT_*` is fixed first. The tube width
may then depend on `j`, and an integer oscillation frequency is chosen after
the tube. This is enough because an operator norm takes a supremum over input
data separately at every time. It does not claim one packet, one frequency,
or one WKB expansion uniform as `j` tends to infinity. Axial rotation after a
meridional circuit is an isometry commuting with the axisymmetric carrier;
the relative return can therefore be iterated while the actual packet center
rotates.

For fixed `j`, shrinking tubes and increasing integer frequencies give a
normalized sequence weakly converging to zero in `X_DA`, with output norms
tending to `lambda_+^j`. Compact operators send that sequence strongly to
zero, proving

    ||S(jT_*)||_ess,X_DA >= lambda_+^j.

The Euclidean symmetry tangent is finite-dimensional, smooth, contained in
the accessible closure by a compact cutoff of each rigid generator near the
carrier support, and invariant under the linearized group. Weak nullness at
input and output makes its finite-dimensional projections vanish. Hence the
same bound holds for the induced quotient norm. This is an essential-norm
argument, not a generator Weyl sequence; no generator approximate spectrum is
claimed.

Consequently,

    limsup_j (jT_*)^-1 log ||S(jT_*)|| >= log(lambda_+)/T_*>0.

The result concerns one exact linear semigroup. It does not give nonlinear
Lyapunov instability, one initial datum with an unbounded orbit, carrier
breakup, a universal carrier obstruction, or any particle or quantum claim.

**Unit B verdict: established after the bounded integer-frequency wording
correction described below.**

## Evidence and oracle audit

| Evidence | Proposition supported | Role | Limit |
|---|---|---|---|
| Baldi/Gavrilov PDFs | Exact carrier, chart, source jets, action/frequency recurrence | Primary exact source | No Kelvin trace or semigroup theorem |
| `verify_trace_metric.py` and first receipt | Full-pressure cross identity, area invariant, return classification | Exact algebraic corroboration | Generic finite-dimensional identities only |
| `verify_pullback_coefficient.py` and second receipt | Metric entries, determinant, cross-operator closure, resonant upper entry | Exact algebraic corroboration | Does not calculate source `C2` |
| `verify_c2_source_recurrence.py` and source receipt | Source recurrence through the trace-determining order | Exact symbolic oracle after bounded correction | Does not establish whole-space packet estimates |
| Shvydkoy--Vishik (2004) | BAS equation and weak-packet mechanism on `T^3` | Source applicability/corroboration | Not an `R^3`, accessible-subspace, or quotient theorem |
| `0039/derivation.md` plus reviewed `0019/0026` operator input | Fixed-time whole-space WKB/Duhamel and accessible weak sequence | Direct analytic proof | Linear finite-time-per-packet quantifier only |
| `verify_accessible_semigroup_lift.py` and receipt | Accessibility inversion, sign, multiplier algebra, per-time frequency choice, exponent | Exact finite-hinge regression | Does not prove the microlocal remainder or essential compactness bridge by itself |

The four original exact executions have zero exits and terminal pass tallies.
They were not rerun because their frozen predicates were unchanged. The
source recurrence verifier did require a changed execution after review found
an order-aware truncation defect; that bounded correction check is recorded
below. No production numerics are involved, and `small-ratio-numerics` does
not bind.

## Finding and minimum correction

The only correction package is bounded and preserves both headline verdicts.

| Finding | Direct evidence | Minimum repair | Effect |
|---|---|---|---|
| The original Unit-A verifier truncated `Lraw` with `mp(z,3)`, discarding cubic spatial terms before substituting `X,Z=O(epsilon)`. Since `H'/sqrt(H)=O(epsilon^-1)` in the swirl-gradient row, those terms contribute to pointwise `C2`. | An order-aware non-mutating replay with `mp(z,4)` changes `C2[1,0]` by `sqrt(2)(3 cos(q)^4/2-13 cos(q)^2/8-9/16)`. It leaves the single-`C2` trace `9*pi^2`, double-`C1` trace `2*pi^2`, and total `11*pi^2` unchanged. | Retain spatial degree three, assert the changed `C2` component, and pin the changed exact execution while preserving the old output as provenance. | Corrects the claimed pointwise source jet; hyperbolicity is unchanged. |
| `exp(i phi/h)` is not globally single-valued on the action torus for arbitrary real `h`. | `phi=n sigma-m beta` is angle-valued with integral periods. | Use `exp(i N phi)`, integer `N`, `h=N^-1`, and choose `N_j` above the same finite threshold. | Pure quantifier/wording repair; WKB and essential-norm bounds are unchanged. |

No further narrowing is supported by contrary evidence. In particular, the
torus scope of the published paper is already honestly separated from the
direct whole-space proof, and the exact linear-only exclusions are present.

## Correction check

The Unit-B correction is complete. Its final hashes are:

- `0039/README.md`:
  `c4f51f772e5eb17c134df23a47304560944693e828594ab761a1e9b654374b94`;
- `0039/derivation.md`:
  `4f78aa4a3f2057c2c883c17a885c13d76a6b95d64f44f4807e99af471cf90845`;
  and
- `0039/review-correction-0038.md`:
  `b1e226025ef3a9fce36ae86fd39b65ce11c20accfcc1d1ac4ce8ec41899319ae`.

The README and derivation use `exp(i N phi)`, `N` a positive integer, in the
injection, remainder, per-`j` threshold, and weak-null construction. The
receipt pins the pre/post hashes. The discrete unbounded sequence supplies
every limiting step formerly written with continuous `h`.

The Unit-A correction is also complete. The source verifier changed from
pre-hash
`194fa4b83f08689bd77dae3cc0da6f0a26b12ffe72dd570e716367754a5d95c8`
to post-hash
`74535827251d9cb3edc23e2a9c016984f7cc0e7b2df70f16fecb0d16d218ddc2`.
The append-only receipt
`0032/c2-spatial-jet-correction-receipt.md`, SHA-256
`3d51f187262a14c5b13bc957080a77efe4a5b34021220ebb0835cf0a656fa297`,
pins the original stdout and both corrected executions. The repository-
interpreter command has SHA-256
`7b47e0c64dd8a4b5d54b2fbbada624c65508c54dc628fac42947ceb0635112ad`,
its exit is zero, stderr is empty, and its stdout SHA-256
`712618c03a43aa2d931608120af61e7b184a9b9a7102ad604e7fc61a89e0a247`
is byte-identical to the first corrected execution. It prints the restored
pointwise component and all eight checks pass.

The final corrected Unit-A prose hashes are
`844a5a46dd1130e3e2b0ffacf3f6cb1c154ba16dcfde73584599b64694e079c2`
for `derivation.md`,
`107ff915fe03e08cb9040bd1c75f103fc0898a539905908da8bfbdce75922186`
for `result.yaml`, and
`fec9b1ae053501283f06dfa833a18586acac58d9b734f7c789d805dd3a0b112c`
for `validation.md`. No unchanged oracle was replayed.

## Four-axis decision and frontier

The final scientific decision is:

- Verification: exact source recurrence plus direct analytic whole-space
  microlocal/functional argument, with symbolic hinge corroboration.
- Review: established after one bounded correction package.
- Compatibility: compatible active-proposal evidence on the same fixed
  compact Euler carrier.
- Epistemic: exact linear fixed-carrier cocycle and semigroup/essential-norm
  theorem.
- Relationship: `0032` supplies the hyperbolic principal cocycle used by
  `0039`; neither inherits the other's missing bridge.

The strongest supported statement is: there exists one smooth compact
Gavrilov carrier and one sufficiently small derivative-resonant regular shell
whose full-pressure dynamically accessible Kelvin return has
`Delta=88*pi^2 I+O(I^(3/2))>0`; if `lambda_+>1` is its expanding multiplier,
then the exact whole-space linearized Euler group restricted to dynamically
accessible kinetic data, modulo Euclidean symmetry tangents, satisfies
`||S(jT_*)||_ess >= lambda_+^j` for every integer `j>=1` and has growth bound
at least `log(lambda_+)/T_*`.

The next parent achievement is not another qualification of this linear
theorem. For the robust-carrier objective it is a different carrier or
restoring sector with a positive periodic physical metric on the complete
accessible cocycle. Separately, any attempt to turn this adverse linear result
into nonlinear Lyapunov instability needs its own topology, differentiable
Euler flow map, modulation, and nonlinear remainder on a growth-compatible
time scale. Particle and quantum conclusions remain unavailable.
