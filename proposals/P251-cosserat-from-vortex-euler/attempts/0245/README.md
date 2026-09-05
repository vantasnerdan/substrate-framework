# 0245 — resumed independent same-cell periodic geometry review

Resume exactly the historical0237 transaction against immutable0236 candidate C,
with a new non-author reviewer in a separate Herdr process. The parent remains
the full issue200/P251 continuum, with its periodic and compact Euclidean
candidate classes distinct. Baseline main42ef058, research integration a2fcb46,
branch-local release v0.181.0. Frozen source/criteria:0237/README.md; the original
review never executed. Write new evidence only0245 and retain0237 unchanged.

Positive target: actual stationary finite-volume invariant periodic tube,
return-map flux action and nonzero twist on the same C016 field, at its
stated persistence/density scope. Strongest practical oracle is exact source
calculus and the precise primary persistence theorem; existing source/output
receipts may be reused. No expected conclusion, no new review of C016, and no
import of contractible Euclidean knot conclusions into the periodic claim.
The reviewer signs with a unique Herdr name/pane and actual independence.
Root owns any source repair, accepted claim and main integration work.

## Reviewer identity and independence

Reviewer: `herdr geometry-review pane w3:p4`, a fresh separate Codex process,
2026-09-05. I did not author or materially implement 0236, C-CST-016, the
0216 fixed-tag construction, or the current 0241 joint supplier. I received
the frozen transaction and criteria in 0237/0245 without an expected verdict.
My work was limited to reading the immutable sources, independently deriving
the geometry and return-map normalization, checking the precise primary EPS
theorem license, and recording this one substantive review in 0245.

This is the previously unexecuted 0237 transaction, not a second review. The
old 0237 reviewer name belongs to its preregistration; this fresh reviewer
inherits its scientific criteria, not that process's identity or involvement.

## Frozen transaction and positive role

The reviewed object is candidate C in immutable commit
`e88763a07f4e9a1ad67a9ae849bea237b0119c75`: the same C016 periodic field,
its finite-radius axial vortex tube, its flux-action return-map twist, the
local analytic persistence implication, and the exact whole-field density
normalization. The source is unchanged at the current issue200 integration
boundary. The accepted input release is branch-local `v0.181.0`; C-CST-016
supplies the fixed stationary cell and its already reviewed physical
preparation, not the new geometry theorem.

The frozen evidence boundary is:

| Artifact | SHA256 |
| --- | --- |
| `0236/README.md` | `ecd8f7efe8cc2f88cb0a20943316c3344eea2fe35f2485c0af0e32167566a855` |
| `0236/periodic-core-and-density.md` | `3ec7e27792e64bdab873c783e680750275c13b7769c6832931caf2d242a9792a` |
| `0236/verify_periodic_core.py` | `ac34a00446bc0fe673d06c8e4c8c44b1fbc4d857d5fe09b43a7c62212173b385` |
| `0236/core-first.stdout` | `c902a76263274bdaa1675d2fd7598e9bcf9be2713390245d9710ccd896c8d1e5` |
| `0236/core-receipt.md` | `74f41ecaf41221dd15854d7ccc023b6eb8ec9d46459a5358decf9841c0f34fb0` |
| `0216/fixed-tag-harmonics.md` | `41df6bf54baaee168ebe0a85b3e1b6fbbcf0e43ec6ef88921a61a16034402c9c` |
| primary EPS PDF | `309808ccd801552587a04a13f8cea3a85894511c3387a3f4fa8c5eb9e6d5e738` |
| primary EPS extracted text | `854f09f74d8806d5dd694c6f72c85eb8a63c3b37d443d8fc331b94dd41f549a1` |

The primary theorem source is Enciso--Peralta-Salas, *Existence of
knotted vortex tubes in steady Euler flows*, arXiv:1210.6271: Definition
7.5 and Theorem 7.6 give the invariant-circle normal torsion and KAM
persistence statement; equations (7.28)--(7.32) in the proof of Theorem
7.10 give the nearby-flux Moser identification and conjugated return map.
Theorem 7.6 is the paper's explicit planar specialization of de la Llave
et al. Theorem 47. This review uses the statement and reduction printed in
the archived primary EPS paper and does not enlarge it into a global knot
or steady-Euler continuation theorem.

## Strongest supported positive statement

Fix `Omega=1/10` and the orientation-preserving normal-coordinate rotation
`a=-Z`, `b=Y`. Then the accepted C016 field is exactly

```text
psi = cos(b) + Omega^2 cos(a),
u = (psi, sin(b), -Omega^2 sin(a))
```

on the flat three-torus. It is analytic, divergence-free, satisfies
`curl u=-u`, and is an exact stationary Euler field with
`p/rho=-|u|^2/2`. For `E=1+Omega^2-psi`, the identity `u dot grad E=0`
and the positive normal Hessian `diag(Omega^2,1)` give, for every fixed
sufficiently small `E_b>0`, an actual positive-volume invariant solid
torus

```text
T(E_b) = S^1_X x {E<E_b}.
```

Its core `a=b=0` is a nonzero periodic streamline and vortex line with
physical period `2pi/(1+Omega^2)`. Its class is the primitive axial class
in `H_1(T^3)` and is noncontractible. The lift is an unbounded line in
Euclidean space; it is not a compact knot.

On the section `X=0`, let `theta` be the true transverse time angle and
`omega(E)>0` its physical frequency. Since `u_X=1+Omega^2-E`, the exact
first-return map is

```text
P(E,theta) = (E, theta + 2pi r(E)),
r(E) = omega(E)/(1+Omega^2-E).
```

The preserved section form is the actual flux
`Lambda=psi da wedge db`, not unweighted normal area. With

```text
J(E) = (2pi)^(-1) integral_{D(E)} psi da db,
```

coarea and the transverse orbit period give
`J'(E)=(1+Omega^2-E)/omega(E)>0`; hence `Lambda=dJ wedge dtheta` with
the selected positive orientation. Quartic Birkhoff averaging is valid
without a cubic correction because this Hamiltonian is even, and gives

```text
omega_E(0) = -(1+Omega^2)/(8 Omega),
r_J(0) = (-1+6 Omega^2-Omega^4)/(8(1+Omega^2)^3)
       = -235025/2060602.
```

Continuity therefore supplies a fixed positive-radius interval on which
the flux twist is strictly nonzero. Since `r` has interval image there,
one can choose a fixed regular `E_b` in that interval for which `r(E_b)`
is Diophantine. No arithmetic property is attributed to the rational
core rotation `r(0)=10/101`.

The local persistence statement is then exact at the following scope.
For `J_b=J(E_b)`, introduce the flux-polar radius
`s=sqrt(J/J_b)`. Near the selected boundary,

```text
Lambda = 2 J_b s ds wedge dtheta,
P(s,theta) = (s, theta + 2pi r(J_b s^2)).
```

Definition 7.5, with identity boundary conjugacy, therefore gives the
normal torsion

```text
N_P = integral_0^{2pi}
        [partial_s P_theta(1,theta)/(2 J_b)] dtheta
    = 4 pi^2 r_J(J_b),
```

which is nonzero. Thus the computed flux twist is precisely the EPS KAM
nondegeneracy, up to a positive coordinate normalization; it is not a
surrogate area twist. The return map is analytic on a slightly larger
disk, its chosen boundary is analytically conjugate to the Diophantine
rotation, and its flux density is positive.

Consequently, for every requested finite `C^m` closeness there are a
sufficiently high finite `k` and a positive closeness threshold such
that an analytic divergence-free field defined on the slightly larger
solid-torus neighborhood, `C^k`-close enough to this field and still
transverse to the section, has a nearby invariant torus boundary. Its
Poincare map initially preserves its own nearby positive flux form;
equations (7.28)--(7.32) supply the local Moser diffeomorphism that
conjugates it to a map preserving `Lambda`, after which Theorem 7.6
applies. The invariant circle suspends to the invariant torus. This is
a local divergence-free persistence license. It neither constructs a
nearby stationary Euler solution nor imports the paper's global
Euclidean-knot construction.

The core return multipliers are
`exp(plus/minus 2pi i 10/101)`. They have unit modulus but are unequal
to `+1` (and to `-1`), so the return-map implicit theorem separately
persists a nearby fixed point. Preservation of a positive area form and
the strict trace margin preserve ellipticity for sufficiently small
perturbations. This periodic-core statement is independent of the KAM
boundary statement.

Finally, the exact unperturbed field itself supplies positive stationary
density. Uniform translation phase in its full periodic cell, followed
by Haar whole-`O(3)` rotation and the declared time-reversal pairing of
the entire field, lattice, tag, and preparation, is a translation-
stationary isotropic law of exact smooth Euler fields with the complete
ambient fluid and pressure. The tube volume fraction is
`|D(E_b)|/(2pi)^2>0`. One primitive axial core of length `2pi` per cell
has line length per volume `1/(2pi)^2`; hence a positive measured
spin/current coefficient `j0` per axial length gives the isotropic
coefficient density

```text
mean_tag_fraction * j0 / (3(2pi)^2) > 0.
```

Here the factor `1/3` is the vector Haar Gram and positive `j0` is the
literal fixed-tag fourth-harmonic coefficient derived in 0216. The
geometry calculation establishes this density normalization and its
independence from preparation-band narrowing. It does not, by itself,
turn the measured `G/S` coefficient into an unrestricted canonical
inertia; the later action-normalization supplier retains that separate
role.

## Evidence and oracle audit

| Evidence | Proposition established | Role | Limit |
| --- | --- | --- | --- |
| `periodic-core-and-density.md` plus the independent derivation above | Exact field, invariant finite tube, return map, flux action, twist, topology, and law normalization | exact proof | The persistence implication still consumes the primary KAM theorem; the density is for the exact whole-field law |
| EPS Definition 7.5, Theorem 7.6, and equations (7.28)--(7.32) | Same-measure analytic invariant-circle persistence and the nearby-flux conjugation used at field level | primary theorem import and applicability | Local analytic divergence-free perturbations only; no global knot construction or stationary-Euler inverse |
| `verify_periodic_core.py` and `core-first.stdout` | Fifteen exact algebraic anchors, including Euler residuals, weighted flux, Birkhoff coefficient, flux twist, multiplier, energy-density, and line-to-volume normalization | corroborating symbolic regression | The final density row assumes positive `j0`; it does not derive 0216's material coefficient or prove KAM |
| `0216/fixed-tag-harmonics.md` | Nonzero positive fourth-harmonic measured `j0` on a fixed positive tag and fixed regular annulus | exact dependency source | Does not alone identify `j0` with every action-mass notion |

The fifteen-check output is a complete successful first execution and its
hash agrees with the 0236 receipt. The checker derives the load-bearing
coefficient from the field and Hamiltonian rather than inserting only the
final rational. Its explicit wrong-axial-speed comparison exposes the
normalization error that would arise from differentiating `omega` alone.
The theorem hypotheses, topology, and density bridge were audited from the
source prose and primary theorem rather than credited to the check tally.

No new computation was warranted. The only plausible uncovered algebraic
risk was the conversion from `r_J` to the primary theorem's normal torsion;
the exact calculation `N_P=4pi^2 r_J(J_b)` above resolves it directly.
Rerunning the unchanged symbolic script would not add an independent oracle.

## Findings and minimum repair

No false statement, circular or absent load-bearing construction, invalid
dependency use, or broken affected consumer was found in the frozen periodic
supplier boundary. The strongest supported statement therefore needs no
scientific scope reduction and no source repair.

Two precision points govern downstream wording:

- `nonunit core multiplier` means multiplier unequal to the identity
  multiplier `+1`; the multipliers actually have unit modulus. The precise
  wording above should be used downstream. This is not a scientific defect.
- The density result is the positive tube fraction and positive measured
  `G/S` coefficient per physical volume of the exact randomized C016 law.
  Uniform persistence of that density under an arbitrary perturbation, a
  compact-ring density, and a standalone canonical-inertia identification
  are not conclusions of 0236. Restoring any of those stronger statements
  would require, respectively, a uniform marked-law perturbation argument,
  the separate stationary Euclidean-ring ensemble construction, or the
  complete action-normalization bridge.

The 0241 joint-residual table consumes only the same-cell periodic geometry
at this scope and keeps the Euclidean-ring and joint-dynamics obligations
separate. Its dependency use is therefore compatible with this review.

## Precise parent supplier license and decision

Issue200/P251 may import 0236 candidate C as exactly this supplier:

> On the same accepted C016 stationary periodic Euler field, there is a
> fixed positive-volume solid torus around the primitive noncontractible
> axial periodic vortex line whose exact Poincare map preserves physical
> section flux and has nonzero flux-action twist on a fixed Diophantine
> boundary. That boundary, and separately the elliptic periodic core,
> persist under sufficiently small analytic divergence-free perturbations
> on a slightly larger local tube with the finite `C^k`/transversality and
> nearby-flux Moser hypotheses stated above. The exact whole-field
> translation/`O(3)`/time-reversal law is stationary and isotropic and has
> fixed positive tube fraction and positive measured spin/current
> coefficient density `mean_tag_fraction*j0/[3(2pi)^2]`.

This supplier licenses the periodic same-field geometry/density row in the
joint continuum construction. It does not license a contractible compact
Euclidean knot, a stationary positive-density ensemble of compact rings,
the global EPS knot theorem for this axial class, existence of a nearby
stationary Euler family, the acoustic or optical history/action/current
join, or parent completion.

- Verification: exact analytic derivation with symbolic corroboration and
  a primary-theorem applicability audit.
- Review: independently audited and supported at the stated supplier scope.
- Compatibility: compatible extension of branch-local C-CST-016 on the
  same periodic field; no accepted statement is changed by this review.
- Epistemic: established route inside the active P251 campaign, not accepted
  canon merely by this artifact.
- `route_verdict: established as stated`
- `evidence_scope: independently reviewed same-C016 periodic finite tube,
  exact flux-action twist, local analytic persistence, and positive exact-law
  tube/measured-coefficient density`
- Correction check: not needed.
- Parent verdict: active; the compact Euclidean topology/density and complete
  coupled continuum obligations remain distinct.

Signed: `herdr geometry-review pane w3:p4`, fresh independent non-author
reviewer of immutable 0236, 2026-09-05.
