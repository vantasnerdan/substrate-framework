# P248 primary-source and predecessor audit

This audit checks only literature-dependent premises used by the issue-184
paper.  Framework claims are inputs only at their pinned accepted scopes; they
are not evidence that the paper is correct.

## Sciama 1953 and 1964

The 1953 primary paper explicitly presents a tentative flat-spacetime,
Maxwell-type vector theory of inertial induction.  Its cosmic relation follows
only after a homogeneous Hubble model, a superluminal-distance cutoff, and an
approximate numerical factor are assumed.  The paper itself calls the scheme
incomplete and says a tensor potential is needed.  It therefore supports
P248-S21 only as a conditional dimensional motivation, not as a cosmological
identity or a derivation of Newton's constant.

Audited copy of D. W. Sciama, *On the Origin of Inertia*:

- official record: https://doi.org/10.1093/mnras/113.1.34
- SHA-256: `8094902abf3640416d693d20a1d89daa8d4660fef32626a64a2c8f8e46fdfcfe`
- pages: 9
- relevant source locations: summary, sections 1--4, equations (1)--(8)

The 1964 primary review uses a contravariant metric density, a flat-background
wave equation, and a divergence condition, and attributes equivalence with the
Einstein equations to Papapetrou.  This supports the historical core of
P248-S03.  It does not make the density a generally covariant local material
stress or supply the issue paper's optical constitutive map.

- official record: https://doi.org/10.1103/RevModPhys.36.463
- relevant source locations: equations (5)--(6) and the immediately following
  equivalence statement
- acquisition note: the publisher landing page and searchable archival copy
  were inspected; direct automated PDF retrieval was denied by the publisher

Papapetrou 1948 was located bibliographically, but a primary full text was not
available through the issue source or an authoritative open endpoint.  P248
therefore reconstructs the relaxed identity independently rather than treating
the citation as a proof attachment.

## Deser self-coupling theorem

Deser's primary paper derives Einstein dynamics from a free massless spin-two
field plus universal self-coupling in a first-order metric-density/connection
formulation.  It expressly avoids requiring harmonic gauge.  It supports the
conditional self-coupling bridge in P248-S23 once the massless spin-two field,
locality, and consistency hypotheses are declared.  It does **not** establish
that the issue paper's displayed Landau--Lifshitz expression is the unique
quadratic polynomial in its particular variables or that such a mode emerges
from the proposed substrate.

Audited copy of S. Deser, *Self-Interaction and Gauge Invariance*:

- primary author deposit: https://arxiv.org/abs/gr-qc/0411023
- SHA-256: `2b665d21441603c882e38ce8d8ae0462e960c804ecc820a73d1b2dea01cddbaf`
- pages: 9
- relevant source locations: abstract and metric-field derivation, equations
  (1)--(12)

## Bottom-up dependencies [16]--[21]

The accessible Zenodo dependencies were acquired and read at their claimed
load-bearing surfaces:

| Reference | SHA-256 | Scope finding |
| --- | --- | --- |
| [16], 1D sine--Gordon synthesis | `2af3e308a3df25adb0f1122c09c58254a64328b2048b1c74e4530630216dc263` | Contains exact homogeneous boosted-breather formulae, but its broad gravity and equivalence-principle language depends on later constitutive and collective-coordinate bridges. |
| [17], 1D metric/geodesic note | `2d9b1c82f940188c73334850f841e7b2597102fe209ff219c0f813e8039469cd` | The einbein/geodesic identities are valid for a supplied metric; the claimed square-root field reduction and material universality are separate premises. |
| [18], 2+1D lump/geometry note | `48208b52976fdee703f43a74cb25057bc9d426c92586be2753edc81278b9406c` | Its metric inverse, general anisotropic impedance, rigid-profile square-root, and universal equivalence-principle bridges require the repairs independently reproduced below. |
| [20], constitutive dictionary | `a286651a0a6b1dca3414ce55a0769c943bc94c0f6e83f1e0ea19bc96d71516c1` | Supplies the determinant-slaved nine-field optical map and calls it one-to-one; it cannot span a general ten-component foliation-adapted metric without an independent lapse. |

Reference [19], *Continuum Dynamics in 3+1 Dimensions*, has no URL in the issue
paper and no authoritative public record was located.  Its asserted
constitutive action, localized solutions, and equivalence-principle theorem
cannot enter the dependency closure.  Reference [21]'s einbein identities are
already independently available as accepted `C-WLN-001` and `C-WLN-002`, and
P248 replays their exact algebra rather than importing the note's physical
interpretation.

P238 previously audited reference [18].  P248 independently reproduces every
reused algebraic issue: determinant-only matching does not imply directional
anisotropic reflectionlessness; a rigid translated profile in a first-time-
derivative quadratic action cannot generate an all-orders square root; the
published inverse metric fails a direct matrix product; and a wave null cone
does not establish universal massive-lump coupling.  P238 remains provenance
for discovery, not scientific authority for the P248 verdict.

## Accepted dependency boundary

The exact positive completion uses only:

- `C-STG-001` for the accepted mostly-plus Einstein--scalar action and Euler
  equations;
- `C-WLN-001` for the massive supplied-metric worldline identity;
- `C-WLN-002` for the massless null/affine supplied-metric worldline identity;
- new P248 block-map theorems, reviewed before their composition.

`C-GOR-001`, `C-GOR-002`, `C-GRV-001`, `C-GRV-002`, and `C-IGR-001`--`004`
are scope comparators only.  They explicitly do not derive a microscopic
substrate action, a physical value of `G`, or general optical equivalence.
