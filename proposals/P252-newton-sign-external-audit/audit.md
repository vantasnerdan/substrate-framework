# P252 external referee audit: the 2026-09-11 "Newton sign flip" announcement

Object under audit: substrate-framework discussion #186, comment 18406566
(JarekDuda, 2026-09-11T21:15:12Z), report revision 294, announced by email
subject "Finally flipped Newton sign, finite nonzero frequency!".  Companion
context: comments 18399375/18406555 (rev 249 closure), 18379326
(anticommutator proposal), 18388685 (OpenWave R19 pre-registration).
Source snapshot: `sources/comment-thread-snapshot.md` (verbatim jq extraction,
MD5 in `sources/MD5SUMS`).

Framing: per the referee template the report is the last artifact.  Every
checkable claim below carries one independent oracle (derived from explicit
generator matrices, static-source field equations, and direct differentiation -
never from the source's own derivations) and one mutation that was observed to
break.  Modes: `--verify` and `--mutate` of `verify_newton_sign_audit.py`.

## Headline verdict (claim C12)

**The announcement's subject line is contradicted by the report's own
content.**  Revision 294's title is "The gravity sector, reorganised: the sign
was never the problem", and its section 1 states the sign "is correct, in
every sector that could carry a gravity-analogue - and has been since the sign
rule was established".  No sign was flipped; the discovery is that the sign
was never wrong.  Three of the four gravity requirements still fail
(range, universality, strength - the report's own table).  "Finite nonzero
frequency" is qualified in section 7: the clock turned out to be massless and
the surviving finite frequency is omega = K/I with K topologically fixed - the
clock survives by topological protection, not by energetic preference.  The
announcement therefore does not resolve the two issues it names ("wrong Newton
sign and diverging omega"): the sign was a sector-confusion artifact, and the
frequency is the already-known topological-protection structure.

## Claim ledger

| #  | Claim (atomic)                                                                     | Premises                              | Location (comment 18406566) |
|----|------------------------------------------------------------------------------------|----------------------------------------|-----------------------------|
| C1 | so(4) generators all antisymmetric; so(1,3) boost generators symmetric             | eta_00 = -1 (source's signature)       | sec. 2                      |
| C2 | for symmetric inputs {A,B} is symmetric (spin 0+2), [A,B] antisymmetric (spin 1)   | C1                                     | sec. 2                      |
| C3 | signature flips the cross term: boost x boost enters the rotation sector with opposite sign in so(1,3) vs so(4); "flipped by eta, not by hand" | C1                                     | sec. 4                      |
| C4 | sign rule: scalar and tensor exchange attract, vector exchange repels like charges | tree-level exchange, conserved sources | sec. 1 table                |
| C5 | tr({A,B}^2) = tr([A,B]^2) + 4 tr(A^2 B^2) (anticommutator square adds no new invariant) | none                                   | rev 249 sec. 1              |
| C6 | eps^{mu nu rh si} tr(X_mu {X_nu, X_rh}) = 0 identically (no topological current)   | none                                   | rev 249 sec. 1              |
| C7 | 1/d needs a zero-derivative vertex (Einstein-Hilbert structure)                     | C8                                     | sec. 3                      |
| C8 | vertex derivative count 0/1/2 gives static potentials 1/d, 1/d^3, 1/d^5            | multipole expansion                    | sec. 3 table                |
| C9 | inc M = 2 G^(1)[M]; M.inc M is the linearized EH action; winding sector makes its core energy divergent (2mp = 2 < 3) | linearized geometry                    | sec. 3                      |
| C10| Gamma_mu = O^T d_mu O is teleparallel: zero curvature, torsion carries everything  | O in SO(1,3)                           | sec. 4                      |
| C11| E_EM/M = 4.87e-4 (He), 2.44e-3 (Fe), 4.44e-3 (U); spread 9.1x vs Eotvos bound 1e-15 | liquid-drop Coulomb energy             | sec. 6                      |
| C12| headline: a sign flip happened and a finite nonzero frequency was won               | -                                      | email subject vs report     |
| C13| "ratio -1.000 at five separations" 48^3 lattice check of V_complete = -<J, K^-1 J> | their run bundle                       | sec. 1                      |

## Oracle table

| #  | Oracle (independent route)                                                                                                   | Mode       | Result |
|----|------------------------------------------------------------------------------------------------------------------------------|------------|--------|
| O1 | explicit so(4)/so(1,3) generator matrices; transpose properties; full commutator tables                                       | exact      | PASS   |
| O2 | spin decomposition: scalar commutes with rotations; traceless-symmetric part obeys the rank-2 transformation law with coefficients read off [J,K] | exact      | PASS   |
| O3 | signed structure triples in the eps convention: s([K,K]) = -s([J,J]) in so(1,3), = +s([J,J]) in so(4); Maurer-Cartan curvature of an explicit O(x,z) field vanishes identically | exact      | PASS   |
| O4 | static-source interaction energies from the Lagrangians (scalar -, vector +, spin-2 -); Green kernels exact off origin + FT numeric (quadosc, rel ~ 1e-16); route B: linearized field equation Newton limit | exact+numeric | PASS   |
| O5 | trace identity on general symbolic 4x4 entries and integer matrices n=2,3,5                                                    | exact      | PASS   |
| O6 | eps-contraction exact vanishing; commutator contrast generically nonzero                                                       | exact      | PASS   |
| O7 | dipole-dipole tail formula and trace-vanishing; quadrupole-quadrupole homogeneous of degree -5                                 | exact      | PASS   |
| O8 | linearized Einstein tensor validated by pure-gauge annihilation and the linearized Bianchi identity; Derrick scaling R^{3-2p} by exact change of variables | exact      | PASS   |
| O9 | inc identity from the comment alone; superseded by B10 (addendum): +2 G^(1) Euclid / -2 G^(1) Lorentz exact, not eps-fixable | exact      | SIGNATURE-CARRIED |
| O10| liquid-drop Coulomb fractions (Z^2 convention, mpmath 30 digits)                                                                | numeric    | PASS   |

## Mutation table (observed breaking)

| #   | Mutation                                                                  | Observed |
|-----|---------------------------------------------------------------------------|----------|
| M1  | wrong eta_00 sign breaks boost symmetry                                    | BROKE    |
| M2  | demanding so(4)'s same-sign relation in so(1,3) fails                      | BROKE    |
| M3  | mixed-symmetry inputs break the {A,B}-symmetry claim                       | BROKE    |
| M4  | trace-identity coefficient 4 -> 2 leaves a nonzero remainder               | BROKE    |
| M5  | symmetric contraction tensor in (nu,rh) breaks the eps-vanishing           | BROKE    |
| M6  | scalar coupling sign flip turns attraction into repulsion                  | BROKE    |
| M7  | antisymmetrized projector kills the Fierz-Pauli 00-00 contraction          | BROKE    |
| M8  | Newton-route sign flip turns the potential positive                        | BROKE    |
| M9  | quadrupole tail is not homogeneous of degree -3                            | BROKE    |
| M10 | Derrick exponent 3-p instead of 3-2p fails for p=2                          | BROKE    |
| M11 | a physical (non pure-gauge) perturbation does not annihilate G^(1)          | BROKE    |
| M12 | Lorentzian K's in the Euclidean slot remove the flip                        | BROKE    |
| M13 | Z(Z-1) Coulomb convention breaks the quoted He-4 number                     | BROKE    |

Full ledger: `attempts/0001/full_run.log` - "ALL 41 CHECKS PASS; ALL 15
MUTATIONS BREAK", exit 0 (main audit 37/13; addendum block B10 adds
4 checks and 2 mutations).

## Debt ledger (the flaws are the product)

| Debt | Introduced by | Disposition |
|------|---------------|-------------|
| D2 the "inc M = 2 G^(1)[M]" identity is asserted without an index convention in the comment; three independent double-eps arrangements failed to reproduce it against a pure-gauge-validated G^(1) | the comment (details live in rev 294 sections) | SHARPENED by B10 (addendum): the canonical double-eps operator gives +2 G^(1) in Euclid and -2 G^(1) in Lorentz exactly - the split is signature-carried, not eps-fixable (B10/C12c), and not re-slot-fixable (B10/M1). The Lorentzian +2 pin needs the round291/rev-294 operator definition; fetch Zenodo 22714918 before any framework use of the EH-analogue term |
| D3 the winding-divergence exponent ("2mp = 2 < 3") is a counting convention not present in the comment; the generic R^{3-2p} law is verified, the specific 2mp claim is not | the comment | PINNED at comment level by the disclosure (m = 1, p = 1 -> 2 < 3), consistent with the verified R^{3-2p} law at m = 1 (B10/C12d); rev-294 definitions still required for framework use |
| D4 the M5-internal mode identifications (dual photon = vector, amplitude/split and clock/twist = scalar, symmetric boost bilinear = spin-2) are source identifications, not audited here | the comment | recorded; any framework use requires its own claim work |
| D5 the 48^3 lattice sign check is not reproducible without the run bundle (Zenodo 22714918, scripts round267-round285) | the comment | provenance_only; not audited by proxy (AP-8/AP-14 gate) |
| D6 the flip is established at the constraint/Maurer-Cartan level only; the report's own caveat says it is "conjectural for the energy" - the Hamiltonian F = [dM, dM] still uses the commutator | the comment | the pre-registered energy-level test is OpenWave R19-1 with NEWTON_SIGN_REVERSED / CANDIDATE_REFUTED gates; gravity from M5 remains pending that run |
| D7 "finite nonzero frequency" rests on the clock being massless and topologically protected (omega = K/I) | the comment, sec. 7 | sharpened by the disclosure (source-asserted, pending rev-294 text): rev 294 withdrew its own sections 263-266 gap derivation in section 268; omega = K/I holds only at the split vacuum; at the degenerate vacuum the rotation is a stabiliser and there is no mode at all.  Consistent with this framework's certified P249/P250 clock structure; no new clock mechanism |

## Findings

1. **C1, C2, C3, C10 - VERIFIED (exact).**  The sector algebra is right:
   boost generators are symmetric, the symmetric bracket carries spin 0+2,
   the commutator spin 1, and the Lorentzian signature does flip the boost x
   boost structure relative to so(4) (signed triples: s(K,K) = -s(J,J) in
   so(1,3), s(K,K) = +s(J,J) in so(4); Maurer-Cartan curvature of an explicit
   field vanishes identically).  This is form-level: the report's own caveat
   that the flip is "conjectural for the energy" stands, and the pre-registered
   energy-level test is OpenWave R19-1 (debt D6).
2. **C4 - VERIFIED (exact + numeric, two routes).**  Scalar attract, vector
   repel, spin-2 attract, with the Fourier kernel verified exactly off origin
   and numerically to ~1e-16; the spin-2 attraction independently reproduced
   through the linearized field-equation route.
3. **C5, C6, C8 - VERIFIED (exact).**  The trace identity, the
   eps-contraction vanishing (with a genuinely nonzero commutator contrast),
   and the 1/d, 1/d^3, 1/d^5 multipole hierarchy all hold as stated.
4. **C9 - PARTIAL, sharpened by the addendum.**  The linearized Einstein
   tensor side is validated (pure-gauge annihilation, linearized Bianchi).
   Block B10 (addendum) proves the inc identity's sign is
   signature-carried: +2 G^(1) in Euclid, -2 G^(1) in Lorentz, exactly;
   the source's Lorentzian +2 pin needs its own operator definition
   (D2).  The generic Derrick scaling law R^{3-2p} is verified exactly,
   and the source's 2mp counting reduces to it at m = 1 (B10/C12d, D3).
5. **C11 - VERIFIED (numeric).**  The three EM mass fractions reproduce to
   better than 1% in the Z^2 liquid-drop convention and the 9.1x spread is
   right; the Z(Z-1) variant is excluded by the quoted numbers (mutation M13).
6. **C12 - REFUTED AS STATED.**  Nothing was flipped.  The honest status of
   "gravity from M5" after this report: the sign was never the obstruction;
   range (Goldstone derivative coupling, 1/d^5 - and their own section 10
   flags that this rests on power counting, not a computed propagator),
   universality (any partial stress-tensor coupling is composition-dependent
   and excluded far beyond Eotvos precision), and strength all remain failed,
   and the energy-level static sign is an open, pre-registered test (R19-1).

## Relation to framework canon

No registry claims are proposed or changed.  P236 (two-clock GEM Newton) and
P245 (nonlinear self-gravity) scopes are adjacent context and remain
unchanged.  If the framework later engages the EH-analogue term
(M . inc M) or the boost-bilinear coupling, that work must first pin the
rev-294 conventions (D2/D3) and run under its own preregistration.

## Addendum: disclosure received during the audit window (2026-09-12)

During the PR review window the auditee posted a pre-audit disclosure on
issue #211 (comment 5643293338, JarekDuda), agreeing with the headline
verdict and adding facts about rev 294 that the audited comment does not
carry.  Rev 294's full text was NOT attached in a form this audit could
fetch (no attachment links on the comment; no copy in the mailbox), so
everything below that cites the disclosure is recorded source-asserted,
pending rev-294 text.  Processing of the disclosure:

1. **C12 verdict confirmed and sharpened.**  Per the disclosure, rev 294's
   section 267 concludes nothing was flipped because nothing was wrong
   ("any scalar attracts; that is a property of scalar exchange, not of
   M5"), and section 266.5 states in so many words "Newton: not repaired".
   The disclosure also states that sections 263-266's clock-gap derivation
   (omega = m, delta = alpha^2) was WITHDRAWN in section 268: the potential
   depends only on trace invariants, so a uniform rotation costs zero
   potential energy; the clock is a massless Goldstone at the split vacuum
   and has NO mode at the degenerate vacuum used since section 2.  The
   omega = K/I statement is therefore valid only at the split vacuum.  This
   audit's verdict wording ("topological protection") is accordingly
   tightened: the announcement's "finite nonzero frequency" is a
   split-vacuum statement whose derivation-by-gap was retracted inside the
   report itself.
2. **B6 scope caveat (accepted).**  The verified multipole hierarchy (O7)
   is the general vertex-count ladder.  The report's sections 256/261
   extended the pattern to vertex structures it did not compute
   (S^gg.T - two derivatives on one leg against a stress tensor), and the
   report's own sections 256.6/261.5 flag this as power counting, not a
   computed propagator.  This audit does NOT validate that extension; the
   requested source-contracted-propagator calculation (the handoff's G-02)
   would settle it and is formally requested - but it is BLOCKED on the
   rev-294 vertex definitions (blocked-on-artifact, see below).
3. **Hamiltonian cross-term oracle requested (accepted, blocked).**  The
   disclosure correctly observes that section 254 verifies [Gamma_mu,
   Gamma_nu] (constraint level) while the Hamiltonian uses F_{mu nu} =
   [d_mu M, d_nu M] with R-decomposition shape factors (Lambda_i -
   Lambda_j).  An oracle on that object is worth having; it is BLOCKED on
   the rev-294 R-decomposition definition.
4. **inc-sign pin (B10, new oracle).**  The disclosure pins
   inc(h) = +2 G^(1)[h] with the standard-sign linearized Ricci
   ("almost certainly a Ricci convention rather than an error").  New
   oracle block B10 settles what can be settled from the comment alone:
   the natural double-eps operator inc_{mu nu}[h] = eps_{mu abc} eps_{nu rst}
   eta^{ar} d^b d^s h^{ct} equals +2 G^(1) EXACTLY in the Euclidean
   signature and -2 G^(1) EXACTLY in Lorentzian (-,+,+,+) with the
   mostly-minus Ricci (B10/C12a, C12b).  The +2/-2 split is therefore
   signature-carried - the same flip mechanism verified in B1/B8, now in
   the double-dual identity.  It is NOT fixable by the eps index convention
   (a double-eps product is invariant under a global eps flip,
   B10/C12c), and re-slotting the operator destroys the identity
   (B10/M1).  The source's Lorentzian +2 pin consequently requires the
   round291/rev-294 operator definition; both parallel efforts are
   internally consistent and no error is derivable from the comment alone
   (D2 stays open, sharpened).
5. **2mp counting pin (D3, resolved at comment level).**  The disclosure
   states the counting (m = 1, p = 1, 2mp = 2 < 3, divergent), which
   reduces to this audit's verified generic Derrick law R^{3-2p} at m = 1
   (B10/C12d).  D3 disposition upgraded from "unpinned" to "pinned at
   comment level; rev-294 definitions still required for framework use".
6. **Section 262.5 datum recorded.**  Per the disclosure: four of five
   parallel efforts found something rev 294 had wrong or incomplete; rev
   294 found errors in none of theirs.  Recorded as context on the
   report's reliability, source-asserted.
7. **Composition caveat (accepted).**  The C11 numbers are a uniform-sphere
   model ignoring exchange, surface diffuseness, and internal nucleon EM
   energy - order-of-magnitude support for the universality exclusion
   (twelve orders of margin), not a precision prediction.

### Blocked-on-artifact (formally requested, cannot proceed from INBOX-reachable material)

- G-02 source-contracted propagator for a two-derivative vertex against a
  stress tensor (settles sections 256/261) - needs the rev-294 vertex
  definitions.
- Hamiltonian cross-term flip oracle (F_{mu nu} shape factors) - needs the
  section 254 R-decomposition.

Both are one message away: the rev-294 bundle, if re-shared on issue #211
or pushed to the repository, unblocks them and this audit extends by one
slice under the same oracle discipline.

Final tally with B10: ALL 41 CHECKS PASS; ALL 15 MUTATIONS BREAK
(`attempts/0001/full_run.log`, exit 0).
