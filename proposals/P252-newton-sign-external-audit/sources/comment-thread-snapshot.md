# Source snapshot: discussion #186 comments under audit (P252)

Provenance: fetched 2026-09-11 via GitHub GraphQL
(`repository(vantasnerdan/substrate-framework){discussion(number:186){comments}}`,
first:100, replies first:50), raw JSON retained off-tree at /tmp/d186_fresh.json.
Bodies below are extracted verbatim by jq; nothing hand-edited.

----- comment 18379326 by JarekDuda at 2026-09-10T03:18:05Z -----

Just talked separately with Opus about reversing Newton sign, and suggested completely new looking approach: replacing commutator with mixed on - in time changing to anticommutator:

commutator for the spatial/EM sector, anticommutator for the time/gravity sector

On the literature: the mixed structure exists, and the closest match is Savvidy's tensor gauge fields. Field strengths "antisymmetric in their first two indices and totally symmetric with respect to the rest" — exactly your structure. And notably his construction was built to "eliminate all negative norm states," yielding "two symmetric polarizations of helicity-two massless tensor gauge bosons and antisymmetric polarization of helicity-zero." The ghost problem you're facing is the one that literature exists to solve.

The clean statement is for the plain brackets, without ξ. For symmetric A, B: [A,B]ᵀ = −[A,B] — antisymmetric, an so(4) element, spin-1. And {A,B}ᵀ = {A,B} — symmetric, decomposing into trace (spin-0) plus traceless symmetric (spin-2). That part is solid.

With the ξ insertion the symbolic test returns "mixed", and that's a real complication rather than noise. A ξ B is not symmetric in the naive sense — the correct statement involves index position, M_αβ versus M_α^β, which your §V already flags as needing care. So the spin decomposition has to be done with indices raised consistently, and my one-line claim skipped that.

What survives, and I think it's the useful part: your F_μν = [∂_μM, ∂_νM] uses the bracket that produces odd spin, and odd spin repels like charges. The trace and traceless-symmetric parts of M — the spin-0 and spin-2 content, both attractive — never enter the field strength at all. The commutator discards exactly what gravity needs.

So the proposal is: keep [·,·] for the rotation sector where like-charge repulsion is wanted, and use {·,·} for the boost/time sector where attraction is. Two brackets, two spins, two signs — each where it belongs, and no sign reversed by hand, hence no ghost.

----- comment 18388685 by xrodz at 2026-09-10T15:41:50Z -----

# R19 planned: the mixed bracket checked at the form level before any run, and the one question that decides what it can change

## TL;DR

1. The proposal is read (the comment here and the list message). Before any relaxation we checked what it is algebraically (the script pinned below, 10 checks that can fail, audited by a second agent with its own code, 8 of 8 confirmed and one qualified): the anticommutator `G_mu nu = A_mu eta A_nu + A_nu eta A_mu` is symmetric in the derivative pair as well as the internal pair, so it is not a 2-form field strength and not Savvidy's structure ("antisymmetric in the first two indices, symmetric in the rest"); its square is a genuinely new operator, `GG = 4 PP - I1` with `PP` the norm of the product `A_mu eta A_nu`, outside everything R1 and R14 screened (rank 8 to 9 against the frozen basis).
2. "The time / gravity sector" has three readings and they differ in what they can change. Keyed on the time DERIVATIVE index, the action is identical to the certified one on every static configuration, so no static Newton read can move and it acts only on the clock. Keyed on the internal BOOST block, the Coulomb sector is unchanged identically (the G1 gate holds by construction) and the boost block becomes the M5.21.16 variant-A flip plus four times a new product norm. Applied everywhere, it changes the Coulomb sector too.
3. Every reading is quartic in omega on the rigid clock (the certified action is quadratic), so `J = 2 kin omega` and `E_J = E_stat + J^2 / 4 kin` have to be restated before any fixed-J number is quoted.
4. R19 as planned: R19-0 the registry entrants and R1's certificate re-run with the new columns (a coefficient window or an exact infeasibility, per reading); R19-1 the static Newton read on the boost-dressed pair (the R3 instrument) under the boost-block reading; R19-2 the clock under the quartic action, with the two-clock cross inertia under the time-index reading.
5. What we need from you before the go: which reading is meant (and whether the `(0, 0)` corner, `G_00 = 2 A_0 eta A_0`, belongs to the gravity sector), and if your own runs of the proposal exist, the script they ran so that the two sides compare the same object. We hold the run until then.

## 1. What the check settles

Notation: jets `A_mu = d_mu M`, `P_mu nu = A_mu eta A_nu`, the certified `F = P - P^T`, the anticommutator `G = P + P^T`, `<X, Y>_eta = sum_ab eta_a eta_b X_ab Y_ab`, and the densities with the FULL derivative contraction `(1/2) sum_{mu nu} eta^mu eta^nu <X_mu nu, X_mu nu>_eta` (the `sum_{mu < nu}` form is a complete contraction only for objects antisymmetric in the derivative pair; on `G` it is not boost invariant).

| Check | Result |
| --- | --- |
| C1 symmetries | `F` antisymmetric in both pairs; `G` symmetric in both; `F eta` is in `so(1, 3)` |
| C1 the "mixed" symbolic result | index placement: `(A eta)(B eta) + (B eta)(A eta)` is not symmetric, the covariant `A eta B + B eta A` is |
| C2 identity | `G^2 + F^2 = 2 (P^2 + P^T2)` entrywise, so `GG + I1 = 4 PP` to 2e-16; the `sum_{mu < nu}` half-sum of `G` is not even rotation invariant and is exactly quadratic in omega |
| C3 span | `GG` outside `{I1..I6, C6a, C6b}`: rank 8 to 9; exactly `GG = 2 J - I1` with `J = sum eta^mu eta^nu tr(eta A_mu eta A_nu eta A_nu eta A_mu)` the new direction |
| C4 invariance | `GG` Lorentz invariant to 6e-16; the fixed-frame boost-block split is rotation invariant, not boost invariant (0.96) |
| C5 time-index reading | equals `I1` exactly on every static configuration (`G_0mu = 0`) |
| C6 boost-block reading | equals `I1` exactly on the Coulomb sector (time row zero) |
| C7 boost-block reading | on static fields with a time row: `B = I1 - 2 I(F_t) + 4 I(P_t)` exactly, `I(F_t) = -sum F_ij[0, k]^2` the boost block of `I1` (a signature identity); on the unit sphere of jets the density's infimum is -0.5 for `I1` and -1.0 for the boost-block reading: indefinite in the same kind, twice as deep |
| C8 omega degree | `I1` quadratic; `GG`, both splits quartic, no odd powers; the quartic comes from the diagonal `G_00` block |

Two remarks the check supports. The spin statement is exact for the internal content: a jet proportional to `eta` drops out of `F` and not of `G`, while traceless-symmetric jets do enter `F` (so "the trace and traceless-symmetric parts never enter the field strength" is half right). And every term here is quartic in the jets, so it has no quadratic fluctuation about the constant vacuum (your § 33.1): the static interaction the record measures is between textures at nonlinear order (the boost dressing), not a linear exchange whose spin fixes the sign. Whether the even-spin content sets the sign of a texture-texture force is what R19-1 measures. One more fact from the audit: the space of covariant bilinears in `(A_mu, A_nu)` antisymmetric in the derivative pair and symmetric in the internal pair is one-dimensional, `tr(A_mu eta) A_nu - tr(A_nu eta) A_mu`, and it vanishes on eta-traceless jets; a Savvidy-type field strength on this field needs an epsilon, a fixed vector, or the field's `u`.

## 2. The three readings

| Reading | Static Newton reads (R0, R3, R11) | Coulomb sector | What it can change |
| --- | --- | --- | --- |
| (T) the time derivative index | identical | identical | the clock only (omega sector) |
| (B) the internal boost block, fixed frame | changes | identical | the boost dressing; not boost invariant |
| (B-u) the boost block on the field's timelike eigenvector `u` | changes | identical (`u = e_0` there) | covariant, field dependent, the case R1 left open |
| (A) everywhere | changes | changes | everything, a control |

## 3. R19 as planned

| Stage | What runs | Pre-registered outcomes |
| --- | --- | --- |
| R19-0 | the four objects enter the registry (sympy and numpy, complex-step gradient, a mutant that reddens; the C2, C5, C6 identities as gates); R1's channel certificate with the new columns: `H2` PSD, the `omega^4` sign, the boost-block sign, the lattice bump descent, at `c` in {0.25, 0.5, 0.75, 1} for (B) and (B-u) and at the literal (A), (T) | a window: `CANDIDATE` for R19-1; an exact infeasibility: `CANDIDATE_REFUTED` at the form level; (T) needs no static certificate |
| R19-1 | the R3 boost-dressed pair and the R11 same-sign instrument under (B) and (B-u), n32 L48 and n48 L72, d in {12, 18, 24, 30}: `E_int(d)`, the far-field fit, the block split, the dressing amplitude trend, the single dressed hedgehog's existence | `dE_int / dd > 0` with the dressing surviving: `NEWTON_SIGN_REVERSED` on that instrument; repulsive: `CANDIDATE_REFUTED`; unbounded or melted: `CANDIDATE_REFUTED` (no object) |
| R19-2 | the omega polynomial on the R10 core and the R12 ring, the fixed-J read with the quartic Legendre bridge (`J = 2 C omega + 4 D omega^3`), the box ladder L in {48, 72, 96}, the two-clock cross inertia (R3.iii) under (T) at d 12 and 24 | `D > 0` and `C`, `D` saturating: a localized fixed-J clock; `C ~ L`: `CANDIDATE_REFUTED`; `D < 0`: unbounded in omega |
| R19-3 | the R18 carry-overs if the day allows (pinned pairs at d 6 and 9, the converged 3D cores from the tube profiles) | reported as in R18 |

## 4. The list message to Faber

Of the three ideas there, the first is this proposal. An imaginary mass needs a term to carry it before it can run, and the near-constant-frequency reading needs a localized clock, which the record does not have (R4, R7, R10). Both are recorded as author-gated, not as candidates.

## 5. Code

| Item | Pinned |
| --- | --- |
| the pre-registration check (10 lines that can fail) and its audit (a second agent, own code) | [`m5_32_r19_precheck_audit.py`](https://github.com/openwave-labs/openwave/blob/23459f8a/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r19_precheck_audit.py), [`m5_32_r19_precheck.py`](https://github.com/openwave-labs/openwave/blob/23459f8a/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r19_precheck.py), its output [`m5_32_r19_precheck.json`](https://github.com/openwave-labs/openwave/blob/23459f8a/openwave/xperiments/m5_liquid_crystal/research/data/m5_32_r19_precheck.json) |
| the R19 proposal, every gate stated before the go | [candidate ledger § 6.8](https://github.com/openwave-labs/openwave/blob/23459f8a/openwave/xperiments/m5_liquid_crystal/research/findings/m5_32_candidate_ledger.md) |

The R18 open items (the pinned pairs at d 6 and 9, the converged 3D cores, the `I = int b^2` and phase-stiffness reads on a condensed core) stay queued behind your reply to the R18 post.


----- comment 18399375 by JarekDuda at 2026-09-11T10:24:26Z -----

# Time treatments and Newton: a structural closure, and where the one live question is

Report revision 249. Scripts `round238_exotic_time.py` … `round243_self_contract.py` in the bundle; §§ are that report's.

## TL;DR

1. **The obstruction to Newton was never the sign.** Three routes to the attractive sign exist in M5 and all three are short of a *protected massless mode coupled to something continuous*, not short of a sign.
2. **A time treatment can only touch one step of the closure** — that mass is the Noether charge of ℝ. Searching the groups rather than listing treatments: two choices break it, and both hit the next step.
3. **The route closes at the field content, not at the time.** In a static configuration `F_ij = [∂_i M, ∂_j M]` lives in the spatial 3×3 block and `g` lives in the (0,0) slot. They do not meet through η, through `M⁻¹`, or through any contraction. **No definition of "time" changes which slot `g` occupies.**
4. The covariant term that *would* couple them is the one Novarax tested (their § 12.3) and found **unbounded below under boosts**. We found the same instability independently in the boost sector, with **V supplying no restoring force**.
5. **One gravitational question in the report is still live**, and it is not about time.

---

## 1. Three routes to the attractive sign, all short of the same thing

| proposal | what it gives | what it lacks |
|---|---|---|
| **anticommutator** for the time sector | `{A,B}` is symmetric (spin 0 + 2) and **even spin attracts** — correct, and `[A,B]` does discard it | `tr({A,B}²) = tr([A,B]²) + 4tr(A²B²)`, already in the R1 basis; and `ε^{μνρσ}tr(∂M{∂M,∂M}) ≡ 0` identically (symmetric pair against antisymmetric ε) — **no topological current, no protected massless mode** |
| **imaginary axis** | M5's `η₀₀ = −1` already *is* an imaginary axis in the precise sense (an imaginary *axis* gives a negative real coefficient; an imaginary *coefficient* is a quarter-turn) — the time sector attracts for free | nothing: it was always there, and it never helped |
| **imaginary coupling** | for a *vector* mediator `g → ig` turns `+g²/r` into `−g²/r` — attractive | costs unitarity (`H = pσ_x + iμσ_z` is anti-Hermitian, `H² = (p²−μ²)I`, low-`p` modes grow); and applied to the dual photon it would break Coulomb repulsion |

**M5 has three routes to the attractive sign and none is short of a sign.** Every one is short of a protected massless mode coupled to something that is not quantised.

## 2. What a time treatment can touch

The closure chain is: **(A)** a medium's sources couple to strain, because uniform shifts cost nothing; **(B)** so the force is 1/d³ unless the charge is a surface integral; **(C)** a medium's surface charges are homotopy classes, hence quantised; **(D)** mass is the Noether charge of time translation, whose group is ℝ, hence continuous; **(E)** a continuous surface charge needs a 1-form symmetry.

**Time appears in exactly one step — (D).** So a time treatment is a choice of the group `G_t`, and (D) fails iff the choice makes the charge discrete or a surface term:

| `G_t` | the charge | verdict |
|---|---|---|
| ℝ | continuous Noether | (D) holds |
| U(1), compact time period T | quantised in 2π/T, still a **volume** charge | (C) blocks it; and the unit is < 10⁻⁴² eV |
| ℤ, discrete/Floquet | quasi-energy | not energy — gravity couples to `T^{μν}` |
| **Diff(ℝ), reparametrisation** | **a boundary term** | **(D) fails** — then (E) asks for a dynamical lapse |
| timeless (Wheeler–DeWitt) | a boundary term | same |
| ℝ × ℝ, two times (Bars) | needs Sp(2,ℝ) | that *is* (E): it adds a gauge structure |
| ℝ with a winding | π₁(RP²) = ℤ₂ | quantised, and unrelated to mass |
| **the clock's phase θ** | **K ∈ 2ℤ, quantised** | **(D) fails** — then (B): `K = ∫f²∂₀θ d³x` is a volume charge, 1/d³ |
| thermal/modular (Connes–Rovelli) | the modular Hamiltonian | needs a quantum state; M5 is classical |

Two break (D). Both hit the next step.

## 3. Why the survivor closes

Reparametrisation needs a dynamical lapse, and M5's only candidate is the time eigenvalue `g`.

It looked promising for exactly one section: on the locked frame `[diag(ġ,0,0,0), ∂_i M] = 0` identically, so **`g` has no kinetic term** — which is what a lapse looks like. Then the same algebra for spatial derivatives:

> **`F_ij = [S_i, S_j]` exactly.** No derivative of `g` — time or space — appears in the quartic. `g` appears only in `V`, and `dV/dg = 0` is algebraic with no coupling to anything else. **`g = 8` everywhere, source or no source.**

A lapse multiplies a constraint involving the *other* fields' spatial gradients. `g` multiplies nothing. A field with neither kinetic nor gradient term is not a gauge field — it is a number.

**And it is not a property of η.** Self-contracting instead, `tr(M⁻¹ F M⁻¹ F)` is `g`-independent too. The reason is structural: in a static configuration every derivative is spatial, so `F_ij` lives in the spatial block, and `g` lives in the (0,0) slot. **They do not meet through any contraction.**

The only two ways they can:
- **time derivatives** — dynamics, not a static potential;
- **a tilted frame** — `∂_i M` acquires (0,k) entries when the time eigenvector varies in space.

We computed the second. A localised radial boost χ(r) is **physical, not gauge** (§ 112's "pure gauge" was a *uniform* vacuum boost): it costs +0.2 at χ₀ = 0.01 with exponent ~2. With χ ≠ 0, `g` does acquire gradient stiffness, scaling as **χ^2.96** — so `g` couples to the medium's **velocity field**. That is a fluid's rest-energy structure, and it vanishes when nothing moves.

## 4. And the term that would fix it destroys the floor

GR's lapse couples through `N × (spatial gradient invariant)`. The covariant M5 analogue is `(u^μ M_μν u^ν) × I_spatial(u)` — which is **exactly the covariant clock completion f(y) in Novarax's § 12.3**, found there to be **unbounded below under boosts**, `E ≤ C − ck⁴ → −∞`, with a positive-definite velocity Hessian.

Independently, in the boost sector: beyond χ₀ ≈ 0.03–0.08 our quartic falls steeply (−79 at χ₀ = 0.08 against a base of −19), with the threshold moving to smaller χ₀ as the grid refines. And **V cannot help** — it is built from the spectrum of `N = Mη`, a boost acts as `N → BNB⁻¹` which preserves that spectrum, so **V is invariant under any boost, local or uniform, and supplies zero restoring force.**

That is now four independent sightings of one disease: § 99's runaway, § 22's hyperbolicity conditions, Novarax's unbounded admissible path, and the boost-sector sign. **Whether the certified action has an energy floor at all is upstream of every prediction in the report.**

## 5. The statement

> **The route through time is closed at the field content, not at the time.** Every treatment is a choice of `G_t`; only two break (D); both then need `g` to carry a field; `g` meets the static gradients through no contraction; and the term that would couple them destroys the energy floor.

This would survive any time treatment proposed next, because it is a property of `F = [∂M, ∂M]` with a locked time row, not of what "time" is taken to mean.

**Caveats.** The reduction to "a choice of `G_t`" covers the chain as written — a treatment that changed what *static* means would attack § 3 above rather than (D). And the block-structure argument uses the frame lock, which is an **exterior** finding: inside cores the frame tilts and `g` does meet the gradients, at short range, in a sector whose own stability is unresolved. The closure is of the **long-range** field.

## 6. Where the one live gravitational question is

Not in the time sector. In the matter sector, and it bypasses the medium theorem entirely because **a fermion bilinear is not a strain**:

1. **The condensed cores can be fermions.** π₄(RP²) = ℤ₂; a *bare* hedgehog is rotation-invariant at every angle to 10⁻¹⁶ (so the FR loop is constant and carries no sign), but the **index-4 condensate** breaks that invariance. The Finkelstein–Rubinstein sign costs nothing classical and θ = π is T-invariant.
2. **A Fermi point needs nodal pairing**, and the inter-core interaction is **dipolar** — § 187's `U = −2P₂(cos θ)/d³` for polarised cores. For identical spin-polarised fermions Pauli requires odd `l`, and `⟨1|U|1⟩ = −0.533/d³`: **attractive in p-wave**. That is Baranov et al.'s mechanism.
3. **The open step:** bulk 3D dipolar gases prefer `p_z`, which has a **line** node. Volovik needs `p_x + ip_y`, with **point** nodes — a chiral state that breaks time reversal. **M5's cores carry a clock, which § 190 identified as angular momentum along the director: T-odd.**

**Does a gas of clocked, spin-½, dipolar cores pair chirally?** It is a BCS calculation with a ³He literature and a binary answer. The honest prior is against it.

**Computed since drafting — and the prior held.** On a 7200-direction grid the dipolar kernel gives pairing eigenvalues `p_z = +32.12` (attractive) and `p_x + ip_y = −15.52` (**repulsive**); Baranov's result reproduced. `p_z` has a **line** node — codimension 2, no Fermi point. The chiral state would need a T-odd coupling worth **148 % of the entire dipolar interaction**, starting from the wrong sign, and the clock is not that. Every attractive channel found is `m = 0`.

**So the route closes at step 3.5.** The strength cancels in the ratio, so the ordering is independent of `f`; only a differently *shaped* interaction could change it. That was the last gravitational question in the report with a method and a computable answer, and the answer is no.

*Rev. 249 attached. §§ 218–222 the sign proposals; §§ 228–230 the time treatments; § 229 the boost sector; § 227 the pairing route. Prior parallel work: OpenWave R14–R18, the M5–teleparallel chain C1–C18 (reviewed in § 221), and the Novarax dossier (Braun AI Wiktorowicz, Rafał Wiktorowicz, Marcin Rosa).*


----- comment 18406555 by JarekDuda at 2026-09-11T21:13:26Z -----

# Time treatments and Newton: a structural closure, and where the one live question is

Report revision 249. Scripts `round238_exotic_time.py` … `round243_self_contract.py` in the bundle; §§ are that report's.

## TL;DR

1. **The obstruction to Newton was never the sign.** Three routes to the attractive sign exist in M5 and all three are short of a *protected massless mode coupled to something continuous*, not short of a sign.
2. **A time treatment can only touch one step of the closure** — that mass is the Noether charge of ℝ. Searching the groups rather than listing treatments: two choices break it, and both hit the next step.
3. **The route closes at the field content, not at the time.** In a static configuration `F_ij = [∂_i M, ∂_j M]` lives in the spatial 3×3 block and `g` lives in the (0,0) slot. They do not meet through η, through `M⁻¹`, or through any contraction. **No definition of "time" changes which slot `g` occupies.**
4. The covariant term that *would* couple them is the one Novarax tested (their § 12.3) and found **unbounded below under boosts**. We found the same instability independently in the boost sector, with **V supplying no restoring force**.
5. **One gravitational question in the report is still live**, and it is not about time.

---

## 1. Three routes to the attractive sign, all short of the same thing

| proposal | what it gives | what it lacks |
|---|---|---|
| **anticommutator** for the time sector | `{A,B}` is symmetric (spin 0 + 2) and **even spin attracts** — correct, and `[A,B]` does discard it | `tr({A,B}²) = tr([A,B]²) + 4tr(A²B²)`, already in the R1 basis; and `ε^{μνρσ}tr(∂M{∂M,∂M}) ≡ 0` identically (symmetric pair against antisymmetric ε) — **no topological current, no protected massless mode** |
| **imaginary axis** | M5's `η₀₀ = −1` already *is* an imaginary axis in the precise sense (an imaginary *axis* gives a negative real coefficient; an imaginary *coefficient* is a quarter-turn) — the time sector attracts for free | nothing: it was always there, and it never helped |
| **imaginary coupling** | for a *vector* mediator `g → ig` turns `+g²/r` into `−g²/r` — attractive | costs unitarity (`H = pσ_x + iμσ_z` is anti-Hermitian, `H² = (p²−μ²)I`, low-`p` modes grow); and applied to the dual photon it would break Coulomb repulsion |

**M5 has three routes to the attractive sign and none is short of a sign.** Every one is short of a protected massless mode coupled to something that is not quantised.

## 2. What a time treatment can touch

The closure chain is: **(A)** a medium's sources couple to strain, because uniform shifts cost nothing; **(B)** so the force is 1/d³ unless the charge is a surface integral; **(C)** a medium's surface charges are homotopy classes, hence quantised; **(D)** mass is the Noether charge of time translation, whose group is ℝ, hence continuous; **(E)** a continuous surface charge needs a 1-form symmetry.

**Time appears in exactly one step — (D).** So a time treatment is a choice of the group `G_t`, and (D) fails iff the choice makes the charge discrete or a surface term:

| `G_t` | the charge | verdict |
|---|---|---|
| ℝ | continuous Noether | (D) holds |
| U(1), compact time period T | quantised in 2π/T, still a **volume** charge | (C) blocks it; and the unit is < 10⁻⁴² eV |
| ℤ, discrete/Floquet | quasi-energy | not energy — gravity couples to `T^{μν}` |
| **Diff(ℝ), reparametrisation** | **a boundary term** | **(D) fails** — then (E) asks for a dynamical lapse |
| timeless (Wheeler–DeWitt) | a boundary term | same |
| ℝ × ℝ, two times (Bars) | needs Sp(2,ℝ) | that *is* (E): it adds a gauge structure |
| ℝ with a winding | π₁(RP²) = ℤ₂ | quantised, and unrelated to mass |
| **the clock's phase θ** | **K ∈ 2ℤ, quantised** | **(D) fails** — then (B): `K = ∫f²∂₀θ d³x` is a volume charge, 1/d³ |
| thermal/modular (Connes–Rovelli) | the modular Hamiltonian | needs a quantum state; M5 is classical |

Two break (D). Both hit the next step.

## 3. Why the survivor closes

Reparametrisation needs a dynamical lapse, and M5's only candidate is the time eigenvalue `g`.

It looked promising for exactly one section: on the locked frame `[diag(ġ,0,0,0), ∂_i M] = 0` identically, so **`g` has no kinetic term** — which is what a lapse looks like. Then the same algebra for spatial derivatives:

> **`F_ij = [S_i, S_j]` exactly.** No derivative of `g` — time or space — appears in the quartic. `g` appears only in `V`, and `dV/dg = 0` is algebraic with no coupling to anything else. **`g = 8` everywhere, source or no source.**

A lapse multiplies a constraint involving the *other* fields' spatial gradients. `g` multiplies nothing. A field with neither kinetic nor gradient term is not a gauge field — it is a number.

**And it is not a property of η.** Self-contracting instead, `tr(M⁻¹ F M⁻¹ F)` is `g`-independent too. The reason is structural: in a static configuration every derivative is spatial, so `F_ij` lives in the spatial block, and `g` lives in the (0,0) slot. **They do not meet through any contraction.**

The only two ways they can:
- **time derivatives** — dynamics, not a static potential;
- **a tilted frame** — `∂_i M` acquires (0,k) entries when the time eigenvector varies in space.

We computed the second. A localised radial boost χ(r) is **physical, not gauge** (§ 112's "pure gauge" was a *uniform* vacuum boost): it costs +0.2 at χ₀ = 0.01 with exponent ~2. With χ ≠ 0, `g` does acquire gradient stiffness, scaling as **χ^2.96** — so `g` couples to the medium's **velocity field**. That is a fluid's rest-energy structure, and it vanishes when nothing moves.

## 4. And the term that would fix it destroys the floor

GR's lapse couples through `N × (spatial gradient invariant)`. The covariant M5 analogue is `(u^μ M_μν u^ν) × I_spatial(u)` — which is **exactly the covariant clock completion f(y) in Novarax's § 12.3**, found there to be **unbounded below under boosts**, `E ≤ C − ck⁴ → −∞`, with a positive-definite velocity Hessian.

Independently, in the boost sector: beyond χ₀ ≈ 0.03–0.08 our quartic falls steeply (−79 at χ₀ = 0.08 against a base of −19), with the threshold moving to smaller χ₀ as the grid refines. And **V cannot help** — it is built from the spectrum of `N = Mη`, a boost acts as `N → BNB⁻¹` which preserves that spectrum, so **V is invariant under any boost, local or uniform, and supplies zero restoring force.**

That is now four independent sightings of one disease: § 99's runaway, § 22's hyperbolicity conditions, Novarax's unbounded admissible path, and the boost-sector sign. **Whether the certified action has an energy floor at all is upstream of every prediction in the report.**

## 5. The statement

> **The route through time is closed at the field content, not at the time.** Every treatment is a choice of `G_t`; only two break (D); both then need `g` to carry a field; `g` meets the static gradients through no contraction; and the term that would couple them destroys the energy floor.

This would survive any time treatment proposed next, because it is a property of `F = [∂M, ∂M]` with a locked time row, not of what "time" is taken to mean.

**Caveats.** The reduction to "a choice of `G_t`" covers the chain as written — a treatment that changed what *static* means would attack § 3 above rather than (D). And the block-structure argument uses the frame lock, which is an **exterior** finding: inside cores the frame tilts and `g` does meet the gradients, at short range, in a sector whose own stability is unresolved. The closure is of the **long-range** field.

## 6. Where the one live gravitational question is

Not in the time sector. In the matter sector, and it bypasses the medium theorem entirely because **a fermion bilinear is not a strain**:

1. **The condensed cores can be fermions.** π₄(RP²) = ℤ₂; a *bare* hedgehog is rotation-invariant at every angle to 10⁻¹⁶ (so the FR loop is constant and carries no sign), but the **index-4 condensate** breaks that invariance. The Finkelstein–Rubinstein sign costs nothing classical and θ = π is T-invariant.
2. **A Fermi point needs nodal pairing**, and the inter-core interaction is **dipolar** — § 187's `U = −2P₂(cos θ)/d³` for polarised cores. For identical spin-polarised fermions Pauli requires odd `l`, and `⟨1|U|1⟩ = −0.533/d³`: **attractive in p-wave**. That is Baranov et al.'s mechanism.
3. **The open step:** bulk 3D dipolar gases prefer `p_z`, which has a **line** node. Volovik needs `p_x + ip_y`, with **point** nodes — a chiral state that breaks time reversal. **M5's cores carry a clock, which § 190 identified as angular momentum along the director: T-odd.**

**Does a gas of clocked, spin-½, dipolar cores pair chirally?** It is a BCS calculation with a ³He literature and a binary answer. The honest prior is against it.

**Computed since drafting — and the prior held.** On a 7200-direction grid the dipolar kernel gives pairing eigenvalues `p_z = +32.12` (attractive) and `p_x + ip_y = −15.52` (**repulsive**); Baranov's result reproduced. `p_z` has a **line** node — codimension 2, no Fermi point. The chiral state would need a T-odd coupling worth **148 % of the entire dipolar interaction**, starting from the wrong sign, and the clock is not that. Every attractive channel found is `m = 0`.

**So the route closes at step 3.5.** The strength cancels in the ratio, so the ordering is independent of `f`; only a differently *shaped* interaction could change it. That was the last gravitational question in the report with a method and a computable answer, and the answer is no.

*Rev. 249 attached. §§ 218–222 the sign proposals; §§ 228–230 the time treatments; § 229 the boost sector; § 227 the pairing route. Prior parallel work: OpenWave R14–R18, the M5–teleparallel chain C1–C18 (reviewed in § 221), and the Novarax dossier (Braun AI Wiktorowicz, Rafał Wiktorowicz, Marcin Rosa).*


----- comment 18406566 by JarekDuda at 2026-09-11T21:15:12Z -----

Working original Lagrangian (to flip Newton sign) in https://arxiv.org/pdf/2108.07896
Updated M5 report: https://zenodo.org/records/22714918


# The gravity sector, reorganised: the sign was never the problem

Report revision 294. Scripts `round267`–`round285` in the bundle; §§ are that report's. Companion to the electron-sector note. This one records a reorganisation rather than a result: **forty-three sections asked how to flip a sign that was already correct**, and identifying that changes what the open question is.

---

## TL;DR

1. **The Newton sign is correct, in every sector that could carry a gravity-analogue** — and has been since the sign rule was established. The sector where like sources repel is electromagnetism, where repulsion is right.
2. **The anticommutator proposal was tested in the wrong sector.** The original suggestion — commutator for rotations, anticommutator for boosts — was stated correctly and then computed for the *spatial* sector. Redone in the boost sector it does give spin 0+2 and attraction.
3. **And M5's connection is teleparallel by construction**, so the EM and GEM sectors are coupled by group structure rather than by any Lagrangian term — with the Lorentzian signature flipping the cross term relative to the Euclidean case the paper's eq. (9) states.
4. **The flip exists, couples to the conserved stress tensor, and gives 1/d⁵**, because the boost mode is a Goldstone and Goldstones couple derivatively.
5. **The one structure that changes a power is the axion term**, and it reaches 1/d between charges.
6. **A quadratic GEM–EM coupling does reach neutral matter** — the neutrality objection I raised was wrong — and it couples to the EM *fraction* of mass, which varies 9× across the periodic table against an Eötvös bound of 10⁻¹⁵.
7. **Cosmological structure was tested four ways and reaches the energy floor, not Newton.**

---

## 1. The sign, sector by sector

The rule, verified on a 48³ lattice (`V_complete = −⟨J, K⁻¹J⟩`, ratio −1.000 at five separations): **scalar and tensor exchange attract; vector exchange repels.**

| sector | spin | like sources | verdict |
|---|---|---|---|
| the dual photon | 1 | repel | **correct — it is electromagnetism** |
| the amplitude / split | 0 | **attract** | **correct** |
| the clock / twist | 0 | **attract** | **correct** |
| the symmetric boost bilinear | 2 | **attract** | **correct** |
| `R^gg = Γ̃ × Γ̃` — the paper's | 1 | repel | wrong for gravity |

> **The sign is correct everywhere it matters. The only wrong one is the paper's own `R^gg`, and that is a choice of bracket rather than a sign.**

| requirement | status |
|---|---|
| **the sign** | **correct** |
| the range | wrong |
| universality | wrong — 7.4 % measured across four objects |
| the strength | wrong — 4·10⁴² unless separately suppressed |

**Three of four fail and the sign is not one of them.** Seven proposals — the anticommutator, the imaginary axis, the imaginary coupling, the NUT charge, the symmetric bilinear, the flipped Hamiltonian term — each found a way to get a sign that was never in question. **And "the sign is correct" is weaker than it sounds: any scalar attracts. That is a property of scalar exchange, not of M5.**

---

## 2. The anticommutator, tested in the right sector

The proposal as originally stated was: *keep `[·,·]` for the rotation sector where like-charge repulsion is wanted, and use `{·,·}` for the boost/time sector where attraction is.* **That is right, and the test I ran computed the spatial sector instead.**

Redone where it belongs: SO(1,3)'s boost generators are **symmetric** where SO(4)'s are antisymmetric, and for two symmetric generators `{A,B}` is symmetric — **spin 0 + 2** — while `[A,B]` is antisymmetric — spin 1. The paper's `R^gg = Γ̃ × Γ̃` is the cross product, so it keeps spin 1 and discards spin 2.

**In the spatial sector the antisymmetric part is the physical one**, because its charge is the π₂ winding. **In the boost sector there is no such winding** — symmetric generators have no cross-product charge to be topological about. **So the antisymmetric choice there is an analogy with electromagnetism, not a derivation, and the symmetric bilinear is equally available and gives the right spin.**

---

## 3. What the flip actually buys

**Every piece is present.** M5 has a conserved stress tensor by Noether's theorem — already used for the von Laue condition — and the symmetric bilinear `S^gg` is spin-2 and attractive. The vertex `S·T` can be written.

**And the mode is massless, for the reason that dooms it.** The boosts are broken by the vacuum, so `Γ̃` is a Goldstone — and a Goldstone couples derivatively. With `S ~ {Γ̃, Γ̃} ~ (∂M)²`, that is two derivatives at the vertex:

| derivatives at the vertex | 0 | 1 | **2** |
|---|---|---|---|
| static potential | **1/d** | 1/d³ | **1/d⁵** |

> **An attractive spin-2 mode coupled to the stress tensor. Right spin, right source, right sign, wrong range by four powers.**

To reach 1/d the vertex needs *zero* derivatives — a field coupling to T directly, with two derivatives in its kinetic term. **That is the Einstein–Hilbert structure, and M5 can build it**: `inc M = 2G⁽¹⁾[M]` means M's incompatibility *is* the linearised Einstein tensor, so `∫M·inc M` is the linearised EH action, built from M alone, vanishing on the vacuum, needing no new field.

**And it has m = 1, p = 1, so 2mp = 2 < 3: linearly divergent energy on a winding texture.**

| action | potential | core energy |
|---|---|---|
| the quartic alone | linear, confining | finite |
| **+ M·inc M** | **1/r — Newton** | **divergent** |
| + (inc M)² — the Stelle R² analogue | linear + Yukawa | finite, **and untested** |

> **M5 can have quadratic gravity. What it cannot have is the Einstein–Hilbert term that supplies the 1/r — and that is a theorem about winding, not a choice of coefficient.**

---

## 4. Coupling by orthogonality

`Γ_μ = Oᵀ∂_μO` is a **Weitzenböck connection** — pure gauge, zero Lorentz curvature, torsion carrying everything. **M5's O-field is teleparallel by construction**, and the EM and GEM sectors are tied by the Maurer–Cartan equation rather than by any interaction term. That is the "coupling by orthogonality" exactly.

**And the signature flips which combination:**

| generator | boosts | spatial block of `[Γ_μ, Γ_ν]` |
|---|---|---|
| SO(4) — the paper's eq. (6) | antisymmetric | `R^ee + R^gg` |
| **SO(1,3) — eq. (39)** | **symmetric** | **`R^ee − R^gg`** |

The paper states eq. (9) for the Euclidean generator and then replaces it at eq. (39). **Figure 6's "tendency for opposite curvatures" is the Euclidean statement; with boosts they have the same sign.** Every sign proposal asked how to flip a coupling by hand; **this one is flipped by η.**

*Caveat:* this verifies the Maurer–Cartan commutator, not the Hamiltonian, which uses `F = [∂M, ∂M]` with shape factors. The flip is established for the constraint and conjectural for the energy.

---

## 5. The one structure that changes a power

A Goldstone's shift symmetry forces the coupling through ∂θ. **The exception is a coupling to a topological density**, because the shift then changes the action by `c × Q` with Q an integer — the axion structure. M5 has every ingredient, and

> `S_θ = κ∫θ ρ_top d⁴x` gives `□θ = −κρ_top`, hence a Coulomb tail, hence **U = −κ²Q₁Q₂/4πd: inverse-square, attractive, long-range.**

Coulomb agrees with QED to ~10⁻¹², so κ/e < 10⁻⁶ — **a bound sitting at 4·10³⁰ times gravity's strength.** The channel could be gravity-strength and entirely invisible.

**It is proportional to Q₁Q₂.** I concluded from this that it vanishes for neutral bodies and that neutrality was the obstruction. **That was wrong** — see the next section.

---

## 6. The GEM–EM coupling does reach neutral matter

Matter is built of protons and neutrons, which carry plenty of charge even when the net is zero. **And the paper's GEM–EM coupling is `R^ee·R^gg` — quadratic in the EM curvature — so it couples to EM *energy density*, which is positive-definite and does not cancel when the charges do.**

So the question becomes what fraction of mass is EM energy:

| nucleus | E_EM/M |
|---|---|
| helium | 4.87·10⁻⁴ |
| iron | 2.44·10⁻³ |
| uranium | **4.44·10⁻³** |

**A spread of 9.1 across the periodic table**, giving an Eötvös parameter of order one between elements, against MICROSCOPE's 10⁻¹⁵. **Excluded by fifteen orders.**

**Which is the same failure as the NUT charge in a new channel**: that one tracks baryon number and misses the binding energy; this one tracks EM energy. Two couplings, one exclusion, because both see a *subset*.

> **Gravity couples to the total stress tensor. Any coupling to a part of it is composition-dependent at the level of that part's variation — ~10⁻³ for every part in the periodic table. Any partial coupling is excluded by ~10¹² unless the part it tracks *is* the whole thing to a part in 10¹⁵.**

**The obstruction is not neutrality. It is that every channel sees a fraction.** Our own δ-coupling was measured at a 7.4 % spread across four objects — 7·10¹³ times the bound.

---

## 7. Cosmological structure, tested four ways

**Nearly constant twist frequency.** For a massive phase field the uniform mode sits at the gap, so ω = m — the de Broglie relation, not an extra hypothesis. *This is now qualified:* the clock turned out to be massless (see the electron note), so "nearly constant frequency" means ω = K/I with K topologically fixed, which is structure the report already had.

**Least action over an eon.** For a rotor, `S(N) = 2π²IN²/T` grows as N², so free extremisation gives **N = 0** — a global variational principle *suppresses* the clock rather than producing it. A theta term fixes K = −θI/T, but |θ| ≤ π and I/T ~ 10⁻³⁹, so K ≤ 9·10⁻³⁹ against the ½ the electron needs.

> **Over an eon, dynamics wants K = 0 and topology forces K = ½. The clock survives because it is topologically protected, not because it is energetically favoured.** Which generalises: continuous charges relax to zero under a global principle, quantised topological ones cannot, and **the surviving structure at late times is the topologically protected part.**

**Extremising over the eon's length too.** ∂S/∂T = K²/2I > 0, so free extremisation drives **T → 0**. What opposes the collapse is the Friedmann equation, with **T ∝ GM** — so the eon's duration is a gravitational quantity, M5 has no G, and its longest intrinsic time scale is 10⁻²² s against an observed 4·10¹⁷ s. **Forty orders.**

**Compact space and the energy floor.** A closed universe does bound the unbounded-below witness (its own scope note says "open-space or increasing-box") — but the bound is attained by the negative phase *filling the universe*, so the certified vacuum is metastable at best. And a periodic time boundary condition **forbids the two *dynamical* runaways and neither static one.**

> **The cyclic structure's real content is the exclusion of runaway solutions — which matters, because two of the four sightings of the energy-floor problem are dynamical.** It does not reach Newton.

**A black hole as an inward-boost hedgehog.** This is exact rather than analogical: `v = tanh χ(r)` gives Painlevé–Gullstrand, and the outward boost is its time reverse — a white hole. M5 would have three horizons, one per branch. **But the boost sector destabilises at 1.85 % of c** (converged across resolutions) against the 71–100 % a horizon needs, **so M5's collapse is a phase change rather than a horizon** — and the Schwarzschild profile p = ½ sits exactly on the convergence boundary.

---

## 8. What the literature already said

Volovik's ³He-A result: the induced Einstein–Hilbert term **is** produced and is **swamped by non-diffeomorphism-invariant terms** from the superfluid's own Lagrangian. And Barceló–Liberati–Visser: Sakharov assumes the spacetime has **no prior dynamics**, while a medium's prior dynamics is *what makes* the effective metric.

> **The medium's own equations of motion are both the source of the effective metric and the obstruction to its being Einsteinian. They cannot be dropped — they make the metric. They cannot be kept — they are not diffeomorphism invariant.**

**M5's quartic is a contaminating term in exactly that sense.** Our closures are instances of a known obstruction; what is added here is a convergence theorem making it quantitative, a quantisation argument independent of diffeomorphism invariance, and an explicit name for the covariant term. **That is a narrower contribution than the report has been implying.**

---

## 9. The open question, restated

Not *how do we flip the sign* — nothing needed flipping. The question is:

> **Why is the attractive sector gapped or derivatively coupled, and what would carry a source that a neutral body has in proportion to *all* of its energy?**

Everything in §§145–275 stops at the second half. **Setting the transverse splitting to zero would restore the exact symmetry and remove the sector the model was built to have — that is the trade, in the model's own variables.**

---

## 10. What is worth someone else's time

- **`(inc M)²`, the Stelle R² analogue.** Admissible, convergent, built from M alone via the incompatibility identity, and **completely untested** — including its ghost content, which is generic for four-derivative terms. It is the one new term the survey of admissible terms left standing.
- **The source-contracted propagator** for the boost sector. Every force law in the last forty sections rests on derivative counting at the vertex rather than on a computed propagator. This is the handoff's G-02 and it is the only calculation that could still overturn a negative result.
- **The energy floor.** Four independent sightings — the coupled solve's runaway, the boost instability at 1.85 % of c, the unbounded-below witness, and the hyperbolicity conditions. **Nobody has computed the decay rate of the certified vacuum into the negative phase**, which is a Coleman bounce and decides whether the theory has a usable ground state.

*Rev. 294 attached. §§254–256 the teleparallel structure and the flip; §267 the sign; §§245, 261 the admissible terms; §§269–270 the axion term and neutral matter; §§236, 249–252 the cosmological tests; §247 the literature. Parallel efforts: OpenWave (R14–R18), the M5–teleparallel chain (C1–C18), Novarax (Braun AI Wiktorowicz, Rafał Wiktorowicz, Marcin Rosa), the 11 September technical handoff, and Mikulski's twelve reports.*


